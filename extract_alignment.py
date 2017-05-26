#!/usr/bin/python
import sys
import time
import pickle
import os
import cPickle
import alignment
import hypergraph
from fragment_hypergraph import FragmentHGNode, FragmentHGEdge
import smatch
from smatch import get_amr_line
import amr_graph
from amr_graph import *
import logger
import re
from rule import Rule
from rule import retrieve_edges
from collections import defaultdict
from pos_processor import readPOSs

import gflags
FLAGS = gflags.FLAGS

#gflags.DEFINE_string(
#    'fragment_nonterminal',
#    'X',
#    'Nonterminal used for phrase forest.')
#gflags.DEFINE_bool(
#    'delete_unaligned',
#    False,
#    'Delete unaligned words in phrase decomposition forest.')
#
#FRAGMENT_NT = '[%s]' % FLAGS.fragment_nonterminal
def is_num(s):
    regex = re.compile('[0-9]+([^0-9\s]*)')
    match = regex.match(s)
    return match and len(match.group(1)) == 0

def filter_vars(line):
    fields = line.split('|||')
    parts = fields[2].split('/')
    for i in xrange(len(parts)):
        try:
            while parts[i][-1] in '0123456789':
                parts[i] = parts[i][:-1]
        except:
            print line
            print parts
    fields[2] = '/'.join(parts)
    return '|||'.join(fields)

def is_decades(s):
    regex = re.compile('[0-9]+([^0-9\s]*)')
    match = regex.match(s)
    return match and len(match.group(0)) == 5 and match.group(1) == 's'

#For each word, transform it into word sequence
def get_lemma_seq(wrd_seq, infl_lemma_map, trans_lemma_map, der_lemma_map):
    seqs = []
    words = wrd_seq.strip().split()
    for wrd in words:
        if wrd in infl_lemma_map:
            for lemma in infl_lemma_map[wrd]:
                seqs.append(lemma)
                break
        elif wrd in trans_lemma_map:
            for lemma in trans_lemma_map[wrd]:
                seqs.append(lemma)
                break
        elif wrd in der_lemma_map:
            for lemma in der_lemma_map[wrd]:
                seqs.append(lemma)
                break
        else:
            seqs.append(wrd)
    return ' '.join(seqs)

def get_wordseq_part(line):
    return line.strip().split('|||')[1].strip()

def replace_wordseq_part(line, new_seq):
    try:
        parts = line.split('##')
        parts[1] = ' %s ' % new_seq
    except:
        print 'what'
        #print new_seq
    return '##'.join(parts)

def get_seq_mapping(file, to_low = True):
    seq_map = {}
    has_A1 = set()
    has_ARG = set()
    with open(file, 'r') as f:
        for line in f:
            if line.strip() == '':
                continue
            wrd_seq = get_wordseq_part(line)
            if to_low:
                wrd_seq = wrd_seq.lower()
            if all_stop_wrds(wrd_seq):
                continue
            seq_map.setdefault(wrd_seq, set()).add(line.strip().replace('|||', '##'))

            if int(line[2]) == 1 and not (line[3] >= '0' and line[3] <= '9'):
                has_A1.add(wrd_seq)
            if 'ARG' in line and '-of' not in line:
                has_ARG.add(wrd_seq)
    #assert 'it' not in seq_map
    return (seq_map, has_A1, has_ARG)

def all_stop_wrds(wrd_seq):
    parts = wrd_seq.split()
    all_stop = True
    for part in parts:
        if part not in stop_words:
            all_stop = False
            break
    return all_stop



def build_lemma_map(seq_map, has_A1, has_ARG, infl_lemma_map, trans_lemma_map, der_lemma_map):
    lemma_to_rules = defaultdict(set)
    lemma_A1 = set()
    lemma_ARG = set()
    for wrd_seq in seq_map:
        lemma_seq = get_lemma_seq(wrd_seq, infl_lemma_map, trans_lemma_map, der_lemma_map)
        if wrd_seq == lemma_seq:
            continue
        if all_stop_wrds(lemma_seq):
            continue
        lemma_to_rules[lemma_seq] = seq_map[wrd_seq]
        if wrd_seq in has_A1:
            lemma_A1.add(lemma_seq)
        if wrd_seq in has_ARG:
            lemma_ARG.add(lemma_seq)
    return (lemma_to_rules, lemma_A1, lemma_ARG)

def retrieve_NN(word, der_lemma_map, lemma_to_rules, train_dict):
    #assert (word not in lemma_to_rules) and (word not in train_dict), 'Weird, why still unaligned:%s' % word
    if word not in der_lemma_map:
        return None

    der_lem = list(der_lemma_map[word])[0]
    if der_lem in lemma_to_rules or der_lem in train_dict:
        return der_lem
    return None

def retrieve_NNS(word, infl_lemma_map, der_lemma_map, lemma_to_rules, train_dict):
    #assert (word not in lemma_to_rules) and (word not in train_dict), 'Weird, why still unaligned:%s' % word
    if word not in infl_lemma_map:
        return None
    infl_lem = list(infl_lemma_map[word])[0]
    if infl_lem in lemma_to_rules or infl_lem in train_dict:
        return infl_lem
    if infl_lem not in der_lemma_map:
        return infl_lem
    der_lem = list(der_lemma_map[infl_lem])[0]
    if der_lem in lemma_to_rules or der_lem in train_dict:
        return der_lem
    return infl_lem

def to_written_rule(start, end, rule_str, str_side = None):
    return '%d-%d####%s' % (start, end, rule_str if str_side is None else replace_wordseq_part(rule_str, str_side))

def is_cap(s):
    return s.lower() != s

def check_entity(toks, n_tok, unaligned_toks):
    entity = toks[n_tok]
    start = n_tok
    end = n_tok+1
    #is_entity = True
    for i in xrange(n_tok-1, 0, -1):
        if is_cap(toks[i]):
            if toks[i] not in unaligned_toks: #Already been aligned
                return None
            elif toks[i] == '-':
                return None
            elif toks[i] == '\'':
                return None
            else:
                entity = '%s %s' % (toks[i], entity)
                start = i
        elif toks[i] == '-':
            return None
        elif toks[i] == '\'':
            return None
        else:
            break
    for i in xrange(n_tok+1, len(toks)):
        if is_cap(toks[i]):
            if toks[i] not in unaligned_toks: #Already been aligned
                return None

            elif toks[i] == '-':
                return None

            elif toks[i] == '\'':
                return None

            else:
                entity = '%s %s' % (entity, toks[i])
                end = i+1
        elif toks[i] == '-':
            return None
        elif toks[i] == '\'':
            return None
        else:
            return (entity, start, end)
    return None

def build_entity(Entity):
    s, start, end = Entity
    operands = s.split()
    graph_str = '(. :n/name'
    for i in xrange(len(operands)):
        graph_str += (' :op%d (. :"%s" )' % (i+1, operands[i]))
    graph_str += ')'
    return '%d-%d####[A1-1] ## %s ## %s' % (start, end, s, graph_str)

#Each rule in the format like:
# [A1-1] ## United States ## (. :c/country  :name (. :n/name  :op1 (. :"States" ) :op2 (. :"United" )))
#Jeff's aligner does not deal with the order well
def get_entity_name(rule_str):
    op_regex = re.compile(':op[0-9]+ \(\. :([^\(\)]*) \)')
    position = 0
    ops = []
    old_op_str = ''
    while position < len(rule_str):
        match = op_regex.search(rule_str, position)
        if not match:
            break
        curr_op = match.group(1)
        old_op_str += ' %s' % match.group(0)
        assert curr_op[0] == '"' and curr_op[-1] == '"'
        if curr_op[0] == '"' and curr_op[-1] == '"':
            curr_op = curr_op[1:-1]
        ops.append(curr_op)
        position = match.end()
    return (ops, old_op_str.strip())

def get_new_op_str(op_toks):
    result = ''
    for i in xrange(len(op_toks)):
        result += ' :op%d (. :"%s" )' % (i+1, op_toks[i])
    return result.strip()

def alike(str1, str2):
    min_len = min(len(str1), len(str2))
    if min_len > 4:
        return str1[:4].lower() == str2[:4].lower()
    if len(str1) > len(str2):
        return str1[:min_len].lower() == str2.lower()
    return str2[:min_len].lower() == str1.lower()

def single_alike(str1, str2):
    return str1[:2].lower() == str2[:2].lower()

def same_seq(seq1, seq2):
    for (s1, s2) in zip(seq1, seq2):
        if not alike(s1, s2):
            return False
    return True

def alignm(str1, str2):
    return '%s: %s' % (str1, str2)

def shortened_toks(toks):
    s_toks = []
    i = 0
    while i < len(toks):
        if i + 1 < len(toks) and toks[i+1] == '\'s':
            s_toks.append(toks[i]+'\'s')
            i += 2
            continue
        elif i + 2 < len(toks) and toks[i+1] == '-':
            s_toks.append(''.join(toks[i:i+3]))
            i += 3
            continue
        else:
            s_toks.append(toks[i])
            i += 1
    return s_toks

def reversed(tmp_entities, start, end, toks):
    #if end - start != len(entities):
    #    print ' '.join(entities)
    #    print ' '.join(toks[start:end])
    #    print ' '.join(toks[start:start+len(entities)])

    entities = [w for w in tmp_entities if w != '-']
    #First assume it's in the correct order
    #Find the first entity from start, then left to right
    #ltr = True
    if 'of Iran' in ' '.join(entities):
        #print alignm(' '.join(toks[start:end]), ' '.join(tmp_entities))
        return (start, end, entities)

    elif 'Iran of' in ' '.join(entities):
        tmp_entities.reverse()
        #print alignm(' '.join(toks[start:end]), ' '.join(tmp_entities))
        return (start, end, shortened_toks(tmp_entities))

    #Usually in this case toks will have to concatenate things
    elif '\'s' in ' '.join(entities) or '-' in ' '.join(entities):
        #correct order
        s_toks = shortened_toks(toks[start:])
        if s_toks[0] == entities[0]:
            #print alignm(' '.join(toks[start:end]), ' '.join(tmp_entities))
            return (start, end, shortened_toks(tmp_entities))
        else:
            tmp_entities.reverse()
            entities.reverse()
            entities = shortened_toks(entities)
            new_toks = shortened_toks(toks)
            (new_start, new_end) = subseq(entities, new_toks, 0)
            try:
                assert new_start != -1,'%s#%s#%s' % (' '.join(entities), ' '.join(toks[start:end]), ' '.join(new_toks))
                #assert s_toks[0] == entities[0], '%s#%s#%s' % (' '.join(entities), ' '.join(toks[start:end]), s_toks[0])
                #print alignm(' '.join(toks[start:end]), ' '.join(tmp_entities))
            except:
                entities[1] = entities[0] + '-' + entities[1]
                entities = entities[1:]
                (new_start, new_end) = subseq(entities, new_toks, 0)
                assert new_start != -1,'%s#%s#%s' % (' '.join(entities), ' '.join(toks[start:end]), ' '.join(new_toks))

            new_str = orig_str(new_toks[new_start:new_end])
            or_str = ' '.join(new_toks[new_start:new_end])

            length = len(new_str.split())
            whole_str = ' '.join(toks)
            try:
                s_index = whole_str.index(new_str)
            except:
                #print whole_str
                #print new_str
                #sys.exit(-1)
                try:
                    s_index = whole_str.index(or_str)
                    length = len(or_str.split())
                except:
                    t_toks = new_str.split()
                    t_toks[-3] = t_toks[-3] + t_toks[-2] + t_toks[-1]
                    t_toks = t_toks[:-2]
                    try:
                        s_index = whole_str.index(' '.join(t_toks))
                        length = len(t_toks)
                    except:
                        print ' '.join(t_toks)
                        print whole_str
                        sys.exit(-1)
                    #print or_str
                    #print whole_str

            new_start = len(whole_str[:s_index].split())
            #print alignm(' '.join(toks[new_start:new_start+length]), ' '.join(tmp_entities))

            #print alignm(' '.join(new_toks[new_start:new_end]), ' '.join(tmp_entities))
            return (new_start, new_start+length, shortened_toks(tmp_entities))

    pos = start
    #while pos < len(toks) and toks[pos] != entities[0]:
    while pos < len(toks) and not alike(toks[pos], entities[0]):
        pos += 1
    try:
        assert pos != len(toks), '%s (position: %d-%d) is not found in sentence: %s' % (' '.join(entities), start, end, '#'.join(toks))
    except:
        assert len(entities) == 1
        assert single_alike(toks[start], entities[0])
        return (start, start + 1, shortened_toks(tmp_entities))
        #print alignm(toks[start], ' '.join(tmp_entities))

    if len(entities) == 1:
        return (start, start + 1, shortened_toks(tmp_entities))
        #print alignm(toks[start], ' '.join(tmp_entities))

    search_start = pos
    pos += 1
    n_sign = 0
    while pos < len(toks) and pos - search_start- n_sign < len(entities):
        if toks[pos] == '-':
            pos +=1
            n_sign += 1
            continue

        if not alike(toks[pos], entities[pos-search_start-n_sign]):
            break
        pos += 1

    if pos-search_start-n_sign == len(entities):
        #print alignm(' '.join(toks[search_start:pos]), ' '.join(tmp_entities))
        return (search_start, pos, shortened_toks(tmp_entities))

    #Should be right to left
    entities.reverse()
    tmp_entities.reverse()
    #try:
    assert search_start-len(entities)+1 >= 0, '%s\n%s'% (' '.join(entities), ' '.join(toks))
    #except:
    #    if 'Islamic' in ' '.join(entities) and 'Iran' in ' '.join(entities):
    #        print ' '.join(entities)
    #        print ' '.join(toks[start:end])
    #        entities.reverse()
    #        return (start, start + len(entities))
    #assert ' '.join(entities) == ' '.join(toks[search_start-len(entities)+1: search_start+1])
    (new_start, new_end) = subseq(entities, toks, start)
    try:
        assert new_start != -1
    except:
        print ' '.join(tmp_entities), ' '.join(toks), start
        sys.exit(-1)
    return (new_start, new_end, shortened_toks(tmp_entities))
    #print alignm(' '.join(toks[new_start:new_end]), ' '.join(tmp_entities))
    #try:
    #    #if 'Iran of' in ' '.join(entities) and 'Islamic' in ' '.join(entities):
    #    #    print ' '.join(entities)
    #    #    print ' '.join(toks[start:end])
    #    #    entities.reverse()
    #    #    return (start, start + len(entities))
    #    assert same_seq(entities, toks[search_start-len(entities)+1: search_start+1]), '%s\n%s\n%s' % (' '.join(entities), ' '.join(toks[search_start-len(entities)+1: search_start+1]), ' '.join(toks))
    #except:
    #    tmp_start = subseq(entities, toks, start)
    #    try:
    #        assert tmp_start != -1,'%s\n%s'% (' '.join(entities), ' '.join(toks))
    #    except:
    #        assert '-' in ' '.join(entities) or ' ' in toks[start:start+len(entities)] or '-' in ' '.join(toks[search_start: search_start+len(entities)]),  '%s\n%s'% (' '.join(entities), ' '.join(toks))

    #        print 'here:', toks[start:start+len(entities)]
    #        return (start, start + len(entities))
    #    return (tmp_start, tmp_start+ len(entities))
    #return (search_start-len(entities)+1, search_start+1)

#to verify seq1 is subsequence of seq2 from start position of seq2
def subseq(seq1, seq2, start):
    #i = 0
    tmp = start
    while start < len(seq2):
        while start < len(seq2):
            if seq2[start] != seq1[0]:
                start += 1
            else:
                break
        if start == len(seq2):
            if tmp == 0:
                return (-1, -1)
            return subseq(seq1, seq2, 0)
        tmp_seq = seq2[start:]
        tmp_seq = [s for s in tmp_seq if s != '-']
        if ' '.join(seq1) == ' '.join(tmp_seq[:len(seq1)]):
            end_pos = start
            i = 0
            while i < len(seq1):
                if seq2[end_pos] == '-':
                    end_pos += 1
                elif seq2[end_pos] == seq1[i]:
                    end_pos += 1
                    i += 1
            return (start, end_pos)
        start += 1

    if tmp == 0:
        return (-1, -1)
    return subseq(seq1, seq2, 0)

#Return the original sequence separated by white space
def orig_str(new_toks):
    ret = ''
    for tok in new_toks:
        if tok.endswith('\'s'):
            ret += ' %s' % tok[:-2]
            ret += ' \'s'
        elif '-' in tok:
            parts = tok.split('-')
            assert len(parts) == 2
            ret += ' %s' % parts[0]
            ret += ' -'
            ret += ' %s' % parts[1]
        else:
            ret += ' %s' % tok
    return ret.strip()
#Before this, should filter with stop words file
#Situation 1: multiple ";"
#Situation 2: some aligned and some not
def extract_aligned_rules(amr_file, sent_file, align_file, train_dict_file, person_dict_file, grammar_file, infl_lemma_file, trans_lemma_file, der_lemma_file, pos_file):
    f1 = open(amr_file, 'r')
    f2 = open(align_file, 'r')
    f3 = open(sent_file, 'r')

    global stop_words
    stop_words = set()
    with open('../run-decoder/stop_words', 'r') as f:
        for line in f:
            stop_words.add(line.strip())

    pos_tag_seqs = readPOSs(pos_file)

    infl_lemma_map = initialize_lemma(infl_lemma_file)
    trans_lemma_map = initialize_lemma(trans_lemma_file)
    der_lemma_map = initialize_lemma(der_lemma_file)

    (train_dict, train_A1, has_ARG_train) = get_seq_mapping(train_dict_file)

    (lemma_to_rules, lemma_A1, lemma_ARG) = build_lemma_map(train_dict, train_A1, has_ARG_train, infl_lemma_map, trans_lemma_map, der_lemma_map)
    #assert 'authorizing' in lemma_to_rules

    for s in stop_words:
        assert s not in train_dict
        assert s not in lemma_to_rules

    (person_dict, _, _) = get_seq_mapping(person_dict_file, False)

    gram_f = open(grammar_file, 'w')

    amr_line = get_amr_line(f1)
    alignment_line = f2.readline()
    sent = f3.readline()

    unaligned_f = open('unaligned_words', 'w')
    sent_num = 0

    concept_voc_set = set()

    has_arg = set()
    while amr_line and amr_line != '':
        assert alignment_line != '', 'The following amr graph does not have a corresponding alignment:\n%s' % amr_line

        amr_graph = AMRGraph(amr_line)
        sent = sent.strip()
        #print sent
        amr_graph.set_sentence(sent.split())

        alignments = alignment_line.strip().split()
        alignments = [tuple(align.split('|')) for align in alignments]
        alignments = sorted(alignments, key=lambda align: int(align[0].split('-')[0]))
        aligned_fragments = {}

        curr_start = 0
        aligned_toks = set()

        toks = sent.split()
        n_toks = len(toks)

        word_pos_seq = pos_tag_seqs[sent_num]
        try:
            assert len(word_pos_seq) == len(toks)
        except:
            print sent
            print word_pos_seq

        frag_str_set = set()
        aligned_strs = []
        for align in alignments:
            s_side = align[0]
            f_side = align[1]
            s_start = s_side.split('-')[0]
            s_start = (int)(s_start)
            s_end = s_side.split('-')[1]
            s_end = (int)(s_end)
            #print (s_start, s_end)

            curr_start = s_end
            fragment = amr_graph.retrieve_fragment(f_side)

            fragment.set_span(s_start, s_end)
            aligned_fragments[(s_start, s_end)] = fragment
            if s_end - s_start == 1 and toks[s_start] in has_ARG_train:
                aligned_toks |= set(xrange(s_start, s_end))
                continue
            else:

                new_rule = Rule()
                new_rule.lhs = intern('[A1-1]')
                new_rule.f = fragment.str_list()
                visited_index = set()
                fragment.init_new_hgraph(new_rule.e, fragment.root, {}, visited_index, None, None, None, {}, {})
                frag_str = filter_vars(str(new_rule)).replace('|||', '##')

                if ':name' in frag_str:
                    assert ':op' in frag_str, frag_str
                    #print frag_str
                    (entities, old_str) = get_entity_name(frag_str)
                    #old_span = '%d-%d'
                    (s_start, s_end, op_toks) = reversed(entities, s_start, s_end, toks)
                    fragment.set_span(s_start, s_end)
                    new_rule.f = fragment.str_list()
                    frag_str = filter_vars(str(new_rule)).replace('|||', '##')

                    new_str = get_new_op_str(op_toks)
                    if old_str != new_str:
                        print '######'
                        print ' '.join(toks[s_start:s_end])
                        print frag_str
                        frag_str = frag_str.replace(old_str, new_str)
                        print frag_str

                aligned_toks |= set(xrange(s_start, s_end))
                frag_str_set.add(frag_str)
                edges = retrieve_edges(str(new_rule))
                if len(edges) <= 1 and ('/' not in edges[0]):
                    assert len(new_rule.f) == 1
                    if is_num(edges[0]) or '@@@@' in edges[0]:
                        if is_decades(new_rule.f[0]):
                            if '@@@@' in edges[0]:
                                edges[0] = edges[0].split('@@@@')[1]
                            aligned_strs.append('%d-%d####[A1-1] ## %s ## (. :d/date-entity  :decade (. :%s ))' % (fragment.start, fragment.end, new_rule.f[0], edges[0]))
                        else:
                            if '@@@@' in edges[0]:
                                edges[0] = edges[0].split('@@@@')[1]
                            aligned_strs.append('%d-%d####[A1-0] ## %s ## (. :quant (. :%s ))' % (fragment.start, fragment.end, new_rule.f[0], edges[0]))
                    else:
                        aligned_strs.append('%d-%d####[A1-1] ## %s ## (. :n/%s )' % (fragment.start, fragment.end, new_rule.f[0], edges[0]))
                    continue

                aligned_strs.append('%d-%d####%s' % (fragment.start, fragment.end, frag_str))


        #print aligned_fragments
        visited = set()
        num_comma = 0
        first_comma = -1
        #assert 'Sarkozy' in person_dict, 'What?'
        for start in xrange(n_toks):
            if start in visited:
                continue
            for length in xrange(n_toks+1, 0, -1):
                end = start + length
                if toks[start] == ';':
                    num_comma += 1
                    if first_comma == -1:
                        first_comma = start
                    break

                if length == 1 and toks[start] == ';':
                    num_comma += 1
                    if first_comma == -1:
                        first_comma = start
                    break

                if end > n_toks:
                    continue

                seq_str = ' '.join(toks[start:end])

                if length < 5:
                    no_align = True
                    #if length == 1 and toks[start] == 'Sarkozy':
                        #print 'get her'
                        #print toks
                    #assert start not in aligned_toks, 'what'
                    for n_pos in xrange(start, end):
                        if n_pos in aligned_toks:
                            no_align = False
                            break
                    if no_align and seq_str in person_dict:
                        for rule_str in person_dict[seq_str]:

                            aligned_strs.append('%d-%d####%s' % (start, end, rule_str))
                            aligned_toks |= set(xrange(start, end))


                fake = 0

                get_out = False
                if start == 0:
                    seq_str = seq_str.lower()

                lemma_str = get_lemma_seq(seq_str, infl_lemma_map, trans_lemma_map, der_lemma_map)

                if seq_str in train_dict:
                    assert seq_str not in stop_words
                    #aligned_toks |= set(xrange(start, end))
                    for rule_str in train_dict[seq_str]:

                        edges = retrieve_edges(rule_str)
                        if len(edges) <= 1 and '/' not in edges[0]:
                            continue

                        if len(edges) <=1 and seq_str in has_ARG_train:
                            continue

                        if seq_str in train_A1 and ':name' in rule_str and int(rule_str[2]) > 2:
                            continue

                        #if frag_str and rule_str == frag_str:
                        #    continue
                        if rule_str in frag_str_set:
                            continue
                        value = int(rule_str[2])
                        aligned_strs.append('%d-%d####%s' % (start, end, replace_wordseq_part(rule_str, ' '.join(toks[start:end]))))
                        if end - start > 1 and end-start >= value:
                            if ';' not in rule_str:
                                for i in xrange(start, end):
                                    if toks[i] in has_ARG_train:
                                        for tok_rule_str in train_dict[toks[i]]:
                                            edges = retrieve_edges(tok_rule_str)
                                            if len(edges) <=1:
                                                continue

                                            aligned_strs.append('%d-%d####%s' % (i, i+1, tok_rule_str))
                                visited |= set(xrange(start, end))
                                get_out = True
                                break

                        if value <= 2 and end - start > 1:
                            if ';' not in rule_str:
                                for i in xrange(start, end):
                                    if toks[i] in has_ARG_train:
                                        for tok_rule_str in train_dict[toks[i]]:
                                            edges = retrieve_edges(tok_rule_str)

                                            if len(edges) <=1:
                                                continue
                                            aligned_strs.append('%d-%d####%s' % (i, i+1, tok_rule_str))
                                visited |= set(xrange(start, end))
                                get_out = True
                                break
                    #if get_out:
                    #    break
                elif lemma_str in train_dict:
                    assert lemma_str not in stop_words
                    for rule_str in train_dict[lemma_str]:

                        edges = retrieve_edges(rule_str)
                        if len(edges) <= 1 and '/' not in edges[0]:
                            continue

                        if len(edges) <=1 and lemma_str in has_ARG_train:
                            continue

                        if lemma_str in train_A1 and ':name' in rule_str and int(rule_str[2]) > 2:
                            continue

                        if rule_str in frag_str_set:
                            continue
                        value = int(rule_str[2])
                        aligned_strs.append('%d-%d####%s' % (start, end, replace_wordseq_part(rule_str, ' '.join(toks[start:end]))))
                        if end - start == 1:
                            aligned_toks |= set(xrange(start, end))

                        if end - start > 1 and end-start >= value:
                            if ';' not in rule_str:
                                for i in xrange(start, end):
                                    if toks[i] in has_ARG_train:
                                        for tok_rule_str in train_dict[toks[i]]:
                                            edges = retrieve_edges(tok_rule_str)
                                            if len(edges) <=1:
                                                continue

                                            aligned_strs.append('%d-%d####%s' % (i, i+1, tok_rule_str))
                                        #if seq_str in train_A1 and ':name' in rule_str and int(rule_str[2]) > 2:
                                        #    continue
                                visited |= set(xrange(start, end))
                                get_out = True
                                break

                        if value <= 2 and end - start > 1:
                            if ';' not in rule_str:
                                for i in xrange(start, end):
                                    if toks[i] in has_ARG_train:
                                        for tok_rule_str in train_dict[toks[i]]:
                                            edges = retrieve_edges(tok_rule_str)

                                            if len(edges) <=1:
                                                continue
                                            aligned_strs.append('%d-%d####%s' % (i, i+1, tok_rule_str))
                                visited |= set(xrange(start, end))
                                get_out = True
                                break
                elif lemma_str in lemma_to_rules:

                    assert lemma_str not in stop_words
                    for rule_str in lemma_to_rules[lemma_str]:

                        edges = retrieve_edges(rule_str)
                        if len(edges) <= 1 and '/' not in edges[0]:
                            continue

                        if len(edges) <=1 and lemma_str in lemma_ARG:
                            continue

                        if lemma_str in lemma_A1 and ':name' in rule_str and int(rule_str[2]) > 2:
                            continue

                        if rule_str in frag_str_set:
                            continue
                        value = int(rule_str[2])
                        aligned_strs.append('%d-%d####%s' % (start, end, replace_wordseq_part(rule_str, ' '.join(toks[start:end]))))

                        if end - start == 1:
                            aligned_toks |= set(xrange(start, end))

                        if end - start > 1 and end-start >= value:
                            if ';' not in rule_str:
                                for i in xrange(start, end):
                                    lemma_i = get_lemma_seq(toks[i], infl_lemma_map, trans_lemma_map, der_lemma_map)
                                    if lemma_i in lemma_ARG:
                                        for tok_rule_str in lemma_to_rules[lemma_i]:
                                            edges = retrieve_edges(tok_rule_str)
                                            if len(edges) <=1:
                                                continue

                                            aligned_strs.append('%d-%d####%s' % (i, i+1, tok_rule_str))
                                visited |= set(xrange(start, end))
                                get_out = True
                                break

                        if value <= 2 and end - start > 1:
                            if ';' not in rule_str:
                                for i in xrange(start, end):
                                    lemma_i = get_lemma_seq(toks[i], infl_lemma_map, trans_lemma_map, der_lemma_map)
                                    if lemma_i in lemma_ARG:
                                        for tok_rule_str in lemma_to_rules[lemma_i]:
                                            edges = retrieve_edges(tok_rule_str)

                                            if len(edges) <=1:
                                                continue
                                            aligned_strs.append('%d-%d####%s' % (i, i+1, tok_rule_str))
                                visited |= set(xrange(start, end))
                                get_out = True
                                break

        if num_comma > 0:
            minor = 'jfdkajk'
            comma_symbol = 'jfkajkjljl'
            higher = 'jfakjljlll'
            if num_comma <= 2:
                minor = 'A%d' % (num_comma+1)
                comma_symbol = 'A%d' % (num_comma+2)
                higher = 'A%d' % (num_comma+3)
            else:
                #minor = 'A4'
                comma_symbol = 'A5'


        sent_num += 1
        #print aligned_toks
        unaligned_toks = set(xrange(n_toks)) - aligned_toks
        unaligned_nouns = ''
        processed = set()
        for tok in unaligned_toks:
            if tok in processed:
                continue
            if word_pos_seq[tok][1][:2] == 'NN':
                w = word_pos_seq[tok][0]
                assert w == toks[tok]
                if tok == 0:
                    w = w.lower()
                if (word_pos_seq[tok][1] == 'NN'):
                    if w.lower() == w and w != '-':
                        lemma = retrieve_NN(w, der_lemma_map, lemma_to_rules, train_dict)
                        if lemma: #Have found something
                            if lemma in train_dict:
                                for rule_str in train_dict[lemma]:
                                    aligned_strs.append(to_written_rule(tok, tok+1, rule_str, toks[tok]))
                                    processed.add(tok)
                                continue
                            elif lemma in lemma_to_rules:
                                for rule_str in lemma_to_rules[lemma]:
                                    aligned_strs.append(to_written_rule(tok, tok+1, rule_str, toks[tok]))
                                    processed.add(tok)
                                continue
                            elif (not w.endswith('ion')) and (not w.endswith('ist')) and w != '-' and (not w.endswith('er')):
                                aligned_strs.append('%d-%d####[A1-1] ## %s ## (. :%s/%s )' % (tok, tok+1, toks[tok], lemma[0], lemma))
                                processed.add(tok)
                                continue
                        elif (not w.endswith('ion')) and (not w.endswith('ist')) and w != '-' and (not w.endswith('er')):
                            aligned_strs.append('%d-%d####[A1-1] ## %s ## (. :%s/%s )' % (tok, tok+1, toks[tok], toks[tok][0], toks[tok]))
                            processed.add(tok)
                            continue

                    elif w.lower() != w:
                        valid = True
                        if tok > 0 and is_cap(toks[tok-1]):
                            valid = False
                        if tok < len(toks)-1 and is_cap(toks[tok+1]):
                            valid = False
                        if valid:
                            lemma = toks[tok].lower()
                            aligned_strs.append('%d-%d####[A1-1] ## %s ## (. :%s/%s )' % (tok, tok+1, toks[tok], lemma[0], lemma))
                            processed.add(tok)
                            continue

                if (word_pos_seq[tok][1] == 'NNS'):
                    nn_form = ''
                    if w in infl_lemma_map:
                        nn_form = list(infl_lemma_map[w])[0]
                    elif w.endswith('s'):
                        nn_form = w[:-1]
                    if nn_form != '':
                        if w.lower() == w and w != '-':
                            lemma = retrieve_NNS(w, infl_lemma_map, der_lemma_map, lemma_to_rules, train_dict)
                            if lemma: #Have found something
                                if lemma in train_dict:
                                    for rule_str in train_dict[lemma]:
                                        aligned_strs.append(to_written_rule(tok, tok+1, rule_str, toks[tok]))
                                        processed.add(tok)
                                    continue
                                elif lemma in lemma_to_rules:
                                    for rule_str in lemma_to_rules[lemma]:
                                        aligned_strs.append(to_written_rule(tok, tok+1, rule_str, toks[tok]))
                                        processed.add(tok)
                                    continue
                                elif (not w.endswith('ion')) and (not w.endswith('ist')) and w != '-' and (not w.endswith('er')):
                                    aligned_strs.append('%d-%d####[A1-1] ## %s ## (. :%s/%s )' % (tok, tok+1, toks[tok], lemma[0], lemma))
                                    processed.add(tok)
                                    continue
                        elif w != '-': #Upper case
                            assert w.endswith('s')
                            valid = True
                            if tok > 0 and is_cap(toks[tok-1]):
                                valid = False
                            if tok < len(toks)-1 and is_cap(toks[tok+1]):
                                valid = False

                            if valid:
                                graph_str = '(. :n/name :op1 (. :"%s" ))' % w[:-1]
                                aligned_strs.append('%d-%d####[A1-1] ## %s ## %s' % (tok, tok+1, toks[tok], graph_str))
                                processed.add(tok)
                                continue

                elif (word_pos_seq[tok][1][:3] == 'NNP'):
                    tmp_tok = ''
                    if word_pos_seq[tok][1] == 'NNPS' and toks[tok].endswith('s'):
                        #tmp_tok = toks[tok][:-1]
                        w = w[:-1]
                    if w.lower() == w:
                        lemma = retrieve_NNS(w, infl_lemma_map, der_lemma_map, lemma_to_rules, train_dict)
                        if lemma: #Have found something
                            if lemma in train_dict:
                                for rule_str in train_dict[lemma]:
                                    aligned_strs.append(to_written_rule(tok, tok+1, rule_str, toks[tok]))
                                    processed.add(tok)
                                continue
                            elif lemma in lemma_to_rules:
                                for rule_str in lemma_to_rules[lemma]:
                                    aligned_strs.append(to_written_rule(tok, tok+1, rule_str, toks[tok]))
                                    processed.add(tok)
                                continue
                            elif (not w.endswith('ion')) and (not w.endswith('ist')) and w != '-' and (not w.endswith('er')):
                                aligned_strs.append('%d-%d####[A1-1] ## %s ## (. :%s/%s )' % (tok, tok+1, toks[tok], lemma[0], lemma))
                                processed.add(tok)
                                continue
                    else: #Check if this capital word is an entity
                        if word_pos_seq[tok][1] == 'NNP':
                            Entity = check_entity(toks, tok, unaligned_toks)
                            if Entity:
                                aligned_strs.append(build_entity(Entity))
                                processed |= set(range(Entity[1], Entity[2]))
                                continue
                        else:
                            #assert w.endswith('s')
                            valid = True
                            if tok > 0 and is_cap(toks[tok-1]):
                                valid = False
                            if tok < len(toks)-1 and is_cap(toks[tok+1]):
                                valid = False

                            if valid:
                                graph_str = '(. :n/name :op1 (. :"%s" ))' % w
                                aligned_strs.append('%d-%d####[A1-1] ## %s ## %s' % (tok, tok+1, toks[tok], graph_str))
                                processed.add(tok)
                                continue

                unaligned_nouns += '%s/%s ' % (word_pos_seq[tok][0], word_pos_seq[tok][1])


            unaligned_f.write(toks[tok] + '\n')
            #try:
            #    assert word_pos_seq[tok][0] == toks[tok]
            #except:
            #    print sent
            #    print word_pos_seq
            #    print toks[tok]
            #    sys.exit(-1)

        #print unaligned_nouns
        unaligned_toks -= processed

        gram_f.write('%s ||| %s ||| %s\n' % (sent, ' '.join([str(k) for k in unaligned_toks]), '++'.join(aligned_strs)))

        #break
        amr_line = get_amr_line(f1)
        alignment_line = f2.readline()
        sent = f3.readline()

    f1.close()
    f2.close()
    f3.close()
    unaligned_f.close()
    gram_f.close()

#if __name__ == '__main__':
#    amr_file = sys.argv[1]
#    sentence_file = sys.argv[2]
#    alignment_file = sys.argv[3]
#    train_dict_file = sys.argv[4]
#    person_dict_file = sys.argv[5]
#    grammar_file = sys.argv[6]
#    #pos_file = sys.argv[7]
#
#    extract_aligned_rules(amr_file, sentence_file, alignment_file, train_dict_file, person_dict_file, grammar_file, sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
