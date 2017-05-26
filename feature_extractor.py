#!/usr/bin/python
import sys
import time
import pickle
import os
import cPickle
import alignment
import hypergraph
from fragment_hypergraph import FragmentHGNode, FragmentHGEdge
from amr_graph import *
from amr_utils import *
import logger
import gflags
from HRGSample import *
from rule import Rule, retrieve_edges
import argparse
from re_utils import *
from collections import defaultdict
from filter_stop_words import *
from lemma_util import initialize_lemma
from feature import *

FLAGS = gflags.FLAGS

gflags.DEFINE_string(
    'fragment_nonterminal',
    'X',
    'Nonterminal used for phrase forest.')
gflags.DEFINE_bool(
    'delete_unaligned',
    False,
    'Delete unaligned words in phrase decomposition forest.')
gflags.DEFINE_bool(
    'href',
    False,
    'Delete unaligned words in phrase decomposition forest.')
gflags.DEFINE_integer(
    'max_type',
    7,
    'Set the maximum attachment nodes each nontermial can have.')

FRAGMENT_NT = '[%s]' % FLAGS.fragment_nonterminal
p_rule_f = open('poisoned_rule.gr', 'w')
q_rule_f = open('another_poisoned_rule.gr', 'w')
#unalign_f = open('unaligned_info', 'w')
def filter_with_maxtype(curr_node):
    root_index = curr_node.frag.root
    ext_set = curr_node.frag.ext_set
    nonterm_type = len(ext_set) if root_index in ext_set else (len(ext_set) + 1)
    if nonterm_type > FLAGS.max_type:
        curr_node.set_nosample(True)

#Enlarge a chart with a set of items
def enlarge_chart(prev_chart, new_items):
    for node1 in new_items:
        flag = False
        for node2 in prev_chart:
            if node1.frag == node2.frag:
                for edge in node1.incoming:
                    node2.add_incoming(edge)
                    flag = True

        if not flag: #node1 has never appeared in the previous chart
            prev_chart.add(node1)

#Make sure this node has constant edge
def is_predicate(curr_node):
    concept_label = curr_node
    if len(concept_label.split('/')) < 2:
        return False
    label = concept_label.split('/')[-1]
    if '-' not in label:
        return False
    number = label.split('-')[1].strip()
    return number.isdigit()

#Add one item to a chart
def add_one_item(prev_chart, item):
    flag = False
    for node in prev_chart:
        if node.frag == item.frag:
            for edge in item.incoming:
                node.add_incoming(edge)
                flag = True
    if not flag:
        prev_chart.add(item)

#To verify if a chart item has covered the graph
#i.e. covered all nodes and all edges
def is_goal_item(chart_item):
    fragment = chart_item.frag
    nodes = fragment.nodes
    edges = fragment.edges
    return len(edges) == edges.count()
    #return (len(nodes) == nodes.count()) and (len(edges) == edges.count())

def initialize_edge_alignment(aligned_fragments, edge_alignment):
    for frag in aligned_fragments:
        edge_alignment |= frag.edges

#This method output all the unaligned node information of the current AMR graph
def output_all_unaligned_nodes(edge_alignment, amr_graph):
    un_seq = []
    un_nodes = []
    for i in xrange(len(amr_graph.nodes)):
        curr_node = amr_graph.nodes[i]
        c_edge_index = curr_node.c_edge
        if edge_alignment[c_edge_index] == 0: #Found a concept that is not aligned
            un_seq.append(str(curr_node))
            un_nodes.append(curr_node)
    #print >> unalign_f, ' '.join(un_seq)
    return un_nodes

def output_all_unaligned_edges(edge_alignment, amr_graph):
    for i in xrange(len(edge_alignment)):
        if edge_alignment[i] == 0:
            un_seq.append(str(amr_graph.edges[i]))
    #print >> unalign_f, ' '.join(un_seq)

#Build one node that would cover unaligned edges going out of a node
def enlarge_frag(curr_frag, amr_graph, edge_alignment):
    curr_node_index = curr_frag.root
    curr_graph_node = amr_graph.nodes[curr_node_index]
    logger.writeln(str(curr_frag))

    if edge_alignment[curr_graph_node.c_edge] == 0:
        frag_label = str(curr_frag).replace(' ', '')
        return (curr_frag, False, frag_label)

    #To remember the unaligned relation going out of each entity
    root_arcs = []
    head_arcs = []
    visited = set()

    is_pred = False #Use for deciding the category of the root node
    is_op = False

    #Dealing with parent edges, ARGx-of
    if len(curr_graph_node.p_edges) > 0:
        for curr_edge_index in curr_graph_node.p_edges:
            curr_edge = amr_graph.edges[curr_edge_index]
            edge_label = curr_edge.label

            if edge_alignment[curr_edge_index] == 1: #This edge has already been aligned
                if curr_frag.edges[curr_edge_index] == 1 and (edge_label[:3] == 'ARG' and 'of' in edge_label):
                    logger.writeln("what the hell is this")
                    logger.writeln(str(curr_frag))
                    is_pred = True
                continue

            #tail_node_index = curr_edge.tail

            #Our intuition: ARGs and ops goes with the root
            if (edge_label[:3] == 'ARG' and 'of' in edge_label):
                is_pred = True
                head_arcs.append((curr_edge_index, curr_edge.head))

    if len(curr_graph_node.v_edges) > 0:
        for curr_edge_index in curr_graph_node.v_edges:
            curr_edge = amr_graph.edges[curr_edge_index]
            edge_label = curr_edge.label

            if edge_alignment[curr_edge_index] == 1: #This edge has already been aligned
                if curr_frag.edges[curr_edge_index] == 1 and is_root_arc(edge_label): #Special case, there is already args attached
                    if 'ARG' in edge_label:
                        is_pred = True
                        logger.writeln("what the hell is this")
                        logger.writeln(str(curr_frag))
                    else:
                        is_op = True
                continue

            #Our intuition: ARGs and ops goes with the root
            if is_root_arc(edge_label):
                if 'ARG' in edge_label:
                    is_pred = True
                else:
                    assert 'op' in edge_label
                    is_op = True
                root_arcs.append((curr_edge_index, curr_edge.tail))

    n_nodes = len(amr_graph.nodes)
    n_edges = len(amr_graph.edges)
    frag = AMRFragment(n_edges, n_nodes, amr_graph)
    if len(root_arcs) > 0 or len(head_arcs) > 0:
        #assert len(head_arcs) < 2, str(amr_graph)
        frag.set_root(curr_node_index)
        for rel_index, tail_index in root_arcs:
            edge_alignment[rel_index] = 1
            frag.set_edge(rel_index)
            frag.set_node(tail_index)

        if len(head_arcs) > 0:
            (rel_index, head_index) = head_arcs[0]
        #for rel_index, head_index in head_arcs:
            edge_alignment[rel_index] = 1
            frag.set_edge(rel_index)
            frag.set_root(head_index)

        frag.build_ext_list()
        frag.build_ext_set()
        frag_label = str(frag).replace(' ', '')

        new_frag = combine_fragments(curr_frag, frag)
        return (new_frag, is_pred, frag_label)
    else:
        return (curr_frag, is_pred, str(curr_frag).replace(' ', ''))

def build_bimap(tok2frags):
    frag2map = defaultdict(set)
    index2frags = defaultdict(set)
    for index in tok2frags:
        for frag in tok2frags[index]:
            index2frags[index].add(frag)
            frag2map[frag].add(index)
            #matched_list = extract_patterns(str(frag), '~e\.[0-9]+(,[0-9]+)*')
            #matched_indexes = parse_indexes(matched_list)
            #for matched_index in matched_indexes:
            #    frag2map[frag].add(matched_index)
    return (index2frags, frag2map)

#Here we try to make the tok to fragment mapping one to one
def rebuild_fragment_map(tok2frags):
    (index2frags, frag2map) = build_bimap(tok2frags)
    for index in tok2frags:
        if len(tok2frags[index]) > 1:
            new_frag_list = []
            min_frag = None
            min_length = 100
            for frag in tok2frags[index]:
                index_set = frag2map[frag]
                assert index in index_set
                if len(index_set) > 1:
                    if len(index_set) < min_length:
                        min_length = len(index_set)
                        min_frag = frag
                    index_set.remove(index)
                else:
                    new_frag_list.append(frag)
            if len(new_frag_list) == 0:
                assert min_frag is not None
                new_frag_list.append(min_frag)
            tok2frags[index] = new_frag_list
    return tok2frags

def extract_fragments(alignments, amr_graph):
    #alignments = s2g_alignment.strip().split()
    tok2frags = defaultdict(list)

    num_nodes = len(amr_graph.nodes)
    num_edges = len(amr_graph.edges)

    op_toks = []
    role_toks = []
    for curr_align in reversed(alignments):
        curr_tok = curr_align.split('-')[0]
        curr_frag = curr_align.split('-')[1]

        span_start = int(curr_tok)
        span_end = span_start + 1

        (index_type, index) = amr_graph.get_concept_relation(curr_frag)
        frag = AMRFragment(num_edges, num_nodes, amr_graph)
        if index_type == 'c':
            frag.set_root(index)
            curr_node = amr_graph.nodes[index]

            #Extract ops for entities
            if len(curr_node.p_edges) == 1:
                par_edge = amr_graph.edges[curr_node.p_edges[0]]
                if 'op' in par_edge.label:
                    op_toks.append((span_start, curr_node.c_edge))

            if curr_node.is_entity():
                role_toks.append((span_start, curr_node.c_edge))

            frag.set_edge(curr_node.c_edge)

        else:
            frag.set_edge(index)
            curr_edge = amr_graph.edges[index]
            frag.set_root(curr_edge.head)
            frag.set_node(curr_edge.tail)

            if 'ARG' in curr_edge.label:
                continue

        frag.build_ext_list()
        frag.build_ext_set()

        tok2frags[span_start].append(frag)

    for index in tok2frags:
        if len(tok2frags[index]) > 1:
            tok2frags[index] = connect_adjacent(tok2frags[index], logger)

    tok2frags = rebuild_fragment_map(tok2frags)
    for index in tok2frags:
        for frag in tok2frags[index]:
            frag.set_span(index, index+1)

    return (op_toks, role_toks, tok2frags)

#Verify this fragment contains only one edge and return it
def unique_edge(frag):
    assert frag.edges.count() == 1, 'Not unify edge fragment found'
    amr_graph = frag.graph
    root_node = amr_graph.nodes[frag.root]
    for edge_index in root_node.v_edges:
        if frag.edges[edge_index] == 1:
            return edge_index
    assert True, 'This is impossible'
    return None

def construct_forest(args):
    amr_graphs = None
    if args.amr_file:
        amr_graphs = load_amr_graphs(args.amr_file)

    if args.parallel:
        global rule_f
        run_nodes = args.nodes.split('+')
        if amr_graphs: #This means the amr graphs are not split yet
            for (i, curr_node) in enumerate(run_nodes):
                curr_graph_file = os.path.join(args.data_dir, 'graph_%d' % i)
                #curr_graph_file = os.path.join(args.data_dir, 'new_graph_%d' % args.sub_id)
                curr_f = open(curr_graph_file, 'wb')
                curr_start = args.sent_per_node * i
                curr_end = args.sent_per_node * (i+1)
                cPickle.dump(amr_graphs[curr_start:curr_end], curr_f)
                curr_f.close()

            if args.dump_graph:
                print 'A total of %d AMR graphs' % len(amr_graphs)
                return

        os.system('rm -rf %s' % args.save_dir)
        os.mkdir(args.save_dir)

        for (i, curr_node) in enumerate(run_nodes):
            cmd = 'python %s --data_dir %s --save_dir %s --sub_id %d --lemma %s --stop %s' % (sys.argv[0], args.data_dir, args.save_dir, i, args.lemma, args.stop)
            print 'start to launch program in %s' % curr_node
            alignment_log_file = os.path.join(args.save_dir, 'alignment_log_%d' % i)

            os.system(r'ssh %s "cd %s; nohup %s >& %s" &' % (curr_node, os.getcwd(), cmd, alignment_log_file))

    else:
        sys.setrecursionlimit(sys.getrecursionlimit() * 30)

        logger.file = open(os.path.join(args.save_dir, 'logger_%d' % args.sub_id), 'w')
        amr_graph_file = os.path.join(args.data_dir, 'graph_%d' % args.sub_id)
        #amr_graph_file = os.path.join(args.data_dir, 'new_graph_%d' % args.sub_id)
        alignment_file = os.path.join(args.data_dir, 'align_%d' % args.sub_id)
        sent_file = os.path.join(args.data_dir, 'sent_%d' % args.sub_id)
        lemma_file = os.path.join(args.data_dir, 'lemma_%d' % args.sub_id)
        pos_file = os.path.join(args.data_dir, 'pos_%d' % args.sub_id)

        feature_file = os.path.join(args.save_dir, 'feature_%d' % args.sub_id)
        used_tok_file = os.path.join(args.save_dir, 'tok_%d' % args.sub_id)
        used_lemma_file = os.path.join(args.save_dir, 'lemma_%d' % args.sub_id)
        used_pos_file = os.path.join(args.save_dir, 'pos_%d' % args.sub_id)

        predicate_words = set()
        predicate_lemmas = set()
        predicate_labels = set()

        entity_words = set()
        entity_lemmas = set()
        entity_map = defaultdict(set)

        non_pred_map = defaultdict(set)
        non_pred_lemma_map = defaultdict(set)

        unknown_set = set()

        word2predicate = {}  #To a defaultdict
        ent2frag = {}   #
        arg2frag = {}
        label2frag = {}

        amr_f = open(amr_graph_file, 'rb')
        amr_graphs = cPickle.load(amr_f)
        amr_f.close()

        s2g_map = []
        with open(alignment_file, 'r') as align_f:
            for line in align_f:
                s2g_map.append(line.strip())
            align_f.close()

        sents = []
        with open(sent_file, 'r') as sent_f:
            for line in sent_f:
                sents.append(line.strip())
            sent_f.close()

        tokenized = []
        with open(lemma_file, 'r') as lemma_f:
            for line in lemma_f:
                tokenized.append(line.strip())
            lemma_f.close()

        poss = []
        with open(pos_file, 'r') as pos_f:
            for line in pos_f:
                poss.append(line.strip())
            pos_f.close()

        assert len(amr_graphs) == len(s2g_map) and len(amr_graphs) == len(sents) and len(sents) == len(poss), '%d %d %d' % (len(amr_graphs), len(s2g_map), len(sents))

        constructed_forests = []
        sent_indexes = []

        lemma_map = initialize_lemma(args.lemma)

        global print_sign
        stop_words = set([line.strip() for line in open(args.stop, 'r')])

        num_self_cycle = 0
        used_sents = 0

        feature_f = open(feature_file, 'w')
        tok_f = open(used_tok_file, 'w')
        lemma_f = open(used_lemma_file, 'w')
        pos_f = open(used_pos_file, 'w')

        label_set = set()
        for (sent_index, (sent, lemma_sent, tags, s2g_alignment, amr_graph)) in enumerate(zip(sents, tokenized, poss, s2g_map, amr_graphs)):

            logger.writeln('Sentence #%d' % sent_index)
            logger.writeln(sent)
            logger.writeln(str(amr_graph))
            logger.writeln(s2g_alignment)

            print_sign = False

            pos_seq = tags.split()

            edge_alignment = bitarray(len(amr_graph.edges))
            if edge_alignment.count() != 0:
                edge_alignment ^= edge_alignment

            assert edge_alignment.count() == 0
            if amr_graph.check_self_cycle():
                num_self_cycle += 1
                logger.writeln('self cycle detected')

            if s2g_alignment == '':
                logger.writeln('totally unaligned')
                continue

            toks = sent.split()
            lemmas = lemma_sent.split()
            amr_graph.set_sentence(toks)
            aligned_fragments = []

            reentrancies = {}  #Map multiple spans as reentrancies, keeping only one as original, others as connections

            has_multiple = False
            no_alignment = False

            aligned_set = set()

            all_frags = []

            (opt_toks, role_toks, aligned_fragments) = extract_fragments(s2g_alignment, amr_graph)
            logger.writeln(str(opt_toks))
            logger.writeln(str(role_toks))

            if not aligned_fragments:
                logger.writeln('wrong alignments')
                continue

            temp_aligned = set(aligned_fragments.keys())
            aligned_fragments = sorted(aligned_fragments.items(), key=lambda frag: frag[0])

            temp_unaligned = set(xrange(len(pos_seq))) - temp_aligned

            ####Extract entities#####
            for (frag, entity_tag) in amr_graph.extract_entities():
                if len(opt_toks) == 0:
                    logger.writeln("No alignment for the entity found")
                    no_alignment = True

                (frag_start, frag_end, multiple) = extract_entity_spans(frag, opt_toks, role_toks, temp_unaligned)
                if frag_start is None:
                    logger.writeln("No alignment found")
                    logger.writeln(str(frag))

                    no_alignment = True
                    continue

                if multiple:
                    has_multiple = True
                    logger.writeln("Multiple found here!")

                frag.set_span(frag_start, frag_end)

                new_aligned = set(xrange(frag_start, frag_end))
                if len(new_aligned & aligned_set) != 0:
                    has_multiple = True
                    continue

                aligned_set |= new_aligned
                all_frags.append((frag, 'ent+%s' % entity_tag, False, True))

                if (edge_alignment & frag.edges).count() != 0:
                    has_multiple = True

                edge_alignment |= frag.edges
                logger.writeln('Exracted entities:')
                logger.writeln(' '.join(toks[frag_start:frag_end]))
                logger.writeln(str(frag))

            if no_alignment:
                continue

            one2many = False
            #####Extra other alignments######
            logger.writeln('Aligned fragments:')
            for (index, frag_list) in aligned_fragments:
                if index in aligned_set:
                    continue

                assert len(frag_list) > 0
                non_conflict = 0
                non_conflict_list = []
                for frag in frag_list:
                    if (edge_alignment & frag.edges).count() == 0:
                        non_conflict += 1
                        non_conflict_list.append(frag)

                if non_conflict != 1:
                    one2many = True

                used_frag = None
                if non_conflict == 0:
                    used_frag = frag_list[0]
                else:
                    used_frag = non_conflict_list[0]

                edge_alignment |= used_frag.edges
                (used_frag, is_pred, frag_label) = enlarge_frag(used_frag, amr_graph, edge_alignment)
                #if is_pred:
                #    pred_label = str(amr_graph.nodes[used_frag.root])
                #    if tok2
                #    word2predicate[toks[index]] = pred_label
                curr_tok = toks[index]
                curr_pos = pos_seq[index]
                curr_frag_str = str(used_frag).replace(' ', '')

                all_frags.append((used_frag, frag_label, is_pred, False))
                aligned_set.add(index)

            logger.writeln("%d aligned edges out of %d total" % (edge_alignment.count(), len(edge_alignment)))
            used_sents += 1

            assert len(toks) == len(pos_seq)

            unaligned_toks = [(i, tok) for (i, tok) in enumerate(toks) if i not in aligned_set]

            (aligned, unaligned) = amr_graph.recall_unaligned_concepts(edge_alignment, unaligned_toks, lemma_map, stop_words)

            for frag in aligned:
                aligned_set.add(frag[0].start)

            all_frags += aligned

            un_set = set([i for i in xrange(len(toks)) if i not in aligned_set])
            for index in un_set:
                n_nodes = len(amr_graph.nodes)
                n_edges = len(amr_graph.edges)
                frag = AMRFragment(n_edges, n_nodes, amr_graph)
                frag.set_span(index, index+1)
                all_frags.append((frag, 'UNKNOWN', False, False))
                unknown_set.add(toks[index])

            aligned_fragments = sorted(all_frags, key=lambda frag: (frag[0].start, frag[0].end))

            for (frag, frag_label, is_pred, is_ent) in aligned_fragments:
                if frag_label != 'UNKNOWN':
                    frag.build_ext_list()
                    frag.build_ext_set()

                if is_ent:
                    entity_str = '_'.join(toks[frag.start:frag.end])
                    entity_words.add(entity_str)
                    entity_map[entity_str].add(frag_label)

                    ent2frag[entity_str] = frag

                elif frag_label != 'UNKNOWN':
                    assert frag.end - frag.start == 1
                    curr_tok = toks[frag.start]
                    curr_lem = lemmas[frag.start]
                    if is_pred:
                        if curr_tok not in stop_words:
                            predicate_words.add(curr_tok)

                        if curr_lem not in stop_words:
                            predicate_lemmas.add(curr_lem)

                        predicate_labels.add(frag_label)
                        pred_label = str(amr_graph.nodes[frag.root])
                        assert '/' in pred_label, pred_label

                        if curr_tok not in word2predicate:
                            word2predicate[curr_tok] = defaultdict(int)

                        word2predicate[curr_tok][pred_label] += 1

                        if curr_tok != curr_lem:
                            if curr_lem not in word2predicate:
                                word2predicate[curr_lem] = defaultdict(int)
                            word2predicate[curr_lem][pred_label] += 1

                    else:
                        non_pred_map[curr_tok].add(frag_label)
                        non_pred_lemma_map[curr_lem].add(frag_label)

                        label2frag[frag_label] = frag

                fs = []
                f_start = frag.start
                f_end = frag.end-1
                #Word feature
                fs += extract_span(toks, f_start, f_end, 3, 'word')
                fs += extract_bigram(toks, f_start, f_end, 3, 'word')
                fs += extract_curr(toks, f_start, f_end, 'word')
                if is_ent:
                    fs += extract_seq_feat(toks, f_start, f_end, 'word')

                #Lemma feature
                fs += extract_span(lemmas, f_start, f_end, 3, 'lemma')
                fs += extract_bigram(lemmas, f_start, f_end, 3, 'lemma')
                fs += extract_curr(lemmas, f_start, f_end, 'lemma')

                #Pos tag feature
                fs += extract_span(pos_seq, f_start, f_end, 3, 'POS')
                fs += extract_bigram(pos_seq, f_start, f_end, 3, 'POS')
                fs += extract_curr(pos_seq, f_start, f_end, 'POS')

                #Length of span feature
                fs.append('Length=%d' % (frag.end - frag.start))

                #Suffix feature
                if not is_ent and f_start == f_end:
                    fs += suffix(toks[f_start])

                print >>feature_f, '%s %d-%d %s %s %s' % (frag_label, frag.start, frag.end, '1' if is_pred else '0', '1' if is_ent else '0', ' '.join(fs))

                label_set.add(frag_label)
            print >>feature_f, ''
            print >>tok_f, ' '.join(toks)
            print >>lemma_f, ' '.join(lemmas)
            print >>pos_f, ' '.join(pos_seq)

            unaligned_words = set(range(len(toks))) - aligned_set

            un_seq = []
            for pos in unaligned_words:
                un_seq.append(toks[pos])
            logger.writeln("Unaligned toks: %s" % ' '.join(un_seq))

        feature_f.close()
        tok_f.close()
        lemma_f.close()
        pos_f.close()

        statistics = []
        statistics.append(predicate_words)
        statistics.append(predicate_lemmas)
        statistics.append(predicate_labels)

        statistics.append(entity_words)
        statistics.append(entity_lemmas)
        statistics.append(entity_map)

        statistics.append(non_pred_map)
        statistics.append(non_pred_lemma_map)
        statistics.append(unknown_set)

        stat_file = os.path.join(args.save_dir, 'stat_%d' % args.sub_id)
        stat_f = open(stat_file, 'wb')
        cPickle.dump(statistics, stat_f)
        stat_f.close()

        word_pred_file = os.path.join(args.save_dir, 'word2pred_%d' % args.sub_id)
        stat_f = open(word_pred_file, 'wb')
        cPickle.dump(word2predicate, stat_f)
        stat_f.close()

        ent_frag_file = os.path.join(args.save_dir, 'ent2frag_%d' % args.sub_id)
        stat_f = open(ent_frag_file, 'wb')
        cPickle.dump(ent2frag, stat_f)
        stat_f.close()

        label_frag_file = os.path.join(args.save_dir, 'label2frag_%d' % args.sub_id)
        stat_f = open(label_frag_file, 'wb')
        cPickle.dump(label2frag, stat_f)
        stat_f.close()

        logger.writeln('A total of %d unknown words' % len(unknown_set))
        logger.writeln('A total of: %d labels' % len(label_set))
        #logger.writeln(str(label_set))
        logger.writeln('total used sents is: %d' % used_sents)
        logger.writeln('total dump forests: %d' % len(constructed_forests))
        logger.writeln('finished')

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()

    argparser.add_argument("--amr_file", type=str, help="the original AMR graph files", required=False)
    argparser.add_argument("--stop", type=str, help="stop words file", required=False)
    argparser.add_argument("--lemma", type=str, help="lemma file", required=False)
    argparser.add_argument("--dump_graph", action="store_true", help="if only to dump graph object")
    argparser.add_argument("--preprocess", action="store_true", help="if needed to preprocess the AMR graphs into dumped AMR graph objects")
    argparser.add_argument("--parallel", action="store_true", help="if to run multiple process to run the forest construction")
    argparser.add_argument("--dump_rules", action="store_true", help="if to dump lexical rules")
    argparser.add_argument("--refine", action="store_true", help="if to refine the nonterminals")
    argparser.add_argument("--data_dir", type=str, help="the data directory for dumped AMR graph objects, alignment and tokenized sentences")
    argparser.add_argument("--save_dir", type=str, help="the output directory for saving the constructed forest")
    argparser.add_argument("--nodes", type=str, help="nodes for running processes")
    argparser.add_argument("--sent_per_node", type=int, help="number of sentences for each node")
    argparser.add_argument("--sub_id", type=int, help="if a subprocess, the subprocess id")

    args = argparser.parse_args()
    construct_forest(args)
