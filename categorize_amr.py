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
from rule import Rule
import argparse
from re_utils import *
from collections import defaultdict
from filter_stop_words import *
from lemma_util import initialize_lemma

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
def build_one_node(curr_frag, curr_start, curr_end, amr_graph, edge_alignment, refine=False):
    curr_node_index = curr_frag.root
    curr_graph_node = amr_graph.nodes[curr_node_index]

    if edge_alignment[curr_graph_node.c_edge] == 0:
        new_node = FragmentHGNode(FRAGMENT_NT, curr_start, curr_end, curr_frag)
        return new_node

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
                    #logger.writeln("what the hell is this")
                    #logger.writeln(str(curr_frag))
                    is_pred = True
                continue

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
                    else:
                        is_op = True
                continue

            tail_node_index = curr_edge.tail

            #Our intuition: ARGs and ops goes with the root
            if is_root_arc(edge_label):
                if 'ARG' in edge_label:
                    is_pred = True
                else:
                    assert 'op' in edge_label
                    is_op = True
                root_arcs.append((curr_edge_index, tail_node_index))

    unaligned_node = None
    if refine:
        init_ext_frag(curr_frag, is_pred, is_op) #Initialize the current fragment

    if len(root_arcs) > 0 or len(head_arcs) > 0:
        n_nodes = len(amr_graph.nodes)
        n_edges = len(amr_graph.edges)
        frag = AMRFragment(n_edges, n_nodes, amr_graph)
        frag.set_root(curr_node_index)

        for rel_index, tail_index in root_arcs:
            edge_alignment[rel_index] = 1
            frag.set_edge(rel_index)
            frag.set_node(tail_index)

        if head_arcs:
            (rel_index, head_index) = head_arcs[0]
            edge_alignment[rel_index] = 1
            frag.set_edge(rel_index)
            frag.set_root(head_index)

        if refine:
            init_ext_frag(frag, is_pred, is_op)

        frag.build_ext_list()
        frag.build_ext_set()
        new_frag = combine_fragments(curr_frag, frag, refine)
        assert new_frag, 'Weird combination found'

        new_node = FragmentHGNode(FRAGMENT_NT, curr_start, curr_end, new_frag)

    else: #Should be either an entity or a single concept
        new_node = FragmentHGNode(FRAGMENT_NT, curr_start, curr_end, curr_frag)

    s = Sample(hypergraph.Hypergraph(new_node), 0)
    new_node.cut = 1
    new_rule, _ = s.extract_one_rule(new_node, None, new_node.frag.ext_list, refine)
    rule_str = '%s ||| %s\n' % (filter_vars(new_rule.dumped_format()), context_str(new_node.frag, amr_graph))
    rule_f.write(rule_str)
    fields = rule_str.split(' ||| ')
    fields[1] = ' '.join(amr_graph.lems[new_node.frag.start: new_node.frag.end])
    lem_rule_str = ' ||| '.join(fields)
    lemma_rule_f.write(lem_rule_str)
    return new_node

def build_bimap(tok2frags):
    frag2map = defaultdict(set)
    index2frags = defaultdict(set)
    for index in tok2frags:
        for frag in tok2frags[index]:
            index2frags[index].add(frag)
            frag2map[frag].add(index)
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
    #assert frag.edges.count() == 1, 'Not unify edge fragment found'
    amr_graph = frag.graph
    edge_list = []
    n_edges = len(frag.edges)
    for i in xrange(n_edges):
        if frag.edges[i] == 1:
            edge_list.append(i)
    assert len(edge_list) == frag.edges.count()
    return tuple(edge_list)
    #root_node = amr_graph.nodes[frag.root]
    #for edge_index in root_node.v_edges:
    #    if frag.edges[edge_index] == 1:
    #        return edge_index
    #assert True, 'This is impossible'
    #return None

class AMR_stats(object):
    def __init__(self):
        self.num_reentrancy = 0
        self.num_predicates = defaultdict(int)
        self.num_nonpredicate_vals = defaultdict(int)
        self.num_consts = defaultdict(int)
        self.num_entities = defaultdict(int)

    def update(self, local_re, local_pre, local_non, local_con, local_ent):
        self.num_reentrancy += local_re
        for s in local_pre:
            self.num_predicates[s] += local_pre[s]

        for s in local_non:
            self.num_nonpredicate_vals[s] += local_non[s]

        for s in local_con:
            self.num_consts[s] += local_con[s]

        for s in local_ent:
            self.num_entities[s] += local_ent[s]

    def dump2dir(self, dir):
        def dump_file(f, dict):
            sorted_dict = sorted(dict.items(), key=lambda k:(-k[1], k[0]))
            for (item, count) in sorted_dict:
                print >>f, '%s %d' % (item, count)
            f.close()

        pred_f = open(os.path.join(dir, 'pred'), 'w')
        non_pred_f = open(os.path.join(dir, 'non_pred_val'), 'w')
        const_f = open(os.path.join(dir, 'const'), 'w')
        entity_f = open(os.path.join(dir, 'entities'), 'w')

        dump_file(pred_f, self.num_predicates)
        dump_file(non_pred_f, self.num_nonpredicate_vals)
        dump_file(const_f, self.num_consts)
        dump_file(entity_f, self.num_entities)

    def __str__(self):
        s = ''
        s += 'Total number of reentrancies: %d\n' % self.num_reentrancy
        s += 'Total number of predicates: %d\n' % len(self.num_predicates)
        s += 'Total number of non predicates variables: %d\n' % len(self.num_nonpredicate_vals)
        s += 'Total number of constants: %d\n' % len(self.num_consts)
        s += 'Total number of entities: %d\n' % len(self.num_entities)
        return s


def linearize_amr(args):
    logger.file = open(os.path.join(args.run_dir, 'logger'), 'w')

    #amr_file = os.path.join(args.data_dir, 'amr')
    amr_file = os.path.join(args.data_dir, 'aligned_amr_nosharp')
    alignment_file = os.path.join(args.data_dir, 'alignment')
    sent_file = os.path.join(args.data_dir, 'sentence')
    tok_file = os.path.join(args.data_dir, 'token')
    #lemma_file = os.path.join(args.data_dir, 'lemma')
    pos_file = os.path.join(args.data_dir, 'pos')

    amr_graphs = load_amr_graphs(amr_file)
    alignments = [line.strip().split() for line in open(alignment_file, 'r')]
    sents = [line.strip().split() for line in open(sent_file, 'r')]
    toks = [line.strip().split() for line in open(tok_file, 'r')]
    #lemmas = [line.strip().split() for line in open(lemma_file, 'r')]
    poss = [line.strip().split() for line in open(pos_file, 'r')]

    assert len(amr_graphs) == len(alignments) and len(amr_graphs) == len(sents) and len(amr_graphs) == len(toks) and len(amr_graphs) == len(poss), '%d %d %d %d %d' % (len(amr_graphs), len(alignments), len(sents), len(toks), len(poss))
    #assert len(amr_graphs) == len(alignments) and len(amr_graphs) == len(sents) and len(amr_graphs) == len(toks) and len(amr_graphs) == len(lemmas) and len(amr_graphs) == len(poss), '%d %d %d %d %d %d' % (len(amr_graphs), len(alignments), len(sents), len(toks), len(lemmas), len(poss))

    #lemma_map = initialize_lemma(args.lemma)
    num_self_cycle = 0
    used_sents = 0

    amr_statistics = AMR_stats()

    for (sent_index, (sent_seq, tok_seq, pos_seq, alignment_seq, amr_graph)) in enumerate(zip(sents, toks, poss, alignments, amr_graphs)):

        logger.writeln('Sentence #%d' % (sent_index+1))
        logger.writeln(str(amr_graph))

        #if sent_index > 100:
        #    break

        edge_alignment = bitarray(len(amr_graph.edges))
        if edge_alignment.count() != 0:
            edge_alignment ^= edge_alignment
        assert edge_alignment.count() == 0

        has_cycle = False
        if amr_graph.check_self_cycle():
            num_self_cycle += 1
            has_cycle = True
            #logger.writeln('self cycle detected')

        amr_graph.set_sentence(toks)
        #amr_graph.set_lemmas(lemma_seq)
        amr_graph.set_poss(pos_seq)

        aligned_fragments = []
        reentrancies = {}  #Map multiple spans as reentrancies, keeping only one as original, others as connections

        has_multiple = False
        no_alignment = False

        aligned_set = set()

        all_frags = []

        (opt_toks, role_toks, aligned_fragments) = extract_fragments(alignment_seq, amr_graph)
        ##logger.writeln(str(opt_toks))
        ##logger.writeln(str(role_toks))

        #if not aligned_fragments:
        #    logger.writeln('wrong alignments')
        #    continue

        #temp_aligned = set(aligned_fragments.keys())
        #aligned_fragments = sorted(aligned_fragments.items(), key=lambda frag: frag[0])

        #temp_unaligned = set(xrange(len(pos_seq))) - temp_aligned

        (entity_frags, root2entityfrag, root2entitynames) = amr_graph.extract_all_entities()

        logger.writeln(amr_graph.collapsed_form(root2entityfrag, root2entitynames))
        (entity_nums, predicate_nums, variable_nums, const_nums, reentrancy_nums) = amr_graph.statistics(root2entityfrag, root2entitynames)

        amr_statistics.update(reentrancy_nums, predicate_nums, variable_nums, const_nums, entity_nums)
        #logger.writeln(stats[0])
        #logger.writeln(stats[1])
        #logger.writeln(stats[2])
        #logger.writeln(stats[3])
        #logger.writeln(stats[4])

        ####Extract entities#####
        #for (frag, frag_label) in amr_graph.extract_entities():
        #    if len(opt_toks) == 0:
        #        logger.writeln("No alignment for the entity found")
        #        #no_alignment = True

        #    (frag_start, frag_end, multiple) = extract_entity_spans(frag, opt_toks, role_toks, temp_unaligned)
        #    if frag_start is None:
        #        logger.writeln("No alignment found")
        #        logger.writeln(str(frag))

        #        no_alignment = True
        #        continue

        #    if multiple:
        #        has_multiple = True
        #        logger.writeln("Multiple found here!")

        #    frag.set_span(frag_start, frag_end)

        #    amr_graph.collapse_entities(frag, frag_label)

        #    new_aligned = set(xrange(frag_start, frag_end))
        #    if len(new_aligned & aligned_set) != 0:
        #        print str(amr_graph)
        #        print str(frag)
        #        has_multiple = True
        #        break
        #        #continue

        #    aligned_set |= new_aligned
        #    all_frags.append(frag)

        #    if (edge_alignment & frag.edges).count() != 0:
        #        has_multiple = True

        #    edge_alignment |= frag.edges

        #if no_alignment:
        #    continue

        #one2many = False
        ######Extra other alignments######
        #logger.writeln('Aligned fragments:')
        #for (index, frag_list) in aligned_fragments:
        #    if index in aligned_set:
        #        continue

        #    assert len(frag_list) > 0
        #    non_conflict = 0
        #    non_conflict_list = []
        #    for frag in frag_list:
        #        if (edge_alignment & frag.edges).count() == 0:
        #            non_conflict += 1
        #            non_conflict_list.append(frag)

        #    if non_conflict != 1:
        #        one2many = True

        #    used_frag = None
        #    if non_conflict == 0:
        #        used_frag = frag_list[0]
        #    else:
        #        used_frag = non_conflict_list[0]

        #    edge_alignment |= used_frag.edges
        #    all_frags.append(used_frag)

        #    aligned_set.add(index)

        #logger.writeln("%d aligned edges out of %d total" % (edge_alignment.count(), len(edge_alignment)))
        #used_sents += 1

        #assert len(toks) == len(pos_seq)

        #unaligned_toks = [(i, tok) for (i, tok) in enumerate(toks) if i not in aligned_set]
        #(aligned, unaligned) = amr_graph.recall_unaligned_concepts(edge_alignment, unaligned_toks, lemma_map, stop_words)
        #aligned = [x for (x, y, z, k) in aligned]

        #all_frags += aligned
    amr_statistics.dump2dir(args.run_dir)
    print str(amr_statistics)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()

    argparser.add_argument("--amr_file", type=str, help="the original AMR graph files", required=False)
    argparser.add_argument("--stop", type=str, help="stop words file", required=False)
    argparser.add_argument("--lemma", type=str, help="lemma file", required=False)
    argparser.add_argument("--data_dir", type=str, help="the data directory for dumped AMR graph objects, alignment and tokenized sentences")
    argparser.add_argument("--run_dir", type=str, help="the output directory for saving the constructed forest")

    args = argparser.parse_args()
    linearize_amr(args)
