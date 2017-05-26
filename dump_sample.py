#!/usr/bin/python
from __future__ import print_function
from __future__ import division

import sys
import time
import pickle
import os
import random
import cPickle
from math import log, exp, factorial, lgamma
import copy
import socket

import logprob
import logger
import gflags
from common import INF, ZERO
from fragment_forest import *
from collections import deque
from HRGSample import *

FLAGS = gflags.FLAGS

PHRASE_NT = 'X'

gflags.DEFINE_string(
    'base',
    'poisson',
    'Base distribution')
gflags.DEFINE_float(
    'alpha',
    5.0,
    'Concentration parameter in Dirichlet process.')
gflags.DEFINE_float(
    'discount',
    0.5,
    'Discount parameter in Pitman-yor process.')
gflags.DEFINE_integer(
    'maxcut',
    7,
    'no cut sampling when target or source span larger than maxcut.')
gflags.DEFINE_boolean(
    'sample_cut_only',
    False,
    'Sample at only nodes already cut (wrong implementation).')
gflags.DEFINE_integer(
    'sample_level',
    None,
    'Sample only nodes with level <= sample_level."')
gflags.DEFINE_integer(
    'level_inc',
    None,
    'Sample each level for #level_inc# iterations."')
gflags.DEFINE_boolean(
    'double_count',
    False,
    'Use double counting."')
gflags.DEFINE_boolean(
    'variable_alpha',
    False,
    'concentration parameter is different for each length')
gflags.DEFINE_boolean(
    'correct_edge_sampling',
    False,
    '"correct" path sampling. Number of incoming nodes is considered')
gflags.DEFINE_boolean(
    'lhs_conditional',
    False,
    'normalized over lhs.')
gflags.DEFINE_boolean(
    'seed_random',
    False,
    'seed random number generation.')
gflags.DEFINE_boolean(
    'sample_cut',
    True,
    'sample cut point (set to no to sample only minimal rules)')
gflags.DEFINE_boolean(
    'sample_edge',
    True,
    'sample edge switch')
gflags.DEFINE_integer(
    'splits',
    2,
    'two-way nt split by default, but you can change it.')
gflags.DEFINE_string(
    'model',
    'PY',
    'model. choose from DP|PY')
gflags.DEFINE_integer(
    'split_iter',
    10,
    'do symbol split every #split_iter iterations')
gflags.DEFINE_boolean(
    'split_at_iter_0',
    True,
    'allow split to happen at iter 0')
gflags.DEFINE_boolean(
    'type',
    False,
    'use type-based sampling')
gflags.DEFINE_boolean(
    'refine',
    False,
    'use symbol refinement')
gflags.DEFINE_boolean(
    'check_index',
    False,
    'Check cut index for errors.')
gflags.DEFINE_list(
    'nodes',
    'node82,node83,node84,node85',
    'paralleled nodes')
gflags.DEFINE_integer(
    'port',
    12345,
    'port number')
gflags.DEFINE_string(
    'host',
    'cycle2.cs.rochester.edu',
    'host machine')
gflags.DEFINE_boolean(
    'slave',
    False,
    'Slave node')
gflags.DEFINE_string(
    'current',
    'cycle2',
    'current child node')
gflags.DEFINE_integer(
    'currid',
    0,
    'current child id')


#base = poisson
#rule_size_prob = poisson
def update_sampler(samples, g_sampler):
    for sample in timed(samples):
        for n, rule in sample.composed_rules_under(sample.hg.root):
            g_sampler.count(rule)

'''
xpeng: init the sampler with the current settings of the samples
rule, total rule and rule_size count. Also initiate the type-based indexer
'''
def init_split(samples, split=True):
    #global SAMPLER
    logger.writeln('initialization. split=%s' % split)
    SAMPLER = init_sampler()
    for sample in timed(samples):
        sample.set_sampler(SAMPLER)
        if split:
            for node in sample.hg.nodes:
                node.pnt = node.nt
                node.nt = random.choice(child_symbols(node.pnt))
        total_count = 0
        #print type(sample.hg.root)
        for n, rule in sample.composed_rules_under(sample.hg.root):
            total_count += 1
            SAMPLER.count(rule)
        logger.writeln('Total number of rules: %d' % total_count)
    #for rule, c in SAMPLER.counts.items():
    #    logger.writeln('Rule hash: %d, Rule counts: %d' % (hash(rule), c))
    return SAMPLER

def child_symbols(nt):
    "Return range of symbol indices given parent symbol index"
    return range(nt*FLAGS.splits+1, (nt+1)*FLAGS.splits+1)

def parent_symbol(nt):
    return (nt - 1)//2

'''
xpeng: return the non-terminal # of the node
1.if is symbol-refined, also return the node.nt attribute
2.else just PHRASE_NT
'''
def get_nt(node):
    if FLAGS.refine:
        return '[%s-%s]' % (PHRASE_NT, node.nt)
    else:
        return '[%s]' % PHRASE_NT

def init_sampler():
    if FLAGS.lhs_conditional:
        sampler = NTSampler()
    else:
        sampler = NPSampler()

    if FLAGS.model == 'PY':
        NPSampler.likelihood = NPSampler.rule_size_likelihood
        NPSampler.posterior = NPSampler.pitman_yor_posterior_rule_size
    elif FLAGS.model == 'DP':
        NPSampler.likelihood = NPSampler.dp_likelihood
        NPSampler.posterior = NPSampler.simple_dirichlet_posterior
    else:
        assert False, 'unsupported model'
    return sampler

class TreeFile():
    def __init__(self, filename):
        self.f = open(filename, 'w')
        self.i = 1

    def dump(self, sample):
        self.f.write('# %s\n' % self.i)
        self.f.write('%s\n' % sample.tree_str())
        self.i += 1

    def close(self):
        self.f.close()

def dump_trees(samples, filename):
    file_handle = open(filename, 'w')
    for sample in samples:
        sample.dump_hrg_rules(file_handle)
    file_handle.close()

    #logger.writeln('dump trees')
    #treefile = TreeFile(filename)
    #for s in timed(samples):
    #    # call this before dumping rules for each sample!
    #    LEXICAL_WEIGHTER.compute_lexical_weights(s.a)
    #    treefile.dump(s)
    #treefile.close()

def choose_k(n):
    return random.randrange(n+1)

'''
xpeng: find the nearest cut-node ancester
'''
def cut_parent(node):
    #print(node)
    assert hasattr(node, 'parent')
    p = node.parent
    while not p.cut:
        p = p.parent
    return p

def index(node):
    p = node.parent
    n = node
    result = []
    while True:
        for i, c in enumerate(p.incoming[p.edge].tail):
            if c is n:
                result.append(i)
        if p.cut:
            break
        n = p
        p = p.parent
    result.reverse()
    return result

if __name__ == '__main__':
    gflags.DEFINE_integer(
        'interval',
        5000,
        'Print stat every #interval# sentences.')
    gflags.DEFINE_integer(
        'iter',
        1,
        'Number of sampling iterations.')
    gflags.DEFINE_integer(
        'dump_iter',
        1,
        'Dump trees every #dump_iter# iterations.')
    gflags.DEFINE_string(
        'dump',
        'dump',
        'Dump directory.')
    gflags.DEFINE_string(
        'data',
        'data',
        'Data directory.')
    gflags.DEFINE_string(
        'forest_dir',
        'forest_dir',
        'Forest file directory.')
    gflags.DEFINE_string(
        'prefix',
        'forest',
        'Separate forest file prefix.')
    gflags.DEFINE_string(
        'file_indexes',
        '0',
        'Different forest file indexes, separted by +.')

    try:
        argv = FLAGS(sys.argv)  # parse flags
    except gflags.FlagsError as e:
        print('%s\nUsage: %s ARGS\n%s' % (e, sys.argv[0], FLAGS))
        sys.exit(1)

    file_ids = FLAGS.file_indexes.split('+')
    sys.setrecursionlimit(sys.getrecursionlimit() * 30)
    for id in file_ids:
        hg_file = os.path.join(FLAGS.data, 'forest_%s' % id)
        f = open(hg_file, 'rb')
        hypergraphs = cPickle.load(f)
        f.close()

        sent_no_file = os.path.join(FLAGS.data, 'used_sent_%s' % id)
        f = open(sent_no_file, 'rb')
        sent_nos = cPickle.load(f)
        f.close()

        samples = []
        for (hg, sent_num) in zip(hypergraphs, sent_nos):
            samples.append(Sample(hg, sent_num))

        sample_file = os.path.join(FLAGS.data, 'sample_%s' % id)
        f = open(sample_file, 'wb')
        cPickle.dump(samples, f)
        f.close()
        samples[:] = []

