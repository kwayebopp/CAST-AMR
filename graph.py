import amr_graph_0
import math
#from graph_utils import *
from nltk.model import *
import dill
from collections import OrderedDict

print 'building LMs'
with open('bigram.lm', 'rb') as f:
    bg_model = dill.load(f)
with open('trigram.lm', 'rb') as f:    
    tg_model =  dill.load(f)
with open('4gram.lm', 'rb') as f:    
    qg_model = dill.load(f)
print 'LMs complete'

sent_con_map = OrderedDict()
sent_node_map = OrderedDict()
sent_rule_map = OrderedDict()
sent_amr_map = OrderedDict()
sent_graph_map = OrderedDict()

dev_rules = './devout/dumped_rules_0'
#dev_cons = './dev/amr.lin'
dev_amrs = './dev/amr-aligned'

#train_rules = './trainout/dumped_rules_0'
#train_cons = './train/amr.lin'
#train_amrs = '.train/amr-aligned'  #make nosharp file

test_rules = './testout/dumped_rules_0'
#test_cons = './test/amr.lin'
test_amrs = './test/amr-aligned' #make nosharp files

rules = 'YOUR RULE FILE GOES HERE'
amrs = 'YOUR AMR FILE GOES HERE'


#does the rule have an instance of concept?
def rule_has_concept(rule, concept):
    return concept in rule.split('|||')[2] #is concept in AMR frag of rule?

#generate list of nodes <==> (concept, rule) pairs
def node_factory(sentence):
    nodes = []
    for rule in sent_rule_map[sentence]: 
        for concept in sent_con_map[sentence]:
            if rule_has_concept(rule, concept):#match concepts and rules
                n = Node(concept, rule)
                if n not in nodes: 
                    nodes.append(n)
    return nodes

def ordered_concepts(rule):
    ordered_nodes = []

    frags = rule.split('|||')
    frag = ''
    if len(frags) >= 3:
        frag = frags[2]
    parts = frag.split()

    #build list of nodes in order
    for idx, part in enumerate(parts):
        if ':op' in part:
            idx = idx + 1
            continue
        elif '"' in part:
            part = part.strip(':').strip('"')
            if "_" in part:
                if '#' in part:
                    wordz = part.split('#')
                    for word in wordz:
                        ordered_nodes += word
            elif '#' in part:
                wordz = part.split('#')
                ordered_nodes += wordz  
            else:
                ordered_nodes.append(part)      
        elif '/' in part:
            ordered_nodes.append(part.strip(':').split('/')[1])

    return ordered_nodes

#build dict mapping sents to concepts
def build_sent_rule_map(rule_file):
    global sent_rule_map
    sent = None
    with open(rule_file, 'r') as rules:
        for line in rules:
            if "#Sentence: " in line:
                sent = line.replace("#Sentence: ", "").strip('\n') #grab sentence
                sent_rule_map[sent] = [] #initialize list of values
            elif len(line) <= 1:
                continue
            else: sent_rule_map[sent].append(line) #append rules to value list
        return sent_rule_map

# builds dictionary mapping sentences to concepts 
# (where concepts are lemmatized/verbalized words excluding stopwords)
def build_sent_con_map():
    global sent_con_map
    for sent in sent_rule_map.keys():
        sent_con_map[sent] = []
        for rule in sent_rule_map[sent]:
            o_n = ordered_concepts(rule)
            for concept in o_n:
                if concept not in sent_con_map[sent]:
                    sent_con_map[sent].append(concept)
    return sent_con_map 

def build_sent_amr_map(aligned_amr_file): #use nosharp
    
    with open(aligned_amr_file, 'r') as data:
        global sent_amr_map  
        sent = ""
        amr = ''
        lines = data.readlines()
        for idx in range(0, len(lines)):
            if lines[idx][0] == '#':
                if '# ::tok ' in lines[idx]:
                    sent = lines[idx].replace('# ::tok ', "").strip()
                else: continue
            elif len(lines[idx].strip()) != 0:
                amr += lines[idx].strip('\n').strip()
            else:
                if sent not in sent_amr_map.keys():
                    sent_amr_map[sent] = amr
                    amr = ''
                   
        return sent_amr_map
                        

#map sentences to nodes
def build_sent_node_map():
    global sent_node_map
    for sentence in sent_rule_map.keys():
        sent_node_map[sentence] = node_factory(sentence) #map sentence to all the nodes it generates
    return sent_node_map

#build maps
def build_maps(rule_file, amr_file):
    global sent_con_map, sent_rule_map, sent_node_map, sent_amr_map
    sent_rule_map = build_sent_rule_map(rule_file)
    sent_con_map = build_sent_con_map()
    sent_node_map = build_sent_node_map()
    sent_amr_map = build_sent_amr_map(amr_file)

class Node:
    def __init__(self, concept, r):
        self.concept = concept
        self.rule = r

    def __str__(self): return "(" + self.concept + ", " + self.rule + ")"
    def __repr__(self): return "(" + self.concept + ", " + self.rule + ")"

class Graph:
    def get_path_length(self, n_i, n_j):
        #print sent_amr_map[self.sentence]
        amr = amr_graph_0.AMRGraph(sent_amr_map[self.sentence])
        i = None
        j = None
        for node in amr.nodes:
            if str(node) == n_i.concept:
                i = node
            elif str(node) == n_j.concept:
                j = node
        
        return len(amr.get_path_v2(i, j)[0])

    def traveling_cost(self, n_i, n_j, model=None):            
        trans_i = []
        trans_j = []

        if n_i.concept == 'i':
            trans_i = ['I']
        elif n_i.concept == '<s>':
            trans_i = [n_i.concept]
        else:
            trans_i = n_i.rule.split('|||')[1].strip().split()
  
        if n_j.concept == 'i':
            trans_j = ['I']
        elif n_j.concept == '</s>':
            trans_j = [n_j.concept]
        else:
            try: 
                trans_j = n_j.rule.split('|||')[1].strip().split()
            except IndexError:
                print n_j.rule.split('|||')

        

        if trans_i == '<s>':
            model = bg_model
            try:
                return math.log(1.0/model.score(trans_j[0], trans_i))
            except ZeroDivisionError:
                return float(1000)

        elif trans_j == '</s>':
            model = bg_model
            try:
                return math.log(1.0/model.score(trans_j, trans_i[len(trans_i) - 2 : len(trans_i) - 1]))
            except ZeroDivisionError:
                return float(1000)

        else:
            try:
                word_count = len(trans_i) + len(trans_j)
            except:
                print n_i.rule.split('|||')
                print n_j.rule.split('|||')

            if len(trans_i) == 1:
                model = bg_model
                try:
                    return math.log(1.0/model.score(trans_j[0], trans_i) * (math.exp(word_count * self.get_path_length(n_i, n_j))))
                except ZeroDivisionError as e:
                    return float(1000)
    
            elif len(trans_i) == 2:
                model = tg_model
                try:
                    return math.log(1.0/model.score(trans_j[0], trans_i) * (math.exp(word_count * self.get_path_length(n_i, n_j))))

                except ZeroDivisionError as e:
                    return float(1000)            
                
            elif len(trans_i) >= 3:
                model = qg_model
                try:
                    return math.log(1.0/model.score(trans_j[0], trans_i[len(trans_i) - 3:]) * (math.exp(word_count * self.get_path_length(n_i, n_j))))
                except ZeroDivisionError as e:
                    return float(1000)
            
            else: return 1000.0

    def build_cost_matrix(self):
        n_start = Node("<s>", "<s>")
        n_end = Node("</s>", "</s>")
        self.nodes.append(n_start)
        self.nodes.append(n_end)
        self.update_graph_size()

        mat = [[float(1000) for x in range(self.size)] for y in range(self.size)]
        for node in self.nodes:
            x = ordered_concepts(node.rule)
            frag_first = ''
            if len(x) != 0:
                frag_first = x[0]
            if node.concept == frag_first:
                mat[self.nodes.index(n_start)][self.nodes.index(node)] = self.traveling_cost(n_start, node)
                #print str(mat[self.nodes.index(n_start)][self.nodes.index(node)]) + ' @ ' + str(self.nodes.index(node)) + ', ' + str(self.nodes.index(n_start)) 

            else:
                mat[self.nodes.index(n_start)][self.nodes.index(node)] = float(1000)

            frag_last = ''
            if len(x) >= 1:
                frag_last = x[len(x) - 1] #clean up string and save
            if node.concept == frag_last:
                mat[self.nodes.index(node)][self.nodes.index(n_end)] = self.traveling_cost(node, n_end)
                #print str( mat[self.nodes.index(node)][self.nodes.index(n_end)]) + ' @ ' + str(self.nodes.index(node)) + ', ' + str(self.nodes.index(n_end))
 
            else:
                mat[self.nodes.index(node)][self.nodes.index(n_end)] = float(1000)


        for node1 in self.nodes:
            for node2 in self.nodes:
                if node1 != node2:
                    o_n = ordered_concepts(node1.rule)
                    o_n_2 = ordered_concepts(node2.rule)
                    if node1.rule == node2.rule:
                       # print 'same rule'
                        if node1.concept in o_n:
                           # print 'node1.concept in o_n'
                            c_idx = o_n.index(node1.concept)
                            if c_idx < len(o_n) - 1 and o_n[c_idx + 1] == node2.concept:
                              #  print 'node2.concept is next'
                                mat[self.nodes.index(node1)][self.nodes.index(node2)] = 0.0001    
 
                    elif set(o_n).intersection(set(o_n_2)) == set():
                        #print 'null intersection'
                        if len(o_n) != 0 and node1.concept == o_n[len(o_n) - 1]:
                            #print 'node1.concept is frag.last'
                            if len(o_n_2) != 0 and node2.concept == o_n_2[0]:
                                #print 'node2.concept  is frag.first'
                                mat[self.nodes.index(node1)][self.nodes.index(node2)] = self.traveling_cost(node1, node2)
                                #print str(mat[self.nodes.index(node1)][self.nodes.index(node2)]) + ' @ '  + str(self.nodes.index(node1)) + ', ' + str(self.nodes.index(node2))

                    else: mat[self.nodes.index(node1)][self.nodes.index(node2)] = float(1000)
        
        mat[self.nodes.index(n_start)][self.nodes.index(n_end)] = float(1000)
        mat[self.nodes.index(n_end)][self.nodes.index(n_start)] = 0.0001                    
        self.cost_matrix = mat
        return self.cost_matrix

    def __init__(self, sentence=None): #sentence must be from rule file
        if sentence is not None:
            self.sentence = sentence
            self.rules = sent_rule_map[sentence]
            self.concepts = sent_con_map[sentence]
            self.groups = len(self.concepts)
            self.nodes = node_factory(sentence)
            self.size = len(self.nodes)
            self.AMR = sent_amr_map[sentence]
            self.cost_matrix = self.build_cost_matrix()
        
        else:
            self.sentence = None
            self.rules = None
            self.rules = None
            self.concepts = None
            self.groups = None
            self.nodes = None
            self.size = None
            self.cost_matrix = None
        

    def build_edges(self):
        for node_i in self.nodes:
            for node_j in self.nodes:
                if node_i != node_j:
                    e = Edge(node_i, node_j, self.cost_matrix[node_i][node_j])
                    node_i.edges_out.append(e)
                    node_j.edges_in.append(e)

    def update_graph_size(self):
        self.size = len(self.nodes)