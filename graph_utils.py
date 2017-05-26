from collections import OrderedDict
import graph
import pickle

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
                n = graph.Node(concept, rule)
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

def build_graphs(rule_file, amr_file):
    global sent_amr_map, sent_con_map, sent_rule_map, sent_node_map, sent_graph_map
    build_maps(rule_file, amr_file)
    sr = sent_rule_map
    sc = sent_con_map
    sn = sent_node_map
    sa = sent_amr_map

    for s in sr.keys():
        g = graph.Graph()
        g.sentence = s
        g.rules = sr[s]
        g.concepts = sc[s]
        g.groups = len(g.concepts)
        g.nodes = sn[s]
        g.cost_matrix = g.build_cost_matrix()
        g.amr = sa[s]
        sent_graph_map[g.sentence] = g 
        
# def clean_sense_tag(con):
#     con = con.split('-')
#     try:
#         if int(con[len(con) - 1]) >= 0 and int(con[len(con) - 1]) != 91: #check for sense tag
#             con = ' '.join(con[: len(con) - 1]) #drop it
#             return con
#         elif int(con[len(con) - 1]) == 91:
#             return 'code-91'
#     except:
#             con = ' '.join(con)   
#             return con

