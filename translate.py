from antcolony import *
from antgraph import *
from graph import *
sent_graph_map = OrderedDict()
sent_trans_map = OrderedDict()

def build_sent_graph_map(rules, amrs):
    global sent_graph_map
    
    build_maps(rules, amrs)
    
    for sent in sent_rule_map:
        g = Graph(sent)
        ag = AntGraph(g.size, g.nodes, g.cost_matrix)
        sent_graph_map[sent] = ag

def test_graphs():
    global sent_trans_map
    for sent in sent_graph_map:
        ac = AntColony(sent_graph_map[sent], 10, 3500)
        trans = ac.start()
        sent_trans_map[sent] = trans
    
    return sent_trans_map
        