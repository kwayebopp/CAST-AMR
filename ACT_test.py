from graph import *
from antcolony import *
from antgraph import *

import random

print "building maps"
build_maps(test_rules, test_amrs)
print "maps completed"

def BaselineTest():
    with open('ACT_q0_experiment', 'a') as outfile: 
        outfile.write('#------------BASELINE TEST------------#\n')
    
    sents = [] #['The only thing that surprises me is how rapidly this is happening .', 'The reality', 'It \'s the same old problem .', 'It \'s ok as long as there is something to eat :loveliness :', 'Viegas stated that the limitation was necessary .']
    while len(sents) < 100:
        sent = sent_rule_map.keys()[random.randint(0, len(sent_rule_map.keys()) - 1)]
        while len(sent.split()) > 30:
            sent = sent_rule_map.keys()[random.randint(0, len(sent_rule_map.keys()) - 1)]
            if len(sent.split()) > 0 and len(sent.split()) < 30 and sent.split() not in sents:
                sents.append(sent)
                break
    for sent in sents:
        print 'building AGTSP'  
        g = Graph(sent)
        print 'AGTSP built'
        print 'Building ant graph'
        ag = AntGraph(g)
        print 'ant graph complete'

        ac = AntColony(ag, 10, 30)
        print "running ants"
        with open('ACT_q0_experiment', 'a') as outfile: 
            outfile.write('%s ' % (sent))
        ac.start()
        print "translation acquired"

       # with open('ACT_q0_experiment', 'a') as outfile: 
       #     outfile.write('\n')
        
        


def main():
    BaselineTest()



if  __name__ =='__main__':
    main()