CAST-AMR: The Cooperative Ant System Translator for Abstract Meaning Representations
=================

CAST-AMR generates English text from [Abstract Meaning Representations](http://amr.isi.edu/) using algorithms inspired by the foraging behavior of ants.


#Building

Provided are a few data files that you can use to generate rules:
    
* `.\dev`   ---   dev data directory containing AMRs, alignments, etc. \n
* `.\devout` ---  contains rules generated from dev set
* `.\test`  ---   test data directory
* `.\testout` --- contains rules generated from dev set

If you'd like to generate your own rules, consult the files in dev or test to see the data you need. You can generate AMRs and alignments of your own with J. Flanigan's [JAMR](https://github.com/jflanigan/jamr) parser/generator. 


Once you have the necessary files, you open and edit:

    ./run.sh 

and set inputDir= and outputDir= to your desired input and output directories. Save, close, then run the script, and your rule data should be generated 

With rule generation completed, you can build the graph structures that ants can start working on! Simply modify the code at the top of graph_utils.py if you want to use your own rules. Otherwise, you can use the built-in dev and test data.

#Running
CAST-AMR works very well in Python's interactive shell.

First we import dependencies:

    from antcolony import *
    from antgraph import AntGraph
    from graph import * 
    build_maps(rules, amrs)  #use your rules/amrs, dev_rules/amrs or test_rules/amrs

You may be interrupted by a message saying "building LMs", just wait for it to finish .

To build a graph from the 3rd sentence in our list:

    sent = sent_rule_map.keys()[2] #choose 3rd sentence in data
    g = Graph(sent) #takes a sentence in your data
    ag = AntGraph(g) #takes graph object

To translate a sentence:    
    
    ac = AntColony(ag, 10, 50) #takes antgraph object, # of ants, # of iterations
    ac.start()

The ant colony print the best tour data of each iteration until it reaches the iteration cutoff specified by the user. From empirical recommend between 50 and 100 iterations, erring on the lower side--convergence seems to happen pretty quickly! When system is done, you should see best tour information, a translation, and a BLEU evaluation score of the translation from 0 to 1. 

To translate another sentence, update sent and proceed as normal:
    sent = sent_rule_map.keys()[0] #choose 1st sentence in data
    g = Graph(sent) 
    ag = AntGraph(g) 
    ac = AntColony(ag, 10, 50) 
    ac.start()

To exit the CAST-AMR and the interactive shell, press ctrl + z.
    

Acknowledgements
=================== 
We would like to thank: 

Xiaochang Peng at the University of Rochester for providing us with the code and data to make this project happen
Trev Lovett (github: trevlovett), for his ACS implementation here: https://github.com/trevlovett/Python-Ant-Colony-TSP-Solver
Zak Kincaid, my adviser, for his good questions and insights which helped shape this project
Terrace F. Club
Mom