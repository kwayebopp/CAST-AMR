from ant import Ant
#from threading import Lock, Condition

import random
import sys
import bleu_score

class AntColony:
    def __init__(self, graph, num_ants, num_iterations):
        self.graph = graph
        self.num_ants = num_ants
        self.num_iterations = num_iterations
        self.Alpha = 0.1
        self.dead_ends = []
        

        # condition var
        #self.cv = Condition()

        self.reset()

    def reset(self):
        self.best_path_cost = sys.maxint
        self.best_path_vec = None
        self.best_path_mat  = None
        self.last_best_path_iteration = 0
        self.dead_ends = []

    def start(self):
        self.ants = self.create_ants()
        self.iter_counter = 0

        while self.iter_counter < self.num_iterations:
            self.iteration()

           # self.cv.acquire()
            #print 'colony.cv.acquire in ac.start'
            # wait until update calls notify()

            #print 'colony.cv.waiting'
            #self.cv.wait()
            #print 'colony.cv done waiting'


            #lock = self.graph.lock
            #lock.acquire()
            #print 'graph lock acquired: ac.start'
            self.global_updating_rule()
            #lock.release()
            #print 'graph lock released: ac.start'

            #self.cv.release()
            #print 'colony.cv.release in ac.start'
        

    # one iteration involves spawning a number of ant threads
    def iteration(self):
        self.avg_path_cost = 0
        self.ant_counter = 0
        self.iter_counter += 1
        print "iter_counter = %s" % (self.iter_counter,)
        for ant in self.ants:
            #print "starting ant = %s" % (ant.ID)
            ant.run()

    def num_ants(self):
        return len(self.ants)

    def num_iterations(self):
        return self.num_iterations

    def iteration_counter(self):
        return self.iter_counter

    # called by individual ants
    def update(self, ant):
        #lock = Lock()
        #lock.acquire()
        #print 'lock acquired: ac.update'

        #print "Update called by %s" % (ant.ID,)
        self.ant_counter += 1

        self.avg_path_cost += ant.path_cost

        # book-keeping
        if ant.path_cost < self.best_path_cost:
            self.best_path_cost = ant.path_cost
            self.best_path_mat = ant.path_mat
            self.best_path_vec = ant.path_vec
            self.last_best_path_iteration = self.iter_counter

        if self.ant_counter == len(self.ants):
            self.avg_path_cost /= len(self.ants)
            print "Best: %s, %s, %s, %s" % (self.best_path_vec, self.best_path_cost, self.iter_counter, self.avg_path_cost)
            #frags = []
            self.bleu_tester()
           #self.cv.acquire()
            #self.cv.notifyAll()
            #self.cv.release()
        #outfile.close()
        #lock.release()
        #print 'lock released: ac.update'

    def bleu_tester(self):
        trans = []
        t = ""
        for concept in self.graph.concepts:  
            for item in self.best_path_vec:
                if concept == self.graph.nodes[item].concept and self.graph.nodes[item].rule.split('|||')[1].strip() not in trans:
                    trans.append(self.graph.nodes[item].rule.split('|||')[1].strip())
                    break
            
        t = " ".join(trans)
        print t
        print bleu_score.sentence_bleu([t], self.graph.sentence)
        print trans
        if self.iter_counter == self.num_iterations:
            outfile = open("ACT_q0_experiment", "a")
            outfile.write("%d\t%f\t%f\t%s\n" % (self.iter_counter, self.best_path_cost, bleu_score.sentence_bleu([t], self.graph.sentence), t))
            outfile.close()       
      

    def done(self):
        return self.iter_counter == self.num_iterations

    # assign each ant a random start-node
    def create_ants(self):
        self.reset()
        ants = []
        for i in range(0, self.num_ants):
            ant = Ant(i, self.graph.num_nodes - 2, self)
            ants.append(ant)
        
        return ants

    # changes the tau matrix based on evaporation/deposition 
    def global_updating_rule(self):
        evaporation = 0
        deposition = 0

        for r in range(0, self.graph.num_nodes):
            for s in range(0, self.graph.num_nodes):
                if r != s:
                    delt_tau = self.best_path_mat[r][s] / self.best_path_cost
                    evaporation = (1 - self.Alpha) * self.graph.tau(r, s)
                    deposition = self.Alpha * delt_tau

                    self.graph.update_tau(r, s, evaporation + deposition)
