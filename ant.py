import math
from numpy import random
import sys
#from threading import *
from graph import Node

class Ant:
    def __init__(self, ID, start_node, colony, beta=None, q0=None, rho=None):
        #Thread.__init__(self)
        self.ID = ID
        self.start_node = start_node
        self.colony = colony
        self.concepts = []

        self.curr_node = self.start_node
        self.graph = self.colony.graph
        self.concepts.append(self.colony.graph.nodes[start_node].concept)
        self.path_vec = []
        self.path_vec.append(self.start_node)
        self.path_cost = 0

        # same meaning as in standard equations
        if beta == None:
            self.Beta = 2
        #self.Q0 = 1  # Q0 = 1 works just fine for 10 city case (no explore)
        if q0 == None:
            self.Q0 = .9
        if rho == None:    
            self.Rho = .1

        # store the nodes remaining to be explored here
        self.nodes_to_visit = {}

        for i in range(0, self.graph.num_nodes):
            if i != self.start_node and i not in self.colony.dead_ends:
                self.nodes_to_visit[i] = i    

        # create n X n matrix 0'd out to start
        self.path_mat = []

        for i in range(0, self.graph.num_nodes):
            self.path_mat.append([0]*self.graph.num_nodes)

    # overide Thread's run()
    def run(self):
        graph = self.colony.graph
        while not self.end():
            # we need exclusive access to the graph
            #graph.lock.acquire()
            #print 'graph lock acquired: ant.run'

            new_node = self.state_transition_rule(self.curr_node) 

            self.path_cost += graph.delta(self.curr_node, new_node)

            self.path_vec.append(new_node)
            self.path_mat[self.curr_node][new_node] = 1  #adjacency matrix representing path

            #print "Ant %s : %s, %s" % (self.ID, self.path_vec, self.path_cost,)
            
            self.local_updating_rule(self.curr_node, new_node) 
            #graph.lock.release()
            #print 'graph lock released: ant.run'

            self.curr_node = new_node

        # don't forget to close the tour
        self.path_cost += graph.delta(self.path_vec[-1], self.path_vec[0])

        # send our results to the colony
        self.colony.update(self)
        #print "Ant thread %s terminating." % (self.ID,)


        # allows thread to be restarted (calls Thread.__init__)
        self.__init__(self.ID, self.start_node, self.colony)

    def end(self):
        return not self.nodes_to_visit

    # described in report -- determines next node to visit after curr_node
    def state_transition_rule(self, curr_node):
        graph = self.colony.graph
        q = random.random()
        max_node = -1

        if q < self.Q0:
            #print "Exploitation"
            max_val = -1
            val = None

            for node in self.nodes_to_visit:
                if graph.tau(curr_node, node) == 0:
                    raise Exception("tau = 0")

                val = graph.tau(curr_node, node) * math.pow(graph.etha(curr_node, node), self.Beta)
                if val > max_val:
                    max_val = val
                    max_node = node
        else:
            #print "Exploration"
            sum = 0
            node = -1

            for node in self.nodes_to_visit:
                if graph.tau(curr_node, node) == 0:
                    raise Exception("tau = 0")
                sum +=  graph.tau(curr_node, node) * math.pow(graph.etha(curr_node, node), self.Beta)
            #if sum == 0:
                #raise Exception("sum = 0")

            #avg = sum / len(self.nodes_to_visit)

            #print "avg = %s" % (avg,)

            p = {}
            for node in self.nodes_to_visit.values():
                    p[node] = graph.tau(curr_node, node) * math.pow(graph.etha(curr_node, node), self.Beta) / sum
                   
            max_node = random.choice(p.keys(), p=p.values())

            if max_node == -1:
                max_node = node
        
        if max_node < 0:
            raise Exception("max_node < 0")

        del self.nodes_to_visit[max_node]
        self.concepts.append(graph.nodes[max_node].concept)

        for i in self.nodes_to_visit.values():
            if graph.nodes[i].concept in self.concepts or i in self.colony.dead_ends:
                del self.nodes_to_visit[i]
    
        
        return max_node

    # phermone update rule for indiv ants
    def local_updating_rule(self, curr_node, next_node):
        graph = self.colony.graph
        val = max(graph.tau0, (1 - self.Rho) * graph.tau(curr_node, next_node) + (self.Rho * graph.tau0))
        graph.update_tau(curr_node, next_node, val)

