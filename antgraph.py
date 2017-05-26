#from threading import Lock

class AntGraph:
    def __init__(self, agtsp_graph, tau_mat=None):
        # print len(delta_mat)
        # if len(delta_mat) != num_nodes:
        #     raise Exception("len(delta) != num_nodes")

        self.num_nodes = agtsp_graph.size
        self.delta_mat = agtsp_graph.cost_matrix 
        self.nodes = agtsp_graph.nodes
        self.num_concepts = len(agtsp_graph.concepts)
        self.concepts = agtsp_graph.concepts
        self.sentence = agtsp_graph.sentence
        #self.lock = Lock()

        # tau mat contains the amount of phermone at node x,y
        if tau_mat is None:
            self.tau_mat = []
            for i in range(0, self.num_nodes):
                self.tau_mat.append([0]*self.num_nodes)
        
        self.reset_tau()

    def delta(self, r, s):
        return self.delta_mat[r][s]

    def tau(self, r, s):
        return self.tau_mat[r][s]

    # 1 / delta = eta or etha 
    def etha(self, r, s):
        if self.delta_mat[r][s] == None:
            print r + " " + s
            print self.num_nodes
        try:
            return 1.0 / self.delta(r, s)
        except ZeroDivisionError:
            return float(1000)

    # inner locks most likely not necessary
    def update_tau(self, r, s, val):
        #lock = Lock()
        #lock.acquire()
        #print 'graph lock acquired: ag.update_tau'
        self.tau_mat[r][s] = val
        #lock.release()
        #print 'graph lock released: ag.update_tau'

    def reset_tau(self):
        #lock = Lock()
        #lock.acquire()
        #print 'graph lock acquired: ag.reset_tau'
        avg = self.average_delta()

        # initial tau 
        self.tau0 = 1.0/ (avg * self.num_concepts) #rough estimate of tour length = num_concepts * avg

        print "Average = %s" % (avg,)
        print "Tau0 = %s" % (self.tau0)

        for r in range(0, self.num_nodes):
            for s in range(0, self.num_nodes):
                self.tau_mat[r][s] = self.tau0
        #lock.release()
        #print 'graph lock released: ag.reset_tau'

    # average delta in delta matrix
    def average_delta(self):
        return self.average(self.delta_mat)

    # average tau in tau matrix
    def average_tau(self):
        return self.average(self.tau_mat)

    # average val of a matrix
    def average(self, matrix):
        sum = 0
        ct = 0
        for r in range(0, self.num_nodes):
            for s in range(0, self.num_nodes):
                if  r != s and matrix[r][s] != float(1000):
                    ct += 1
                    #print "r: "+ str(r) + " s: " + str(s) + " val: " + str(matrix[r][s])
                    sum += matrix[r][s]

        avg = sum / (ct)
        return avg

