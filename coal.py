"""
Date: 5/29/19

@author: tuk32868
"""
import random
import numpy as np
from collections import defaultdict
import math
import classes 

class coal_methods(classes.tree):
    def discrete_coal(self):
        """
            create tree based on discrete model
        """
        nodelist = self.notcoal_branches()
        list2 = list(nodelist)
        random.shuffle(list2) # random shuffle of nodelist to decrease bias
        while len(list2) >= 2: # if nodelist has more than one node in it
            x = list2[0] #
            list2.remove(x)
            templist = list(list2)
            for y in templist:
                prob = random.random()
                if prob < (1/self.popsize):
                    ## set new branch
                    nextlabel = self.n
                    self[nextlabel] = classes.branch(id =nextlabel, child = [x, y], par = -1, t = 0.0)
                    self.n +=1
                    self.p -=1
                    ##update old branches
                    self[x].update(par = nextlabel)
                    self[y].update(par = nextlabel)
                    ## get rid of a branch from the possible combinations once its coalescenced
                    list2.remove(y)
                    break
        self += 1
        self.totalgen +=1

        #returns 1 if just single root (coalescenced)
        if len(nodelist) == 1:
            keylist = self.keys()
            self.rootid = max(keylist)
            return 0
        else:
            return 1

    def continuous_coal(self, mutrate):
        """
        create tree based on continuous model
        """
        ## get nodes that have not coalescenced
        nodelist = self.notcoal_branches()

        # pull the time to next event from exponential distribution
        #print(self.p)
        t = np.random.exponential(scale = (4*self.popsize)/(self.p * (self.p -1)))


        ## add time to tree
        self += t
        self.totalgen += t

        ## pull two random non coalescenced nodes from list
        node1, node2 = random.sample(nodelist, 2)

        ## create new branch
        nextlabel = self.n
        self[nextlabel] = classes.branch(id =nextlabel, child = [node1, node2], par = -1, t = 0.0, totalpos = self[node1].totalpos)
        self.n +=1
        self.p -=1

        ##update old branches
        self[node1].update(par = nextlabel)
        #num1 = self[node1].mutate(mutrate)
        self[node2].update(par = nextlabel)
        #num2 = self[node2].mutate(mutrate)
        #self.totalmut += (num1 + num2)

        if self.p == 1:
            keylist = self.keys()
            self.rootid = max(keylist)
            for x in keylist:
                num = self[x].mutate(mutrate)
                self.totalmut += num
            return 0
        else:
            return 1
