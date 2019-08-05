
"""
    Implements a bifurcating ultrametric rooted tree

    two main classes:

    edge class
        an edge on a tree
    tree class
        a subclass of dictionary
        tree is always ultrametric (all external edges come up to the same time point)
        tree is always rooted

        growrandom(self,addtime):
            add time to all external edges, and then split one random external edge into two

        getedgetotaltime(self,e):
            return the time at the bottom of edge e
            have to sum all the descendant edges

        newick(self):
            return a newick string for self

        sumbranchlengths(self):
            return a simple sum of branch lengths

        checktree(self):
            does several checks of integrity
            kind of slow

        gettreestats(self):
            print length and edge count summaries

"""
import os
import sys
import random
import math
import matplotlib.pyplot as plt


class edge:
    """
        an edge on the tree
    """

    def __init__(self,id=-1, up = None, dn = -1, t = 0.0):
        """
            up is a list of the descendants, can be empty
            dn is the down ancestor
            t is the length of the edge  (not the time at the bottom)
        """
        self.id = id
        if up != None:
            self.up = list(up)
        else:
            self.up = []
        self.dn = dn
        self.t = max(0.0,t)

    def update(self,id = None, upedge = None, up = None, dn = None, t = None):
        """
            for updating any of the basic variables belonging to an edge
        """
        if id != None:
            self.id = id
        if upedge != None: ## add the upedge to the list up
            assert upedge not in self.up
            self.up.append(upedge)
        if up != None:
            self.up = up
        if dn != None:
            self.dn = dn
        if t != None:
            self.t = max(0.0,t)

    def isexternal(self):
        """
            True if no daughter edges
        """
        return self.up == []

#    def __str__(self):
#        """
#            create a string with information for using print()
#        """
#        s =  "id %d nup %d "%(self.id,self.nup)
#        if self.up == []:
#            s += "[]"
#        else:
#            s += '['
#            for i in range(self.nup):
#                if i < self.nup-1:
#                    s += "%d,"%self.up[i]
#                else:
#                    s += "%d"%self.up[i]
#            s += ']'
#        s += " dn %d"%self.dn
#        s += " time %.3f"%self.t
#        return s

class tree(dict):
    """
        a tree class, is a subclass (inherits from) the dictionary type,
        contains a phylogenetic tree
        positive integers are used as labels
    """

    def __init__(self):
        """
            initialize the dictionary with a single root edge
        """
        self.tlen = 0.0
        self.tmrca = 0.0
        label = 0
        self[label] = edge(id = label,up = [],dn = -1,t = 0.0)
        self.nextlabel = label + 1
        self.rootid = label
        self.numexternal = 1 # number of external edges


    def __add__(self, addtime):
        """
        overload addition +=  add addtime to all extern edges
        """
        for key in self.keys():
            if self[key].isexternal():
                self[key].t += addtime
        self.tlen += self.numexternal * addtime ## increment the length of the tree
        self.tmrca += addtime
        return self


    def growrandom(self, rand, addtime = None):
        """
            if rand is true,
                add random amount of time to tree from 0 to 1
            if rand is false, 
                use fixed amount of time 
        
            if tree has more than one edge
                add time to all external edge lengths and then

            if only one edge,  this is the root and do not add time
            randomly pick an external edge to split into two daughter edges
        """
        def add_daughters ():
            """
                pick an external edge at random and make two daughters
            """
            ##print(self)
            while True:
                x = random.choice(list(self.keys()))
                if  self[x].isexternal():
                    break
            self[self.nextlabel] = edge(id = self.nextlabel,up = [], dn = x,t = 0.0)

            tempup = [self.nextlabel]
            self.nextlabel += 1
            self[self.nextlabel] = edge(id = self.nextlabel,up = [], dn = x,t = 0.0)
            tempup.append(self.nextlabel)
            self.nextlabel += 1
            self[x].update(up = tempup)
            self.numexternal += 1
            return

        if rand == False:
            if len(self) > 1:
                self += addtime;
            add_daughters()
            return
        if rand == True:
            random.seed()
            addtime = random.random()
            if len(self) > 1:
                self += addtime;
            add_daughters()
            return 
    def extinction(self):
        ##randomly select edge to be deleted
        ##make sure its not the root
        list_keys = list(self.keys())
        nonroot = list_keys[1:]
        deled = random.choice(nonroot)
        ##create counter to count number of deletions
        counter = 0
        ##print("deleted edge: " + str(deled))
        ##delete edge and its descendants
        if self[deled].isexternal != True:
            ##create set to list all descendants
            dels = set()
            ##for each edge, if external map back to root
            for x in self.keys():
                if self[x].isexternal():
                    path = [x]
                    ##while ancestor isn't root add it to the path list
                    anc = self[x].dn
                    while anc > -1:
                        path.append(anc)
                        anc = self[anc].dn
                    ##if deleted edge is in path list, if descendant add to delete list
                    if deled in path:
                        for item in path:
                            if path.index(item) <= path.index(deled):
                                dels.add(item)  
            ##create list from set                     
            dell = list(dels)
        else: 
            ##if external, just need to add the edge to delete list 
            dell.append(deled)
        ##join sister and down edge 
        ##in list of descendants of mother edge (to one being deleted) is the edge the first or 2nd
        if self[self[deled].dn].up[0] == deled:
            ##set sister to second if deleted edge is first
            sister = self[self[deled].dn].up[1]
        else:
            ##set sister to first if deleted edge is second 
            sister = self[self[deled].dn].up[0]
        ##set mother edge/ancestor to the ancestor (dn) of sister edge
        ancestor = self[sister].dn
        ##if ancestor is root, set sister as new root
        if self[ancestor].dn == -1:
            self[sister].update(dn = -1, t = 0)
            ##print(self[ancestor])
            ##deleter ancestor, add one to counter 
            del(self[ancestor])
            counter +=1
        else:
            ##if ancestor is not root, update sister with grandpa as new dn and update time 
            self[sister].update(dn = self[ancestor].dn, t = self[sister].t + self[ancestor].t)
            ##update grandfather to have sister as new descendant and remove deleted ancestor 
            if self[self[ancestor].dn].up[0] == ancestor:
                ##print(self[self[ancestor].dn].up[0])
                del(self[self[ancestor].dn].up[0])
                self[self[ancestor].dn].update(upedge = sister)
            else:
                ##print(self[self[ancestor].dn].up[1])
                del(self[self[ancestor].dn].up[1])
                self[self[ancestor].dn].update(upedge = sister)
            ##print(self[ancestor].id)
            ##delete ancestor and add to counter 
            del(self[ancestor])
            counter +=1
        ##for each value in delete list, delete that edge and add 1 to counter 
        for x in dell:
            ##print(self[x])
            del(self[x])
            counter+=1
        ##counter is divided by 2 because grow random as two for every call 
        counter = counter/2
        ##return int of counter to call that number of grow random 
        return int(counter) 

            
    def getedgetotaltime(self,e):
        """
            return the time at the bottom of edge e
            have to sum all the descendant edges
        """
        tempt = self[e].t
        while self[e].up != []:
            e = self[e].up[0]
            tempt += self[e].t
        return tempt

    def newick(self):
        """
            return a newick string for self

        """

        def innernewick(self,tempedge):
            """
                recursive
            """
            if self[tempedge].isexternal():
                return str(tempedge)
            else:
                tempnodeinfolist = []
                for i in range(2):
                    tempupid  = self[tempedge].up[i]
                    tempnodestr = innernewick(self,tempupid)
                    temptime = self[self[tempedge].up[i]].t
                    tempnodeinfolist.append(tempnodestr + ":%.5f"%(temptime))
                tempbuildnode = "(" + tempnodeinfolist[0] + ',' + tempnodeinfolist[1] + ')' + str(tempedge)
                return tempbuildnode
            return

        newickstr = innernewick(self,self.rootid) + ";"
        return newickstr

    def sumbranchlengths(self):
        """
            return a simple sum of branch lengths

        """
        sumt = 0.0
        for e in self:
            sumt += self[e].t
        return sumt

    def checktree(self):
        """
            does several checks of integrity
            kind of slow
        """
        sumt = 0.0
        countexternal = 0

        for e in self:
            if e != self.rootid:
                if (e == self[self[e].dn].up[0] or e == self[self[e].dn].up[1]) is False: # check that an up edge and its dn ege identify each other
                    print (self[e])
                    raise Exception(" up dn failure in checktree:  edge %d"%(e))
                if self[e].up != [] and self[e].up[0] == self[e].up[1]:  # check that both up edges are different
                    raise Exception (" up failure for edge %d:  up[0] %d up[1] %d"%(e,self[e].up[0],self[e].up[1]))
            else:  # its the root and should not have a dn edge or a time
                if self[e].dn != -1:
                    raise Exception (" rootid: %d  has a down edge"%(self.rootid))
                if self[e].t != 0.0:
                    raise Exception (" rootid: %d  has a t value"%(self.rootid))
        for e in self:   # loop over all edges,  sum length, and check that each edge's path to the root has length = tmrca
            sumt += self[e].t
            if self[e].isexternal():
                countexternal += 1
                d = e
                tmrcacheck = 0.0
                while self[d].dn != -1:
                    tmrcacheck += self[d].t
                    d = self[d].dn
                if math.isclose(tmrcacheck,self.tmrca) is False:
                    raise Exception("tmrca check failed: tmrcacheck %.4f  self.tmrca %.4f"%(tmrcacheck,self.tmrca))
        if math.isclose(sumt,self.tlen) is False:
            raise Exception("tlen check failed: sumt %.4f  self.tlen %.4f"%(sumt,self.tlen))
        if countexternal != self.numexternal:
            raise Exception("external count failed: count %d  self.numexternal %d"%(countexternal,self.numexternal))

    def gettreestats(self):
        """
            print length and edge count summaries
        """
        sumexlen = 0.0
        sumintlen = 0.0
        exc = intc = 0
        for e in self:
            if self[e].isexternal():
                exc += 1
                sumexlen += self[e].t
            else:
                intc += 1
                sumintlen += self[e].t
        try:
            if intc > 0:
                info = "Length : %.3f  TMRCA: %.3f #edges: %d  #internal: %d (mean length: %.4f)  #external: %d (mean length: %.4f)"%(self.tlen,self.tmrca,len(self),intc,sumintlen/intc,exc,sumexlen/exc)
            else:
                assert exc == 1
                info = "Length : %.3f  TMRCA: %.3f #edges: %d  #internal: %d #external: %d (mean length: %.4f)"%(self.tlen,self.tmrca,len(self),intc,exc,sumexlen/exc)
        except Exception:
            info = "problem getting stats"
        return info

    def __str__(self):
        s =  self.gettreestats() + "\n"
        s += self.newick() + "\n"
        return s
    
    def simultree(self, exedge, random, time = None):
        """
            creates a tree with a certain number of external edges (exedge)
            with either random or non-random growth (random)
            if non-random growth (time codes for fixed growth amount)
        """
        ##if random is true, randomly add time for the number of edges (exedge)
        if random == True: 
            for i in range(exedge):
               self.growrandom(random)
        else:
            ##not random add number of edges for specified time 
            for i in range(exedge):
                self.growrandom(random, time)
        ##print tree 
        print(self)
        
    
    def burnin(self, cycles):
        """ 
        random extinction of one or more edges, addition of edges for number lost, for cycles 
        """
        ##for number of cycles perform pruning and growth cycles 
        for x in range(cycles):
            k = self.extinction()
            for x in range(k):
                self.growrandom(True)
        print("End of burning in stage")
        
    def sampling(self, cycles):
        """
        random extinction and growth for x cycles, every 10 cycles find tree length to be returned as list
        """
        ##for number of cycles perform pruning and growth cycles as in burnin 
        counter = 0
        length = []
        for x in range(cycles):
            k = self.extinction()
            counter+=1
            for x in range(k):
                self.growrandom(True)
            ##when counter is 9 (10 total) find total tree length by adding up all edges length, reset counter
            if counter == 9:
                leng = 0 
                for x in self.keys():
                    leng += self[x].t
                length.append(leng)
                ##print(self.tlen)
                counter = 0
        ##print(self)
        return length
##random seed        
random.seed(11)
##create tree, simulate tree with 20 external edges, burnin cycle 100, sampling cyle 10,000
t = tree()
t.simultree(20, True)        
t.burnin(100)
leng = t.sampling(10000)

##create histogram with x and y labels and title 
plt.hist(leng)
plt.title("Distribution of Tree Lengths")
plt.xlabel("Tree Lengths")
plt.ylabel("Frequency")

##output as png file 
fig = plt.gcf()
fig.savefig("Distribution.png")
plt.close(fig)
   
        