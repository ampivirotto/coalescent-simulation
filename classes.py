# -*- coding: utf-8 -*-
"""
Date: 5/29/19

@author: tuk32868
"""
import random
import numpy as np
from collections import defaultdict
import math

class branch:
    """
    structure to contain branch of tree
    """
    def __init__(self,id=-1, child = None, par = -1, t = 0.0, mut = 0, totalpos = 0):
        """
            child is the descendant nodes/branches
            par is the parent node/branch
                if par is -1 branches have not coalescenced or its at the root
            t is the length of the edge  (not the time at the bottom)
            mut is the number of mutations on the branch
        """
        self.id = id
        if child != None:
            self.child = list(child)
        else:
            self.child = []
        self.par = par
        self.t = max(0.0,t)
        self.mut = mut
        self.totalpos = max(0.0, totalpos)

    def update(self,id = None, child = None, par = None, t = None, mut = None, totalpos = None):
        """
            update a branch
        """
        if id != None:
            self.id = id
        if par != None:
            self.par = par
        if child != None:
            self.child = child
        if t != None:
            self.t = t
        if mut != None:
            self.mut += mut
        if totalpos != None:
            self.totalpos = totalpos

    def hasnotcoalescenced(self):
        """
            returns True if branch has not coalescenced
        """
        return self.par == -1

    def isexternal(self):
        """
            returns True if branch is external
        """
        return self.child == []

    def mutate(self, mutrate):
        """
            add mutations to branch based on Poisson distribution (Length of branch * mutation)
        """
        nummutate = np.random.poisson((self.t * mutrate))
        #print(nummutate)
        self.update(mut = nummutate)
        return nummutate


class tree(dict):
    """
    data container for a tree
    """
    def __init__(self, popsize, samplesize, seed):
        self.tlen = 0.0 ## total length of tree
        self.tmrca = 0.0 ## time of most recent common ancestor (number of gen to coalescent)
        self.totalgen = 0 # total number of generations
        for x in range(0, samplesize):
            self[x] = branch(id=x, par = -1, t = 0.0) ## initialize sample branches
        self.rootid = -1 # set root id, it'll be updated when the tree finishes
        self.n = samplesize # number of edges
        self.p = samplesize # number of edges with no parent
        self.popsize = popsize # initialize population size
        self.samplesize = samplesize
        self.totalmut = 0 # total number of mutations across tree
        self.interval = None
        random.seed(seed) ## initialize random seed

    def __add__(self, addtime):
        """
        overload addition +=  add addtime to all extern edges
        """
        for key in self.keys():
            if self[key].hasnotcoalescenced():
                self[key].t += addtime ## add time to each branch that has not coalescenced
                self[key].totalpos += addtime
        self.tlen += self.p * addtime ## increment the length of the tree
        self.tmrca += addtime ## add time to the
        return self

    def notcoal_branches(self):
        """
            return list of nodes that have not coalescenced
        """
        ##identify the branches that do not have parents
        nodelist = []
        for x in self.keys():
            if self[x].hasnotcoalescenced():
                nodelist.append(x)
        ## return these nodes as a list
        return nodelist

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
                    self[nextlabel] = branch(id =nextlabel, child = [x, y], par = -1, t = 0.0)
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
        self[nextlabel] = branch(id =nextlabel, child = [node1, node2], par = -1, t = 0.0, totalpos = self[node1].totalpos)
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

    def coal_events_dict_creation(self):
        """
        create dictionary of coalescent events
        """
        ## create empty dictionary
        dict1 = defaultdict(list)

        ## cycle through each branch
        for x in self.keys():
            ## for just the internal branches - thus branches created by coalescent event
            if not self[x].isexternal():
                ## find one of the children of selected branch
                child1 = self[x].child[0]

                ## check to make sure code is working correctly
                if child1 == None:
                    break

                ## the time of coalescent event is the total time of child1
                time_event = self[child1].totalpos
                if round(time_event) == 0:
                        print("possible error")

                ## cycle through branches again to find branches that fall within time event
                for x in self.keys():

                    # set end time of each branch
                    end = self[x].totalpos

                    # set start time of each branch
                    start = end - self[x].t

                    # add the branch to the dict entry if the start time falls before the coalescent event and
                    #end time is equal to or after the coal event
                    if round(time_event, 9) > round(start, 9) and round(time_event, 9) <= round(end, 9):
                        dict1[time_event].append(x)
        if dict1 == None:
            raise Exception("Error: Coalescent Dictionary is empty.")
        return dict1

    def find_k(self, c_event_dict, pos):
        """
        identify the number of branches (k) at a certain position of the tree
        """
        # set k to 0
        k = 0
        # set temp to 0
        temp = 0
        # set the event to 0
        event = 0

        # cycle through sorted keys so that you go lowest to highest
        for x in sorted(c_event_dict.keys()):
            # check to see if passed pos is higher than temp and lower than x len from dict
            if pos > temp and pos< x:
                #if it is get the number of branches k from len of list
                k = len(c_event_dict[x])
                # set event to be x time
                event = x
                # return both values
                return k, event
            #if not, set temp to now be x and cycle to next x value
            else:
                temp = x

    def remove_prior_events(self, coaldict, pos):
        update_dict = defaultdict(list)
        for x in coaldict.keys():
            if x >=pos:
                update_dict[x] = coaldict[x]
        return update_dict

    def recomb_branch(self, pos):
        """
        find location on tree where selected branch will recombine
        """
        ## create dictionary of all coalescent events
        c_event_dict = self.coal_events_dict_creation()

        if not c_event_dict:
            raise Exception("ERROR: Coalescence dictionary is empty.")

        ## remove events in coalescent dictionary prior to the position
        c_event_dict = self.remove_prior_events(c_event_dict, pos)

        ##
        counter = 0

        ## identify the earliest point that could be selected
        startpos = pos

        ## while there are still items in dictionary, do these things
        while c_event_dict:
            ## find number of branches time of next coalescent event near selected branch position
            k, event = self.find_k(c_event_dict, pos)

            ## choose time to event from exponential distribution
            time_to_event = np.random.exponential((2*self.popsize)/k) + startpos

            ## if this time is less then next coal event
            if time_to_event < event and time_to_event >= startpos:

                ## randomly select location

                ## find length between the starting position and the next coal event
                from_pos_to_event = event - startpos

                ## find the total tree distance by taking
                #this distnace times the number of branches at this position in the tree
                total_dist = k * from_pos_to_event

                ## pick  a random point along this total
                random_point = random.random() * total_dist

                ## figure out what branch this point falls on
                which_branch = 0
                xlabel = None
                ## loop through the possible branches at this point
                for x in c_event_dict[event]:
                    ## for each branch add its length from starting point to next coal event to total
                    which_branch += from_pos_to_event
                    ## if total is larger than random point, then point falls on branch
                    if which_branch > random_point:
                        ## set label for recombining branch
                        xlabel = x
                        break
                ## find position along this branch
                distance_from_end = random_point % from_pos_to_event
                xpos = self[xlabel].totalpos - distance_from_end

                ## if position is less then 1 print warning
                if xpos < 1:
                    print("warning")
                ## return branch and position
                return xlabel, xpos
            else:
                ## if it doesn't fall between start and event
                ## delete event and set new start as event
                del c_event_dict[event]
                startpos = event

        ## if dictionary is empty, you're at the root
        if not c_event_dict:
            ## if at root k branches = 1
            k = 1

            ## check to see if startpos is tmrca
            if startpos != self.tmrca:
                raise Exception("Choosing position along root for recombination error.")

            ## find time to event, exponentially picked plus startpos which is last event
            time_to_event = np.random.exponential((2*self.popsize)/k) + startpos

            ## return recombine branch (root) and new position
            return self.rootid, time_to_event

    def select_branch(self):
        """
        select branch to undergo sequential markove coalescent process
        """
        ## randomly select position where branch snip will occur
        random_len = random.random() * self.tlen

        ## find which branch and position this random selection is occuring on
        label = 0
        pos = 0
        total = 0
        ## loop through all branches
        for x in self.keys():
            ## set temp value as total before addition
            temp = total
            ## add next branches length to the total
            total += self[x].t
            ## if the total is more than random length, random point must occur on the branch
            if total > random_len:
                ## set the branch as label
                label = x
                ## find position by finding start by taking random point minus previous total (total of all previous branches)
                distance_from_start = random_len - temp
                ## find distance from end by taking length of branch minus start
                distance_from_end = self[x].t - distance_from_start
                ## find total pos of snip by taking totalpos of branch minus distance from end
                pos = self[x].totalpos - distance_from_end
                break

        ## check to make sure, the position is positive value
        if pos < 0:
            raise Exception("ERROR: negative branch length has occured.")

        ## call recombination to find branch and position of recombination
        xbranch, xpos = self.recomb_branch(pos)

        ## check to make sure recombining point is after snip position
        if xpos < pos:
            raise Exception("Error: Recombining point earlier in tree than snip position")

        ## how is the tree updated
        ##
        if label == xbranch:
            #self.check_totalpos()
            print("same")
        ## if the parent of the selected branch has a parent of -1, if the parent is the root (max(self))
        elif self[self[label].par].par == -1:
            self.root_update(label, pos, xbranch, xpos)
        elif self[label].par == xbranch:
            self.parent_update(label, xbranch, xpos, pos)
            #self.check_totalpos()
        else:
            self.nonpar_update(label, xbranch, xpos, pos)

        error = self.check_allbranches()
        if error == 0:
            print("ERROR")
        new_sum = self.sumbranchlengths()
        self.tlen = new_sum
        #self.check_totalpos()

    def root_update(self, label, pos, xlabel, xpos):
        """
        update that occurs if the root is in the mix
        """

        ## find sister branch
        temp = list(self[self[label].par].child)
        temp.remove(label)
        sister = temp[0]

        ## if the recombining branch is root
        if xlabel == self.rootid:
            ## find new sis length
            new_sis_tot_len = xpos - self[sister].totalpos + self[sister].t
            ## udpate sister
            self[sister].update(totalpos = xpos, t = new_sis_tot_len)
            ## find new selected branch length
            update_branch_len = xpos - self[label].totalpos + self[label].t
            ## update selected branch
            self[label].update(totalpos = xpos, t = update_branch_len)
            ## update the root to new totalpos
            self[self[label].par].update(totalpos = xpos)
            ## update the tmrca
            self.tmrca = xpos

        ## if the recombining branch is the sister
        elif xlabel == sister:
            ## find length from reocmbining position to end of branch
            temp_len = self[sister].totalpos - xpos
            ## find new sister branch length using the length from the end
            new_sis_len = self[sister].t - temp_len
            ## update the sister branch
            self[sister].update(totalpos = xpos, t = new_sis_len)
            ## find new selected branch length
            update_branch_len = self[label].t - temp_len
            ## update selected branch
            self[label].update(totalpos = xpos, t = update_branch_len)
            ## update the root
            self[self[label].par].update(totalpos = xpos)
            ## update the tmrca
            self.tmrca = xpos

        ## if recombining branch is not root or sister
        else:
            ## calculate new lengths
            new_empty_len = self[xlabel].totalpos - xpos
            new_recombine_len = self[xlabel].t - new_empty_len
            if self[label].isexternal():
                new_label_len = xpos
            else:
                distance_from_start = pos - (self[label].totalpos - self[label].t)
                new_label_len = xpos - pos + distance_from_start
            new_empty_total = self[xlabel].totalpos

            ## figure out new relationships
            par_xlabel = self[xlabel].par
            child_par_xlabel = list(self[par_xlabel].child)
            child_par_xlabel.remove(xlabel)
            sis_xlabel = child_par_xlabel[0]

            new_root_children = list(self[sister].child)
            if xlabel in new_root_children:
                new_root_children.remove(xlabel)
                new_root_children = [new_root_children[0], sister]
                new_sis_parent = self.rootid
                self[sis_xlabel].update(par = new_sis_parent)

            else:
                new_sis_parent = self[xlabel].par
                self[par_xlabel].update(child = [sister, sis_xlabel])


            ## calc new tmrca
            new_child_root = new_root_children[0]
            new_tmrca = self[new_child_root].totalpos
            self.tmrca = new_tmrca


            ## update new branches
            self[xlabel].update(t = new_recombine_len, totalpos = xpos, par = sister)
            self[label].update(par = sister, t = new_label_len, totalpos = xpos)
            self[sister].update(t = new_empty_len, totalpos = new_empty_total, par = new_sis_parent, child = [label, xlabel])
            self[self.rootid].update(totalpos = new_tmrca, child = new_root_children)



        ## check to make sure the root is behaving correctly
        self.root_check()

    def parent_update(self, label, xlabel, xpos, pos):
        """
        if the recombining branch is the parent of the selected branch
        """
        ## if position of recombining is the same position of selection, nothing happens
        if xpos == pos:
            #self.check_totalpos()
            print("same - wrong update tho")
        ## if the position is any other on the parent branch
        else:
            ## calculate new parent time
            new_par_time = self[xlabel].totalpos - xpos

            ## calculate new sister time
            children = list(self[xlabel].child)
            children.remove(label)
            sister = children[0]
            new_sis_time = (self[xlabel].t - new_par_time) + self[sister].t

            ## calculate new selected branch time
            new_branch_time = (self[xlabel].t - new_par_time) + self[label].t

            ## update all branches
            self[xlabel].update(t = new_par_time)
            self[sister].update(t = new_sis_time, totalpos = xpos)
            self[label].update(t = new_branch_time, totalpos = xpos)


    def nonpar_update(self, label, xlabel, xpos, pos):
        """
        update the branches when the snip position branch and the recombining branch are not parent/child
        """
        ## find parent of branch to be pruned
        parent_branch = self[label].par

        ## find children of parent branch to identify sister
        children_of_par = list(self[parent_branch].child)

        ## find sister of branch to be pruned
        children_of_par.remove(label)
        sister = children_of_par[0]

        ## update time of sister
        new_t = self[parent_branch].t + self[sister].t

        ## update sister branch to new parent  and new time
        self[sister].update(par = self[parent_branch].par, totalpos = self[parent_branch].totalpos, t = new_t)

        ## update the grandparent to new child
        temp = list(self[self[parent_branch].par].child)
        temp.remove(parent_branch)
        sis_of_par = temp[0]
        self[self[parent_branch].par].update(child = [sis_of_par, sister])

        ## remove parent branch
        del self[parent_branch]
        emptylabel = parent_branch

        ############ RECOMBINATION ################################

        if xlabel == self.rootid:
            self.recomb_at_root(label, xpos, emptylabel)
        else:
            ## find new lengths of branch x and new branch
            temp = self[xlabel].totalpos - xpos
            new_x_length = self[xlabel].t - temp

            ## update parent of xlabel
            gpar = self[xlabel].par
            if gpar == -1:
                print("warning")
            children_of_gpar = list(self[gpar].child)
            children_of_gpar.remove(xlabel)
            self[gpar].update(child = [emptylabel, children_of_gpar[0]])

            ## create new sister branch
            self[emptylabel] = branch(par = self[xlabel].par, child = [xlabel, label], t = temp, totalpos = self[xlabel].totalpos)

            ## update x branch with new children
            self[xlabel].update(par = emptylabel, t = new_x_length, totalpos = xpos)


            if label == None or emptylabel == None:
                raise Exception("ERROR: either selected branch or recombining branch is set to None.")

            ## update b branch with new parent, and new times
            selected_branch_time = xpos - self[label].totalpos + self[label].t
            self[label].update(par = emptylabel, t = selected_branch_time, totalpos = xpos)

    def recomb_at_root(self, label, xpos, emptylabel):
        """
        if you're recombining past the root, particular set of updates for recombination
        """
        ## find children of root
        old_root = self.rootid
        children_of_root = self[old_root].child

        ## update old root with new time
        new_time_old_root = xpos - self[old_root].totalpos
        self[old_root].update(par=emptylabel, t = new_time_old_root, totalpos = xpos)

        ## update label
        new_time_selected_branch = xpos - self[label].totalpos + self[label].t
        self[label].update(par = emptylabel, t= new_time_selected_branch, totalpos = xpos)

        ## make emptylabel branch as new root
        self[emptylabel] = branch(child = [label, old_root], par = -1, totalpos = xpos)

        self.rootid = emptylabel

        self.tmrca = xpos


    def find_interval(self, recombrate, start_interval, length):
        """
        find interval along chr for tree
        """
        ## find time to recombination event as distance along chromosome
        time_to_recombination_event = np.random.exponential(scale = self.tlen * recombrate)

        ## find new end time as length to event plus start
        new_end = time_to_recombination_event + start_interval

        ## if end is longer than total length set end as total length
        if new_end > length:
            new_end = length

        ## set interval
        self.interval = (start_interval, new_end)

        ## return end to use as start of next interval
        return new_end

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
                    tempupid  = self[tempedge].child[i]
                    tempnodestr = innernewick(self,tempupid)
                    temptime = self[self[tempedge].child[i]].t
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


    def gettreestats(self):
        """
            print length and edge count summaries
        """
        info = "Length : %.3f  TMRCA: %.3f  number edges: %d"%(self.tlen,self.tmrca, self.n)
        return info


    def __str__(self):
        s = self.gettreestats() + "\n"
        s += self.newick()
        return s

##############################CHECK #############################################
    def check_totalpos(self):
        """
        check totalpos for each branch to make sure it matches the total length when loop through all children to leaf
        """
        ## loop through all branches
        for x in self.keys():
            ## add up the sum of each branch plus its line of descendants to tip (leaves)
            sum_len = self[x].t
            temp = x
            ## loop through descendants
            while True:
                children = self[temp].child
                if len(children) > 0:
                    temp = children[0]
                    sum_len += self[temp].t
                else:
                    break
            if math.isclose(sum_len, self[x].totalpos):
                continue
            else:
                raise Exception("error: Branch {} has totalpos of {} but it should be {}".format(x, self[x].totalpos, sum_len))



    def length_check(self):
        """
        check to see if each leaf to root length is the same
        """
        lengthlist = []
        def leaf_to_root(self, x, totalleng):
            totalleng = totalleng + self[x].t
            x = self[x].par
            if x == -1:
                lengthlist.append(round(totalleng, 6))
            else:
                leaf_to_root(self, x, totalleng)

        for x in self.keys():
            totalleng = 0
            if self[x].isexternal():
                leaf_to_root(self, x, totalleng)
        if lengthlist.count(lengthlist[0]) == len(lengthlist):
            return False, lengthlist
        else:
            return True, lengthlist

    def treecomp(self):
        """
        check tree consistency
            1. all parents have 2 children
            2. all children have 1 parent
        """
        ## loop through all branches
        for x in self.keys():
            ## find parent of each branch
            par = self[x].par
            ## if the parent is root, you're finished with tree
            if par == -1:
                ## return error = False, and None
                return False, None
            ## find children of parent
            children = self[par].child
            ## if children length is not 2 for each parent
            if not len(children) == 2:
                ## return error is true, and the parent with error
                return True, par
            ## if the branch isn't in the children
            if not x in children:
                ## return error is true, and the parent branch with error
                return True, par

    def check(self, branchpaths = False, isTree = False):
        """
        overall check to check two things:
            1. all root to tip lengths are equal
            2. there are no tree inconsistencies
        """
        ## root to leaf check
        if branchpaths == True:
            error, llist = self.length_check()
            if error == True:
                raise Exception("All leaf-to-root lengths were not equal.  The lengths were: {}".format(llist))

        ## tree consistency check
        if isTree == True:
            error, errorbranch = self.treecomp()
            if error == True:
                raise Exception("Tree had inconsistency.  Error occured at {} branch.".format(errorbranch))

    def root_check(self):
        """
        check to make sure root id totalpos matches the tmrca
        """
        if self[self.rootid].totalpos != self.tmrca:
            raise Exception("ERROR: root position does not match TMRCA time")
        if self[self.rootid].par != -1:
            raise Exception("ERROR: rootid wrong")

    def check_allbranches(self):
        for x in self.keys():
            if self[x].t < 0:
                return 0

    def check_lengths_chr(treelist):
        lenglist = []
        for x in treelist:
            lenglist.append(x.tlen)

        totalleng = 0
        for x in lenglist:
            totalleng += x

        avglen = totalleng / len(lenglist)