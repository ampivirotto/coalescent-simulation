#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      tuk32868
#
# Created:     31/07/2019
# Copyright:   (c) tuk32868 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import math
import classes

class check(classes.tree):


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

    def check_lengths_chr(self, treelist):
        lenglist = []
        for x in treelist:
            lenglist.append(x.tlen)

        totalleng = 0
        for x in lenglist:
            totalleng += x

        avglen = totalleng / len(lenglist)
        return avglen


