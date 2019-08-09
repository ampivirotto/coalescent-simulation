# -*- coding: utf-8 -*-
"""

@author: tuk32868
"""

import classes
import copy

N = 1000
ss = 30
mutrate = 0.0001
seed = 14
length = 20
recombrate = 0.0001
output_name = "trees_output.txt"
f = open(output_name, "w")
f.write("Interval\tTree\tLength\n")

treelist = []

t = classes.tree(N, ss, seed)

while True:
    finish = t.continuous_coal()
    if finish == 0:
        break

t.check(True, True)
t.check_totalpos()

## smc initialize
start_interval = 0
while True:
    end_int = t.find_interval(recombrate, start_interval, length)
    t.mutate(mutrate)
    t.check(True, True)
    t.check_totalpos()
    f.write("{}\t{}\t{}\n".format(t.interval, t.newick(), t.tlen))
    copied_tree = copy.deepcopy(t)
    treelist.append(copied_tree)
    if end_int >= length:
        break
    else:
        start_interval = end_int
        t.select_branch()

f.close()

def return_mutes(tree):
    mutationdict = {}
    blist = []
    def rec_mute(tree, x):
        children = tree[x].child
        for y in children:
            blist.append(y)
            if not tree[x].isexternal():
                rec_mute(tree, y)
        return blist
    for x in tree.keys():
        mutation = tree[x].mut
        for y in mutation:
            if tree[x].isexternal():
                mutationdict[y] = x
            else:
                newlist = rec_mute(tree, x)
                mutationdict[y] = newlist
                blist = []
    return mutationdict

m = open("mutations.txt", "w")
for x in treelist:
    mutdict = return_mutes(x)
    for key in sorted(mutdict.keys()):
        m.write("{}\t{}\n".format(key, mutdict[key]))

m.close()


