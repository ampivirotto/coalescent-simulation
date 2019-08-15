# -*- coding: utf-8 -*-
"""

@author: tuk32868
"""

import classes
import copy

N = 1000
ss = 10
mutrate = 1e-3
seed = 14
length = 200
recombrate = 1e-3
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
inter_interval = 0
while True:
    end_int = t.find_interval(recombrate, inter_interval, length, start_interval)
    t.mutate(mutrate)
    t.check(True, True)
    t.check_totalpos()
    copied_tree = copy.deepcopy(t)
    if end_int >= length:
        f.write("{}\t{}\t{}\n".format(copied_tree.interval, copied_tree.newick(), copied_tree.tlen))
        treelist.append(copied_tree)
        break
    else:
        same = t.select_branch()
        if same == True:
            inter_interval = end_int
        else:
            inter_interval = end_int
            start_interval = end_int
            f.write("{}\t{}\t{}\n".format(copied_tree.interval, copied_tree.newick(), copied_tree.tlen))
            treelist.append(copied_tree)


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
                mutationdict[y] = [x]
            else:
                newlist = rec_mute(tree, x)
                mutationdict[y] = newlist
                blist = []
    return mutationdict

m = open("mutations.txt", "w")
header = "Position\t"
for br in range(0,ss):
    header += str(br) + "\t"
m.write(header + "\n")
for x in treelist:
    mutdict = return_mutes(x)
    for key in sorted(mutdict.keys()):
        gt = ""
        derived = mutdict[key]
        for branch in x.keys():
            if x[branch].isexternal():
                if branch in derived:
                    gt += "1\t"
                else:
                    gt += "0\t"
        m.write("{}\t{}\n".format(key, gt))

m.close()

avglen, avgtmrca = t.check_lengths_chr(treelist)

actlen = 0
for i in range(1, ss-1):
    actlen += 1/i
actlen = actlen * 4 *N
print(actlen, avglen)

acttmrca = 4 * N * (1 - (1/ss))
print(acttmrca, avgtmrca)


