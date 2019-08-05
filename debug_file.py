# -*- coding: utf-8 -*-
"""

@author: tuk32868
"""

import classes

N = 100
ss = 3
mutrate = 1e-7
seed = 14
length = 100
recombrate = 0.01
output_name = "trees_output.txt"
f = open(output_name, "w")
f.write("Interval\tTree\tLength\n")

treelist = []

t = classes.tree(N, ss, seed)

while True:
    finish = t.continuous_coal(mutrate)
    if finish == 0:
        break

t.check(True, True)
t.check_totalpos()

## smc initialize
start_interval = 0
while True:
    end_int = t.find_interval(recombrate, start_interval, length)
    t.check(True, True)
    t.check_totalpos()
    f.write("{}\t{}\t{}\n".format(t.interval, t.newick(), t.tlen))
    treelist.append(t.tlen)
    if end_int >= length:
        break
    else:
        start_interval = end_int
        t.select_branch()

f.close()

totalleng = 0
for x in treelist:
    totalleng += x

avglen = totalleng/len(treelist)

print(avglen)
