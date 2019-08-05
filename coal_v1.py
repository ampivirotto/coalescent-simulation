#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      ampiv
#
# Created:     26/05/2019
# Copyright:   (c) ampiv 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import sys
import classes 


def main(popsize, samplesize, mutrate, seed):
    t = classes.tree(popsize, samplesize, seed)
    while True:
        finish = t.generation(mutrate)
        if finish == 0:
            #print(t)
            return t.tmrca, t.tlen

if __name__ == '__main__':
    ## input
    ## N = population size, ss = sample size, mutrate = mutation rate
    #N = sys.argv[1]
    N = 1000
    #ss = sys.argv[2]
    ss = 10
    #mutrate = sys.argv[3]
    mutrate = 0.0001
    
    genlist = []
    lenglist = []
    for x in range(1000, 2000):
        seed = x
        gen, leng = main(N, ss, mutrate, seed)
        lenglist.append(leng)
        genlist.append(gen)
    
    avggen = sum(genlist)/len(genlist)
    avglen = sum(lenglist)/len(lenglist)
    print(avglen)
    print(avggen)