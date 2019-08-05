#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      ampiv
#
# Created:     15/06/2019
# Copyright:   (c) ampiv 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import sys
import classes


def main(popsize, samplesize, mutrate, seed, length, recombrate):
    
    ## initialize list to contain trees 
    treelist = [] 
    
    ## initialize tree 
    t = classes.tree(popsize, samplesize, seed)
    
    ## initialize continuous coalescent until single root
    while True:
        finish = t.continuous_coal(mutrate)
        if finish == 0:
            break
    
    ## smc initialize 
    start_interval = 0
    while True:
        end_int = t.find_interval(recombrate, start_interval, length)
        treelist.append(t)
        if end_int >= length:
            break
        else:
            start_interval = end_int + 1
            t.select_branch()          
    return treelist   
            
    

if __name__ == '__main__':
    ## input
    ## N = population size, ss = sample size, mutrate = mutation rate
    #N = sys.argv[1]
    N = 1000
    #ss = sys.argv[2]
    ss = 10
    #mutrate = sys.argv[3]
    mutrate = 1e-7
    #seed = sys.argv[4]
    seed = 89
    #length = sys.argv[5]
    length = 1e4
    #recombrate = sys.argv[6]
    recombrate = 1e-7
    
    tree_list = main(N, ss, mutrate, seed, length, recombrate)
    print(tree_list)

