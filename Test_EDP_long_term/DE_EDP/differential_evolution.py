# Tutorial - https://nathanrooy.github.io/posts/2017-08-27/simple-differential-evolution-with-python/
#
#1) Initialize a random population of individuals throughout the search space.
#2) while iter <= max num of generations
#    3) cycle through each individual in the population   
#        3.A) perform mutation            
#        3.B) perform recombination ("crossover" in GA lingo)            
#        3.C) perform selection           
#    4) if stopping criterion has been met:
#            exit and return best individual            
#        else:
#            iter = iter + 1
#            go back to step #3

#------------------------------------------------------------------------------+
#
#   Nathan A. Rooy
#   A simple, bare bones, implementation of differential evolution with Python
#   August, 2017
#
#------------------------------------------------------------------------------+

#--- IMPORT DEPENDENCIES ------------------------------------------------------+

import random
import numpy as np
from depend.cost import viralmodelfit
from depend.bounds import ensure_bounds
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import sys
import os

#--- MAIN ---------------------------------------------------------------------+


#--- CONSTANTS ----------------------------------------------------------------+

cost_func = viralmodelfit                                  # Cost function
bounds = [(0.990,0.999),(0.900,0.960),(0.200,0.400),(0.1,0.3),(15,24),(5.0,9.9),(9.7,14)]  # Bounds [(x1_min, x1_max), (x2_min, x2_max),...]
popsize = 50                                               # Population size, must be >= 4
mutate = 0.5                                               # Mutation factor [0,2]
recombination = 0.7                                        # Recombination rate [0,1]
maxiter = 10                                               # Max number of generations (maxiter)

#Vetor com todos os pacientes
patients = [ ]
# experimental data
#PATB06 = [6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378]#, 2.3404, 2.3345, 2.2355, 2.0492, 2.1173]
#PATB07 = [4.2227, 4.0733, 3.0599, 2.2122, 1.7924, 1.6435, 1.7160, 1.6435, 1.2304]#, 1.0792, 1.0000, 1.0000, 1.0000, 1.0000]
#PATB08 = [5.6546, 5.6486, 4.6203, 3.6876, 3.1526, 2.6355, 2.5302, 2.6294, 2.2788]#, 2.0719, 1.9031, 1.7924, 1.6335, 1.5682]
#PATB09 = [6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385]#, 4.1949, 3.0931, 2.1790, 1.5682, 1.5051]
#PATB16 = [6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529]#, 3.0120, 3.0302, 2.7528, 2.3838, 2.1818, 1.9243]
PATB17 = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128]#, 2.8432, 2.7474, 2.7016, 2.3541, 2.0453, 1.4914]


#PATC05 = [6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339]#, 2.609594, 2.527630, 2.743510, 2.694605, 2.227887, 1.863323]
#PATC06 = [6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317]#, 2.484300, 2.509203, 2.369216, 1.949390, 1.623249, 1.556303]
#PATC09 = [6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341]#, 2.089905, 2.025306, 1.698970, 1.278754, 1.342423]
#PATC10 = [5.583842, 5.455846, 4.754914, 3.447468, 2.749736, 2.311754, 2.158362, 1.944483, 1.707570]#, 1.322219, 1.322219, 1.000000, 1.000000, 1.000000]

'''
PAT8  = [5.64, 5.31, 4.23, 3.36, 3.14, 2.86, 2.75, 2.50, 2.32, 1.56]
PAT42 = [5.65, 5.00, 3.98, 3.84, 2.94, 2.82, 2.87, 2.53, 2.31, 2.61]
PAT68 = [7.15, 7.02, 6.19, 5.50, 4.96, 4.29, 4.11, 3.75, 3.68, 3.35]
PAT69 = [6.14, 5.87, 4.73, 4.17, 3.55, 3.14, 2.87, 2.60, 2.55, 2.58]
PAT83 = [5.45, 5.38, 4.73, 4.00, 3.39, 2.89, 2.68, 2.72, 2.97, 1.93]
'''
'''
PAT = sys.argv[1]

if(sys.argv[1] == "PAT8"):
    patients = PAT8
if(sys.argv[1] == "PAT42"):
    patients = PAT42
if(sys.argv[1] == "PAT68"):
    patients = PAT68
if(sys.argv[1] == "PAT69"):
    patients = PAT69
if(sys.argv[1] == "PAT83"):
    patients = PAT83
'''
'''
patients.append(PAT8)
patients.append(PAT42)
patients.append(PAT68)
patients.append(PAT69)
patients.append(PAT83)
'''

#patients.append(PATB06)
#patients.append(PATB07)
#patients.append(PATB08)
#patients.append(PATB09)
#patients.append(PATB16)
patients.append(PATB17)
#patients.append(PATC05)
#patients.append(PATC06)
#patients.append(PATC09)
#patients.append(PATC10)

#t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2 ]

t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99]#, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03] #PATB17
    
#--- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    
    saida = open("D:\Faculdade\IC\GitHub\PrmFitting\Test_EDO_EDP\DE_EDP\outputs/relatorio.txt", "a")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: '+ str(popsize)+ '\nNumber of generations: '+ str(maxiter)+ '\n')

    best_solves = []
    pat_cont = 1
    for pat in patients:
        print(pat_cont, " Patient")
        saida.writelines(str(pat_cont) + " Patient\n\n")
        
        sol_pat = differential_evolution(cost_func, bounds, args=(pat,10**pat[0]), maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination)
        #sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
        print(sol_pat.x, "\n", sol_pat.fun)
        saida.writelines('\nCusto do melhor conjunto de parametros: '+ str(sol_pat.fun) +'\n\n')
        saida.writelines('\nO melhor conjunto de parametros: '+str(sol_pat.x) +'\n\n')
        plt.clf()
        #Plot experimental        
        plt.plot(t_exp, pat, 'ro')
        #Plot da solucao com os melhores parametros
        cost_func(sol_pat.x, pat, (10**pat[0]))
        plt.savefig("pat_"+str(pat_cont)+"NGem_"+str(maxiter)+"NPop_"+str(popsize)+".png")
        
        pat_cont += 1

    
    saida.close()
    '''
    cwd = os.getcwd()
    print(cwd)
    saida = open(cwd+"/HCV_DE_scipy/docs/" + str(PAT) + ".txt", "a+")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: '+ str(popsize)+ '\nNumber of generations: '+ str(maxiter)+ '\n')

    best_solves = []
    
    sol_pat = differential_evolution(cost_func, bounds, args=(patients,10**patients[0]), maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination)
    #sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
    print(sol_pat.x, "\n", sol_pat.fun)
    saida.writelines('\nCusto do melhor conjunto de parametros: '+ str(sol_pat.fun) +'\n\n')
    saida.writelines('\nO melhor conjunto de parametros: '+str(sol_pat.x) +'\n\n')
    plt.clf()
    #Plot experimental        
    plt.plot(t_exp, patients, 'ro')
    #Plot da solucao com os melhores parametros
    cost_func(sol_pat.x, patients, (10**patients[0]))
    plt.savefig(PAT + ".png")
    
    saida.close()
    '''
#--- END ----------------------------------------------------------------------+

