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

#--- MAIN ---------------------------------------------------------------------+


#--- CONSTANTS ----------------------------------------------------------------+

cost_func = viralmodelfit                                  # Cost function
bounds = [(0.09,0.99),(0.09,0.99),(0.09,0.99),(0.01,1.2),(25,35),(0.9,2.0),(6,12)]  # Bounds [(x1_min, x1_max), (x2_min, x2_max),...]
popsize = 10                                               # Population size, must be >= 4
mutate = 0.5                                               # Mutation factor [0,2]
recombination = 0.7                                        # Recombination rate [0,1]
maxiter = 1                                              # Max number of generations (maxiter)

#Vetor com todos os pacientes
patients = [ ]

PATB09 = [6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385, 4.1949, 3.0931, 2.1790, 1.5682, 1.5051]

patients.append(PATB09)

t_exp = np.array([0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90, 6.94, 10.96, 14.95, 20.96, 28.03]) #PATB09
    
#--- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    saida = open("docs/relatorio.txt", "a+")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: '+ str(popsize)+ '\nNumber of generations: '+ str(maxiter)+ '\n')

    best_solves = []
    pat_cont = 1
    for pat in patients:
        print(pat_cont, " Patient")
        saida.writelines(str(pat_cont) + " Patient\n\n")
        
        sol_pat = differential_evolution(cost_func, bounds, args=(pat,10**pat[0]), maxiter=5, popsize=10, mutation=0.5, recombination=0.7)
        #sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
        print(sol_pat.x, "\n", sol_pat.fun)
        saida.writelines('\nCusto do melhor conjunto de parametros: '+ str(sol_pat.fun) +'\n\n')
        saida.writelines('\nO melhor conjunto de parametros: '+str(sol_pat.x) +'\n\n')
        plt.clf()
        #Plot experimental        
        plt.plot(t_exp, pat, 'ro')
        #Plot da solucao com a media dos parametros
        cost_func(sol_pat.x, pat, (10**pat[0]))
        plt.savefig("docs/figura"+str(pat_cont)+".png")
        
        pat_cont += 1

    
    saida.close()
#--- END ----------------------------------------------------------------------+

