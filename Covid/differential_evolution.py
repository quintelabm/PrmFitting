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
bounds = [(10.0, 30.0), (0.4, 1.80), (4.0, 25.0)]          # Bounds [(x1_min, x1_max), (x2_min, x2_max),...]
popsize = 10                                               # Population size, must be >= 4
mutate = 0.5                                               # Mutation factor [0,2]
recombination = 0.7                                        # Recombination rate [0,1]
maxiter = 1                                                # Max number of generations (maxiter)

#Vetor com todos os pacientes
patients = []

PATCOV2 = [4.923722693632887, 7.575302146580392, 6.801792298849835, 0.5711651663789257]

patients.append(PATCOV2)

t_exp = [4.0, 5.0, 6.0, 8.0]
    
#--- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    saida = open("/home/matheus/Documents/PrmFitting/HCV_DE_scipy/docs/relatorio.txt", "a")
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
        cost_func(sol_pat.x, pat, 10**pat[0])
        plt.savefig("pat_"+str(pat_cont)+"NGem_"+str(maxiter)+"NPop_"+str(popsize)+".png")
        
        pat_cont += 1

    
    saida.close()
#--- END ----------------------------------------------------------------------+
