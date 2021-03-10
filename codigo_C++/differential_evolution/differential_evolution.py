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
bounds = [(0.1,0.8),(0.8,0.99),(1,6),(1,4),(0.9,0.999),(0.01,0.8)]  # Bounds [(x1_min, x1_max), (x2_min, x2_max),...]
popsize = 10                                               # Population size, must be >= 4
mutate = 0.5                                               # Mutation factor [0,2]
recombination = 0.7                                        # Recombination rate [0,1]
maxiter = 5                                                # Max number of generations (maxiter)

#Vetor com todos os pacientes
patients = [ ]

PAT8  = [5.64, 5.31, 4.23, 3.36, 3.14, 2.86, 2.75, 2.50, 2.32, 1.56]
PAT42 = [5.65, 5.00, 3.98, 3.84, 2.94, 2.82, 2.87, 2.53, 2.31, 2.61]
PAT68 = [7.15, 7.02, 6.19, 5.50, 4.96, 4.29, 4.11, 3.75, 3.68, 3.35]
PAT69 = [6.14, 5.87, 4.73, 4.17, 3.55, 3.14, 2.87, 2.60, 2.55, 2.58]
PAT83 = [5.45, 5.38, 4.73, 4.00, 3.39, 2.89, 2.68, 2.72, 2.97, 1.93]

patients.append(PAT8)
patients.append(PAT42)
patients.append(PAT68)
patients.append(PAT69)
patients.append(PAT83)

t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2 ]
    
#--- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    os.system("make")
    cwd = os.getcwd()
    print(cwd)
    saida = open(cwd+"/relatorio_DE.txt", "a")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: '+ str(popsize)+ '\nNumber of generations: '+ str(maxiter)+ '\n')

    best_solves = []
    pat_cont = 1
    for pat in patients:
        print(pat_cont, " Patient")
        saida.writelines(str(pat_cont) + " Patient\n\n")
        sol_pat = differential_evolution(cost_func, bounds, args=(pat,10**pat[0], pat_cont), maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination)
        #sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
        print(sol_pat.x, "\n", sol_pat.fun)
        saida.writelines('\nCusto do melhor conjunto de parametros: '+ str(sol_pat.fun) +'\n\n')
        saida.writelines('\nO melhor conjunto de parametros: '+str(sol_pat.x) +'\n\n')
        plt.clf()
        #Plot experimental        
        plt.plot(t_exp, pat, 'ro')
        if pat_cont==1:
            plt.title("PAT8")
        elif pat_cont==2:
            plt.title("PAT42")
        elif pat_cont==3:
            plt.plot("PAT68")
        elif pat_cont==4:
            plt.plot("PAT69")
        else:
            plt.plot("PAT83")

        plt.xlabel("dias")
        plt.ylabel("Carga viral $log_{10}$")
        #Plot da solucao com os melhores parametros
        cost_func(sol_pat.x, pat, (10**pat[0]), pat_cont)
        plt.savefig("pat_"+str(pat_cont)+"NGem_"+str(maxiter)+"NPop_"+str(popsize)+".png")
        
        pat_cont += 1

    
    saida.close()
#--- END ----------------------------------------------------------------------+

