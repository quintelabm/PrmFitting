# --- IMPORT DEPENDENCIES ------------------------------------------------------+

import random
import numpy as np
from depend.cost import viralmodelfit
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, NonlinearConstraint
import os
import time
# --- MAIN ---------------------------------------------------------------------+


# --- CONSTANTS ----------------------------------------------------------------+

cost_func = viralmodelfit                                  #Cost function
bounds = [(0.1,0.99),(0.1,0.99),(0.1,0.99),(20,60),(0.1,10),(0.01,2),(0.1,2),(1,15), (1,2), (1,2), (10,25)] 
array_param = "epsilon_r, epsilon_alpha, epsilon_s, alpha, r, delta, mu_c, rho, theta, sigma, c" 
popsize = 50                                               #Population size
# Mutation factor [0,2]
mutate = 0.7
# Recombination rate [0,1]
recombination = 0.5
# Max number of generations (maxiter)
maxiter = 100

t_exp = [
        np.array([0.0, 0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90, 6.94, 10.96, 14.95, 20.96, 28.03]) #PATB09
        ]
patients = []

PATB09 = [6.3541, 6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385, 4.1949, 3.0931, 2.1790, 1.5682, 1.5051]
patients.append(PATB09)

# --- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    
    cwd = os.getcwd()
    saida = open(cwd+"/relatorio_DE.txt", "a")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: ' + str(popsize) +
                     '\nNumber of generations: ' + str(maxiter) + '\n')
    saida.writelines(array_param + '\n')

    saida.close()
    pat_cont = 1
    for pat in patients:
        saida = open(cwd+"/relatorio_DE.txt", "a")
        os.system("make clean")
        os.system("make")
        print("2 Patient")
        saida.writelines("2 Patient\n\n")
        tic = time.perf_counter()
        sol_pat = differential_evolution(cost_func, bounds, args=(
            pat, 10**pat[0], pat_cont-1, t_exp), maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination)#, tol=0.1)
        toc = time.perf_counter()
        # sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
        saida.writelines(f"Tempo de execução: {toc - tic:0.4f} segundos")
        saida.writelines(
            '\nCusto do melhor conjunto de parametros: ' + str(sol_pat.fun) + '\n\n')
        saida.writelines(
            '\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
        plt.clf()
        # Plot experimental
        plt.plot(t_exp[pat_cont-1], pat, 'ro')
        
        plt.title("PATB09")

        plt.xlabel("dias")
        plt.ylabel("Carga viral $log_{10}$")
        plt.legend()
        # Plot da solucao com os melhores parametros
        cost_func(sol_pat.x, pat, (10**pat[0]), pat_cont-1, t_exp)
        print('\nCusto do melhor conjunto de parametros: ' +
              str(sol_pat.fun) + '\n\n')
        print('\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
        plt.savefig("../figs/pat_2NGem_" +
                    str(maxiter)+"NPop_"+str(popsize)+".png")

        pat_cont += 1
        saida.close()

    
# --- END ----------------------------------------------------------------------+
