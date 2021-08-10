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
popsize = 3                                               #Population size
# Mutation factor [0,2]
mutate = 0.7
# Recombination rate [0,1]
recombination = 0.5
# Max number of generations (maxiter)
maxiter = 1

t_exp = [
        np.array([0, 0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85, 13.85, 20.86, 27.96]) #PATB06
        ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50,3.00,6.01,7.01]) #PATB15
        ,np.array([0,0.04,0.09,0.17,0.33,0.50,1.00,1.50,3.00,4.03,6.99,9.99]) #PATB11
        ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50,3.00,3.95,6.98,9.98]) #PATB10
        ,np.array([0,0.04,0.08,0.17,0.33,0.51,1.00,1.50,2.87,3.83,8.84]) #PATC04
        ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50,2.94,5.93,7.93,9.94,13.97,20.96]) #PATC16
        ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50]) #PATC15
        ,np.array([0,0.04,0.09,0.17,0.34,0.50,1.00,1.50,2.91,5.96,8.96]) #PATC17
        ,np.array([0,0.04,0.08,0.17,0.34,0.50,1.00,1.50,2.97,5.99,7.97]) #PATC14
        ]
patients = []
PATB06 = [6.2584, 6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378, 2.3404, 2.3345, 2.2355, 2.0492, 2.1173]
PATB15 = [5.3108, 5.2514, 5.2002, 4.0835, 3.2017, 2.9759, 2.1038, 2.1959, 1.3222, 1.0000, 1.1139]
PATB11 = [6.0881,6.1478,6.0223,5.2421,3.8581,3.4216,2.8808,2.7168,2.2718,2.0170,1.6721,1.3222]
PATB10 = [6.3840,6.4587,6.4587,6.1718,5.4932,4.9461,3.5412,3.1535,2.4669,2.3945,1.9294,1.6990]
PATC04 = [6.7323,6.4158,6.5225,5.6844,4.4105,3.8132,3.1732,2.7987,2.5611,2.3962,1.9031]
PATC16 = [6.1903,6.0958,6.0654,5.4680,4.5812,4.2582,3.7706,3.5206,3.0637,2.6212,2.4728,2.1790,1.7782,1.1461]
PATC15 = [3.2524, 2.5145, 2.4116, 2.4014, 1.9243, 1.8633, 1.5052, 1.5441]
PATC17 = [6.0135,5.9922,5.9374,5.3949,3.7766,3.1000,2.6335,2.6304,2.0934,1.4771,1.1761]
PATC14 = [5.9922,5.8673,5.8764,5.6997,5.0688,4.5903,4.0052,3.3286,2.6702,1.8865,1.5185]

patients.append(PATB06)
patients.append(PATB15)
patients.append(PATB11)
patients.append(PATB10)
patients.append(PATC04)
patients.append(PATC16)
patients.append(PATC15)
patients.append(PATC17)
patients.append(PATC14)

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
        print(pat_cont, " Patient")
        saida.writelines(str(pat_cont) + " Patient\n\n")
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
        if pat_cont == 1:
            plt.title("PATB06")
        elif pat_cont == 2:
            plt.title("PATB15")
        elif pat_cont == 3:
            plt.title("PATB11")
        elif pat_cont == 4:
            plt.title("PATB10")
        elif pat_cont == 5:
            plt.title("PATC04")
        elif pat_cont == 6:
            plt.title("PATC16")
        elif pat_cont == 7:
            plt.title("PATC15")
        elif pat_cont == 8:
            plt.title("PATC17")
        elif pat_cont == 9:
            plt.title("PATC14")

        plt.xlabel("dias")
        plt.ylabel("Carga viral $log_{10}$")
        plt.legend()
        # Plot da solucao com os melhores parametros
        cost_func(sol_pat.x, pat, (10**pat[0]), pat_cont-1, t_exp)
        print('\nCusto do melhor conjunto de parametros: ' +
              str(sol_pat.fun) + '\n\n')
        print('\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
        plt.savefig("../figs/pat_"+str(pat_cont)+"NGem_" +
                    str(maxiter)+"NPop_"+str(popsize)+".png")

        pat_cont += 1
        saida.close()

    
# --- END ----------------------------------------------------------------------+
