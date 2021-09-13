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
popsize = 30                                               #Population size
# Mutation factor [0,2]
mutate = 0.7
# Recombination rate [0,1]
recombination = 0.5
# Max number of generations (maxiter)
maxiter = 30

t_exp = [np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.51, 1.00, 1.50, 2.91, 4.92, 6.91, 10.87, 13.86, 20.87, 25.92, 39.96]) #PATB08
        ,np.array([0.0, 0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90, 6.94, 10.96, 14.95, 20.96, 28.03, 42.00]) #PATB09
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01, 7.00, 9.99, 14.00, 20.98, 30.01, 44.01]) #PATB16
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85, 13.85, 20.86, 27.96, 41.97, 53.92]) #PATB06
        ,np.array([0.0, 0.04, 0.09, 0.17, 0.33, 0.50, 1.00, 1.50, 3.00, 4.03, 6.99, 9.99, 16.07, 21.01]) #PATB11
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 3.00, 3.95, 6.98, 9.98, 15.97, 20.99]) #PATB10
        ]
patients = []
PATB08 = [5.7831, 5.6546, 5.6486, 4.6203, 3.6876, 3.1526, 2.6355, 2.5302, 2.6294, 2.2788, 2.0719, 1.9031, 1.7924, 1.6335, 1.5682, 1.0000]
PATB09 = [6.3541, 6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385, 4.1949, 3.0931, 2.1790, 1.5682, 1.5051, 1.0000]
PATB16 = [6.2943, 6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120, 3.0302, 2.7528, 2.3838, 2.1818, 1.9243, 1.7160]
PATB06 = [6.2584, 6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378, 2.3404, 2.3345, 2.2355, 2.0492, 2.1173, 1.3424, 1.7559]
PATB11 = [6.0881, 6.1478, 6.0223, 5.2421, 3.8581, 3.4216, 2.8808, 2.7168, 2.2718, 2.0170, 1.6721, 1.3222, 1, 1]
PATB10 = [6.3840, 6.4587, 6.4587, 6.1718, 5.4932, 4.9461, 3.5412, 3.1535, 2.4669, 2.3945, 1.9294, 1.6990, 1, 1]

#*********** PATIENTS C***************

# t_exp = [np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94, 5.99, 7.01, 9.91, 15.00, 22.00, 28.02, 44.01]) #PATC06
#         ,np.array([0,0.04,0.08,0.17,0.33,0.51,1.00,1.50,2.87,3.83,8.84,15.87,22.85]) #PATC04
#         ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50,2.92,7.07]) #PATC15
#         ,np.array([0,0.04,0.09,0.17,0.34,0.50,1.00,1.50,2.91,5.96,8.96,12.96,14.93,20.96]) #PATC17
#         ,np.array([0,0.04,0.08,0.17,0.34,0.50,1.00,1.50,2.97,5.99,7.97,10.06,14.06,21.02]) #PATC14
# ]

# PATC06 = [6.808956, 6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317, 2.484300, 2.509203, 2.369216, 1.949390, 1.623249, 1.556303, 1.000000]
# PATC04 = [6.7323,6.4158,6.5225,5.6844,4.4105,3.8132,3.1732,2.7987,2.5611,2.3962,1.9031, 1, 1]
# PATC15 = [3.2524, 2.5145, 2.4116, 2.4014, 1.9243, 1.8633, 1.5052, 1.5441, 1, 1]
# PATC17 = [6.0135,5.9922,5.9374,5.3949,3.7766,3.1000,2.6335,2.6304,2.0934,1.4771,1.1761, 1, 1, 1]
# PATC14 = [5.9922,5.8673,5.8764,5.6997,5.0688,4.5903,4.0052,3.3286,2.6702,1.8865,1.5185, 1, 1, 1]

patients.append(PATB08)
patients.append(PATB09)
patients.append(PATB16)
patients.append(PATB06)
patients.append(PATB11)
patients.append(PATB10)

pat_name = []
pat_name.append("PATB08")
pat_name.append("PATB09")
pat_name.append("PATB16")
pat_name.append("PATB06")
pat_name.append("PATB11")
pat_name.append("PATB10")

# pat_name.append("PATC06")
# pat_name.append("PATC04")
# pat_name.append("PATC15")
# pat_name.append("PATC17")
# pat_name.append("PATC14")
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
        plt.title(str(pat_name[pat_cont-1]))

        plt.xlabel("dias")
        plt.ylabel("Carga viral $log_{10}$")
        plt.legend()
        # Plot da solucao com os melhores parametros
        cost_func(sol_pat.x, pat, (10**pat[0]), pat_cont-1, t_exp)
        print('\nCusto do melhor conjunto de parametros: ' +
              str(sol_pat.fun) + '\n\n')
        print('\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
        plt.savefig("../figs/pat_"+str(pat_name[pat_cont-1])+"NGem_" +
                    str(maxiter)+"NPop_"+str(popsize)+".png")

        pat_cont += 1
        saida.close()

    
# --- END ----------------------------------------------------------------------+
