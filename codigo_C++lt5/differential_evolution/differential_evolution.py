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
popsize = 5                                               #Population size
# Mutation factor [0,2]
mutate = 0.7
# Recombination rate [0,1]
recombination = 0.5
# Max number of generations (maxiter)
maxiter = 2

t_exp = [np.array([0.0, 0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.91, 6.92, 10.94, 13.96, 20.91, 27.98]) #PATB07
        ,np.array([0.0, 0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90, 6.94, 10.96, 14.95, 20.96, 28.03]) #PATB09
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.51, 1.00, 1.50, 2.91, 4.92, 6.91, 10.87, 13.86, 20.87]) #PATB08
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01, 7.00, 9.99, 14.00, 20.98, 30.01]) #PATB16
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03]) #PATB17
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85, 13.85, 20.86, 27.96]) #PATB06
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 2.86, 6.91, 8.86, 9.88, 13.90, 21.87, 30.87]) #PATC05
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94, 5.99, 7.01, 9.91, 15.00, 22.00, 28.02]) #PATC06
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94, 7.94, 9.95, 14.94, 24.03, 30.94]) #PATC09
        ,np.array([0.0, 0.04, 0.09, 0.18, 0.33, 0.50, 1.00, 1.50, 2.96, 3.98, 7.00, 9.98, 13.97, 21.96]) #PATC10
        ]
patients = []
PATB07 = [4.1510, 4.2227, 4.0733, 3.0599, 2.2122, 1.7924, 1.6435, 1.7160, 1.6435, 1.2304, 1.0792, 1.0000, 1.0000, 1.0000, 1.0000]
PATB09 = [6.3541, 6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385, 4.1949, 3.0931, 2.1790, 1.5682, 1.5051]
PATB08 = [5.7831, 5.6546, 5.6486, 4.6203, 3.6876, 3.1526, 2.6355, 2.5302, 2.6294, 2.2788, 2.0719, 1.9031, 1.7924, 1.6335]
PATB16 = [6.2943, 6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120, 3.0302, 2.7528, 2.3838, 2.1818, 1.9243]
PATB17 = [6.1927, 6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432, 2.7474, 2.7016, 2.3541, 2.0453, 1.4914]
PATB06 = [6.2584, 6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378, 2.3404, 2.3345, 2.2355, 2.0492, 2.1173]
PATC05 = [6.208568, 6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339, 2.609594, 2.527630, 2.743510, 2.694605, 2.227887, 1.863323]
PATC06 = [6.808956, 6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317, 2.484300, 2.509203, 2.369216, 1.949390, 1.623249, 1.556303]
PATC09 = [6.296968, 6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341, 2.089905, 2.025306, 1.698970, 1.278754, 1.342423]
PATC10 = [5.547272, 5.583842, 5.455846, 4.754914, 3.447468, 2.749736, 2.311754, 2.158362, 1.944483, 1.707570, 1.322219, 1.322219, 1.000000, 1.000000]

patients.append(PATB07)
patients.append(PATB09)
patients.append(PATB08)
patients.append(PATB16)
patients.append(PATB17)
patients.append(PATB06)
patients.append(PATC05)
patients.append(PATC06)
patients.append(PATC09)
patients.append(PATC10)
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
            plt.title("PATB07")
        elif pat_cont == 2:
            plt.title("PATB09")
        elif pat_cont == 3:
            plt.title("PATB08")
        elif pat_cont == 4:
            plt.title("PATB16")
        elif pat_cont == 5:
            plt.title("PATB17")
        elif pat_cont == 6:
            plt.title("PATB06")
        elif pat_cont == 7:
            plt.title("PATC05")
        elif pat_cont == 8:
            plt.title("PATC06")
        elif pat_cont == 9:
            plt.title("PATC09")
        elif pat_cont == 10:
            plt.title("PATC10")

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
