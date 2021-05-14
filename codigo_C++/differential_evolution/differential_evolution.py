#--- IMPORT DEPENDENCIES ------------------------------------------------------+

import random
import numpy as np
from depend.cost import viralmodelfit
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, NonlinearConstraint
import sys
import os
import time
#--- MAIN ---------------------------------------------------------------------+


#--- CONSTANTS ----------------------------------------------------------------+

cost_func = viralmodelfit                                  # Cost function
# Bounds   delta,    mu_t,       r,    mu_c,  epsilon_alpha, epsilon_r,   sigma,      theta,     rho,         alpha
# bounds = [(0.1,0.9),(0.8,0.9),(1,5.2),(1.7,3.5),(0.8,0.999),(0.1,0.5),(1.29,1.31),(1.19,1.21),(8.179,8.181),(29.99,30.01)]  
bounds = [(0.1,0.9),(0.4,0.9),(1,5.8),(1.1,4.5),(0.3,0.999),(0.01,0.8),(1.29,1.31),(1.19,1.21),(8.179,8.181),(29.99,30.01)]  
popsize = 50# 50                                               # Population size, must be >= 4
mutate = 0.7                                               # Mutation factor [0,2]
recombination = 0.5                                        # Recombination rate [0,1]
maxiter = 200#50,100,200,400,600,800,1000                                                # Max number of generations (maxiter)

#Vetor com todos os pacientes
patients = [ ]

PAT8  = [5.64, 5.31, 4.23, 3.36, 3.14, 2.86, 2.75, 2.50, 2.32, 1.56]
PAT42 = [5.65, 5.00, 3.98, 3.84, 2.94, 2.82, 2.87, 2.53, 2.31, 2.61]
PAT68 = [7.15, 7.02, 6.19, 5.50, 4.96, 4.29, 4.11, 3.75, 3.68, 3.35]
PAT69 = [6.14, 5.87, 4.73, 4.17, 3.55, 3.14, 2.87, 2.60, 2.55, 2.58]
PAT83 = [5.45, 5.38, 4.73, 4.00, 3.39, 2.89, 2.68, 2.72, 2.97, 1.93]

# patients.append(PAT8)
# patients.append(PAT42)
# patients.append(PAT68)
# patients.append(PAT69)
patients.append(PAT83)

t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2 ]
    
#--- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":

    def restricoes(param):
        
        mu_t = param[1]
        r     = param[2]
        mu_c = param[3]
        sigma = param[6]
        theta = param[7]
        rho   = param[8]
        alpha = param[9]
        restricao1 = sigma + rho + mu_c - sigma*theta/(theta + rho + mu_t)
        restricao2 = alpha*r - (sigma + rho + mu_c - sigma*theta/(theta + rho + mu_t))*mu_c
        #print(str(restricao1)+"" )
        return [restricao1, restricao2]
    func_restricoes = NonlinearConstraint(restricoes, [0,0], [np.inf,np.inf])
    cwd = os.getcwd()
    saida = open(cwd+"/relatorio_DE.txt", "a")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: '+ str(popsize)+ '\nNumber of generations: '+ str(maxiter)+ '\n')

    best_solves = []
    pat_cont = 1
    for pat in patients:
        
        os.system("make clean")
        os.system("make")
        print(pat_cont, " Patient")
        saida.writelines(str(pat_cont) + " Patient\n\n")
        tic = time.perf_counter()
        sol_pat = differential_evolution(cost_func, bounds, args=(pat,10**pat[0], pat_cont), maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination)#, constraints=(func_restricoes))
        toc = time.perf_counter()
        #sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
        print(sol_pat.x, "\n", sol_pat.fun)
        saida.writelines(f"Tempo de execução: {toc - tic:0.4f} segundos")
        saida.writelines('\nCusto do melhor conjunto de parametros: '+ str(sol_pat.fun) +'\n\n')
        saida.writelines('\nO melhor conjunto de parametros: '+str(sol_pat.x) +'\n\n')
        plt.clf()
        #Plot experimental        
        plt.plot(t_exp, pat, 'ro')
        if pat_cont==1:
            #plt.title("PAT8")
            plt.title("PAT83")
        elif pat_cont==2:
            plt.title("PAT42")
        elif pat_cont==3:
            plt.title("PAT68")
        elif pat_cont==4:
            plt.title("PAT69")
        else:
            plt.title("PAT83")

        plt.xlabel("dias")
        plt.ylabel("Carga viral $log_{10}$")
        #Plot da solucao com os melhores parametros
        cost_func(sol_pat.x, pat, (10**pat[0]), pat_cont)
        print('\nCusto do melhor conjunto de parametros: '+ str(sol_pat.fun) +'\n\n')
        print('\nO melhor conjunto de parametros: '+str(sol_pat.x) +'\n\n')
        plt.savefig("pat_"+str(pat_cont)+"NGem_"+str(maxiter)+"NPop_"+str(popsize)+".png")
        
        pat_cont += 1

    
    saida.close()
#--- END ----------------------------------------------------------------------+

