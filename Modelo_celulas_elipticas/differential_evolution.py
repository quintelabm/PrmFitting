#--- IMPORT DEPENDENCIES ------------------------------------------------------+

import random
import numpy as np
from solver_TIEV import custo, solver
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import sys
import os

#--- CONSTANTS ----------------------------------------------------------------+

cost_func = custo                                  # Cost function
bounds = [(0.01, 1.8),(0.01,0.9999),(1,15),(10,25),(1,10)]  # Bounds [(x1_min, x1_max), (x2_min, x2_max),...]
popsize = 50                                               # Population size, must be >= 4
mutate = 0.5                                               # Mutation factor [0,2]
recombination = 0.7                                        # Recombination rate [0,1]
maxiter = 5                                                # Max number of generations (maxiter)

#Vetor com todos os pacientes
patients = []
t_exp_total = []
# --- for patient B06
PATB06 = [6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378,
            2.3404, 2.3345]#, 2.2355, 2.0492, 2.1173]

t_exp1 = [0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85]#, 13.85, 20.86, 27.96]

# --------------------------------------------------------------------------------------------------------

# --- for patient B16
PATB16 = [6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120,
            3.0302, 2.7528]#, 2.3838, 2.1818, 1.9243]

t_exp2 = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01, 7.00, 9.99]#, 14.00, 20.98, 30.01]

# --------------------------------------------------------------------------------------------------------

# --- for patient B17
PATB17 = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432,
            2.7474, 2.7016]#, 2.3541, 2.0453, 1.4914]

t_exp3 = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98]#, 14.04, 21.03, 30.03]

# --------------------------------------------------------------------------------------------------------

# --- for patient C05
PATC05 = [6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339,
            2.609594, 2.527630, 2.743510]#, 2.694605, 2.227887, 1.863323]

t_exp4 = [0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 2.86, 6.91, 8.86, 9.88]#, 13.90, 21.87, 30.87]

# --------------------------------------------------------------------------------------------------------

# --- for patient C06
PATC06 = [6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317, 2.484300,
            2.509203, 2.369216]#, 1.949390, 1.623249, 1.556303]

t_exp5 = [0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94, 5.99, 7.01, 9.91]#, 15.00, 22.00, 28.02]

# --------------------------------------------------------------------------------------------------------

# --- for patient C09
PATC09 = [6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341,
            2.089905, 2.025306]#, 1.698970, 1.278754, 1.342423]

t_exp6 = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94, 7.94, 9.95]#, 14.94, 24.03, 30.94]

patients.append(PATB06)
patients.append(PATB16)
patients.append(PATB17)
patients.append(PATC05)
patients.append(PATC06)
patients.append(PATC09)

t_exp_total.append(t_exp1)
t_exp_total.append(t_exp2)
t_exp_total.append(t_exp3)
t_exp_total.append(t_exp4)
t_exp_total.append(t_exp5)
t_exp_total.append(t_exp6)

#--- Post processing ----------------------------------------------------------+

def plot_graficos(t_range, solve, t_exp, data_exp,pwd,i):
    
    v = solve[:,3]
    log_v = np.log10(v)
    plt.figure()
    plt.plot(t_range, log_v, '-g', label = "Model")
    plt.plot(t_exp, data_exp, 'or', label = "Data")
    plt.title(f'Carga viral no tempo Paciente {str(i+1)}')
    plt.ylabel("Carga viral ($log_{10}$ UI/mL)")
    plt.xlabel("Tempo (dias)")
    plt.legend()
    #plt.show()
    plt.savefig(pwd + '/Modelo_celulas_elipticas/Results/'+ str(i +1) + 'patiente' + ".png",dpi=300)
    

#--- MAIN ---------------------------------------------------------------------+

if __name__ == "__main__":
    
    days = 10
    pwd = os.getcwd()
    saida = open(f'{pwd}/Modelo_celulas_elipticas/Results/relatorio.txt', "a+")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: '+ str(popsize)+ '\nNumber of generations: '+ str(maxiter)+ '\n')
    for i in range(0, len(patients)):
        
        saida.writelines(str(i +1) + ' patiente')
        sol_pat = differential_evolution(cost_func, bounds, args=(t_exp_total[i], patients[i], days), maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination)
        saida.writelines('\nCusto do melhor conjunto de parametros: '+ str(sol_pat.fun) +'\n\n')
        saida.writelines('\nO melhor conjunto de parametros: '+str(sol_pat.x) +'\n\n')
        #Plot da solucao com os melhores parametros
        #sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
        t_range, sol =solver(sol_pat.x[0], sol_pat.x[1], sol_pat.x[2], sol_pat.x[3], sol_pat.x[4], 10**patients[i][0], days)
        plot_graficos(t_range,sol,t_exp1,patients[i],pwd,i)


    #saida.close()


#--- END ----------------------------------------------------------------------+
