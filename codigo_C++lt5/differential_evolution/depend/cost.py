import numpy as np 
from scipy.spatial import distance
import matplotlib.pyplot as plt
import os

#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0, pat_cont, t_exp):
    
    epsilon_r = poi[0]
    epsilon_alpha = poi[1]

    param_fixados = [[48.93502213,  0.14376176,  1.9840078,   1.90872316,  9.15149142],
    [30.5660933,   1.34154019,  1.04650851,  0.72719804, 14.97288525],
    [39.49677098,  1.00408143,  1.44299392,  0.85239269,  2.77250136],
    [39.3568868,   2.0914808,   1.55277302,  1.57137463,  3.72419731],
    [36.40337668,  2.17761915,  1.51627981,  1.60288674, 14.58361958],
    [24.24964255,  5.98279124,  1.94674456,  1.82706024,  4.03590666],
    [31.34421558,  1.75082034,  0.90815144,  1.39364569,  3.89987979],
    [23.45791309,  5.62054536,  1.76144527,  0.81695743, 13.20626882],
    [41.33574756,  2.72810448,  1.73425789,  1.10120866,  2.11799729],
    [33.09870264,  0.37335046,  0.36376974,  1.23295938, 12.17081886]]

    param_pat_fixados = param_fixados[pat_cont]
    
    with open('parametros_DE.txt', 'w') as filep:
        filep.write(str(V0)+","+str(epsilon_r)+","+str(epsilon_alpha)+","+str(param_pat_fixados[0])+","+str(param_pat_fixados[1])+","+str(param_pat_fixados[2])+","+str(param_pat_fixados[3])+","+str(param_pat_fixados[4]))
    
    os.system("make run")#Executa o modelo C++
    tempoPt = np.empty(0)
    V = np.empty(0)
    with open("saida.txt", 'r') as f:
        lista = [line.split(',')  for line in f]
        for linha in lista:
            tempoPt = np.append(tempoPt, float(linha[0]))
            V = np.append(V, float(linha[1]))
    try:
        # Passa para a base log o resultado
        V_log = np.log10(V)
        
        V_pts = []
        for t in t_exp[pat_cont]:
            V_pts.append(V_log[int(t*100+1)])
        plt.plot(tempoPt, V_log, '-g')
        
        dst = distance.euclidean(V_pts, exp)/len(V_pts)
    except:
        dst = 10000
        pass
    return dst