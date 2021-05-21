import numpy as np 
from scipy.spatial import distance
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import os

#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0, pat_cont):
    delta = poi[0]
    mu_t = poi[1]
    r     = poi[2]
    mu_c = poi[3]
    epsilon_alpha = poi[4]
    epsilon_r = poi[5]
    sigma = poi[6]
    theta = poi[7]
    rho   = poi[8]
    alpha = poi[9]
    
    with open('parametros_DE.txt', 'w') as filep:
        filep.write(str(V0)+","+str(delta)+","+str(mu_t)+","+str(r)+","+str(mu_c)+","+str(epsilon_alpha)+","+str(epsilon_r)+","+str(sigma)+","+str(theta)+","+str(rho)+","+str(alpha))
    
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
        
        t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2 ]
        
        V_pts = [V_log[0], V_log[int(0.083*100)-1], V_log[int(0.167*100)-1] , V_log[int(0.25*100)-1] \
            , V_log[int(0.333*100)-1] , V_log[int(0.5*100)-1] , V_log[int(0.667*100)-1] , \
                V_log[int(1*100)-1] , V_log[int(1.5*100)-1] , V_log[int(2*100)-1]]

        plt.plot(tempoPt, V_log, '-g')
        
        dst = distance.euclidean(V_pts, exp)/len(V_pts)
    except:
        dst = 10000
        pass
    
    return dst

