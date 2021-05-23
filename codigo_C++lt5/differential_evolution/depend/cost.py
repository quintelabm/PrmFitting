import numpy as np 
from scipy.spatial import distance
import matplotlib.pyplot as plt
import os

#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0, pat_cont, t_exp):
    
    alpha = poi[0]
    r     = poi[1]
    delta = poi[2]
    mu_c = poi[3]
    rho   = poi[4]
    epsilon_r = poi[5]
    epsilon_alpha = poi[6]
    
    with open('parametros_DE.txt', 'w') as filep:
        filep.write(str(V0)+","+str(alpha)+","+str(r)+","+str(delta)+","+str(mu_c)+","+str(rho)+","+str(epsilon_r)+","+str(epsilon_alpha))
    
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
        
        dst = distance.euclidean(V_pts, exp[pat_cont])/len(V_pts)
    except:
        dst = 10000
        pass
    return dst