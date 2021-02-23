import numpy as np 
from scipy.spatial import distance
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import os
#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0):
    epsilon_s = poi[0]
    epsilon_alpha = poi[1]
    epsilon_r = poi[2]
    delta = poi[3]
    alpha = poi[4]
    r     = poi[5]
    rho   = poi[6]
    print(delta)
    print(rho)
    with open('parametros.txt', 'w') as filep:
        filep.write(str(epsilon_s)+","+str(epsilon_alpha)+","+str(epsilon_r)+","+str(delta)+","+str(alpha)+","+str(r)+","+str(rho))
    os.system("make run")#Executa o modelo C++

    tempoPt = np.empty(0)
    V = np.empty(0)
    with open("saida.txt", 'r') as f:
        lista = [line.split(' ')  for line in f]
        for linha in lista:
            tempoPt = np.append(tempoPt, float(linha[0]))
            V = np.append(V, float(linha[1]))

    # Passa para a base log o resultado
    V_log = np.log10(V)
    
    t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2 ]
    
    plt.plot(tempoPt, V_log, '-g')
    # fazer interpolacao usando os pontos experimentais
    ius = InterpolatedUnivariateSpline(t_exp, exp)
    
    # aplicar a funcao interpolada nos pontos do metodo numerico do modelo
    yi = ius(tempoPt)
    
    # calcular a distancia euclidiana entre os resultados experimentais interpolados
    # e o resultado do modelo
    
    dst = distance.euclidean(V_log, yi)
    
    return dst

