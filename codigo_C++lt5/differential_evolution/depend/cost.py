import numpy as np 
from scipy.spatial import distance
import matplotlib.pyplot as plt
import os
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0, pat_cont, t_exp):
    epsilon_r = poi[0]
    epsilon_alpha = poi[1]
    epsilon_s = poi[2]
    kappa_t = poi[3]
    kappa_c = poi[4]

    param_fixados = [[50.79045657,  0.28635711,  1.44320029,  1.62737769,  5.54053499,  1.35017632,
  1.08556533, 19.75186589],
    [47.43399187,  1.31253404,  1.34881359,  1.56992606 , 7.17867972 , 1.45425331,
  1.18953957 ,17.68123129],
    [21.89653548 , 0.49846047 , 0.62618439 , 1.40645206, 10.26735378 , 1.69266101,
  1.01005749, 10.28375646],
    [23.77001288,  2.86177985,  1.72165976,  0.83277745 , 3.06722576 , 1.59824232,
  1.68726557 ,16.12291485],
    [57.2214857,   1.07688969 , 1.76686689 , 1.67198125 ,13.26593168 , 1.36026391,
  1.35933788, 15.67626232],
    [39.75636044,  1.26967926 , 1.74260023 , 0.68097692, 13.90182105  ,1.33475382,
  1.41948718 ,16.92846897],
    [24.73436476,  1.90594692,  1.24535908 , 0.95938638 , 5.46734905  ,1.31107016,
  1.55127437 ,22.57429968],
    [24.06566535,  2.25882618 , 1.18414759,  0.90286285 , 7.69214418,  1.30902344,
  1.79985456 ,11.97889563],
    [31.31166921 , 3.24107381 , 1.72588821,  0.61336909 , 1.06576933 , 1.31907584,
  1.09225015 ,13.90359195],
    [29.2164192 ,  0.43273399  ,1.48494351  ,1.28492415, 13.71920283 , 1.99506034,
  1.54772895, 14.45624281]]

    param_pat_fixados = param_fixados[pat_cont]
    
    with open('parametros_DE.txt', 'w') as filep:
        filep.write(str(V0)+","+str(epsilon_r)+","+str(epsilon_alpha)+","+str(epsilon_s)+
        ","+str(kappa_t)+","+str(kappa_c)+","+str(param_pat_fixados[0])+","+str(param_pat_fixados[1])+
        ","+str(param_pat_fixados[2])+","+str(param_pat_fixados[3])+","+str(param_pat_fixados[4])+
        ","+str(param_pat_fixados[5])+","+str(param_pat_fixados[6])+","+str(param_pat_fixados[7]))
    
    os.system("make run")#Executa o modelo C++
    tempoPt = np.empty(0)
    V = np.empty(0)
    V_log = np.empty(0)
    with open("saida.txt", 'r') as f:
        lista = [line.split(',')  for line in f]
        for linha in lista:
            tempoPt = np.append(tempoPt, float(linha[0]))
            V = np.append(V, float(linha[1]))
    try:
      # Passa para a base log o resultado
      V_log = np.log10(V)
    except:
      dst = 10000
      print(dst)
      return dst
    V_pts = []
    for t in t_exp[pat_cont]:
      V_pts.append(V_log[int(t*100+1)])
    plt.plot(tempoPt, V_log, '-g')
    
    plt.plot(t_exp[pat_cont], exp)
    
    dst = distance.euclidean(V_pts, exp)/len(V_pts)
  
    print(dst)
    return dst