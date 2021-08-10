import uncertainpy as un
import chaospy as cp
import numpy as np
import time
import os
import matplotlib.pyplot as plt
####SELECIONE AQUI O PACIENTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pat_select = 1#1#2#...até 9
#pegar pontos experimentais
# t_exp = [np.array([0.0, 0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.91, 6.92, 10.94, 13.96, 20.91, 27.98, 41.85, 56.00, 66.88]) #PATB07
#         ,np.array([0.0, 0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90, 6.94, 10.96, 14.95, 20.96, 28.03, 42.00, 55.97, 66.92]) #PATB09
#         ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.51, 1.00, 1.50, 2.91, 4.92, 6.91, 10.87, 13.86, 20.87, 25.92, 39.96, 53.93, 67.93]) #PATB08
#         ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01, 7.00, 9.99, 14.00, 20.98, 30.01, 44.01, 55.01, 72.02]) #PATB16
#         ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03, 44.03, 55.00, 72.03]) #PATB17
#         ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85, 13.85, 20.86, 27.96, 41.97, 53.92, 67.89]) #PATB06
#         ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 2.86, 6.91, 8.86, 9.88, 13.90, 21.87, 30.87, 41.93, 55.88]) #PATC05
#         ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94, 5.99, 7.01, 9.91, 15.00, 22.00, 28.02, 44.01, 57.01, 70.00]) #PATC06
#         ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94, 7.94, 9.95, 14.94, 24.03, 30.94, 44.91, 58.92, 72.11]) #PATC09
#         ,np.array([0.0, 0.04, 0.09, 0.18, 0.33, 0.50, 1.00, 1.50, 2.96, 3.98, 7.00, 9.98, 13.97, 21.96, 30.02, 41.89, 55.97, 70.01]) #PATC10
#         ]
# patients = []
# PATB07 = [4.1510, 4.2227, 4.0733, 3.0599, 2.2122, 1.7924, 1.6435, 1.7160, 1.6435, 1.2304, 1.0792, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000]
# PATB09 = [6.3541, 6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385, 4.1949, 3.0931, 2.1790, 1.5682, 1.5051, 1.0000, 1.0000, 1.0000]
# PATB08 = [5.7831, 5.6546, 5.6486, 4.6203, 3.6876, 3.1526, 2.6355, 2.5302, 2.6294, 2.2788, 2.0719, 1.9031, 1.7924, 1.6335, 1.5682, 1, 1, 1]
# PATB16 = [6.2943, 6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120, 3.0302, 2.7528, 2.3838, 2.1818, 1.9243, 1.7160, 1, 1]
# PATB17 = [6.1927, 6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432, 2.7474, 2.7016, 2.3541, 2.0453, 1.4914, 1, 1, 1]
# PATB06 = [6.2584, 6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378, 2.3404, 2.3345, 2.2355, 2.0492, 2.1173, 1.3424, 1.7559, 1.1461]
# PATC05 = [6.208568, 6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339, 2.609594, 2.527630, 2.743510, 2.694605, 2.227887, 1.863323, 1.785330, 1.431364]
# PATC06 = [6.808956, 6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317, 2.484300, 2.509203, 2.369216, 1.949390, 1.623249, 1.556303, 1.000000, 1.000000, 1.000000]
# PATC09 = [6.296968, 6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341, 2.089905, 2.025306, 1.698970, 1.278754, 1.342423, 1.000000, 1.000000, 1.000000]
# PATC10 = [5.547272, 5.583842, 5.455846, 4.754914, 3.447468, 2.749736, 2.311754, 2.158362, 1.944483, 1.707570, 1.322219, 1.322219, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000]

# patients.append(PATB07)
# patients.append(PATB09)
# patients.append(PATB08)
# patients.append(PATB16)
# patients.append(PATB17)
# patients.append(PATB06)
# patients.append(PATC05)
# patients.append(PATC06)
# patients.append(PATC09)
# patients.append(PATC10)

t_exp = [
        np.array([0, 0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85, 13.85, 20.86, 27.96]) #PATB06
        ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50,3.00,6.01,7.01,10.02,14.01,20.97,28.01,42.02,56.02,70.02]) #PATB15
        ,np.array([0,0.04,0.09,0.17,0.33,0.50,1.00,1.50,3.00,4.03,6.99,9.99,16.07,21.01,28.02,42.95,56.05,70.02]) #PATB11
        ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50,3.00,3.95,6.98,9.98,15.97,20.99,28.02,42.99,55.97,70.01]) #PATB10
        ,np.array([0,0.04,0.08,0.17,0.33,0.51,1.00,1.50,2.87,3.83,8.84,15.87,22.85,29.83,44.90,58.84,72.85]) #PATC04
        ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50,2.94,5.93,7.93,9.94,13.97,20.96,28.05,42.06,55.92,71.08]) #PATC16
        ,np.array([0,0.04,0.08,0.17,0.33,0.50,1.00,1.50,2.92,7.07,8.04,9.07,13.06,23.08,26.92,42.00,57.01,69.93]) #PATC15
        ,np.array([0,0.04,0.09,0.17,0.34,0.50,1.00,1.50,2.91,5.96,8.96,12.96,14.93,20.96,28.94,43.94,55.90,70.92]) #PATC17
        ,np.array([0,0.04,0.08,0.17,0.34,0.50,1.00,1.50,2.97,5.99,7.97,10.06,14.06,21.02,28.04,44.07,56.03,70.01]) #PATC14
        ]
patients = []
PATB06 = [6.2584, 6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378, 2.3404, 2.3345, 2.2355, 2.0492, 2.1173]
PATB15 = [5.3108, 5.2514, 5.2002, 4.0835, 3.2017, 2.9759, 2.1038, 2.1959, 1.3222, 1.0000, 1.1139, 1, 1, 1, 1, 1, 1, 1]
PATB11 = [6.0881,6.1478,6.0223,5.2421,3.8581,3.4216,2.8808,2.7168,2.2718,2.0170,1.6721,1.3222,1,1,1,1,1,1]
PATB10 = [6.3840,6.4587,6.4587,6.1718,5.4932,4.9461,3.5412,3.1535,2.4669,2.3945,1.9294,1.6990,1,1,1,1,1,1]
PATC04 = [6.7323,6.4158,6.5225,5.6844,4.4105,3.8132,3.1732,2.7987,2.5611,2.3962,1.9031, 0, 0, 0, 0, 0, 0]
PATC16 = [6.1903,6.0958,6.0654,5.4680,4.5812,4.2582,3.7706,3.5206,3.0637,2.6212,2.4728,2.1790,1.7782,1.1461, 0, 0, 0, 0]
PATC15 = [3.2524, 2.5145, 2.4116, 2.4014, 1.9243, 1.8633, 1.5052, 1.5441, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
PATC17 = [6.0135,5.9922,5.9374,5.3949,3.7766,3.1000,2.6335,2.6304,2.0934,1.4771,1.1761, 0, 0, 0, 0, 0, 0, 0]
PATC14 = [5.9922,5.8673,5.8764,5.6997,5.0688,4.5903,4.0052,3.3286,2.6702,1.8865,1.5185, 0, 0, 0, 0, 0, 0, 0]

patients.append(PATB06)
patients.append(PATB15)
patients.append(PATB11)
patients.append(PATB10)
patients.append(PATC04)
patients.append(PATC16)
patients.append(PATC15)
patients.append(PATC17)
patients.append(PATC14)

param = [ 
        [10**4.1510, 0.11195621,  0.83223552 , 0.20835353, 37.99280026  ,1.00997408,  0.32654547,
  0.15854512 ,14.32745201 , 1.94130815  ,1.99091771, 12.81053005],
        [10**6.3541,  0.1873109,   0.38749771 , 0.46654267, 31.95011327,  0.42494102,  0.80153783,
  0.81106208 , 6.82583957 , 1.17112151,  1.59968312, 12.62840383],
        [10**5.7831, 0.13231509,  0.85276596,  0.64494253, 32.87861619,  5.82721283,  0.1521823 
,  0.5340823  ,12.84616902  ,1.03839509  ,1.59181091 ,12.54744969],
        [10**6.2943, 0.22180002 , 0.19188497 , 0.98370764 ,20.37489472 , 3.50278207 , 0.86479336 ,
  1.05523751 ,13.86097434 , 1.75409632 , 1.92212859 ,11.24023732],
        [10**6.1927, 0.52213334 , 0.71553654 , 0.97173176 ,53.96696789 , 3.49846152 , 0.99463084 ,
  1.26739641 ,11.82016039 , 1.2023713  , 1.76838477 ,10.00244143],
        [10**6.2584, 0.29191191 , 0.83220358 , 0.9864726  ,44.4227656  , 2.75136521 , 1.47322094 ,
  0.4375029  , 7.53112131 , 1.53665179 , 1.38476878 ,10.08267574],
        [10**6.208568, 0.55387014 , 0.81193914 , 0.19842489 ,28.29531535 , 6.82039395 , 0.22693666 
,  0.11897197 ,11.33354543  ,1.85385136  ,1.43795092 ,10.00245583], 
        [10**6.808956, 0.22731839 , 0.96644911 , 0.75696487 ,30.85093474 , 4.62074022 , 0.22526165, 
  0.31712111,  9.47433421,  1.24957259,  1.73691962, 10.00837386],
        [10**6.296968, 0.2864288 ,  0.83966836,  0.1217267 , 43.969472  ,  1.44357305,  0.69492007
,  0.2202072  , 5.16750968  ,1.69379056  ,1.32072105 ,10.00338683],
        [10**5.547272, 0.91412008,  0.92584465,  0.95821606, 25.1334769 ,  2.12234684,  0.22813056,
  0.26382177, 10.85085191 , 1.23614772 , 1.89389546 ,11.83417401]
]

V0, epsilon_r, epsilon_alpha, epsilon_s, alpha, r, delta, mu_c, rho, theta, sigma, c = param[pat_select]

with open('parametros_DE.txt', 'w') as filep:
        filep.write(str(V0)+","+str(epsilon_r)+","+str(epsilon_alpha)+","+str(epsilon_s)+
        ","+str(alpha)+","+str(r)+","+str(delta)+
        ","+str(mu_c)+","+str(rho)+","+str(theta)+
        ","+str(sigma)+","+str(c))
cwd = os.getcwd()
os.system("make clean")
os.system("make")
os.system("make run")
tempo = np.empty(0)
viral_load = np.empty(0)
with open(str(cwd)+"/saida.txt", 'r') as f:
    lista = [line.split(',')  for line in f]
    for linha in lista:
        tempo = np.append(tempo, float(linha[0]))
        viral_load = np.append(viral_load, float(linha[1]))

viral_load = np.log10(viral_load)

pat_name = []
pat_name.append("PATB07")
pat_name.append("PATB09")
pat_name.append("PATB08")
pat_name.append("PATB16")
pat_name.append("PATB17")
pat_name.append("PATB06")
pat_name.append("PATC05")
pat_name.append("PATC06")
pat_name.append("PATC09")
pat_name.append("PATC10")
plt.plot(t_exp[pat_select], patients[pat_select], 'or', label="experimental data")
plt.plot(tempo, viral_load, '-b', label="results model C++")
plt.title(pat_name[pat_select])
plt.ylabel("Viral load $log_{10}$")
plt.xlabel("Days")
plt.legend()
plt.savefig(str(cwd)+"/figuras_teste_dos_parametros/"+"patient"+str(pat_select)+".png", dpi=300)
plt.show()