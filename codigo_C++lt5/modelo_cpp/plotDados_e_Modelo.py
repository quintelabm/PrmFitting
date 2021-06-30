import uncertainpy as un
import chaospy as cp
import numpy as np
import time
import os
import matplotlib.pyplot as plt
cwd = os.getcwd()
os.system("make")
tempo = np.empty(0)
viral_load = np.empty(0)
with open(str(cwd)+"/../differential_evolution/depend/saida.txt", 'r') as f:
    lista = [line.split(',')  for line in f]
    for linha in lista:
        tempo = np.append(tempo, float(linha[0]))
        viral_load = np.append(viral_load, float(linha[1]))

viral_load = np.log10(viral_load)

#pegar pontos experimentais
t_exp = [np.array([0.0, 0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.91]) #PATB07
        ,np.array([0.0, 0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90]) #PATB09
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.51, 1.00, 1.50, 2.91, 4.92]) #PATB08
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01]) #PATB16
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03]) #PATB17
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89]) #PATB06
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 2.86, 6.91]) #PATC05
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94, 5.99]) #PATC06
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94, 7.94]) #PATC09
        ,np.array([0.0, 0.04, 0.09, 0.18, 0.33, 0.50, 1.00, 1.50, 2.96, 3.98, 7.00]) #PATC10
        ]
patients = []
PATB07 = [4.1510, 4.2227, 4.0733, 3.0599, 2.2122, 1.7924, 1.6435, 1.7160, 1.6435, 1.2304]
PATB09 = [6.3541, 6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385]
PATB08 = [5.7831, 5.6546, 5.6486, 4.6203, 3.6876, 3.1526, 2.6355, 2.5302, 2.6294, 2.2788]
PATB16 = [6.2943, 6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120]
PATB17 = [6.1927, 6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432]
PATB06 = [6.2584, 6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378]
PATC05 = [6.208568, 6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339, 2.609594]
PATC06 = [6.808956, 6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317, 2.484300]
PATC09 = [6.296968, 6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341, 2.089905]
PATC10 = [5.547272, 5.583842, 5.455846, 4.754914, 3.447468, 2.749736, 2.311754, 2.158362, 1.944483, 1.707570, 1.322219]

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

####SELECIONE AQUI O PACIENTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pat_select = 0#1#2#...

plt.plot(t_exp[pat_select], patients[pat_select], 'or', label="dados experimentais")
plt.plot(tempo, viral_load, '-b', label="resultado do modelo em C++")
plt.ylabel("Viral load $log_{10}$")
plt.xlabel("Days")
plt.legend()
plt.savefig(str(cwd)+"/figs/"+str(pat_select)+".png", dpi=300)
plt.show()