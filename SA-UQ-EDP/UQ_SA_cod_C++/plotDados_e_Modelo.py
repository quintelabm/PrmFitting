import uncertainpy as un
import chaospy as cp
import numpy as np
import time
import os
import matplotlib.pyplot as plt

tempo = np.empty(0)
viral_load = np.empty(0)
with open("saida.txt", 'r') as f:
    lista = [line.split(',')  for line in f]
    for linha in lista:
        tempo = np.append(tempo, float(linha[0]))
        viral_load = np.append(viral_load, float(linha[1]))

viral_load = np.log10(viral_load)

tempo_exp = np.empty(0)
viral_load_exp = np.empty(0)
with open("data2day_pats.csv", 'r') as f:
    f.readline()
    lista = [line.split(',')  for line in f]
    for linha in lista:
        tempo_exp = np.append(tempo_exp, float(linha[0]))
        viral_load_exp = np.append(viral_load_exp, float(linha[1]))

plt.plot(tempo_exp, viral_load_exp, 'or', label="dados experimentais do pat8")
plt.plot(tempo, viral_load, '-b', label="resultado do modelo em C++")
plt.ylabel("Viral load $log_{10}$")
plt.xlabel("Days")
plt.legend()
plt.show()