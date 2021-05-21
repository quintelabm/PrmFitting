import uncertainpy as un
import chaospy as cp
import numpy as np
import time
import os
import matplotlib.pyplot as plt
cwd = os.getcwd()
tempo = np.empty(0)
viral_load = np.empty(0)
with open(str(cwd)+"/saida.txt", 'r') as f:
    lista = [line.split(',')  for line in f]
    for linha in lista:
        tempo = np.append(tempo, float(linha[0]))
        viral_load = np.append(viral_load, float(linha[1]))

viral_load = np.log10(viral_load)

tempo_exp = np.empty(0)
viral_load_exp = np.empty(0)
num_pat = 'nome do paciente'
with open(str(cwd)+"/parametros.txt", 'r') as f:
    vet = f.readline()
    vet = vet.split(',')
    num_pat = vet[0]

with open(str(cwd)+"/../data/data2day_pats.csv", 'r') as f:
    f.readline()
    lista = [line.split(',')  for line in f]
    if num_pat=='PAT8':
        for linha in lista:
            tempo_exp = np.append(tempo_exp, float(linha[0]))
            viral_load_exp = np.append(viral_load_exp, float(linha[1]))
    elif num_pat=='PAT42':
        for linha in lista:
            tempo_exp = np.append(tempo_exp, float(linha[0]))
            viral_load_exp = np.append(viral_load_exp, float(linha[2]))
    elif num_pat=='PAT68':
        for linha in lista:
            tempo_exp = np.append(tempo_exp, float(linha[0]))
            viral_load_exp = np.append(viral_load_exp, float(linha[3]))
    elif num_pat=='PAT69':
        for linha in lista:
            tempo_exp = np.append(tempo_exp, float(linha[0]))
            viral_load_exp = np.append(viral_load_exp, float(linha[4]))
    elif num_pat=='PAT83':
        for linha in lista:
            tempo_exp = np.append(tempo_exp, float(linha[0]))
            viral_load_exp = np.append(viral_load_exp, float(linha[5]))

plt.title(str(num_pat))
plt.plot(tempo_exp, viral_load_exp, 'or', label="dados experimentais")
plt.plot(tempo, viral_load, '-b', label="resultado do modelo em C++")
plt.ylabel("Viral load $log_{10}$")
plt.xlabel("Days")
plt.legend()
plt.savefig(str(cwd)+"/figs/"+str(num_pat)+".png", dpi=300)
plt.show()