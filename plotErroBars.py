
import matplotlib.pyplot as plt
 

#Tempo das medicoes e da discretizacao do modelo
x = [0.04, 0.09, 0.18, 0.33, 0.50, 1.00, 1.50, 2.96, 3.98, 7.00, 9.98, 13.97, 21.96, 30.02]

#Media do resultado experimental
y = [6.47991433, 6.42897983, 5.857277, 4.63766183, 4.08108333, 3.38275467,
     3.18124267, 2.88121,    2.66053917, 2.54078967, 2.487822,   2.21939417,
     1.90103167, 1.7158415 ]

#Desvio padrao dos dados
yerr = [0.182129952, 0.182630388, 0.337284052,	0.390495914, 0.307469046, 0.18942229,
        0.200947972, 0.239112485, 0.222385492, 0.324597051, 0.294209462, 0.350716883,
        0.372253077, 0.297196377]

#plt.ylim(0,120)
#plt.errorbar(x, y, yerr = yerr, marker='s', fmt='.k', uplims=True, lolims=True)
#plt.errorbar(x, y, yerr = yerr, marker='s', fmt='.k')
plt.errorbar(x, y, yerr = yerr, fmt='.k')

plt.show()