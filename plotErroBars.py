
import matplotlib.pyplot as plt
 

#Tempo das medicoes e da discretizacao do modelo
x = [0,0.04,0.08,0.17,0.33,0.5,1,1.5,3,10,15,24,30];
#resultado do modelo
yy = [10,20,30,40,50,60,70,80,90,100]
#Media do resultado experimental
y = [6.34331,6.47991,6.42897,5.857277,4.63766,4.08108 , 3.38275,3.18124,2.88121,2.487822,2.21939,1.90103,1.7158415]
#Desvio padrao dos dados
yerr = [0.23213,0.18213,0.18262,0.33728 ,0.39049,0.30748, 0.18942,0.20093,0.23913,0.29420,0.35071 ,0.37224 ,0.29720]
#plt.ylim(0,120)
#plt.errorbar(x, y, yerr = yerr, marker='s', fmt='.k', uplims=True, lolims=True)
#plt.errorbar(x, y, yerr = yerr, marker='s', fmt='.k')
plt.errorbar(x, y, yerr = yerr, fmt='.k')
#plt.plot(x,yy)
