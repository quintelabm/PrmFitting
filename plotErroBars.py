
import matplotlib.pyplot as plt
 

#Tempo das medicoes e da discretizacao do modelo
x = [1,2,3,4,5,6,7,8,9,10];
#resultado do modelo
yy = [10,20,30,40,50,60,70,80,90,100]
#Media do resultado experimental
y = [5,25,25,45,50,55,75,80,85,95]
#Desvio padrao dos dados
yerr = [15,15,15,15,15,15,15,15,15,15]
plt.ylim(0,120)
#plt.errorbar(x, y, yerr = yerr, marker='s', fmt='.k', uplims=True, lolims=True)
plt.errorbar(x, y, yerr = yerr, marker='s', fmt='.k')

plt.plot(x,yy)
