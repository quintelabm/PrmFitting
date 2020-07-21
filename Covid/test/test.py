import numpy as np 
from scipy.spatial import distance
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import sys

#variaveis de inicializacao
tempoFim = 10
deltaT   = 0.01 #passo no tempo
tempoNpts = int(tempoFim/deltaT) + 1 
tempoPt = np.linspace(4, tempoFim, tempoNpts)
tempoCont = 0


# define uma funcao feval que recebe o nome da equacao a ser avaliada como 
# string e retorna a funcao a ser avaliada
#def feval(funcName, *args):
 #   return eval(funcName)(*args)

#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0):
    
    def func(x, y):
        ## parametros do sistema de 3 equacoes descrito abaixo
        
        #PATCOV
        s       = 1.3*10**5
        beta    = 5*10**-8
        d       = 0.01
        epsilon = 0.0
        #Parametros ajustados
        c       = poi[0]
        delta   = poi[1]
        p       = poi[2]
            
        ## inicializa com zeros
        dy = np.zeros(3)
        
        ## equacoes: y[0] = T, y[1] = I, y[2] = V
        dy[0] = s - beta*y[2]*y[0] - d*y[0]
        dy[1] = beta*y[2]*y[0] - delta*y[1]
        dy[2] = (1 - epsilon)*p*y[1] - c*y[2]

        return dy

    V = np.zeros(tempoNpts)
    V[0] = V0*400
        
    x = tempoPt[0]

    T0  = 2.9168*10**6
    I0 = 8.7186*10**5
    y = [T0, I0, V0]

    for i in range(tempoNpts):
        
        k1 = func(x, y)

        yp2 = y + k1*(deltaT/2)

        k2 = func(x+deltaT/2, yp2)

        yp3 = y + k2*(deltaT/2)

        k3 = func(x+deltaT/2, yp3)

        yp4 = y + k3*deltaT

        k4 = func(x+deltaT, yp4)

        for j in range(3):
            y[j] = y[j] + (deltaT/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

        x = x + deltaT
   
    
    
    t_exp = [4, 5, 6, 8]


    exp[0] = 10**exp[0]
    exp[1] = 10**exp[1]
    exp[2] = 10**exp[2]
    exp[3] = 10**exp[3]


    plt.xlim([0, 10])
    plt.ylim([0, exp[1]])
    plt.plot(t_exp, exp, 'or')
    plt.plot(tempoPt, V, '-g')
    # fazer interpolacao usando os pontos experimentais
    ius = InterpolatedUnivariateSpline(t_exp, exp)
    
    # aplicar a funcao interpolada nos pontos do metodo numerico do modelo
    yi = ius(tempoPt)
    
    # calcular a distancia euclidiana entre os resultados experimentais interpolados
    # e o resultado do modelo
    
    dst = distance.euclidean(V, yi)
    print(exp)
    return dst

PATCOV2 = []

if(sys.argv[1] == "PATCOV2"):

    PATCOV2 = [4.923722693632887, 7.575302146580392, 6.801792298849835, 0.5711651663789257]
print(PATCOV2)

dist = viralmodelfit([20, 0.5, 60],PATCOV2,10**PATCOV2[0])
print(dist)

plt.show()