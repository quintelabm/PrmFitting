import numpy as np 
from scipy.spatial import distance
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt


#variaveis de inicializacao
tempoFim = 2
deltaT   = 0.01 #passo no tempo
tempoNpts = int(tempoFim/deltaT) + 1 
tempoPt = np.linspace(0, tempoFim, tempoNpts)
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
        epsilon = 0.996
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
    V[0] = V0
        
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
        

    #Se algum valor de V for negativo, defina um custo alto para o conjunto de parametros
    for pos in V:
        if not (pos>=0):
            return 1000000

    # Passa para a base log o resultado
    V_log = np.log10(V)
    
    t_exp = [4.0, 5.0, 6.0, 8.0]
    
    plt.plot(tempoPt, V_log, '-g')
    # fazer interpolacao usando os pontos experimentais
    ius = InterpolatedUnivariateSpline(t_exp, exp)
    
    # aplicar a funcao interpolada nos pontos do metodo numerico do modelo
    yi = ius(tempoPt)
    
    # calcular a distancia euclidiana entre os resultados experimentais interpolados
    # e o resultado do modelo
    
    dst = distance.euclidean(V_log, yi)
    
    return dst


