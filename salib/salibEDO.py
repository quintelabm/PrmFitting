from SALib.sample import saltelli

from SALib.analyze import sobol

import matplotlib.pyplot as plt

import seaborn as sns

import numpy as np

sns.set()
#************* MODELO INICIO ************

# passo
h = 0.01
dias = 2
num_pontos = int(dias/h)+1

# Dias simulados
x_range = np.array([0.0, dias])
x_pt = np.linspace(0, dias, num_pontos)

# define uma funcao feval que recebe o nome da equacao a ser avaliada como 
# string e retorna a funcao a ser avaliada
def feval(funcName, *args):
    return eval(funcName)(*args)


# define a resolucao numerica por runge-kutta de quarta ordem
def RK4thOrder(p):
    delta, epsilon, c, p = p
    '''delta = p[0]
    epsilon = p[1]
    c = p[3]
    p = p[2]'''
    # define a funcao que contem as equacoes do sistema a ser avaliado
    def dinamicaIntracelular(x, y):
        ## parametros do sistema de 3 equacoes descrito abaixo
    
        s       = 1.3*10**5
        beta    = 5*10**-8
        d       = 0.01
    
        ## inicializa com zeros
        dy = np.zeros(3)
        
        ## equacoes: y[0] = T, y[1] = I, y[2] = V
        dy[0] = s - beta*y[2]*y[0] - d*y[0]
        dy[1] = beta*y[2]*y[0] - delta*y[1]
        dy[2] = (1 - epsilon)*p*y[1] - c*y[2]
    
        return dy
        
    # condicoes iniciais
    T0  = 2.9168*10**6
    I0 = 8.7186*10**5
    
    #PAT83:
    V0  = 4.9139*10**5
    
    yinit = np.array([T0,I0,V0], dtype='f')
    m = len(yinit)
    n = int((x_range[-1] - x_range[0])/h)
    
    x = x_range[0]
    y = yinit
    
    # Containers for solutions
    xsol = np.empty(0)
    xsol = np.append(xsol, x)

    ysol = np.empty(0)
    ysol = np.append(ysol, y)

    for i in range(n):
        k1 = dinamicaIntracelular(x, y)
        
        yp2 = y + k1*(h/2)

        k2 = dinamicaIntracelular(x+h/2, yp2)

        yp3 = y + k2*(h/2)

        k3 = dinamicaIntracelular(x+h/2, yp3)

        yp4 = y + k3*h

        k4 = dinamicaIntracelular(x+h, yp4)

        for j in range(m):
            y[j] = y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

        x = x + h
        xsol = np.append(xsol, x)

        for r in range(len(y)):
            ysol = np.append(ysol, y[r])  
        
        yV = ysol[2::3]
        #yV = np.log10(yV)
        
    return yV

#************* MODELO FIM ************


delta_m   = 0.14
epsilon_m = 0.99
p_m       = 8.18
c_m       = 22.3    
num_par = 4

# define parameter dictionary
parameters = {"num_vars": num_par,
              "names": ['delta', 'epsilon', 'p', 'c'],
              "bounds": [[0, 1],
                         [epsilon_m*0.9, epsilon_m*1.09],
                         [p_m*0.9, p_m*1.1],
                         [c_m*0.9, c_m*1.1]
                         ]
            }

#Gera parametros aleatorios
samples = 100
param_values = saltelli.sample(parameters, samples, calc_second_order=False )

#Y = np.zeros([param_values.shape[0]])
#As solucoes sao 
#201 -> numero de pontos
dimension = samples*(2+num_par)
solucoes = np.zeros([dimension,num_pontos])


for i, X in enumerate(param_values):
    solucoes[i] = RK4thOrder(X)

solucoes = solucoes.T
Si = np.zeros(num_pontos)
Sres = np.zeros([num_pontos,num_par])
for i in range(1, num_pontos):
    Si = sobol.analyze(parameters, solucoes[i], calc_second_order=False)
    Sres[i] = Si['S1']

plt.plot(x_pt, Sres[:,0])
plt.plot(x_pt, Sres[:,1], '--')
plt.plot(x_pt, Sres[:,2], '-.')
plt.plot(x_pt, Sres[:,3])

plt.legend(parameters["names"])

plt.xlabel("Time (days)")
plt.ylabel("First-order Sobol index")
plt.savefig("sobol_index.png", dpi=300)

plt.show()