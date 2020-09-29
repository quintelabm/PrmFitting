import uncertainpy as un

import chaospy as cp

import numpy as np

from scipy.integrate import odeint


# passo
h = 0.1
days = 5

# Dias simulados
t_range = np.array([0.0, days])

# condicoes iniciais
T0  = 2.9168*10**6
I0 = 8.7186*10**5
#AVERAGE_PAT
V0 = 2.563292*10**6

yinit = np.array([T0,I0,V0], dtype='f')


# define a funcao que contem as equacoes do sistema a ser avaliado
def dinamicaIntracelular(y, t, delta, epsilon, p, c):
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


def solver(delta, epsilon, p, c):
    y_pt = odeint(dinamicaIntracelular,yinit, t_range, args=(delta, epsilon, p, c))
    y_v = y_pt[:,2]
    y_v = np.log10(y_v)
    return t_range, y_v

model = un.Model(
    run = solver,
    labels=["Tempo (dias)",
        "Carga viral (log10)"]
)

delta_m   = 0.07
epsilon_m = 0.999
p_m       = 12.0
c_m       = 19.0

# create distributions
delta_dist=cp.Uniform(delta_m*0.9, delta_m*1.1)
epsilon_dist=cp.Uniform(epsilon_m*0.9, epsilon_m)
p_dist=cp.Uniform(p_m*0.9, p_m*1.1)
c_dist=cp.Uniform(c_m*0.9, c_m*1.1)

# define parameter dictionary
parameters = {"delta": delta_dist,
            "epsilon": epsilon_dist,
            "p": p_dist,
            "c": c_dist
            }

           
# set up UQ
UQ = un.UncertaintyQuantification(
    model=model,
    parameters=parameters
)

#data = UQ.polynomial_chaos()
data = UQ.monte_carlo(nr_samples=100)