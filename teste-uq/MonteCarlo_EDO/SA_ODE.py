from scipy.integrate import odeint

import numpy as np

import seaborn as sns

import chaospy as cp

import uncertainpy as un

sns.set()


def dinamica_Extracelular(y, t, delta, epsilon, p, c):
    # parametros do sistema de 3 equacoes descrito abaixo

    s       = 1.3*10**5
    beta    = 5*10**-8
    d       = 0.01

    # inicializa com zeros
    dy = np.zeros(3)
    
    # equacoes: y[0] = T, y[1] = I, y[2] = V
    dy[0] = s - beta*y[2]*y[0] - d*y[0]
    dy[1] = beta*y[2]*y[0] - delta*y[1]
    dy[2] = (1 - epsilon)*p*y[1] - c*y[2]

    return dy

    
def solver(delta, epsilon, p, c):
    
    # passo
    h = 0.01
    days = 5
    # Dias simulados
    t_range = np.linspace(0, days, int(days/h))
    
    # condicoes iniciais
    T0  = 2.9168*10**6
    I0 = 8.7186*10**5
    #AVERAGE_PAT
    V0 = 2.563292*10**6
    yinit = np.array([T0,I0,V0], dtype='f')
    
    sol = odeint(dinamica_Extracelular, yinit, t_range, args=(delta, epsilon, p, c))
    
    v = sol[:,2]
    log_v = np.log10(v)
    
    return t_range, log_v

if __name__ == "__main__":
    
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
    
    data = UQ.monte_carlo(nr_samples=100)