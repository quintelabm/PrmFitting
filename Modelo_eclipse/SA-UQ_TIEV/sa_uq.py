from scipy.integrate import odeint

import numpy as np

import seaborn as sns

import chaospy as cp

import uncertainpy as un

sns.set()


def dinamica_Extracelular(y, t, beta, delta, epsilon, p, c, k):
    
    # inicializa com zeros
    dy = np.zeros(4)
    
    # equacoes: y[0] = T, y[1] = I, y[2] = E, y[3] = V

    dy[0] = - beta*y[3]*y[0]
    dy[1] = k*y[2] - delta*y[1]
    dy[2] = beta*y[3]*y[0] - k*y[2]
    dy[3] = (1 - epsilon)*p*y[1] - c*y[3]  

    return dy


def solver(delta, epsilon, p, c, k, V0, days):
    
    # passo
    h = 0.1
    # Dias simulados
    t_range = np.linspace(0, days, int(days/h))
    
    # condicoes iniciais
    T0 = 1.85*10**7
    beta = 5*10**-8
    E0 = beta*T0*V0
    I0 = k*E0
        
    yinit = np.array([T0,I0,E0,V0], dtype='f')
    sol = odeint(dinamica_Extracelular, yinit, t_range, args=(beta, delta, epsilon, p, c, k))
    v = sol[:,3]
    log_v = np.log10(v)
    return t_range, log_v
     


if __name__ == "__main__":
    
    model = un.Model(
        run = solver,
        labels=["Tempo (dias)",
            "Carga viral (log10)"]
    )
    
    delta_m   = 0.096
    epsilon_m = 0.999
    p_m       = 3.13
    c_m       = 9.6
    k_m       = 6.04
        
    # create distributions
    #delta_dist   = cp.Bradford(6, 0.01, 1)
    delta_dist   = cp.Beta(1.8, 2, 0.01, 0.2)
    epsilon_dist = cp.GeneralizedHalfLogistic(shape=1, scale=0.0003, shift=0.9996)
    p_dist       = cp.Bradford(6, 3, 5)
    c_dist       = cp.GeneralizedHalfLogistic(shape=1, scale=5, shift=5)
    k_dist       = cp.GeneralizedHalfLogistic(shape=0.8, scale=4, shift=2)
    
    # define parameter dictionary
    parameters = {"delta": delta_dist,
                "epsilon": epsilon_dist,
                "p": p_dist,
                "c": c_dist,
                "k": k_dist,
                "V0": 10**6.47991433,#average pats
                "days": 30 
                }
                   
    # set up UQ
    UQ = un.UncertaintyQuantification(
        model=model,
        parameters=parameters
    )
    
    data = UQ.monte_carlo(nr_samples=500)