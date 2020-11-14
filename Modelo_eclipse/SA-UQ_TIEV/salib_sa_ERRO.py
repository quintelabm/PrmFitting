from SALib.sample import saltelli

from SALib.analyze import sobol

import matplotlib.pyplot as plt

import seaborn as sns

import numpy as np

import chaospy as cp

from scipy.integrate import odeint

sns.set()
#************* MODELO INICIO ************

def dinamica_Extracelular(y, t, beta, delta, epsilon, p, c, k):
    
    # inicializa com zeros
    dy = np.zeros(4)
    
    # equacoes: y[0] = T, y[1] = I, y[2] = E, y[3] = V

    dy[0] = - beta*y[3]*y[0]
    dy[1] = k*y[2] - delta*y[1]
    dy[2] = beta*y[3]*y[0] - k*y[2]
    dy[3] = (1 - epsilon)*p*y[1] - c*y[3]  

    return dy

def solver(param):
    
    delta, epsilon, p, c, k = param
    # Dias simulados

    V0 = 10**6.47991433#average pats
    days = 30
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
    return log_v

#************* MODELO FIM ************

if __name__ == "__main__":
    
    
    # create distributions
    delta_dist   = cp.Bradford(6, 0.01, 1)
    #delta_dist   = cp.Bradford(6, 0.01, 0.5)
    epsilon_dist = cp.GeneralizedHalfLogistic(shape=1, scale=0.0003, shift=0.9996)
    p_dist       = cp.Bradford(6, 3, 5)
    c_dist       = cp.GeneralizedHalfLogistic(shape=1, scale=5, shift=5)
    k_dist       = cp.GeneralizedHalfLogistic(shape=0.8, scale=4, shift=2)
    num_par = 5

    parameters = {"num_vars": num_par,
    "names": ['delta', 'epsilon', 'p', 'c', 'k'],
    "bounds": [delta_dist, epsilon_dist,
                p_dist, c_dist, k_dist]
    }

    #Gera parametros aleatorios
    samples = 100
    param_values = saltelli.sample(parameters, samples, calc_second_order=False )

    # passo
    h = 0.1
    days = 30
    num_pontos = int(days/h)
    t_pontos = np.linspace(0, days, int(days/h))

    dimension = samples*(2+num_par)
    solucoes = np.zeros([dimension,num_pontos])


    for i, X in enumerate(param_values):
        solucoes[i] = solver(X)

    solucoes = solucoes.T
    Si = np.zeros(num_pontos)
    Sres = np.zeros([num_pontos,num_par])
    for i in range(1, num_pontos):
        Si = sobol.analyze(parameters, solucoes[i], calc_second_order=False)
        Sres[i] = Si['S1']

    plt.plot(t_pontos, Sres[:,0])
    plt.plot(t_pontos, Sres[:,1], '--')
    plt.plot(t_pontos, Sres[:,2], '-.')
    plt.plot(t_pontos, Sres[:,3])

    plt.legend(parameters["names"])

    plt.xlabel("Time (days)")
    plt.ylabel("First-order Sobol index")
    plt.savefig("sobol_index.png", dpi=300)

    plt.show()