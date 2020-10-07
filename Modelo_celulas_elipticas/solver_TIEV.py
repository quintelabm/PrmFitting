import matplotlib.pyplot as plt

from scipy.integrate import odeint

import numpy as np

import seaborn as sns

sns.set()


def dinamica_Extracelular(y, t, beta, delta, epsilon, p, c, k):
    
    # inicializa com zeros
    dy = np.zeros(4)
    
    # equacoes: y[0] = T, y[1] = I, y[2] = E, y[3] = V
     
    dy[0] = - beta*y[3]*y[0] 
    dy[1] = k*y[2] - delta*y[1]
    dy[2] = beta*y[3]*y[0] - k*y[2]
    dy[3] = (1 - epsilon)*p*y[1] - c*y[2]

    return dy

def plot(t_pts, solve):
    
    v = solve[:,3]
    log_v = np.log10(v)
    plt.plot(t_pts, log_v, '-g', label='virus')
    plt.title("Carga viral no tempo")
    plt.ylabel("Carga viral (log10)")
    plt.xlabel("Tempo (dias)")
    plt.legend()
    plt.show()
    
    return 0;
    
def solver(beta, delta, epsilon, p, c, k):
    
    # passo
    h = 0.1
    days = 30
    # Dias simulados
    t_range = np.linspace(0, days, int(days/h))
    
    # condicoes iniciais
    T0 = 2.9168*10**6
    #AVERAGE_PAT
    V0 = 10**6.47991433
    E0 = beta*T0*V0
    I0 = k*E0
    
    
    yinit = np.array([T0,I0,E0,V0], dtype='f')
        
    return t_range, odeint(dinamica_Extracelular, yinit, t_range, args=(beta, delta, epsilon, p, c, k))

if __name__ == "__main__":
    
    delta = 0.07
    epsilon = 0.999
    p = 12.0
    c = 19.0
    beta = 5*10**-8
    k = 4
    
    t_range, sol = solver(beta, delta, epsilon, p, c, k)
    
    plot(t_range, sol)