import matplotlib.pyplot as plt

from scipy.integrate import odeint

import numpy as np

import seaborn as sns

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

def plot(t_pts, solve):
    #
    #
    ###Fazer o plot dos gráficos COM LEGENDA E LABELS 
    #
    #
    v = solve[:,2]
    log_v = np.log10(v)
    plt.plot(t_pts, log_v, '-g', label='virus')
    plt.title("Carga viral no tempo")
    plt.ylabel("Carga viral (log10)")
    plt.xlabel("Tempo (dias)")
    plt.legend()
    plt.show()
    
    return 0;
    
def solver():
    
    # passo
    h = 0.1
    days = 30
    # Dias simulados
    t_range = np.linspace(0, days, int(days/h))
    
    # condicoes iniciais
    T0  = 2.9168*10**6
    I0 = 8.7186*10**5
    #AVERAGE_PAT
    V0 = 10**6.47991433
    yinit = np.array([T0,I0,V0], dtype='f')
    
    delta = 0.07
    epsilon = 0.999
    p = 12.0
    c = 19.0
    
    return t_range, odeint(dinamica_Extracelular, yinit, t_range, args=(delta, epsilon, p, c))

if __name__ == "__main__":
    
    t_range, sol = solver()
    
    plot(t_range, sol)