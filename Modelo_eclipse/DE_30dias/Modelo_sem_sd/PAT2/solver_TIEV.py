import matplotlib.pyplot as plt

from scipy.integrate import odeint
from scipy.spatial import distance
from scipy.interpolate import InterpolatedUnivariateSpline

import numpy as np

import seaborn as sns

import os

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

def plot(t_pts, solve):
    
    v = solve[:,3]
    log_v = np.log10(v)
    plt.plot(t_pts, log_v, '-g')
    plt.title("Carga viral no tempo")
    plt.ylabel("Carga viral ($log_{10}$ UI/mL)")
    plt.xlabel("Tempo (dias)")
    plt.legend(["Data", "Model"])

    
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
    return t_range, odeint(dinamica_Extracelular, yinit, t_range, args=(beta, delta, epsilon, p, c, k))

def custo(param_adj, t_exp, data_exp, days):
    
    t_range, sol = solver(param_adj[0], param_adj[1], param_adj[2], param_adj[3], param_adj[4], 10**data_exp[0], days)

    solV = sol[:,3]
    for v in solV:
        if v < 0:
            return 100000
    logsolV = np.log10(solV)
    
    ius = InterpolatedUnivariateSpline(t_exp, data_exp)
    
    yi = ius(t_range)
    
    dst = distance.euclidean(logsolV, yi)
    return dst


if __name__ == "__main__":
    
    delta = 0.117803334
    epsilon = 0.998846228
    p = 1.209706709
    c = 12.29316611
    k = 7.518650959

    beta = 5*10**-8
    days = 30

    # --- for patient B06
    PATB06 = [6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378,
              2.3404, 2.3345, 2.2355, 2.0492, 2.1173]
    t_exp1 = [0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85, 13.85, 20.86, 27.96]
    # --------------------------------------------------------------------------------------------------------
    # --- for patient B16
    PATB16 = [6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120,
              3.0302, 2.7528, 2.3838, 2.1818, 1.9243]
    t_exp2 = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01, 7.00, 9.99, 14.00, 20.98, 30.01]
    # --------------------------------------------------------------------------------------------------------
    # --- for patient B17
    PATB17 = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432,
              2.7474, 2.7016, 2.3541, 2.0453, 1.4914]
    t_exp3 = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03]
    # --------------------------------------------------------------------------------------------------------
    # --- for patient C05
    PATC05 = [6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339,
              2.609594, 2.527630, 2.743510, 2.694605, 2.227887, 1.863323]
    t_exp4 = [0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 2.86, 6.91, 8.86, 9.88, 13.90, 21.87, 30.87]
    # --------------------------------------------------------------------------------------------------------
    # --- for patient C06
    PATC06 = [6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317, 2.484300,
              2.509203, 2.369216, 1.949390, 1.623249, 1.556303]
    t_exp5 = [0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94, 5.99, 7.01, 9.91, 15.00, 22.00, 28.02]
    # --------------------------------------------------------------------------------------------------------
    # --- for patient C09
    PATC09 = [6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341,
              2.089905, 2.025306, 1.698970, 1.278754, 1.342423]
    t_exp6 = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94, 7.94, 9.95, 14.94, 24.03, 30.94]
    # --------------------------------------------------------------------------------------------------------
    # --- average of pats
    AVERAGE_PATS = [6.47991433, 6.42897983, 5.857277, 4.63766183, 4.08108333, 3.38275467, 3.18124267, 2.88121,
                    2.66053917, 2.54078967, 2.487822, 2.21939417, 1.90103167, 1.7158415]

    #plt.plot(t_exp1, PATB06, 'o')
    plt.plot(t_exp2, PATB16, 'o')
    #plt.plot(t_exp3, PATB17, 'o')
    #plt.plot(t_exp4, PATC05, 'o')
    #plt.plot(t_exp5, PATC06, 'o')
    #plt.plot(t_exp6, PATC09, 'o')
    #plt.plot(t_exp6, AVERAGE_PATS, 'o')
    
    t_range, sol = solver(delta, epsilon, p, c, k, 10**PATB16[0], days)

    plot(t_range, sol)

'''
    # Tempo das medicoes e da discretizacao do modelo
    x = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94, 7.94, 9.95, 14.94, 24.03, 30.94]

    # Media do resultado experimental
    y = [6.47991433, 6.42897983, 5.857277, 4.63766183, 4.08108333, 3.38275467,
         3.18124267, 2.88121, 2.66053917, 2.54078967, 2.487822, 2.21939417,
         1.90103167, 1.7158415]

    # Desvio padrao dos dados
    yerr = [0.182129952, 0.182630388, 0.337284052, 0.390495914, 0.307469046, 0.18942229,
            0.200947972, 0.239112485, 0.222385492, 0.324597051, 0.294209462, 0.350716883,
            0.372253077, 0.297196377]

    plt.errorbar(x, y, yerr=yerr, fmt='.k')
    '''

cwd = os.getcwd()
plt.savefig(cwd+"/P2_Average.png", dpi=300)

plt.show()
