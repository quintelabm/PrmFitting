import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()

# define uma funcao feval que recebe o nome da equacao a ser avaliada como 
# string e retorna a funcao a ser avaliada
def feval(funcName, *args):
    return eval(funcName)(*args)


# define a resolucao numerica por runge-kutta de quarta ordem
def RK4thOrder(func, yinit, x_range, h):
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
        k1 = feval(func, x, y)

        yp2 = y + k1*(h/2)

        k2 = feval(func, x+h/2, yp2)

        yp3 = y + k2*(h/2)

        k3 = feval(func, x+h/2, yp3)

        yp4 = y + k3*h

        k4 = feval(func, x+h, yp4)

        for j in range(m):
            y[j] = y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

        x = x + h
        xsol = np.append(xsol, x)

        for r in range(len(y)):
            ysol = np.append(ysol, y[r])  

    return [xsol, ysol]

# define a funcao que contem as equacoes do sistema a ser avaliado
def dinamicaIntracelular(x, y):
    ## parametros do sistema de 6 equacoes descrito acima
    #VALORES DA FIGURA 2 DO MULTISCALE MODEL ALTERADOS!!!!!!
    s       = 1.3e5
    beta    = 5.8e-8
    d       = 0.01

    delta   = 0.07    #0.14
    epsilon = 0.999   #0.996
    p       = 12.0    #8.18
    c       = 19.0   #22.3
    ## inicializa com zeros
    dy = np.zeros(3)

    T = y[0]
    I = y[1]
    V = y[2]

    dy[0] = s - beta*V*T - d*T
    dy[1] = beta*V*T - delta*I
    dy[2] = (1-epsilon)*p*I - c*V

    return dy

# passo
h = 0.1

# Dias simulados
x = np.array([0.0, 31.0])

# condicoes iniciais
T0  = 2.9168*10**6
I0 = 8.7186*10**5

#V0  = 10**6.3780 # PATB06  1
#V0  = 10**6.3541 # PATB16  2
#V0  = 10**6.4885 # PATB17  3

#V0  = 10**6.394490 # PATC05  4
#V0  = 10**6.839431 # PATC06  5
#V0  = 10**6.424965 # PATC09  6

V0  = 10**6.47991433 # MEAN_PATS

yinit = np.array([T0,I0,V0], dtype='f')

# Chama o método de runge-kutta definido com a função e as condições iniciais

[ts, ys] = RK4thOrder('dinamicaIntracelular', yinit, x, h)

# Separa cada variável em um vetor

node = len(yinit)
ys1 = ys[0::node]
ys2 = ys[1::node]
ys3 = ys[2::node]

ys3 = np.log10(ys3)

plt.figure()

#plt.plot(ts, ys1, 'r')
#plt.plot(ts, ys2, 'b')
plt.plot(ts, ys3)
#plt.legend(["Target", "Infectada", "Virus"])
plt.xlim(-1.0, 32.0)

# Tempo Experimentos

#30 dias#

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
MEAN_PATS = [6.47991433, 6.42897983, 5.857277, 4.63766183, 4.08108333, 3.38275467, 3.18124267, 2.88121, 2.66053917,
             2.54078967, 2.487822, 2.21939417, 1.90103167, 1.7158415]

#plt.plot(t_exp1, PATB06, 'o')
#plt.plot(t_exp2, PATB16, 'o')
#plt.plot(t_exp3, PATB17, 'o')

#plt.plot(t_exp4, PATC05, 'o')
#plt.plot(t_exp5, PATC06, 'o')
#plt.plot(t_exp6, PATC09, 'o')

plt.plot(t_exp6, MEAN_PATS, 'o')
plt.legend(["Model", "Average Data"])

plt.xlim(-1.0, 32.0)
#plt.legend(["Model", "PAT 1", "PAT 2", "PAT 3", "PAT 4", "PAT 5", "PAT 6"])
plt.xlabel("Time after beginning of combination therapy (days)")
plt.ylabel("HCV RNA ($log_{10}$ UI/mL)")



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

# plt.ylim(0,120)
# plt.errorbar(x, y, yerr = yerr, marker='s', fmt='.k', uplims=True, lolims=True)
# plt.errorbar(x, y, yerr = yerr, marker='s', fmt='.k')
plt.errorbar(x, y, yerr=yerr, fmt='.k')



#plt.savefig("D:\Results/average_stdev.png", dpi=300)

plt.show()