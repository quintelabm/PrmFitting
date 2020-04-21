import matplotlib.pyplot as plt
import numpy as np

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

    delta   = 0.063
    epsilon = 0.997
    p       = 1.819
    c       = 27.377
    
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
V0  = 10**6.3541
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
plt.plot(ts, ys3, 'g')
#plt.legend(["Target", "Infectada", "Virus"])
plt.xlim(-1.0, x[1])

# Tempo Experimentos

#2 dias#
'''
t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2]

# --- for patient 8
PAT8 = [5.64, 5.31, 4.23, 3.36, 3.14, 2.86, 2.75, 2.50, 2.32, 1.56]
'''

#30 dias#

# --- for patient B07
#PATB07= [4.2227, 4.0733, 3.0599, 2.2122, 1.7924, 1.6435, 1.7160, 1.6435, 1.2304, 1.0792, 1.0000, 1.0000, 1.0000, 1.0000]
#PAT = PATB07

#t_exp = [0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.91, 6.92, 10.94, 13.96, 20.91, 27.98]

#V0 = 4.1510

#--------------------------------------------------------------------------------------------------------

# --- for patient B09
PATB09 = [6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385, 4.1949, 3.0931, 2.1790, 1.5682, 1.5051]
PAT = PATB09

t_exp = [0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90, 6.94, 10.96, 14.95, 20.96, 28.03]

#V0 = 6.3541

#--------------------------------------------------------------------------------------------------------

# --- for patient B08
#PATB08 = [5.6546, 5.6486, 4.6203, 3.6876, 3.1526, 2.6355, 2.5302, 2.6294, 2.2788, 2.0719, 1.9031, 1.7924, 1.6335, 1.5682]
#PAT = PATB08

#t_exp = [0.04, 0.08, 0.17, 0.34, 0.51, 1.00, 1.50, 2.91, 4.92, 6.91, 10.87, 13.86, 20.87, 25.92]

#V0 = 5.7831

#--------------------------------------------------------------------------------------------------------

# --- for patient B16
#PATB16 = [6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120, 3.0302, 2.7528, 2.3838, 2.1818, 1.9243]
#PAT = PATB16

#t_exp = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01, 7.00, 9.99, 14.00, 20.98, 30.01]

#V0 = 6.2943

#--------------------------------------------------------------------------------------------------------

# --- for patient B17
#PATB17 = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432, 2.7474, 2.7016, 2.3541, 2.0453, 1.4914]
#PAT = PATB17

#t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 603, 7.02, 9.98, 14.04, 21.03, 30.03]

#V0 = 6.1927

#--------------------------------------------------------------------------------------------------------

# --- for patient B06
#PATB06 = [6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378, 2.3404, 2.3345, 2.2355, 2.0492, 2.1173]
#PAT = PATB06

#t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85, 13.85, 20.86, 27.96]

#V0 = 6.2584

#--------------------------------------------------------------------------------------------------------

# --- for patient C05
#PATC05 = [6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339, 2.609594, 2.527630, 2.743510, 2.694605, 2.227887, 1.863323]
#PAT = PATC05

#t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 2.86, 6.91, 8.86, 9.88, 13.90, 21.87, 30.87]

#V0 = 6.208568

#--------------------------------------------------------------------------------------------------------

# --- for patient C06
#PATC06 = [6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317, 2.484300, 2.509203, 2.369216, 1.949390, 1.623249, 1.556303]
#PAT = PATC06

#t_exp = [0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94, 5.99, 7.01, 9.91, 15.00, 22.00, 28.02]

#V0 = 6.808956

#--------------------------------------------------------------------------------------------------------

# --- for patient C09
#PATC09 = [6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341, 2.089905, 2.025306, 1.698970, 1.278754, 1.342423]
#PAT = PATC09

#t_exp = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94, 7.94, 9.95, 14.94, 24.03, 30.94]

#V0 = 6.296968

#--------------------------------------------------------------------------------------------------------

# --- for patient C10
#PATC10 = [5.583842, 5.455846, 4.754914, 3.447468, 2.749736, 2.311754, 2.158362, 1.944483, 1.707570, 1.322219, 1.322219, 1.000000, 1.000000, 1.000000]
#PAT = PATC10

#t_exp = [0.04, 0.09, 0.18, 0.33, 0.50, 1.00, 1.50, 2.96, 3.98, 7.00, 9.98, 13.97, 2196, 30.02]

#V0 = 5.547272

plt.plot(t_exp, PAT, 'ro')

plt.show()