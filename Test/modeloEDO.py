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
    s       = 1.3*10**5
    beta    = 5*10**-8
    d       = 0.01
    delta   = 0.51
    epsilon = 0.33
    p       = 0.28
    c       = 4.92
    
    ## inicializa com zeros
    dy = np.zeros(3)
    
    ## equacoes: y[0] = T, y[1] = I, y[2] = V
    dy[0] = s - beta*y[2]*y[0] - d*y[0]
    dy[1] = beta*y[2]*y[0] - delta*y[1]
    dy[2] = (1 - epsilon)*p*y[1] - c*y[2]

    return dy

# passo
h = 0.01

# Dias simulados
x = np.array([0.0, 2.0])

# condicoes iniciais
T0  = 2.9168*10**6
I0 = 8.7186*10**5
V0  = 6.9139*10**5
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
plt.xlim(x[0], x[1])

# Tempo Experimentos
t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2 ]
#t_exp = [0 0.083 0.167 0.25 0.333 0.5 0.667 1 1.5 2 3 6]

# --- for patient 8
PAT8 = [5.64, 5.31, 4.23, 3.36, 3.14, 2.86, 2.75, 2.50, 2.32, 1.56]
#PAT = [5.64 5.31 4.23 3.36 3.14 2.86 2.75 2.50 2.32 1.56 1.53 1.40]

# --- for patient 42
#PAT42 = [5.65 5.00 3.98 3.84 2.94 2.82 2.87 2.53 2.31 2.61];
#PAT = [5.65 5.00 3.98 3.84 2.94 2.82 2.87 2.53 2.31 2.61 2.29 2.18];

# --- for patient 68
#PAT68 =  [7.15 7.02 6.19 5.50 4.96 4.29 4.11 3.75 3.68 3.35];
#PAT =  [7.15 7.02 6.19 5.50 4.96 4.29 4.11 3.75 3.68 3.35 3.07 2.26];

# --- for patient 69
#PAT69 = [6.14 5.87 4.73 4.17 3.55 3.14 2.87 2.60 2.55 2.58];
#PAT = [6.14 5.87 4.73 4.17 3.55 3.14 2.87 2.60 2.55 2.58 2.49 3.57];

# --- for patient 83
#PAT83 = [5.45 5.38 4.73 4.00 3.39 2.89 2.68 2.72 2.97 1.93];
#PAT = [5.45 5.38 4.73 4.00 3.39 2.89 2.68 2.72 2.97 1.93 2.01 1.54];

#V <- c(3.241529, 2.555326, 2.547130, 2.541043) #pat1
PAT1 = [3.241529, 2.555326, 2.547130, 2.541043]
#V <- c(3.522593, 2.920457, 2.547203, 2.615682) #pat2
#V <- c(3.555587, 3.316164, 2.595167, 2.543369) #pat3
#V <- c(4.175440, 2.992317, 2.744784, 2.543393) #pat4

plt.plot(t_exp, PAT8, 'ro')

plt.show()