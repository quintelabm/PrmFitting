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
    s       = 1.3e3
    beta    = 5.8e-8
    d       = 0.1
    epsilon = 0.996

    delta   = 0.59
    p       = 5.0
    c       = 23.0
    
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
t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2]

# --- for patient 8
PAT8 = [5.64, 5.31, 4.23, 3.36, 3.14, 2.86, 2.75, 2.50, 2.32, 1.56]

plt.plot(t_exp, PAT8, 'ro')

plt.show()