import uncertainpy as un

import chaospy as cp

import numpy as np

#************* MODELO INICIO ************

# passo
h = 0.1
days = 2

# Dias simulados
x = np.array([0.0, days])

# define uma funcao feval que recebe o nome da equacao a ser avaliada como 
# string e retorna a funcao a ser avaliada
def feval(funcName, *args):
    return eval(funcName)(*args)


# define a resolucao numerica por runge-kutta de quarta ordem
def RK4thOrder(delta, epsilon, p, c, x_range, h):
    
    # define a funcao que contem as equacoes do sistema a ser avaliado
    def dinamicaIntracelular(x, y):
        ## parametros do sistema de 3 equacoes descrito abaixo
    
        s       = 1.3*10**5
        beta    = 5*10**-8
        d       = 0.01
    
        ## inicializa com zeros
        dy = np.zeros(3)
        
        ## equacoes: y[0] = T, y[1] = I, y[2] = V
        dy[0] = s - beta*y[2]*y[0] - d*y[0]
        dy[1] = beta*y[2]*y[0] - delta*y[1]
        dy[2] = (1 - epsilon)*p*y[1] - c*y[2]
    
        return dy
        
    # condicoes iniciais
    T0  = 2.9168*10**6
    I0 = 8.7186*10**5
    
    #PAT83:
    V0  = 4.9139*10**5
    
    yinit = np.array([T0,I0,V0], dtype='f')
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
        k1 = dinamicaIntracelular(x, y)
        
        yp2 = y + k1*(h/2)

        k2 = dinamicaIntracelular(x+h/2, yp2)

        yp3 = y + k2*(h/2)

        k3 = dinamicaIntracelular(x+h/2, yp3)

        yp4 = y + k3*h

        k4 = dinamicaIntracelular(x+h, yp4)

        for j in range(m):
            y[j] = y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

        x = x + h
        xsol = np.append(xsol, x)

        for r in range(len(y)):
            ysol = np.append(ysol, y[r])  
        
        yV = ysol[2::3]
        yV = np.log10(yV)
        
    return [xsol, yV]



#Para pacientes com HCV
#[ts, ys] = RK4thOrder('dinamicaIntracelular', yinit, x, h)

# Separa cada variável em um vetor


#************* MODELO FIM ************

model = un.Model(
    run = RK4thOrder,
    labels=["Tempo (dias)",
        "Carga viral (log10)"]
)

delta_m   = 0.14
epsilon_m = 0.99
p_m       = 8.18
c_m       = 22.3

# create distributions
delta_dist=cp.Uniform(0.01, 1.8)
epsilon_dist=cp.Uniform(0.8, 0.999)
p_dist=cp.Uniform(5, 10)
c_dist=cp.Uniform(15, 25)

# define parameter dictionary
parameters = {"delta": delta_dist,
            "epsilon": epsilon_dist,
            "p": p_dist,
            "c": c_dist,
            "x_range": x,
            "h": h
            }

           
# set up UQ
UQ = un.UncertaintyQuantification(
    model=model,
    parameters=parameters
)

#data = UQ.polynomial_chaos()
data = UQ.quantify()