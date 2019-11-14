
#--- EXAMPLE COST FUNCTIONS ---------------------------------------------------+
import numpy as np
from scipy.spatial import distance

def func1(x):
    # Sphere function, use any bounds, f(0,...,0)=0
    return sum([x[i]**2 for i in range(len(x))])

def func2(x):
    # Beale's function, use bounds=[(-4.5, 4.5),(-4.5, 4.5)], f(3,0.5)=0.
    term1 = (1.500 - x[0] + x[0]*x[1])**2
    term2 = (2.250 - x[0] + x[0]*x[1]**2)**2
    term3 = (2.625 - x[0] + x[0]*x[1]**3)**2
    return term1 + term2 + term3

    # define uma funcao feval que recebe o nome da equacao a ser avaliada como 
# string e retorna a funcao a ser avaliada
def feval(funcName, *args):
    return eval(funcName)(*args)

# define a resolucao numerica por runge-kutta de quarta ordem
def RK4thOrder(func, yinit, x_range, h, poi):
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
        k1 = feval(func, x, y, poi)

        yp2 = y + k1*(h/2)

        k2 = feval(func, x+h/2, yp2, poi)

        yp3 = y + k2*(h/2)

        k3 = feval(func, x+h/2, yp3, poi)

        yp4 = y + k3*h

        k4 = feval(func, x+h, yp4, poi)

        for j in range(m):
            y[j] = y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

        x = x + h
        xsol = np.append(xsol, x)

        for r in range(len(y)):
            ysol = np.append(ysol, y[r])  

    return [xsol, ysol]


def viralmodel(x,y,poi):
    s = 1.3e3
    d = 0.1
    beta = 5.8e-8

    U = y[0]
    I = y[1]
    V = y[2]

    myDelta = poi[0]
    myEpsilon = poi[1]
    myP = poi[2]
    myC = poi[3]

    ## inicializa com zeros
    dy = np.zeros(3)    

    dy[0] = s - beta*U*V - d*U
    dy[1] = beta*U*V - myDelta*I
    dy[2] = (1 - myEpsilon)*myP*I - myC*V
    return dy

def viralmodelfit(x, poi, exp):   
            
    # passo
    h = 0.01

    # Dias simulados
    x = np.linspace(0, 2.0)

    # condicoes iniciais
    T0  = 2.9168*10**6
    I0 = 8.7186*10**5
    V0  = 6.9139*10**5
    yinit = np.array([T0,I0,V0], dtype='f')

    # Chama o método de runge-kutta definido com a função e as condições iniciais

    [ts, ys] = RK4thOrder('viralmodel', yinit, x, h, poi)
    
    #calculate the euclidean norm
    node = len(yinit)
    ys3 = ys[2::node]
    ys3 = np.log10(ys3)

    # improvisadamente estou pegando diferença entre primeiro e ultimo valor da simulacao
    # seria interessante ver uma forma de pegar nos mesmos pontos 
    V = [ys3[0], ys3[-1]]
    exp = [exp[0], exp[-1]]
    #print('V = ', V, ' exp = ' , exp)
    #return np.linalg.norm(V-exp)
    dst = distance.euclidean(V, exp)

    return dst
