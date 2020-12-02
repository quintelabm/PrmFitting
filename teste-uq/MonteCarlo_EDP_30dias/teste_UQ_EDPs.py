import uncertainpy as un

import chaospy as cp

import numpy as np

#variaveis de inicializacao
ageFim  = 50
deltaA   = 0.1 #passo no age
ageNpts = int(ageFim/deltaA ) + 1
agePt = np.linspace(0, ageFim, ageNpts)
ageCont = 0

tempoFim = 31
deltaT   = 0.01 #passo no tempo
tempoNpts = int(tempoFim/deltaT) + 1
tempoPt = np.linspace(0, tempoFim, tempoNpts)
tempoCont = 0

def calcIntegral(I,Rp,Rt):
    soma = 0.0
    for a in range(0, ageNpts):
        soma = soma + I[a]*Rp[a]*Rt[a]
    return soma/float(ageNpts)


#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(alpha, r, rho, epsilon_s, epsilon_alpha, epsilon_r, delta):


    #Inicializacao

    #parametros
    s     = 1.3*10**5
    d     = 0.01
    beta  = 5*10**-8
    c     = 22.30

    k     = 0.80 # coeficiente da funcao exponencial de atraso na exportacao de RNA positivo
    tau   = 0.50 # tempo de atraso para a exportacao de RNA positivo
    Rmax  = 50.0 # numero maximo de RNA negativo / complexo de replicacao (Rn)
    sigma = 4.0  # taxa de formacao de complexos de replicacao
    mu_t  = 0.8  # 0.8 no artigo # decaimento natural de Rt
    theta = 1.20 # taxa de disponibilidade para traducao
    mu_c  = 2.8  # 0.89 no artigo # decaimento natural de Rc e Rn

    #variaveis

    I = np.zeros((tempoNpts, ageNpts))
    Rp = np.zeros((tempoNpts, ageNpts))
    Rt = np.zeros((tempoNpts, ageNpts))
    Rn = np.zeros((tempoNpts, ageNpts))
    V = np.zeros(tempoNpts)
    T = np.zeros(tempoNpts)

    def diferencasfinitas():

        #Solve

        kappa_t       = 1.0 #poi[3]
        kappa_c       = 1.0 #poi[4]

        for t in range(1, tempoNpts):

            T[t] = (s - d*T[t-1] - beta*V[t-1]*T[t-1])*deltaT + T[t-1]

            V[t] =  deltaT*((1 - epsilon_s)*rho*calcIntegral(I[t-1],Rp[t-1],Rt[t-1])
            - c*V[t-1]) + V[t-1]

            rho1 = 0.0

            Rp[t][0] = 0
            Rn[t][0] = 0
            Rt[t][0] = 1
            I[t][0]  = beta*V[t]*T[t]

            for a in range(1, ageNpts):

                if(float(a*deltaA) < tau):
                    rho1 = 0
                else:
                    rho1 = rho*(1 - np.exp(-k*(float(a*deltaA) - tau)))

                I[t][a] = (-delta*I[t-1][a] - (I[t-1][a] - I[t-1][a-1])/deltaA)*deltaT + I[t-1][a]

                Rn[t][a] = ((1 - epsilon_r)*r*Rp[t-1][a]*(1 - (Rn[t-1][a]/Rmax)) - kappa_c*mu_c*Rn[t-1][a]
                            - (Rn[t-1][a] - Rn[t-1][a-1])/(deltaA ))*deltaT + Rn[t-1][a]

                Rp[t][a] = ((1 - epsilon_alpha)*alpha*Rn[t-1][a] + sigma*Rt[t-1][a] - theta*Rp[t-1][a]
                            - (1- epsilon_s)*rho1*Rp[t-1][a] - kappa_c*mu_c*Rp[t-1][a]
                            - (Rp[t-1][a] - Rp[t-1][a-1])/(deltaA))*deltaT + Rp[t-1][a]

                Rt[t][a] = (theta*Rp[t-1][a] - sigma*Rt[t-1][a] - (1 - epsilon_s)*rho1*Rt[t-1][a] - kappa_t*mu_t*Rt[t-1][a]
                            - (Rt[t-1][a] - Rt[t-1][a-1])/(deltaA))*deltaT + Rt[t-1][a]
        return V



    #condicoes iniciais

    T0 = 1.3*10**5
    #PAT83:
    V0 = 6.4799*10**5

    I0_t0 = beta*T0*V0

    Rt0_t0 = 1
    Rp0_t0 = 0
    Rn0_t0 = 0


    T[0] = T0
    V[0] = V0

    Rt[0][0] = Rt0_t0
    Rp[0][0] = Rp0_t0
    Rn[0][0] = Rn0_t0
    I[0][0] = I0_t0

    rho1 = 0.00

    # Condicao inicial para as EDPs
    for ageCont in range(1, ageNpts):
        I[0][ageCont] = (beta*T[0]*V[0]*np.exp(-deltaA*ageCont))*deltaA + I[0][ageCont-1]

        if ageCont*deltaA < tau:
            rho1 = 0
        else:
            rho1 = (1 - np.exp(-k*((ageCont*deltaA) - tau)))*rho

        Rn[0][ageCont] = (r*Rp[0][ageCont-1] - r*Rp[0][ageCont-1]*(Rn[0][ageCont-1]/Rmax)
        - mu_c*Rn[0][ageCont-1])*deltaA + Rn[0][ageCont-1]

        Rp[0][ageCont] = (alpha*Rn[0][ageCont-1] + sigma*Rt[0][ageCont-1]
        - theta*Rp[0][ageCont-1]- rho1*Rp[0][ageCont-1]
        - mu_c*Rp[0][ageCont-1])*deltaA + Rp[0][ageCont-1]

        Rt[0][ageCont] = (theta*Rp[0][ageCont-1] - sigma*Rt[0][ageCont-1]
        - rho1*Rt[0][ageCont-1] - mu_t*Rt[0][ageCont-1])*deltaA + Rt[0][ageCont-1]

    #Fim das condicoess iniciais----------------------------------------------+

    # Executa o metedo numerico para resolver
    # Precisa passar como parametro as condicoes iniciais
    diferencasfinitas()


    # Passa para a base log o resultado
    V_log = np.log10(V)

    return tempoPt, V_log

##############FIM MODELO################


model = un.Model(
    run = viralmodelfit,
    labels=["Tempo (dias)",
        "Carga viral (log10)"]
)

# create distributions
epsilon_s = cp.GeneralizedHalfLogistic(shape=1, scale=0.009, shift=0.990)
epsilon_alpha = cp.Bradford(2, 0.9, 0.96)
epsilon_r = cp.Normal(0.3,0.03)
delta = cp.Bradford(2, 0.01, 0.3)
alpha = cp.GeneralizedHalfLogistic(shape=1, scale=9, shift=15)
r = cp.GeneralizedHalfLogistic(shape=1.3, scale=6, shift=5)
rho = cp.GeneralizedHalfLogistic(shape=1, scale=3, shift=9.7)

# define parameter dictionary
parameters = {"alpha": alpha,
        "r": r,
        "rho": rho,
        "epsilon_s": epsilon_s,
        "epsilon_alpha": epsilon_alpha,
        "epsilon_r": epsilon_r,
        "delta": delta
            }

# set up UQ
UQ = un.UncertaintyQuantification(
    model=model,
    parameters=parameters
)

data = UQ.monte_carlo(nr_samples=40)