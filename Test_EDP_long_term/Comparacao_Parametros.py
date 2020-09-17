import numpy as np 
import matplotlib.pyplot as plt

# variaveis de inicializacao
ageFim  = 50
deltaA   = 0.1 # passo no age
ageNpts = int(ageFim/deltaA) + 1
agePt = np.linspace(0, ageFim, ageNpts)
ageCont = 0

tempoFim = 31
deltaT   = 0.01 # passo no tempo
tempoNpts = int(tempoFim/deltaT) + 1 
tempoPt = np.linspace(0, tempoFim, tempoNpts)
tempoCont = 0

def calcIntegral(I,Rp,Rt):
    soma = 0.0
    for a in range(0, ageNpts):
        soma = soma + I[a]*Rp[a]*Rt[a]
    return soma/float(ageNpts)


#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0):
    
    #Inicializacao

    #parametros
    s       = 1.3*10**5  # taxa de crescimento das celulas alvo
    d       = 0.01       # taxa de morte das celulas alvo
    beta    = 5*10**-8   # taxa de formacao de complexos de replicacao
    c       = 19         # taxa de eliminacao do virus pelo sistema imune

    kappa_t = 1.0   # fator para aumentar a degradacao de Rt
    kappa_c = 1.0   # fator para aumentar a degradacao de Rp e Rn no complexo de replicacao

    k       = 0.80  # coeficiente da funcao exponencial de atraso na exportacao de RNA positivo
    tau     = 0.50  # tempo de atraso para a exportacao de RNA positivo
    Rmax    = 50.0  # numero maximo de RNA negativo / complexo de replicacao (Rn)
    sigma   = 4.0   # taxa de formacao de complexos de replicacao
    mu_t    = 0.8   # 0.8 no artigo # decaimento natural de Rt
    theta   = 1.20  # taxa de disponibilidade para traducao
    mu_c    = 2.8   # 0.89 no artigo # decaimento natural de Rc e Rn

    # parametros que serao ajustados

    epsilon_s     = poi[0]  # efetividade da terapia em diminuir ou bloquear a exportacao de Rp
    epsilon_alpha = poi[1]  # efetividade da terapia em diminuir ou bloquear a replicacao de Rp
    epsilon_r     = poi[2]  # efetividade da terapia em diminuir ou bloquear a replicacao de Rn
    delta         = poi[3]  # taxa de morte da celula infectada
    alpha         = poi[4]  # taxa de replicacao do Rp
    r             = poi[5]  # taxa de replicacao do Rn / complexo de replicacao (Rn)
    rho           = poi[6]  # taxa de exportacao do Rp

    #variaveis 

    I = np.zeros((tempoNpts, ageNpts))
    Rp = np.zeros((tempoNpts, ageNpts))
    Rt = np.zeros((tempoNpts, ageNpts))
    Rn = np.zeros((tempoNpts, ageNpts))
    V = np.zeros(tempoNpts)
    T = np.zeros(tempoNpts)

    def diferencasfinitas(poi):
        
        #Solve

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
    diferencasfinitas(poi)

    # Se algum valor de V for negativo, defina um custo alto para o conjunto de parametros
    for pos in V:
        if not (pos >= 0):
            return 1000000
    
    # Passa para a base log o resultado

    V_log = np.log10(V)

    t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03] #PATB17

    plt.plot(tempoPt,V_log, 'r-', label = "$ε_s$=0.999, $ε_α$=0.920")


# Copy best result (epsilon_s, epsilon_alpha, epsilon_r, delta, alpha, r, rho)
poi = [0.999, 0.920, 0.300, 0.07, 22.0, 11.1, 12]

#[(0.990, 0.999), (0.900, 0.960), (0.200, 0.400), (0.01, 0.3), (15, 24), (5.0, 11.1), (9.7, 12.4)]  #bounds dos parâmetros

# experimental data

# PATB17
exp = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432, 2.7474, 2.7016, 2.3541, 2.0453, 1.4914]

viralmodelfit(poi, exp, 10**exp[0])


################################################


# recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0):
    # Inicializacao

    # parametros
    s = 1.3 * 10 ** 5  # taxa de crescimento das celulas alvo
    d = 0.01  # taxa de morte das celulas alvo
    beta = 5 * 10 ** -8  # taxa de formacao de complexos de replicacao
    c = 19  # taxa de eliminacao do virus pelo sistema imune

    kappa_t = 1.0  # fator para aumentar a degradacao de Rt
    kappa_c = 1.0  # fator para aumentar a degradacao de Rp e Rn no complexo de replicacao

    k = 0.80  # coeficiente da funcao exponencial de atraso na exportacao de RNA positivo
    tau = 0.50  # tempo de atraso para a exportacao de RNA positivo
    Rmax = 50.0  # numero maximo de RNA negativo / complexo de replicacao (Rn)
    sigma = 4.0  # taxa de formacao de complexos de replicacao
    mu_t = 0.8  # 0.8 no artigo # decaimento natural de Rt
    theta = 1.20  # taxa de disponibilidade para traducao
    mu_c = 2.8  # 0.89 no artigo # decaimento natural de Rc e Rn

    # parametros que serao ajustados

    epsilon_s = poi[0]  # efetividade da terapia em diminuir ou bloquear a exportacao de Rp
    epsilon_alpha = poi[1]  # efetividade da terapia em diminuir ou bloquear a replicacao de Rp
    epsilon_r = poi[2]  # efetividade da terapia em diminuir ou bloquear a replicacao de Rn
    delta = poi[3]  # taxa de morte da celula infectada
    alpha = poi[4]  # taxa de replicacao do Rp
    r = poi[5]  # taxa de replicacao do Rn / complexo de replicacao (Rn)
    rho = poi[6]  # taxa de exportacao do Rp

    # variaveis

    I = np.zeros((tempoNpts, ageNpts))
    Rp = np.zeros((tempoNpts, ageNpts))
    Rt = np.zeros((tempoNpts, ageNpts))
    Rn = np.zeros((tempoNpts, ageNpts))
    V = np.zeros(tempoNpts)
    T = np.zeros(tempoNpts)

    def diferencasfinitas(poi):

        # Solve

        for t in range(1, tempoNpts):

            T[t] = (s - d * T[t - 1] - beta * V[t - 1] * T[t - 1]) * deltaT + T[t - 1]

            V[t] = deltaT * ((1 - epsilon_s) * rho * calcIntegral(I[t - 1], Rp[t - 1], Rt[t - 1])
                             - c * V[t - 1]) + V[t - 1]

            rho1 = 0.0

            Rp[t][0] = 0
            Rn[t][0] = 0
            Rt[t][0] = 1
            I[t][0] = beta * V[t] * T[t]

            for a in range(1, ageNpts):

                if (float(a * deltaA) < tau):
                    rho1 = 0
                else:
                    rho1 = rho * (1 - np.exp(-k * (float(a * deltaA) - tau)))

                I[t][a] = (-delta * I[t - 1][a] - (I[t - 1][a] - I[t - 1][a - 1]) / deltaA) * deltaT + I[t - 1][a]

                Rn[t][a] = ((1 - epsilon_r) * r * Rp[t - 1][a] * (1 - (Rn[t - 1][a] / Rmax)) - kappa_c * mu_c *
                            Rn[t - 1][a]
                            - (Rn[t - 1][a] - Rn[t - 1][a - 1]) / (deltaA)) * deltaT + Rn[t - 1][a]

                Rp[t][a] = ((1 - epsilon_alpha) * alpha * Rn[t - 1][a] + sigma * Rt[t - 1][a] - theta * Rp[t - 1][a]
                            - (1 - epsilon_s) * rho1 * Rp[t - 1][a] - kappa_c * mu_c * Rp[t - 1][a]
                            - (Rp[t - 1][a] - Rp[t - 1][a - 1]) / (deltaA)) * deltaT + Rp[t - 1][a]

                Rt[t][a] = (theta * Rp[t - 1][a] - sigma * Rt[t - 1][a] - (1 - epsilon_s) * rho1 * Rt[t - 1][
                    a] - kappa_t * mu_t * Rt[t - 1][a]
                            - (Rt[t - 1][a] - Rt[t - 1][a - 1]) / (deltaA)) * deltaT + Rt[t - 1][a]
        return V

    # condicoes iniciais

    T0 = 1.3 * 10 ** 5

    I0_t0 = beta * T0 * V0

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
        I[0][ageCont] = (beta * T[0] * V[0] * np.exp(-deltaA * ageCont)) * deltaA + I[0][ageCont - 1]

        if ageCont * deltaA < tau:
            rho1 = 0
        else:
            rho1 = (1 - np.exp(-k * ((ageCont * deltaA) - tau))) * rho

        Rn[0][ageCont] = (r * Rp[0][ageCont - 1] - r * Rp[0][ageCont - 1] * (Rn[0][ageCont - 1] / Rmax)
                          - mu_c * Rn[0][ageCont - 1]) * deltaA + Rn[0][ageCont - 1]

        Rp[0][ageCont] = (alpha * Rn[0][ageCont - 1] + sigma * Rt[0][ageCont - 1]
                          - theta * Rp[0][ageCont - 1] - rho1 * Rp[0][ageCont - 1]
                          - mu_c * Rp[0][ageCont - 1]) * deltaA + Rp[0][ageCont - 1]

        Rt[0][ageCont] = (theta * Rp[0][ageCont - 1] - sigma * Rt[0][ageCont - 1]
                          - rho1 * Rt[0][ageCont - 1] - mu_t * Rt[0][ageCont - 1]) * deltaA + Rt[0][ageCont - 1]

    # Fim das condicoess iniciais----------------------------------------------+

    # Executa o metedo numerico para resolver
    # Precisa passar como parametro as condicoes iniciais
    diferencasfinitas(poi)

    # Se algum valor de V for negativo, defina um custo alto para o conjunto de parametros
    for pos in V:
        if not (pos >= 0):
            return 1000000

    # Passa para a base log o resultado

    V_log = np.log10(V)

    t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03] #PATB17

    plt.plot(tempoPt, V_log, 'b--', label = "$ε_s$=0.999, $ε_α$=0.900")

    plt.xlim(0.0, 31.0)


# Copy best result (epsilon_s, epsilon_alpha, epsilon_r, delta, alpha, r, rho)
poi = [0.999, 0.900, 0.300, 0.07, 22.0, 11.1, 12]

#[(0.990, 0.999), (0.900, 0.960), (0.200, 0.400), (0.01, 0.3), (15, 24), (5.0, 11.1), (9.7, 12.4)]  #bounds dos parâmetros

# experimental data

# PATB17
exp = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432, 2.7474, 2.7016, 2.3541, 2.0453, 1.4914]

viralmodelfit(poi, exp, 10 ** exp[0])


################################################


# recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0):
    # Inicializacao

    # parametros
    s = 1.3 * 10 ** 5  # taxa de crescimento das celulas alvo
    d = 0.01  # taxa de morte das celulas alvo
    beta = 5 * 10 ** -8  # taxa de formacao de complexos de replicacao
    c = 19  # taxa de eliminacao do virus pelo sistema imune

    kappa_t = 1.0  # fator para aumentar a degradacao de Rt
    kappa_c = 1.0  # fator para aumentar a degradacao de Rp e Rn no complexo de replicacao

    k = 0.80  # coeficiente da funcao exponencial de atraso na exportacao de RNA positivo
    tau = 0.50  # tempo de atraso para a exportacao de RNA positivo
    Rmax = 50.0  # numero maximo de RNA negativo / complexo de replicacao (Rn)
    sigma = 4.0  # taxa de formacao de complexos de replicacao
    mu_t = 0.8  # 0.8 no artigo # decaimento natural de Rt
    theta = 1.20  # taxa de disponibilidade para traducao
    mu_c = 2.8  # 0.89 no artigo # decaimento natural de Rc e Rn

    # parametros que serao ajustados

    epsilon_s = poi[0]  # efetividade da terapia em diminuir ou bloquear a exportacao de Rp
    epsilon_alpha = poi[1]  # efetividade da terapia em diminuir ou bloquear a replicacao de Rp
    epsilon_r = poi[2]  # efetividade da terapia em diminuir ou bloquear a replicacao de Rn
    delta = poi[3]  # taxa de morte da celula infectada
    alpha = poi[4]  # taxa de replicacao do Rp
    r = poi[5]  # taxa de replicacao do Rn / complexo de replicacao (Rn)
    rho = poi[6]  # taxa de exportacao do Rp

    # variaveis

    I = np.zeros((tempoNpts, ageNpts))
    Rp = np.zeros((tempoNpts, ageNpts))
    Rt = np.zeros((tempoNpts, ageNpts))
    Rn = np.zeros((tempoNpts, ageNpts))
    V = np.zeros(tempoNpts)
    T = np.zeros(tempoNpts)

    def diferencasfinitas(poi):

        # Solve

        for t in range(1, tempoNpts):

            T[t] = (s - d * T[t - 1] - beta * V[t - 1] * T[t - 1]) * deltaT + T[t - 1]

            V[t] = deltaT * ((1 - epsilon_s) * rho * calcIntegral(I[t - 1], Rp[t - 1], Rt[t - 1])
                             - c * V[t - 1]) + V[t - 1]

            rho1 = 0.0

            Rp[t][0] = 0
            Rn[t][0] = 0
            Rt[t][0] = 1
            I[t][0] = beta * V[t] * T[t]

            for a in range(1, ageNpts):

                if (float(a * deltaA) < tau):
                    rho1 = 0
                else:
                    rho1 = rho * (1 - np.exp(-k * (float(a * deltaA) - tau)))

                I[t][a] = (-delta * I[t - 1][a] - (I[t - 1][a] - I[t - 1][a - 1]) / deltaA) * deltaT + I[t - 1][a]

                Rn[t][a] = ((1 - epsilon_r) * r * Rp[t - 1][a] * (1 - (Rn[t - 1][a] / Rmax)) - kappa_c * mu_c *
                            Rn[t - 1][a]
                            - (Rn[t - 1][a] - Rn[t - 1][a - 1]) / (deltaA)) * deltaT + Rn[t - 1][a]

                Rp[t][a] = ((1 - epsilon_alpha) * alpha * Rn[t - 1][a] + sigma * Rt[t - 1][a] - theta * Rp[t - 1][a]
                            - (1 - epsilon_s) * rho1 * Rp[t - 1][a] - kappa_c * mu_c * Rp[t - 1][a]
                            - (Rp[t - 1][a] - Rp[t - 1][a - 1]) / (deltaA)) * deltaT + Rp[t - 1][a]

                Rt[t][a] = (theta * Rp[t - 1][a] - sigma * Rt[t - 1][a] - (1 - epsilon_s) * rho1 * Rt[t - 1][
                    a] - kappa_t * mu_t * Rt[t - 1][a]
                            - (Rt[t - 1][a] - Rt[t - 1][a - 1]) / (deltaA)) * deltaT + Rt[t - 1][a]
        return V

    # condicoes iniciais

    T0 = 1.3 * 10 ** 5

    I0_t0 = beta * T0 * V0

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
        I[0][ageCont] = (beta * T[0] * V[0] * np.exp(-deltaA * ageCont)) * deltaA + I[0][ageCont - 1]

        if ageCont * deltaA < tau:
            rho1 = 0
        else:
            rho1 = (1 - np.exp(-k * ((ageCont * deltaA) - tau))) * rho

        Rn[0][ageCont] = (r * Rp[0][ageCont - 1] - r * Rp[0][ageCont - 1] * (Rn[0][ageCont - 1] / Rmax)
                          - mu_c * Rn[0][ageCont - 1]) * deltaA + Rn[0][ageCont - 1]

        Rp[0][ageCont] = (alpha * Rn[0][ageCont - 1] + sigma * Rt[0][ageCont - 1]
                          - theta * Rp[0][ageCont - 1] - rho1 * Rp[0][ageCont - 1]
                          - mu_c * Rp[0][ageCont - 1]) * deltaA + Rp[0][ageCont - 1]

        Rt[0][ageCont] = (theta * Rp[0][ageCont - 1] - sigma * Rt[0][ageCont - 1]
                          - rho1 * Rt[0][ageCont - 1] - mu_t * Rt[0][ageCont - 1]) * deltaA + Rt[0][ageCont - 1]

    # Fim das condicoess iniciais----------------------------------------------+

    # Executa o metedo numerico para resolver
    # Precisa passar como parametro as condicoes iniciais
    diferencasfinitas(poi)

    # Se algum valor de V for negativo, defina um custo alto para o conjunto de parametros
    for pos in V:
        if not (pos >= 0):
            return 1000000

    # Passa para a base log o resultado

    V_log = np.log10(V)

    t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03] #PATB17

    #plt.figure()

    plt.plot(tempoPt, V_log, 'g-.', label = "$ε_s$=0.900, $ε_α$=0.920")

    plt.xlim(0.0, 31.0)
    plt.legend()
    plt.show()



# Copy best result (epsilon_s, epsilon_alpha, epsilon_r, delta, alpha, r, rho)
poi = [0.900, 0.920, 0.300, 0.07, 22.0, 11.1, 12]

#[(0.990, 0.999), (0.900, 0.960), (0.200, 0.400), (0.01, 0.3), (15, 24), (5.0, 11.1), (9.7, 12.4)]  #bounds dos parâmetros

# experimental data

# PATB17
exp = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432, 2.7474, 2.7016, 2.3541, 2.0453, 1.4914]

viralmodelfit(poi, exp, 10 ** exp[0])