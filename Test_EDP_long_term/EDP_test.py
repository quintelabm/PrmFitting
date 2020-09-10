import numpy as np 
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.interpolate import InterpolatedUnivariateSpline

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
    sigma   = 4     # taxa de formacao de complexos de replicacao
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

    #t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89]#, 6.88, 8.85, 13.85, 20.86, 27.96] #PATB06
    #t_exp = [0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.91]#, 6.92, 10.94, 13.96, 20.91, 27.98]PATB07
    #t_exp = [0.04, 0.08, 0.17, 0.34, 0.51, 1.00, 1.50, 2.91, 4.92]#, 6.91, 10.87, 13.86, 20.87, 25.92] #PATB08
    #t_exp = [0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90]#, 6.94, 10.96, 14.95, 20.96, 28.03] #PATB09
    #t_exp = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01]#, 6.01, 7.00, 9.99, 14.00, 20.98, 30.01] #PATB16
    t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03] #PATB17


    #t_exp = [0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 2.86]#, 6.91, 8.86, 9.88, 13.90, 21.87, 30.87] #PATC05
    #t_exp = [0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94]#, 5.99, 7.01, 9.91, 15.00, 22.00, 28.02] #PATC06
    #t_exp = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94]#, 7.94, 9.95, 14.94, 24.03, 30.94] #PATC09
    #t_exp = [0.04, 0.09, 0.18, 0.33, 0.50, 1.00, 1.50, 2.96, 3.98]#, 7.00, 9.98, 13.97, 2196, 30.02] #PATC10

    plt.figure()    
    #plt.plot(t_exp, exp, 'r--')
    plt.plot(tempoPt,V_log, 'r-')   #
    #plt.legend(["Modelo EDP"]) #"Exp data PATB17",
    #plt.savefig("PATB17.png")



    ##########################################
    ##########################################

    # --- for patient B06
    PATB06 = [6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378,
              2.3404, 2.3345, 2.2355, 2.0492, 2.1173]

    t_exp1 = [0.04, 0.08, 0.17, 0.33, 0.50, 0.92, 1.50, 2.95, 4.89, 6.88, 8.85, 13.85, 20.86, 27.96]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient B07
    PATB07 = [4.2227, 4.0733, 3.0599, 2.2122, 1.7924, 1.6435, 1.7160, 1.6435, 1.2304,
              1.0792, 1.0000, 1.0000, 1.0000, 1.0000]

    t_exp2 = [0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.91, 6.92, 10.94, 13.96, 20.91, 27.98]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient B08
    PATB08 = [5.6546, 5.6486, 4.6203, 3.6876, 3.1526, 2.6355, 2.5302, 2.6294, 2.2788,
              2.0719, 1.9031, 1.7924, 1.6335, 1.5682]

    t_exp3 = [0.04, 0.08, 0.17, 0.34, 0.51, 1.00, 1.50, 2.91, 4.92, 6.91, 10.87, 13.86, 20.87, 25.92]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient B09
    PATB09 = [6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385,
              4.1949, 3.0931, 2.1790, 1.5682, 1.5051]

    t_exp4 = [0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.90, 6.94, 10.96, 14.95, 20.96, 28.03]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient B16
    PATB16 = [6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120,
              3.0302, 2.7528, 2.3838, 2.1818, 1.9243]

    t_exp5 = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01, 7.00, 9.99, 14.00, 20.98, 30.01]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient B17
    PATB17 = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432,
              2.7474, 2.7016, 2.3541, 2.0453, 1.4914]

    t_exp6 = [0.04, 0.08, 0.17, 0.33, 0.50, 1.01, 1.50, 2.99, 6.03, 7.02, 9.98, 14.04, 21.03, 30.03]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient C05
    PATC05 = [6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339,
              2.609594, 2.527630, 2.743510, 2.694605, 2.227887, 1.863323]

    t_exp7 = [0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 2.86, 6.91, 8.86, 9.88, 13.90, 21.87, 30.87]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient C06
    PATC06 = [6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317, 2.484300,
              2.509203, 2.369216, 1.949390, 1.623249, 1.556303]

    t_exp8 = [0.04, 0.08, 0.17, 0.33, 0.48, 1.00, 1.50, 2.94, 5.99, 7.01, 9.91, 15.00, 22.00, 28.02]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient C09
    PATC09 = [6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341,
              2.089905, 2.025306, 1.698970, 1.278754, 1.342423]

    t_exp9 = [0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 2.96, 3.94, 7.94, 9.95, 14.94, 24.03, 30.94]

    # --------------------------------------------------------------------------------------------------------

    # --- for patient C10
    PATC10 = [5.583842, 5.455846, 4.754914, 3.447468, 2.749736, 2.311754, 2.158362, 1.944483, 1.707570,
              1.322219, 1.322219, 1.000000, 1.000000, 1.000000]

    t_exp10 = [0.04, 0.09, 0.18, 0.33, 0.50, 1.00, 1.50, 2.96, 3.98, 7.00, 9.98, 13.97, 21.96, 30.02]

    plt.plot(t_exp1, PATB06, 'o')
    plt.plot(t_exp5, PATB16, 'o')
    plt.plot(t_exp6, PATB17, 'o')

    plt.plot(t_exp7, PATC05, 'o')
    plt.plot(t_exp8, PATC06, 'o')
    plt.plot(t_exp9, PATC09, 'o')

    # plt.plot(t_exp2, PATB07)
    # plt.plot(t_exp3, PATB08)
    # plt.plot(t_exp10, PATC10)

    # plt.plot(t_exp4, PATB09)

    plt.xlim(0.0, 31.0)
    plt.legend(["Modelo EDP", "PATB06", "PATB16", "PATB17", "PATC05", "PATC06", "PATC09"])
    #plt.show()

    ##########################################
    ##########################################
    


    plt.show()
    
    # fazer interpolacao usando os pontos experimentais
    ius = InterpolatedUnivariateSpline(t_exp, exp)
    
    # aplicar a funcao interpolada nos pontos do metodo numerico do modelo
    yi = ius(tempoPt)
    
    # calcular a distancia euclidiana entre os resultados experimentais interpolados
    # e o resultado do modelo
    
    dst = distance.euclidean(V_log, yi)
    
    return dst


# Copy best result (epsilon_s, epsilon_alpha, epsilon_r, delta, alpha, r, rho)
poi = [0.999, 0.920, 0.300, 0.1, 22.0, 11.1, 12] #[0.999, 0.93753882, 0.2, 0.1, 19.93715449, 9.29497553, 9.71686742]

#[(0.990, 0.999), (0.900, 0.940), (0.200, 0.400), (0.01, 0.3), (15, 24), (5.0, 11.1), (9.7, 12.4)]  #bounds dos parâmetros

# experimental data
# PATB06
#exp = [6.3780, 6.4109, 5.6277, 4.4948, 3.9268, 3.1973, 2.8537, 2.5340, 2.4378]#, 2.3404, 2.3345, 2.2355, 2.0492, 2.1173]
# PATB07
#exp = [4.2227, 4.0733, 3.0599, 2.2122, 1.7924, 1.6435, 1.7160, 1.6435, 1.2304]#, 1.0792, 1.0000, 1.0000, 1.0000, 1.0000]
# PATB08
# exp = [5.6546, 5.6486, 4.6203, 3.6876, 3.1526, 2.6355, 2.5302, 2.6294, 2.2788]#, 2.0719, 1.9031, 1.7924, 1.6335, 1.5682]
# PATB09
#exp = [6.2734, 6.0103, 5.9685, 5.6755, 5.5739, 5.4454, 4.9222, 5.1703, 4.8385]#, 4.1949, 3.0931, 2.1790, 1.5682, 1.5051]
# PATB16
#exp = [6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529]#, 3.0120, 3.0302, 2.7528, 2.3838, 2.1818, 1.9243]
# PATB17
exp = [6.4885, 6.5486, 5.6096, 4.4266, 3.8878, 3.2076, 3.1626, 2.9128, 2.8432, 2.7474, 2.7016, 2.3541, 2.0453, 1.4914]


# PATC05
#exp = [6.394490, 6.254302, 5.833741, 4.727484, 3.889414, 3.261501, 3.182415, 2.990339]#, 2.609594, 2.527630, 2.743510, 2.694605, 2.227887, 1.863323]
# PATC06
#exp = [6.839431, 6.735814, 6.452393, 5.340039, 4.581198, 3.481012, 3.164055, 2.780317]#, 2.484300, 2.509203, 2.369216, 1.949390, 1.623249, 1.556303]
# PATC09
#exp = [6.424965, 6.303063, 6.028728, 4.642148, 4.349588, 3.484015, 3.243286, 2.816904, 2.576341]#, 2.089905, 2.025306, 1.698970, 1.278754, 1.342423]
# PATC10
#exp = [5.583842, 5.455846, 4.754914, 3.447468, 2.749736, 2.311754, 2.158362, 1.944483, 1.707570]#, 1.322219, 1.322219, 1.000000, 1.000000, 1.000000]

viralmodelfit(poi, exp, 10**exp[0])

