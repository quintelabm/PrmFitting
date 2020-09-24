from SALib.sample import saltelli

from SALib.analyze import sobol

import matplotlib.pyplot as plt

import numpy as np
import seaborn as sns

sns.set()

#variaveis de inicializacao
ageFim  = 50
deltaA   = 0.1 #passo no age
ageNpts = int(ageFim/deltaA ) + 1
agePt = np.linspace(0, ageFim, ageNpts)
ageCont = 0

tempoFim = 30
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
def viralmodelfit(param):
    
    s, d, beta, c, k, tau, Rmax, sigma, mu_t, theta, mu_c, kappa_t, kappa_c, alpha, r, rho, epsilon_s, epsilon_alpha, epsilon_r, delta = param
    #variaveis 

    I = np.zeros((tempoNpts, ageNpts))
    Rp = np.zeros((tempoNpts, ageNpts))
    Rt = np.zeros((tempoNpts, ageNpts))
    Rn = np.zeros((tempoNpts, ageNpts))
    V = np.zeros(tempoNpts)
    T = np.zeros(tempoNpts)

    def diferencasfinitas():
        
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
    #PAT83:
    V0 = 2.563292*10**6
    
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
        
        Rn[0][ageCont] = (r*Rp[0][ageCont-1] - r*Rp[0][ageCont-1]*(Rn[0][ageCont-1]/Rmax) - mu_c*Rn[0][ageCont-1])*deltaA + Rn[0][ageCont-1]
    
        Rp[0][ageCont] = (alpha*Rn[0][ageCont-1] + sigma*Rt[0][ageCont-1] - theta*Rp[0][ageCont-1]- rho1*Rp[0][ageCont-1] - mu_c*Rp[0][ageCont-1])*deltaA + Rp[0][ageCont-1]
    
        Rt[0][ageCont] = (theta*Rp[0][ageCont-1] - sigma*Rt[0][ageCont-1] - rho1*Rt[0][ageCont-1] - mu_t*Rt[0][ageCont-1])*deltaA + Rt[0][ageCont-1]
        
    #Fim das condicoess iniciais----------------------------------------------+
    
    # Executa o metedo numerico para resolver
    # Precisa passar como parametro as condicoes iniciais 
    diferencasfinitas()
    

    # Passa para a base log o resultado
    V_log = np.log10(V)

    for v in V_log:
        if v < 0:
            return np.zeros(tempoNpts)

    return V_log

##############FIM MODELO################
    
s = 1.3*10**5
d     = 0.01
beta  = 5*10**-8
c     = 19.0
k     = 0.80 
tau   = 0.50 
Rmax  = 50.0 
sigma = 4.0 
mu_t  = 0.8 
theta = 1.20 
mu_c  = 2.8
kappa_t       = 1.0
kappa_c       = 1.0   
alpha         = 22.0
r             = 11.1
rho           = 12.0
epsilon_s     = 0.999
epsilon_alpha = 0.920     
epsilon_r     = 0.300
delta         = 0.07

parameters = { 'num_vars': 20,
              'names': ['s', 'd', 'beta', 'c', 'k', 'tau', 'Rmax', 'sigma', 'mu_t', 'theta', 'mu_c', 'kappa_t', 'kappa_c', 'alpha', 'r', 'rho', 'epsilon_s', 'epsilon_alpha', 'epsilon_r', 'delta'],
              'bounds': [[s*0.9, s*1.1],
                         [d*0.9, d*1.1],
                         [beta*0.9, beta*1.1],
                         [c*0.9, c*1.1],
                         [k*0.9, k*1.1],
                         [tau*0.9, tau*1.1],
                         [Rmax*0.9, Rmax*1.1],
                         [sigma*0.9, sigma*1.1],
                         [mu_t*0.9, mu_t*1.1],
                         [theta*0.9, theta*1.1],
                         [mu_c*0.9, mu_c*1.1],
                         [kappa_t*0.9, kappa_t*1.1],
                         [kappa_c*0.9, kappa_c*1.1],
                         [alpha*0.9, alpha*1.1],
                         [r*0.9, r*1.1],
                         [rho*0.9, rho*1.1],
                         [epsilon_s*0.9, epsilon_s],
                         [epsilon_alpha*0.9, epsilon_alpha*1.08],
                         [epsilon_r*0.9, epsilon_r*1.1],
                         [delta*0.9, delta*1.1],
                         ]
            }


#Gera parametros aleatorios
n = 1
param_values = saltelli.sample(parameters, n, calc_second_order=False )

#Y = np.zeros([param_values.shape[0]])
ww = n*(22)
solucoes = np.zeros([ww,tempoNpts])


for i, X in enumerate(param_values):
    solucoes[i] = viralmodelfit(X)

solucoes = solucoes.T
Si = np.zeros(tempoNpts)

Sres = np.zeros([tempoNpts,20])
f = open("arqsaida.txt","w")

for i in range(1, tempoNpts):
    Si = sobol.analyze(parameters, solucoes[i], calc_second_order=False)
    Sres[i] = Si['S1']
    f.write(tempoPt[i]+ " "+ str(list(Sres[i])))

plt.xlabel("tempo(dias)")
plt.ylabel("indice")
plt.figure(0)
plt.plot(tempoPt, Sres[:,0])
plt.savefig("/salibEDO/s.png")
plt.figure(1)
plt.plot(tempoPt, Sres[:,1])
plt.savefig("/salibEDO/d.png")
plt.figure(2)
plt.plot(tempoPt, Sres[:,2])
plt.savefig("/salibEDO/beta.png")
plt.figure(3)
plt.plot(tempoPt, Sres[:,3])
plt.savefig("/salibEDO/c.png")
plt.figure(4)
plt.plot(tempoPt, Sres[:,4])
plt.savefig("/salibEDO/k.png")
plt.figure(5)
plt.plot(tempoPt, Sres[:,5])
plt.savefig("/salibEDO/tau.png")
plt.figure(6)
plt.plot(tempoPt, Sres[:,6])
plt.savefig("/salibEDO/Rmax.png")
plt.figure(7)
plt.plot(tempoPt, Sres[:,7])
plt.savefig("/salibEDO/sigma.png")
plt.figure(8)
plt.plot(tempoPt, Sres[:,8])
plt.savefig("/salibEDO/mu_t.png")
plt.figure(9)
plt.plot(tempoPt, Sres[:,9])
plt.savefig("/salibEDO/theta.png")
plt.figure(10)
plt.plot(tempoPt, Sres[:,10])
plt.savefig("/salibEDO/mu_c.png")
plt.figure(11)
plt.plot(tempoPt, Sres[:,11])
plt.savefig("/salibEDO/kappa_t.png")
plt.figure(12)
plt.plot(tempoPt, Sres[:,12])
plt.savefig("/salibEDO/kappa_c.png")
plt.figure(13)
plt.plot(tempoPt, Sres[:,13])
plt.savefig("/salibEDO/alpha.png")
plt.figure(14)
plt.plot(tempoPt, Sres[:,14])
plt.savefig("/salibEDO/r.png")
plt.figure(15)
plt.plot(tempoPt, Sres[:,15])
plt.savefig("/salibEDO/rho.png")
plt.figure(16)
plt.plot(tempoPt, Sres[:,16])
plt.savefig("/salibEDO/epsilon_s.png")
plt.figure(17)
plt.plot(tempoPt, Sres[:,17])
plt.savefig("/salibEDO/epsilon_alpha.png")
plt.figure(18)
plt.plot(tempoPt, Sres[:,18])
plt.savefig("/salibEDO/epsilon_r.png")
plt.figure(19)
plt.plot(tempoPt, Sres[:,19])
plt.savefig("/salibEDO/delta.png")

f.close()
plt.show()
