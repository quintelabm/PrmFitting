import uncertainpy as un
import chaospy as cp
import numpy as np
import scipy.integrate as integrate
import time
from tqdm import tqdm
from tqdm import tqdm_gui
import matplotlib.pyplot as plt

from uqtools import *

class MethodOfLines:

    def __init__(self,N,tmax):
        #variaveis de inicializacao
        self.N = N # number of ages
        self.xmin = 0
        self.xmax  = 50 # number of classes of ages (50)
        #self.agePt = np.linspace(0, ageFim, ageNpts)
        #self.ageCont = 0
        self.tmax = tmax # number of days for simulation (31)
        self.dt   = 0.01 # time step
        #self.tempoPt = np.linspace(0, tempoFim, tempoNpts)
        #self.tempoCont = 0
        self.Virus = []
        self.v0 = 6.4799*10**5
        self.Target = []
        self.t0 = 1.3*10**5
        #self.RNAt = 1 
        self.RNAp = 1
        #self.RNAn = 0
        self.delta = 0.07
        self.beta = 5*10**-8
        self.theta = 1.2
        self.s     = 1.3*10**5
        self.d     = 0.01
        self.beta  = 5*10**-8
        self.c     = 19.0
        self.rho   = 12.0
        self.prm_alpha = 22
        self.mu    = 1
        self.initializeDomain()
        self.initializeVariables()
        self.initializeParams()


    def initializeDomain(self):
        self.dx = (self.xmax-self.xmin)/self.N
        self.x = np.arange(self.xmin-self.dx, self.xmax+(2*self.dx), self.dx)
        self.time = np.arange(0,self.tmax, self.dt)


    def Infected(self,x,args):
        v,t=args
        var = v*t*self.beta*np.exp(-self.delta*x)
        return var


    def Target(self):
        return self.s - self.d*self.target - self.v*self.target*self.beta


    def initializeVariables(self):
        # target cell
        self.t0 = self.t0+(self.s - self.d*self.t0 - self.v0*self.t0*self.beta)*self.dt
        self.target = self.t0
        # infected cell
        self.u0 = self.v0*self.t0*self.beta*np.exp(-self.delta*self.x)
        self.u = self.u0.copy()
        self.unp1 = self.u0.copy()
        # rna
        self.rna_0 = self.RNAp*np.exp(self.x)
        self.rna = self.rna_0.copy()
        self.unp2 = self.rna_0.copy()
        # virus 
        res = integrate.quad(self.Infected,self.xmin,self.xmax,[self.v0,self.target])
        self.v0 = self.v0+(self.rho*res[0]-self.v0*self.c)*self.dt
        self.v = self.v0


    def initializeParams(self):
        self.nsteps = round(self.tmax/self.dt)
        self.alpha = self.dt/(2*self.dx)
        self.MEAN_PATS = [6.47991433, 6.42897983, 5.857277, 4.63766183, 4.08108333, 3.38275467, 3.18124267, 2.88121, 2.66053917]#, 2.54078967, 2.487822, 2.21939417, 1.90103167, 1.7158415]
        self.t_exp = [0.04, 0.09, 0.17, 0.34, 0.50, 1.00, 1.50, 2.98, 4.91]#, 6.92, 10.94, 13.96, 20.91, 27.98]


    def solve_and_plot(self):
        tc = 0

        for i in range(self.nsteps):
            plt.clf()

            # The Lax-Friedrichs scheme
            #for j in range(self.N):
            #    self.unp1[j] = self.u[j]-self.alpha*(self.u[j+1]-self.u[j-1])+(1/2)*(self.u[j+1]-2*self.u[j]+self.u[j-1])
            
            # finite differences
            for j in range(self.N):
                self.unp1 =  self.u[j]+(-self.delta*self.u[j]+self.alpha*self.u*(self.u[j]-self.u[j-1]))
                self.unp2 =  self.rna[j]+(self.prm_alpha*self.rna[j]-(self.mu+self.rho)*self.rna[j]+self.alpha*self.rna*(self.rna[j]-self.rna[j-1]))

            # infected cells
            self.u = self.unp1.copy() 

            # rna
            self.rna = self.unp2.copy() 
        
            # virus
            exported_virus = integrate.quad(self.Infected,self.xmin,self.xmax,[self.v,self.target])
            self.v = self.v+(self.dt*(self.rho*exported_virus[0]-self.c*self.v))
            self.Virus.append(self.v)
            # target cells
            self.target = self.target + self.dt*(self.s - self.d*self.target - self.v*self.target*self.beta)
            self.Target.append(self.target)
            
            #plt.plot(self.x,self.u,'bo-', label= "")
            #plt.ylabel('Infected cells')
            #plt.ylim(0.0,8.0)
            plt.plot(self.x,self.rna,'bo-', label= "")
            plt.xlabel('day')
            plt.ylabel('RNA')
            plt.legend(loc='upper left', prop={'size':13})
            #plt.pause(0.01)
            plt.grid(True)
            tc += self.dt
            

        plt.figure()
        plt.plot(self.time,self.Virus,'bo-', label= "")
        plt.plot(self.t_exp, self.MEAN_PATS, 'r', label='data', linewidth=4)
        plt.xlabel('day')
        plt.ylabel('HVC RNA log_{10}')
        plt.yscale('log')
        plt.grid(True)


        plt.figure()
        plt.plot(self.time,self.Target,'bo-', label= "")
        plt.xlabel('day')
        plt.ylabel('Target Cells')
        plt.yscale('log')
        plt.grid(True)

        

    '''
    def calcIntegral(I,Rp,Rt):
        soma = 0.0
        for a in range(0, ageNpts):
            soma = soma + I[a]*Rp[a]*Rt[a]
        return soma/float(ageNpts)
    '''


    #recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
    def viralmodelfit(prms):

        [alpha, r, rho, epsilon_alpha, epsilon_r, delta] = prms
        #Inicializacao
        #parametros
        s     = 1.3*10**5
        d     = 0.01
        beta  = 5*10**-8
        c     = 19.0

        k     = 0.80 # coeficiente da funcao exponencial de atraso na exportacao de RNA positivo
        tau   = 0.50 # tempo de atraso para a exportacao de RNA positivo
        Rmax  = 50.0 # numero maximo de RNA negativo / complexo de replicacao (Rn)
        sigma = 4.0  # taxa de formacao de complexos de replicacao
        mu_t  = 0.8  # 0.8 no artigo # decaimento natural de Rt
        theta = 1.20 # taxa de disponibilidade para traducao
        mu_c  = 2.8  # 0.89 no artigo # decaimento natural de Rc e Rn
        epsilon_s = 0.998

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
        T0 = 1.3*10**5 #PAT83:
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
            Rn[0][ageCont] = (r*Rp[0][ageCont-1] - r*Rp[0][ageCont-1]*(Rn[0][ageCont-1]/Rmax)- mu_c*Rn[0][ageCont-1])*deltaA + Rn[0][ageCont-1]
            Rp[0][ageCont] = (alpha*Rn[0][ageCont-1]+sigma*Rt[0][ageCont-1]-theta*Rp[0][ageCont-1]-rho1*Rp[0][ageCont-1]- mu_c*Rp[0][ageCont-1])*deltaA + Rp[0][ageCont-1]
            Rt[0][ageCont] = (theta*Rp[0][ageCont-1] - sigma*Rt[0][ageCont-1]- rho1*Rt[0][ageCont-1] - mu_t*Rt[0][ageCont-1])*deltaA + Rt[0][ageCont-1]
        #Fim das condicoess iniciais----------------------------------------------+
        
        # Executa o metodo numerico para resolver
        # Precisa passar como parametro as condicoes iniciais
        V = diferencasfinitas()

        # Passa para a base log o resultado
        V_log = np.log10(V)

        return tempoPt, V_log

    ##############FIM MODELO################

    '''

    model = un.Model(
        run = viralmodelfit,
        labels=["Tempo (dias)",
            "Carga viral (log10)"]
    )

    epsilon_r = 0.3
    epsilon_alpha = 0.9
    delta = 0.07
    alpha = 22
    r = 11.1
    rho = 12.0

    # create pdf distributions
    min_bound = 0.9
    max_bound = 1.1

    pdf_alpha = cp.Uniform(alpha*min_bound,alpha*max_bound)
    pdf_r = cp.Uniform(r*min_bound,r*max_bound)
    pdf_rho = cp.Uniform(rho*min_bound,rho*max_bound)
    pdf_epsilon_alpha = cp.Uniform(epsilon_alpha*min_bound,epsilon_alpha*max_bound)
    pdf_epsilon_r = cp.Uniform(epsilon_r*min_bound,epsilon_r*max_bound)
    pdf_delta = cp.Uniform(delta*min_bound,delta*max_bound)

    #print('Criando distribuicoes normais')
    ##epsilon_r = cp.Normal(0.3,0.1)
    #epsilon_r = cp.Normal(0.3,0.001)
    ##epsilon_s = cp.Normal(0.998,0.001)
    ##epsilon_alpha = cp.Normal(0.92,0.002)
    #epsilon_alpha = cp.Normal(0.92,0.001)
    #delta = cp.Normal(0.07,0.01)
    #alpha = cp.Normal(22,0.001)
    #r = cp.Normal(11.1,0.1)
    #rho = cp.Normal(12.0,1.0)

    print('Distribuicoes criadas')

    dist = cp.J(pdf_alpha, pdf_r, pdf_rho, pdf_epsilon_alpha, pdf_epsilon_r)

    dist = cp.J(dist, pdf_delta)

    npar = len(dist)

    #grau da aproximacao polinomial
    degpc = 2
    ns = 10#3*P(degpc,npar)
    print("number of input parameters %d" % npar)
    print("polynomial degree %d" % degpc)
    print("min number of samples %d" % ns)


    evals_virus = []

    samples = dist.sample(ns,"L")
    k = 0
    lpts = range(ns)
    samp = samples.T
    print("evaluating samples: ")
    for i in tqdm(lpts,bar_format='{l_bar}{bar:20}{r_bar}{bar:-20b}'):
    #for s in samples.T:
        s = samp[k]
        k = k+1
        [sim_time,virus] = viralmodelfit(s)
        evals_virus.append(virus)

    evals_virus = np.array(evals_virus)

    plt.style.use('../estilo/PlotStyle.mplstyle')
    plt.close('all')

    output_path = './output_monte_carlo/'
    ext = '.png'
    itvl = 6
    caso = 'caso.'

    plt.figure('HCV RNA')
    #dadosViremiaLog10.plot.scatter(x='Day',y='Viral_load',color='m',label='Dados experimentais')
    plt.plot(t_exp, MEAN_PATS, 'o', label='data', linewidth=4)
    plt.xlim(0.0,tempoFim)
    #plt.ylim(0.0,8.0)
    plt.xlabel('day')
    plt.ylabel('HVC RNA log_{10}')
    plt.legend(loc='upper left', prop={'size':13})
    plt.grid(True)
    plt.title("Solution Virus")
    mean_virus= plot_confidence_interval(plt, sim_time, evals_virus, 'red', 'Modelo: Virus')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path+caso+'_output_Virus'+ext)
    plt.show()
    '''


def main():
    #sim = LaxFriedrichs(100, 1.5)
    sim = MethodOfLines(50,5)
    sim.solve_and_plot()
    plt.show()
    
    
if __name__ == "__main__":
    main()


# alpha, r, rho, epsilon_alpha, epsilon_r, delta


'''

# define parameter dictionary
parameters = {"alpha": pdf_alpha,
        "r": pdf_r,
        "rho": pdf_rho,
        "epsilon_alpha": pdf_epsilon_alpha,
        "epsilon_r": pdf_epsilon_r,
        "delta": pdf_delta
            }



# set up UQ
UQ = un.UncertaintyQuantification(
    model=model,
    parameters=parameters
)

print('Inicio UQ --- ')
begin = time.perf_counter()
data = UQ.monte_carlo(nr_samples=1000)
end = time.perf_counter()
print(' --- Fim UQ')
tempo = end-begin
with open('tempo.txt', 'w') as file:
    file.write("Tempo de execução: "+str(tempo)+" segundos")

'''