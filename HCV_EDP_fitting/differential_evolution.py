# Tutorial - https://nathanrooy.github.io/posts/2017-08-27/simple-differential-evolution-with-python/
#
#1) Initialize a random population of individuals throughout the search space.
#2) while iter <= max num of generations
#    3) cycle through each individual in the population   
#        3.A) perform mutation            
#        3.B) perform recombination ("crossover" in GA lingo)            
#        3.C) perform selection           
#    4) if stopping criterion has been met:
#            exit and return best individual            
#        else:
#            iter = iter + 1
#            go back to step #3

#------------------------------------------------------------------------------+
#
#   Nathan A. Rooy
#   A simple, bare bones, implementation of differential evolution with Python
#   August, 2017
#
#------------------------------------------------------------------------------+

#--- IMPORT DEPENDENCIES ------------------------------------------------------+

import random
import numpy as np
from depend.cost import viralmodelfit
from depend.bounds import ensure_bounds
import matplotlib.pyplot as plt

#--- MAIN ---------------------------------------------------------------------+

def main(cost_func, bounds, popsize, mutate, recombination, maxiter, PAT, V0):

    #--- INITIALIZE A POPULATION (step #1) ----------------+
    
    population = []
    for i in range(0,popsize):
        indv = []
        for j in range(len(bounds)):
            indv.append(random.uniform(bounds[j][0],bounds[j][1]))
        population.append(indv)
    #melhores solucoes para cada paciente 
    best_sol = []

    #--- SOLVE --------------------------------------------+

    # cycle through each generation (step #2)
    for i in range(1,maxiter+1):
        print('GENERATION:',i)
        saida.writelines('GENERATION:'+str(i))
        gen_scores = [] # score keeping

        # cycle through each individual in the population
        for j in range(0, popsize):

            #--- MUTATION (step #3.A) ---------------------+
            
            # select three random vector index positions [0, popsize), not including current vector (j)
            candidates = list(range(0,popsize))
            #candidates.remove(j)
            del candidates[j]
            random_index = random.sample(candidates, 3)
            
            x_1 = population[random_index[0]]
            x_2 = population[random_index[1]]
            x_3 = population[random_index[2]]

            x_t = population[j]     # target individual

            # subtract x3 from x2, and create a new vector (x_diff)
            x_diff = [x_2_i - x_3_i for x_2_i, x_3_i in zip(x_2, x_3)]

            # multiply x_diff by the mutation factor (F) and add to x_1
            v_donor = [x_1_i + mutate * x_diff_i for x_1_i, x_diff_i in zip(x_1, x_diff)]
            v_donor = ensure_bounds(v_donor, bounds)

            #--- RECOMBINATION (step #3.B) ----------------+

            v_trial = []
            for k in range(len(x_t)):
                crossover = random.random()
                if crossover <= recombination:
                    v_trial.append(v_donor[k])

                else:
                    v_trial.append(x_t[k])
                    
            #--- GREEDY SELECTION (step #3.C) -------------+

            score_trial  = cost_func(v_trial, PAT, V0)
            score_target = cost_func(x_t, PAT, V0)

            if score_trial < score_target:
                population[j] = v_trial
                gen_scores.append(score_trial)
                print('   >',score_trial, v_trial)
                #Pegando todas as solucoes nas 10% ultimas geracoes
                if i>(maxiter+1)*0.9:
                    best_sol.append([score_trial, v_trial])

            else:
                print('   >',score_target, x_t)
                gen_scores.append(score_target)
                #Pegando todas as solucoes nas 10% ultimas geracoes
                if i>(maxiter+1)*0.9:
                    best_sol.append([score_target, x_t])

        #--- SCORE KEEPING --------------------------------+

        gen_avg = sum(gen_scores) / popsize                         # current generation avg. fitness
        gen_best = min(gen_scores)                                  # fitness of best individual
        gen_sol = population[gen_scores.index(min(gen_scores))]     # solution of best individual

        print('      > GENERATION AVERAGE:',gen_avg)
        print('      > GENERATION BEST:',gen_best)
        print('         > BEST SOLUTION:',gen_sol,'\n')
        saida.writelines('> BEST SOLUTION:'+str(gen_sol)+'\n')
    #retorna a melhor solucao e um vetor com todas as 10 % ultimas solucoes 
    return gen_sol, best_sol

#--- CONSTANTS ----------------------------------------------------------------+

cost_func = viralmodelfit                                  # Cost function
bounds = [(0.09,0.99),(0.09,0.99),(0.09,0.99),(0.01,1.2),(25,35),(0.9,2.0),(6,12)]  # Bounds [(x1_min, x1_max), (x2_min, x2_max),...]
popsize = 10                                               # Population size, must be >= 4
mutate = 0.5                                               # Mutation factor [0,2]
recombination = 0.7                                        # Recombination rate [0,1]
maxiter = 10                                               # Max number of generations (maxiter)

#Vetor com todos os pacientes
patients = [ ]

PAT8 = [5.64, 5.31, 4.23, 3.36, 3.14, 2.86, 2.75, 2.50, 2.32, 1.56]
PAT42 = [5.65, 5.00, 3.98, 3.84, 2.94, 2.82, 2.87, 2.53, 2.31, 2.61]
PAT68 =  [7.15, 7.02, 6.19, 5.50, 4.96, 4.29, 4.11, 3.75, 3.68, 3.35]
PAT69 = [6.14, 5.87, 4.73, 4.17, 3.55, 3.14, 2.87, 2.60, 2.55, 2.58]
PAT83 = [5.45, 5.38, 4.73, 4.00, 3.39, 2.89, 2.68, 2.72, 2.97, 1.93]


patients.append(PAT8)
patients.append(PAT42)
patients.append(PAT68)
patients.append(PAT69)
patients.append(PAT83)

t_exp = [0, 0.083, 0.167, 0.25, 0.333, 0.5, 0.667, 1, 1.5, 2 ]
    
#--- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    saida = open("./docs/relatorio.txt", "a")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: '+ str(popsize)+ '\nNumber of generations: '+ str(maxiter)+ '\n')

    best_solves = []
    pat_cont = 1
    for pat in patients:
        print(pat_cont, " Patient")
        saida.writelines(str(pat_cont) + " Patient\n\n")

        sol_pat, vet_best = main(cost_func, bounds, popsize, mutate, recombination, maxiter, pat, (10**pat[0]))
        best_solves.append(sol_pat)        
        
        #Ordena as solucoes em ordem crecente baseada no custo
        vet_best.sort()
        #Seleciona os 10 melhores pares [custo, [parametros]]
        vet_best = vet_best[:10]
        #Transfere para outro vetor apenas os valores dos parametros
        param_best = []
        for p in vet_best:
            param_best.append(p[1])

        average_param = np.zeros(len(bounds))
        for p in param_best:
            average_param[0] = average_param[0] + p[0]
            average_param[1] = average_param[1] + p[1]
            average_param[2] = average_param[2] + p[2]
            average_param[3] = average_param[3] + p[3]
            average_param[4] = average_param[4] + p[4]
            average_param[5] = average_param[5] + p[5]
            average_param[6] = average_param[6] + p[6]

        #Esse vetor tem a media dos 10 melhores valores dos parametros
        average_param = average_param/10
        plt.clf()
        #Plot experimental        
        plt.plot(t_exp,pat, 'ro')
        
        #Plot da solucao com a media dos parametros
        cost_func(average_param, pat, (10**pat[0]))

        plt.savefig("figura"+str(pat_cont)+".png")
        pat_cont = pat_cont + 1

    #Calculando as medias dos parametros para todos os pacientes
    parameters = zip(best_solves[0],best_solves[1],best_solves[2],best_solves[3],best_solves[4])
    parameters = list(parameters)
    
    average = []
    for param in parameters:
        average.append(sum(param)/len(patients))
    print(average)
    saida.writelines('\n\nNa media, os valores dos parametros sao: ' + str(average) )
    saida.close()
#--- END ----------------------------------------------------------------------+

