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


#--- MAIN ---------------------------------------------------------------------+

def main(cost_func, bounds, popsize, mutate, recombination, maxiter, POI, PAT):

    #--- INITIALIZE A POPULATION (step #1) ----------------+
    
    population = []
    for i in range(0,popsize):
        indv = []
        for j in range(len(bounds)):
            indv.append(random.uniform(bounds[j][0],bounds[j][1]))
        population.append(indv)
            
    #--- SOLVE --------------------------------------------+

    # cycle through each generation (step #2)
    for i in range(1,maxiter+1):
        print('GENERATION:',i)

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

            score_trial  = cost_func(v_trial, PAT)
            score_target = cost_func(x_t, PAT)

            if score_trial < score_target:
                population[j] = v_trial
                gen_scores.append(score_trial)
                print('   >',score_trial, v_trial)

            else:
                print('   >',score_target, x_t)
                gen_scores.append(score_target)

        #--- SCORE KEEPING --------------------------------+

        gen_avg = sum(gen_scores) / popsize                         # current generation avg. fitness
        gen_best = min(gen_scores)                                  # fitness of best individual
        gen_sol = population[gen_scores.index(min(gen_scores))]     # solution of best individual

        print('      > GENERATION AVERAGE:',gen_avg)
        print('      > GENERATION BEST:',gen_best)
        print('         > BEST SOLUTION:',gen_sol,'\n')

    return gen_sol

#--- CONSTANTS ----------------------------------------------------------------+

cost_func = viralmodelfit                   # Cost function
bounds = [(0.1,0.6),(5.0,15.0),(15.0,23.0)]  # Bounds [(x1_min, x1_max), (x2_min, x2_max),...]
popsize = 10                                # Population size, must be >= 4
mutate = 0.5                                # Mutation factor [0,2]
recombination = 0.7                         # Recombination rate [0,1]
maxiter = 200                               # Max number of generations (maxiter)

# --- for patient 1 virus no semen
PAT8 = [5.64, 5.31, 4.23, 3.36, 3.14, 2.86, 2.75, 2.50, 2.32, 1.56]
PAT = PAT8

# initial guess for parameters
# delta p c
POI = [0.1, 5.0, 22.3]

#--- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    main(cost_func, bounds, popsize, mutate, recombination, maxiter, POI, PAT)

#--- END ----------------------------------------------------------------------+

