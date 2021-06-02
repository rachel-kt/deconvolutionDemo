# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 15:39:22 2021

@author: mdouaihy
"""

import pygadForDeconvolution
import numpy as np
from sumSignalDroso import sumSignal1_par

last_fitness=0
def optimize_ga1_par(TolFun,target,Nbr_poly_estimate,num_possible_poly,FreqEchSimu,
                     FreqEchImg, TaillePreMarq, TailleSeqMarq, TaillePostMarq, Polym_speed, frame_num, Intensity_for_1_Polym, generations):

    Nbr_simu_DNA = 400 # number of "chromosome" in GA
    Pattern_polys = np.zeros((Nbr_simu_DNA,num_possible_poly))  
    
    # setting the initial population
    for i in range(Nbr_simu_DNA):
        l=np.random.choice(num_possible_poly, int(Nbr_poly_estimate), replace=False)
       # l = random.sample(list(range(num_possible_poly)),Nbr_poly_estimate)
        Pattern_polys[i,l] = 1 # randomly choose poly position

    def fitness_fun(solution, solution_idx):

        find_ind = np.where(solution==1)[0]
        output = sumSignal1_par(find_ind,FreqEchSimu, FreqEchImg, TaillePreMarq, TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)
        fitness = -np.sum((output-target)**2)
    
        return fitness
    
    def callback_generation(ga_instance):
        print("Generation = {generation}".format(generation=ga_instance.generations_completed))
        print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution()[1]))
    
        
    
    ga_instance = pygadForDeconvolution.GA(TolFun=TolFun,
                                num_generations=generations,
                           num_parents_mating=2,
                           fitness_func=fitness_fun,
                           sol_per_pop=Nbr_simu_DNA,
                           num_genes=num_possible_poly,
                           init_range_low=0,
                           init_range_high=2,
                           gene_type= int,
                           mutation_type="random",
                           mutation_by_replacement=True,
                           mutation_percent_genes=0.03,
                           random_mutation_min_val=0.0,
                           random_mutation_max_val=2.0,
                           parent_selection_type = 'rank',
                           initial_population = Pattern_polys,
                           crossover_type='scattered',
                           crossover_probability = 0.8)
                        #   callback_generation=callback_generation)
                      #     save_best_solutions=True)

    ga_instance.run()
    # After the generations complete, some plots are showed that summarize the how the outputs/fitenss values evolve over generations.
#    ga_instance.plot_result()
    
    # Returning the details of the best solution.
   # solution, solution_fitness, solution_idx = ga_instance.best_solution()
    #print("Parameters of the best solution : {solution}".format(solution=solution))
    #print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
    #print("Index of the best solution : {solution_idx}".format(solution_idx=solution_idx))
    solution, solution_fitness, solution_idx = ga_instance.best_solution()
    pattern_fit = solution
 


    return pattern_fit, ga_instance.best_solution_generation, ga_instance.best_solutions_fitness
    
