import random 
from random import Random
from time import time
import inspyred

import os
import csv
import datetime
from robot_3dof import *
import my_tools

K_MAX = 10000

def ga(prng=None, display=False, **kwargs):
    robot = kwargs['robot']
    runs = kwargs['runs']
    
    filename = kwargs['filename']
    fileheaders = kwargs['fileheaders']

    n_elements = len(fileheaders)
    robot_base = robot.get_r()

    points = kwargs['points']
    parameters = {}
    #parameters['points'] = kwargs['points']
    parameters['n_elements'] = n_elements
    parameters['robot'] = robot
    parameters['robot_base'] = robot_base


    #robot = Robot.init_zero()
    #robot.set_r(1)
    #points = my_tools.make_cube_vertices((0,0,0.8), 1, .5)
    
    if prng is None:
        prng = Random()
        prng.seed(time()) 
    ea = inspyred.ec.EvolutionaryComputation(prng)
    ea.observer = inspyred.ec.observers.file_observer
    #ea.observer = my_observer
    ea.selector = inspyred.ec.selectors.tournament_selection
    ea.variator = [inspyred.ec.variators.blend_crossover,
                   inspyred.ec.variators.gaussian_mutation]
    ea.replacer = inspyred.ec.replacers.steady_state_replacement
    ea.terminator = inspyred.ec.terminators.generation_termination
    #ea.terminator = inspyred.ec.terminators.average_fitness_termination
    projdir = os.path.dirname(os.getcwd())

    #stat_file_name = '{0}/robot_ea_statistics.csv'.format(projdir)
    #final_time_name = filename + '_time.csv'
    #final_time = open(final_time_name, 'w')
    #wr_time = csv.writer(final_time, quoting=csv.QUOTE_ALL)
    #wr_time.writerow(["time_taken"])

    final_file_name = filename + '_solutions.csv'
    final_file = open(final_file_name, 'w')
    wr_2 = csv.writer(final_file, quoting=csv.QUOTE_ALL)
    #wr_2.writerow(["a","b","d","e","c-r"])
    wr_2.writerow(fileheaders + ["fitness"] + ["time_taken"])
    
    for i in range(0, runs):
      started = datetime.datetime.now()
      stat_file_name = filename + 'stats_run%d.csv' %(i+1)
      ind_file_name = filename + 'inds_run%d.csv' %(i+1)
      stat_file = open(stat_file_name, 'w')
      #generation number, population size, worst, best, median, average, standard deviation
      ind_file = open(ind_file_name, 'w')
      #generation number, individual number, fitness, string representation of candidate    
      final_pop = ea.evolve(generator = generate_config,
                            evaluator = inspyred.ec.evaluators.parallel_evaluation_mp,
                            mp_evaluator = evaluate_work_population,
                            mp_num_cpus = 4,
                            pop_size = 500,
                            bounder=inspyred.ec.Bounder(0, K_MAX),
                            maximize = True,
                            crossover_rate = 0.9,
                            #num_crossover_points=1,
                            tournament_size = 5,
                            num_selected = 2,
                            max_generations = 1000,
                            mutation_rate = 0.2,
                            statistics_file = stat_file,
                            individuals_file = ind_file,
                            points = points,
                            robot = robot,
                            n_elements = n_elements,
                            robot_base = robot_base,
                            fileheaders = fileheaders
                            )
      stat_file.close()
      ind_file.close()
      best = max(final_pop)
      elapsed = (datetime.datetime.now() - started).total_seconds()
      print "elapsed", elapsed, "seconds"
      wr_2.writerow(best.candidate +[best.fitness] + [elapsed])
      #print best.fitness
      #wr_time.writerow([elapsed])
      if display:
        best = max(final_pop) 
        print('Best Solution: \n{0}'.format(str(best)))
    #inspyred.ec.analysis.allele_plot(ind_file_name)
    final_file.close()
    #final_time.close()
    #if display:
    #    best = max(final_pop) 
    #    print('Best Solution: \n{0}'.format(str(best)))
    return ea

def generate_config(random, args):
  n_elements = args['n_elements']
  robot_base = args['robot_base']
  temp = [random.random() * robot_base for i in range(n_elements)]
  return temp

def evaluate_work_population(candidates, args):
    points = args['points']
    robot = args['robot']
    robot_base = args['robot_base']
    fileheaders = args['fileheaders']
    fitness = []
    for individual in candidates:
      fit = 0
      robot.set_variables(individual, fileheaders)
      for p in points:
        aux = robot.ik(p[0], p[1], p[2])
        if aux[0] == 1 and aux[2] != 'err':
          fit += aux[2]
      fitness.append(fit)
    #print fitness
    return fitness

def my_observer(population, num_generations, num_evaluations, args):
    best = max(population)
    print('{0:6} -- {1} : {2}'.format(num_generations, 
                                      best.fitness, 
                                      str(best.candidate)))
