from operator import itemgetter
from math import *
from copy import deepcopy
import datetime
import random
import math
import shelve
import csv
import inverse_kinematics as ik



####### PROBLEM PARAMETERS ##############

NUM_TEST_POINTS = 10000
RUNS = 1
GENERATIONS = 250
POP_SIZE = 150
selection_method="tournament"
dif_stop_condition = 50 
NEW_INDI_GENE_DIF = 10
NEW_INDI_PERCENTAGE = 0.92
counter=0
chromosome_size = 5
structures_dic = {}

####### GA PARAMETERS ##############
TOURNAMENT_SIZE = 2
XOVER_PROB = 0.4
INDIV_MUTATION = 0.5
ALLELE_MUTATION = 0.6
ELITISM_PERCENTAGE = 0.02
NUM_CUT_POINTS = 2
INDIV_ZERO_MUTATION = 0.5
GEO_XOVER_PROB = 0.5


DEV_ALLE_MUTATION = 0.1
MAX_K_START = 5.0

###################
def main(prng=None, display=False):

  ##### robot initial configuration ######
  parameters = {}
  parameters['loop_length'] = 1
  parameters['a'] = 0
  parameters['b'] = 0
  parameters['d'] = 0
  parameters['e'] = 0
  parameters['c'] = 0
  parameters['r_b'] = 1
  parameters['dz'] = 0

  population = init_population(POP_SIZE)

  csv_file = "ga_paper_my_method_population.csv"
  myfile = open(csv_file, 'wb')
  wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
  csv_file_2 = "ga_paper_my_method_genome.csv"
  myfile_2 = open(csv_file_2, 'wb')
  wr_2 = csv.writer(myfile_2, quoting=csv.QUOTE_ALL)
  wr_2.writerow(["a","b","d","e","c-r","time_taken"])
  wr.writerow(["generation number","population size","worst","best","median","average","standard deviation","run"])
  for i in range(0,RUNS):
    started = datetime.datetime.now()
    print "run: %d" %i

    best_in_generation,best_guy,results = ga(parameters,generations=GENERATIONS, pop_size=POP_SIZE,selection_method=selection_method,population=population,dif_stop_condition= dif_stop_condition,csv_wr=wr,run=i)
    
    print "best_guy: ",best_guy.get_chromosome()
    print "best_fitness: ",
    print best_guy.get_fitness(),
    print " in Generation: ",
    print best_in_generation

    elapsed = (datetime.datetime.now() - started).total_seconds()
    print "elapsed", elapsed, "seconds"
    tbg=best_guy.get_chromosome()
    wr_2.writerow([tbg[0],tbg[1],tbg[2],tbg[3],tbg[4],elapsed])

  myfile_2.close()
  myfile.close()

########## PROBLEM SPECIFIC ###################
def validate_individual(temp):
  if tuple(temp) in structures_dic:
    return structures_dic[tuple(temp)]
  else:
    return False

def evaluate_work_population(population,parameters,method,xt,yt,zt):
    for individual in population:
        temp = individual.get_chromosome()
        f=validate_individual(temp)
        if not f:
          if method == 2:

            if check_arms_length_contrain(temp):
              parameters['a'] = temp[0]
              parameters['b'] = temp[1]
              parameters['d'] = temp[2]
              parameters['e'] = temp[3]
              parameters['c'] = temp[4] + parameters['r_b']
              l=sum(temp[:4])
              x=[x*l for x in xt]
              y=[y*l for y in yt]
              z=[z*(l+parameters['dz']) for z in zt]
              neu = test_configuration(parameters,x,y,z)
              fitness = neu
            else:
              fitness = -1
            individual.set_fitness(fitness)
            structures_dic[tuple(temp)] = fitness
          if method == 1:
              pass
        else:
           fitness = f
           #print "a 2nd hand confffffffffffffffffffffffffffffffffffff"


def check_arms_length_contrain(param):
  if sum(param[:4]) <= 1:
    return 1
  else:
    return 0


###############################################
##generation number, population size, worst, best, median, average, standard deviation
def ga(parameters,**kwargs):
  wr=kwargs["csv_wr"]
  run_id=kwargs["run"]
  generations = kwargs["generations"]
  population = kwargs['population']
  dif_cond = kwargs['dif_stop_condition']
  
  current_generation = 0
  best_in_generation = 0
  best_dif_generation = 0 
  best_guy = Individual()
  results=[]
  x,y,z=create_test_points(parameters)
  #### TESTING POINTS #####
  #------ 
  evaluate_work_population(population,parameters,2,x,y,z)

  while current_generation < generations and best_dif_generation < dif_cond :
    best_dif_generation = current_generation-best_in_generation
    #generation_fitness.append([individual.get_fitness() for individual in population])
    print
    print "Generation: %d" % (current_generation + 1)
    print
    parents = select_parents(population, kwargs["selection_method"])
    new_population = reproduce_parents(parents)

    population = merge_populations(population, new_population,best_dif_generation, 2)
    #print "new pop after merge:", map(lambda x: x.get_id(), population)

    evaluate_work_population(population,parameters,2,x,y,z)
    #guardar media da populacao
    
    population.sort(key=itemgetter("fitness"), reverse=True)
    pop_fit=[individual.get_fitness() for individual in population]
    avg=average(pop_fit)
    var=variance(pop_fit,avg)
    std=standard_deviation(var)
    #generation number, population size, worst, best, median, average, standard deviation
    wr.writerow([current_generation + 1,len(population),population[-1].get_fitness(),population[0].get_fitness(),median(pop_fit),avg,std,run_id])
    
    #print "population:", map(str, population)
    print "bestguy_fit  : ", best_guy.get_fitness()
    print "gen best_fit : ", population[0].get_fitness()
    print "gen n_best : ", population.count(population[0].get_fitness())
    if best_guy.get_fitness() < population[0].get_fitness():
      best_guy=deepcopy(population[0])
      best_in_generation=current_generation
    current_generation += 1
    aux=[population[0].get_fitness()]
    aux.extend(population[0].get_chromosome())
    results.append(aux)
  return best_in_generation,best_guy,results
def average(s): return sum(s) * 1.0 / len(s)
def variance(s,average):return map(lambda x: (x - average)**2, s)
def standard_deviation(variance):return math.sqrt(average(variance))
def median(populalion_fitness):
    sorts = sorted(populalion_fitness)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]
def test_configuration(parameters,x,y,z):
        points_done = 0
        valid_points=0
        S=0
        radius= parameters['a']+parameters['b']+parameters['d']+parameters['e']
        for i in range(NUM_TEST_POINTS):
                aux = ik.inv_rot_robot(x[i],y[i],z[i],parameters)
                if aux[0]==1 and aux[2]!='err':
                        valid_points+=1
                        S=S+aux[2]
        #W = radius**3 * math.pi*valid_points/NUM_TEST_POINTS
        neu = radius**3 *math.pi*S/NUM_TEST_POINTS
        #return W,neu
        return neu

class Individual:
  def __init__(self, chromo=None):
    self.__size = chromosome_size ##__ private parameter
    self.__chromosome = chromo
    if chromo == None:
      temp=[]
      temp.append(random.uniform(0, 1))
      temp.append(random.uniform(0, 1- temp[0]))
      temp.append(random.uniform(0, 1- temp[0]-temp[1]))
      temp.append( 1- temp[0]-temp[1]-temp[2])
      temp.append(random.uniform(0, 1))
      self.__chromosome = temp
    self.__fitness = 0.0
    self.__random = random.random() * 1000

  def get_id(self):
    return "%04d" % self.__random

  def set_fitness(self, fit):
    self.__fitness = fit

  def get_fitness(self):
    return self.__fitness

  def get_chromosome(self):
    return self.__chromosome

  def set_chromosome(self, new_chromo):
    self.__chromosome = new_chromo

  def get_chromosome_size(self):
    return self.__size

  def __getitem__(self, key):
    if key == "fitness":
      return self.__fitness

  def __str__(self):
    t = "[%04d] Fit(%.2f) Chromo %s" % (self.__random, self.__fitness, str(self.__chromosome))
    return t

######SELECTION MECHANISM#############
def roulette_wheel_draft(population):
  pop = population
  total_fitness = sum( [ind.get_fitness() for ind in pop] ) #soma
  mate_pool = []
  for i in range(POP_SIZE):
    value = random.uniform(0,1) # random sampling of an uniform distribution
    index = 0
    total = pop[index].get_fitness() / float(total_fitness)
    while total < value:
      index += 1
      total += pop[index].get_fitness() / float(total_fitness)
    mate_pool.append(pop[index])
  return mate_pool

def tournament_draft(population):
  pop_size = POP_SIZE
  mate_pool = [ tournament(population) for i in range(pop_size) ]
  return mate_pool

def tournament(population):
  pool = random.sample(population, TOURNAMENT_SIZE)
  pool.sort(key=itemgetter("fitness"), reverse=True) #nao vou buscar o get_fitness por causa do sort usar o itemgetter
  return pool[0]

def select_parents(population,method):
  if method == 'tournament':
    return tournament_draft(population)
  if method == 'roulette':
    return roulette_wheel_draft(population)

def n_point_crossover(number_of_points, ind1, ind2):
  off1 = Individual()
  off2 = Individual()
  n_points = []
  for n in range(number_of_points):
    ind = random.randint(0, ind1.get_chromosome_size() - 1)
    while ind in n_points:
      ind = random.randint(0, ind1.get_chromosome_size() - 1)
    n_points.append(ind)
  n_points.sort()
  genomes = swap(n_points, ind1.get_chromosome(),ind2.get_chromosome())
  off1.set_chromosome(genomes[0][:])
  off2.set_chromosome(genomes[1][:])
  return off1, off2

def pseudo_geometric_crossover( ind1, ind2):
  off1 = Individual()
  off2 = Individual()
  genomes = [(ind2.get_chromosome()[i]+ind1.get_chromosome()[i])/2.0 for i in range(len(ind1.get_chromosome()))]
  off1.set_chromosome(genomes)
  #off2.set_chromosome(genomes)

  
  return off1, off2

def swap(n_points, p1, p2):
  #swap between n-points between 2 parents
  #returns the 2 resultant children as a dictionary
  for i in range(0,len(n_points) - 1,2):
    first_point = n_points[i]
    sec_point = n_points[i+1]
#      print first_point
#      print sec_point
    temp = p1[first_point : sec_point]
    p1 = p1[ : first_point] + p2[ first_point : sec_point] + p1[sec_point : ]
    p2 = p2[ : first_point] + temp + p2[sec_point : ]
  return [p1,p2]
######################################

############## GA WORKERS ########################
def init_population(pop_size):
  print "POP_size %d" % pop_size
  temp = [Individual() for i in range(pop_size)]
  return temp
  
def merge_populations(old_pop, new_pop,best_dif_generation,method):
  if method == 1:
    return new_pop
  if method == 2:
    global counter
    old_pop.sort(key=itemgetter("fitness"), reverse=True)
    elit_number=int(ELITISM_PERCENTAGE*POP_SIZE)
    elite = old_pop[ : elit_number  ]
    print "elite:", map(lambda x: x.get_id() + " with " + "%.6f" % (x.get_fitness()), elite)

    if best_dif_generation > NEW_INDI_GENE_DIF and counter > 3:
      new_indi_number=int(NEW_INDI_PERCENTAGE*POP_SIZE)
      new_pop_sample = random.sample(new_pop, POP_SIZE-elit_number-new_indi_number)
      new_indi_pop = [Individual() for i in range(new_indi_number)]
      new_pop_sample.extend(new_indi_pop)
      print "Just Happened an Armageddon!"
      counter=0
    else:
      new_pop_sample = random.sample(new_pop, POP_SIZE -elit_number)
    
    new_pop_sample.extend(elite)
    counter = counter +1
    return new_pop_sample
  if method == 3:
    size_old_pop = POP_SIZE
    merged_populations = old_pop + new_pop
    merged_populations.sort(key=itemgetter("fitness"), reverse=True)
    return merged_populations[ : size_old_pop]

def mutate_individual_work(off):
  aux = off.get_chromosome()
  for i in range(len(aux)):
    if (random.random() < ALLELE_MUTATION):
      while aux[i] < 0:
        aux[i] = random.gauss(aux[i], DEV_ALLE_MUTATION)
  off.set_chromosome(aux)

def mutate_individual_zero(off):
  aux = off.get_chromosome()
  for i in range(len(aux)):
    if (random.random() < ALLELE_MUTATION):
      aux[i] = 0
  off.set_chromosome(aux)

def reproduce_parents(parents):
  new_population = []
  for i in range(0, len(parents) - 1, 2):
    ind1 = parents[i]
    ind2 = parents[i + 1]
    
    if(random.random() < XOVER_PROB):
      #CROSSOVER APPLICATION
      off1, off2 = n_point_crossover(NUM_CUT_POINTS,ind1, ind2)
    elif(random.random() < GEO_XOVER_PROB):
      off1, off2 = pseudo_geometric_crossover(ind1, ind2)
    else:
      off1 = ind1
      off2 = ind2
    for m in range(0,2):
      
      if(random.random() < INDIV_MUTATION):
        if(m == 0):
          mutate_individual_work(off1)
        else:
          mutate_individual_work(off2)

      if(random.random() < INDIV_ZERO_MUTATION):
        if(m == 0):
          mutate_individual_work(off1)
        else:
          mutate_individual_work(off2)       

    new_population.append(off1)
    new_population.append(off2)
  return new_population

def create_test_points(parameters):
  #random.seed(42)
  l= 1
  n_temp=0
  xt=[]
  yt=[]
  zt=[]
  while(n_temp < NUM_TEST_POINTS):
    #x_temp=np.random.rand()*2-1
    x_temp=random.uniform(-l,l) # IN ORDER TO the loop length
    y_temp=random.uniform(-l,l)
    if(x_temp*x_temp+y_temp*y_temp<=l):
      xt.append(x_temp)
      yt.append(y_temp)
      n_temp=n_temp+1
    zt.append(random.uniform(0,l))
  return xt,yt,zt

if __name__ == "__main__":
  main()
