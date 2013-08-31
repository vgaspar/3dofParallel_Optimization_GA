import genetic_algorithm as ga
import my_tools
from robot_3dof import *

import csv

def minmax_cube_cords(cube_vertices):
	x_max = x_min = cube_vertices[0][0]
	y_max = y_min = cube_vertices[0][1]
	z_max = z_min = cube_vertices[0][2]

	for i in range(len(cube_vertices)):
		if cube_vertices[i][0] > x_max:
			x_max = cube_vertices[i][0]
		elif cube_vertices[i][0] < x_min:
			x_min = cube_vertices[i][0]
		if cube_vertices[i][1] > y_max:
			y_max = cube_vertices[i][1]
		elif cube_vertices[i][1] < y_min:
			y_min = cube_vertices[i][1]
		if cube_vertices[i][2] > z_max:
			z_max = cube_vertices[i][2]
		elif cube_vertices[i][2] < z_min:
			z_min = cube_vertices[i][2]
	return x_min,x_max,y_min,y_max,z_min,z_max

def random_point_on_cube(x_min, x_max, y_min, y_max, z_min, z_max):
	x = random.uniform(x_min, x_max)
	y = random.uniform(y_min, y_max) 
	z = random.uniform(z_min, z_max)

	return (x, y, z)

def n_random_points_on_cube(cube_vertices, number_points):
	minmax_cords = minmax_cube_cords(cube_vertices)
	points = []  
	for i in range(number_points):
		points.append(random_point_on_cube(minmax_cords[0], minmax_cords[1],
		minmax_cords[2], minmax_cords[3],
		minmax_cords[4], minmax_cords[5]))
	return points

def save_points_csv(points, filename):
	file_name = filename + '.csv'
	final_file = open(file_name, 'w')
	wr_2 = csv.writer(final_file, quoting=csv.QUOTE_ALL)
	wr_2.writerow(['x','y','z'])
	for p in points:
		wr_2.writerow(p)
	return 1
def save_points_csv_fitness(points, filename, fitness):
	file_name = filename + '.csv'
	final_file = open(file_name, 'w')
	wr_2 = csv.writer(final_file, quoting=csv.QUOTE_ALL)
	wr_2.writerow(['x','y','z','fitness'])
	i=0
	for p in points:
		wr_2.writerow([p[0],p[1],p[2],fitness[i]])
		i+=1
	return 1

robot = Robot.init_zero()
robot.set_r(1)
cube_vertices = my_tools.make_cube_vertices((0,0,0.8), 1, .2)

save_points_csv(cube_vertices, 'vertice')

par_order = ["a","b","c-r"]

filename = 'first_tests'
runs = 30

points = n_random_points_on_cube(cube_vertices, 100000)
save_points_csv(points, 'cube_cloud')
#ga.ga(display = True , robot = robot, points=points, filename=filename, fileheaders=par_order, runs=runs)

filename = 'first_tests_vertices_3'
#filename = 'aa'
ga.ga(display = True , robot = robot, points=cube_vertices, filename=filename, fileheaders=par_order, runs=runs, eval_method=3)

csvfile = open('vertice.csv', 'rb')
reader = csv.reader(csvfile, delimiter=',')
points_v = []
for row in reader:
	try:
		points_v.append([float(i) for i in row])
	except:
		pass
print points_v

csvfile = open('first_tests_vertices_solutions.csv', 'rb')
reader = csv.reader(csvfile, delimiter=',')

header = reader.next()

f_index = 0
for i in range(len(header)):
	if header[i] == 'fitness':
		f_index = i
		break
max_fit = 0
best = []
for row in reader:
	try:
		 if float(row[f_index]) > max_fit:
		 	best = [float(i) for i in row[:f_index]]
	except:
		pass
print header[:f_index]
print best

robot.set_variables(best, header[:f_index])
points = cube_vertices
print len(points)

fit= robot.cloud_condition_number( points)
filename='bb'
save_points_csv_fitness(points, filename, fit)
print fit