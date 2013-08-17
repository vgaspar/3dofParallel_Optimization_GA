import random
from math import *


def gauss_jordan(m, eps = 1.0/(10**10)):
  """Puts given matrix (2D array) into the Reduced Row Echelon Form.
     Returns True if successful, False if 'm' is singular.
     NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
     Written by Jarno Elonen in April 2005, released into Public Domain"""
  (h, w) = (len(m), len(m[0]))
  for y in range(0,h):
    maxrow = y
    for y2 in range(y+1, h):    # Find max pivot
      if abs(m[y2][y]) > abs(m[maxrow][y]):
        maxrow = y2
    (m[y], m[maxrow]) = (m[maxrow], m[y])
    if abs(m[y][y]) <= eps:     # Singular?
      return False
    for y2 in range(y+1, h):    # Eliminate column y
      c = m[y2][y] / m[y][y]
      for x in range(y, w):
        m[y2][x] -= m[y][x] * c
  for y in range(h-1, 0-1, -1): # Backsubstitute
    c  = m[y][y]
    for y2 in range(0,y):
      for x in range(w-1, y-1, -1):
        m[y2][x] -=  m[y][x] * m[y2][y] / c
    m[y][y] /= c
    for x in range(h, w):       # Normalize row y
      m[y][x] /= c
  return True
def solve(M, b):
  """
  solves M*x = b
  return vector x so that M*x = b
  :param M: a matrix in the form of a list of list
  :param b: a vector in the form of a simple list of scalars
  """
  m2 = [row[:]+[right] for row,right in zip(M,b) ]
  return [row[-1] for row in m2] if gauss_jordan(m2) else None

def inv(M):
  """
  return the inv of the matrix M
  """
  #clone the matrix and append the identity matrix
  # [int(i==j) for j in range_M] is nothing but the i(th row of the identity matrix
  m2 = [row[:]+[int(i==j) for j in range(len(M) )] for i,row in enumerate(M) ]
  # extract the appended matrix (kind of m2[m:,...]
  return [row[len(M[0]):] for row in m2] if gauss_jordan(m2) else None

def zeros( s , zero=0):
    """
    return a matrix of size `size`
    :param size: a tuple containing dimensions of the matrix
    :param zero: the value to use to fill the matrix (by default it's zero )
    """
    return [zeros(s[1:] ) for i in range(s[0] ) ] if not len(s) else zero

def normalize(A):
    s = sqrt(sum([ sum([ (e**2) for e in r ]) for r in A ]))
    return s

def sphere_volume(radius):
  return 4/3 * pi * radius**3

def cylinder_volume(radius,h):
  return pi * radius**2 * h

def point_in_sphere(radius, max_phi=pi, max_theta=2*pi, x0=0, y0=0, z0=0):
  r =  random.uniform(0,radius)
  phi = random.uniform(0,max_phi)
  theta = random.uniform(0,max_theta)

  x = x0 + r * cos(theta) * sin(phi)
  y = y0 + r * sin(theta) * sin(phi)
  z = z0 + r * cos(phi)

  return x,y,z

def point_in_cylinder(radius, h, max_phi=2*pi, x0=0, y0=0, z0=0):
  rho = random.uniform(0, radius)
  phi = random.uniform(0, max_phi)
  h_r = random.uniform(0, h)

  x = x0 + rho * cos(phi)
  y = y0 + rho * sin(phi)
  z = z0 + h_r

  return x, y, z

def half_sphere_cloud_points(radius, n_points):
  points = [point_in_sphere(radius, max_theta=pi) for i in range(n_points)]
  vol = sphere_volume(radius)/2
  return points, vol

def make_cube_vertices(center, side, dH):
  d = side/2.0
  a1 = tuple(map(lambda x, y: x + y, center, (d, d, 0)))  
  a2 = tuple(map(lambda x, y: x + y, center, (d, -d, 0)))
  a3 = tuple(map(lambda x, y: x + y, center, (-d, -d, 0)))
  a4 = tuple(map(lambda x, y: x + y, center, (-d, d, 0)))
  a5 = tuple(map(lambda x, y: x + y, a1, (0, 0, dH)))
  a6 = tuple(map(lambda x, y: x + y, a2, (0, 0, dH)))
  a7 = tuple(map(lambda x, y: x + y, a3, (0, 0, dH)))
  a8 = tuple(map(lambda x, y: x + y, a4, (0, 0, dH)))  
  return [a1, a2, a3, a4, a5, a6, a7, a8]

def average(s): return sum(s) * 1.0 / len(s)

def variance(s,average):return map(lambda x: (x - average)**2, s)

def standard_deviation(variance):return sqrt(average(variance))

def median(populalion_fitness):
    sorts = sorted(populalion_fitness)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]
