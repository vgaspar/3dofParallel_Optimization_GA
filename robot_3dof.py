#import numpy as np
from math import *


class Robot:
  def __init__(self,alpha=0, a=0, d=0, b=0, e=0, r=0, c=0, dz=0):
    self.__alpha = alpha
    self.__r = r
    self.__a = a
    self.__d = d
    self.__b = b
    self.__e = e
    self.__c = c
    self.__dz = dz

  @staticmethod
  def init_simetric( a, d, b, e, r, c, dz=0):
    r = Robot(alpha=120, a = a, b = b, d = d, e = e, r = r, c = c ,dz = dz)
    return r

  @staticmethod 
  def init_fast( a, b, r, c, dz=0):
    r = Robot(alpha=120, a = a,b = b, r = r, c = c ,dz = dz)
    return r

  def ik(self,x,y,z,return_mode=0):
        r_b = self.__r
        a = self.__a  
        b = self.__b
        c = self.__c
        d = self.__d
        e = self.__e  
        dz = self.__dz
        alpha_i = self.__alpha
        
        alpha_values=[0, radians(alpha_i), radians(alpha_i*2)]
        teta_1=[]
        j1=[]
        j2=[]
        j3=[]
        for i in range(3):
                alpha=alpha_values[i]
                pu = x * cos(alpha) + y * sin(alpha) - r_b
                pv = -x* sin(alpha) + y * cos(alpha)
                pw = z - dz

                if b==0:
                        if(return_mode==0):
                              return 0,0,'err','err'
                        else:
                              return 0
                aux = pv / b
                if aux < -1 or aux > 1:
                        if(return_mode==0):
                              return 0,0,'err','err'
                        else:
                              return 0
                
                teta_3 = acos(aux)
                temp_1 = pw**2 + pu**2 + 2*c*pu -b**2* sin(teta_3)**2 - 2*b*e* sin(teta_3)-2*b*d* sin(teta_3)-2*d*e+a**2+c**2-d**2-e**2
                
                l0= temp_1-2*a*c-2*a*pu
                l1=-4*a*pw
                l2= temp_1+2*a*c+2*a*pu
                discriminant=l1**2-4*l2*l0


                if discriminant > 0:
                        t1=(-l1- sqrt(discriminant))/(2*l2)
                        teta_1.append( atan2(2*t1,1-t1**2))

                        teta_2= atan2(pw-a*sin(teta_1[i]),pu-a*cos(teta_1[i])+c)
                
                        temp_2=cos(teta_2)*cos(alpha)*sin(teta_3)
                        temp_3=a*sin(teta_1[i]-teta_2)
                        temp_4=temp_3*sin(teta_3)
                        temp_5=cos(teta_3)*sin(alpha)

                        if temp_4 !=0 and temp_3 !=0:
                                j1.append((temp_5-temp_2)/temp_4)
                                j2.append((-1*cos(teta_3)*cos(alpha)-temp_2)/temp_4)
                                j3.append(-sin(teta_2)/temp_3)
                        else:
                                if(return_mode==0):
                                      return 1,teta_1,'err','err'
                                else:
                                      return 1

                else:
                        if(return_mode==0):
                              return 0,0,'err','err'
                        else:
                              return 0
        try:
                A=[j1,j2,j3]
                norm=normalize(A)
                Ai=inv(A)
                normi=normalize(Ai)
                si=1/(norm*normi)
                return 1,teta_1,si,A
        except Exception:
                if(return_mode==0):
                      return 1,teta_1,'err','err'
                else:
                      return 1


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
###################