#import numpy as np
from math import *
from my_tools import *

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
    r = Robot(alpha = 120, a = a, b = b, r = r, c = c, dz = dz)
    return r

  @staticmethod 
  def init_zero():
    r = Robot(alpha = 120)
    return r

  def set_variables(self, values, order):
    if len(values) != len(order):
      return -1
    for i in range(len(order)):
      if order[i] == 'a':
        self.__a = values[i]
      elif order[i] == 'b':
        self.__b = values[i]
      elif order[i] == 'e':
        self.__e = values[i]
      elif order[i] == 'd':
        self.__d = values[i]
      elif order[i] == 'c':
        self.__c = values[i]
      elif order[i] == 'r':
        self.__r = values[i]
      elif order[i] == 'c-r':
        self.__c = self.__r + values[i]
      elif order[i] == 'alpha':
        self.__alpha = values[i]
      elif order[i] == 'dz':
        self.__dz = values[i]
    return 1
    
  def set_alpha(self, alpha):
    self.__alpha = alpha

  def set_a(self, a):
    self.__a = a

  def set_b(self, b):
    self.__b = b

  def set_r(self, r):
    self.__r = r

  def set_c(self, c):
    self.__c = c

  def set_d(self, d):
    self.__d = d

  def set_e(self, e):
    self.__e = e

  def set_dz(self, dz):
    self.__dz = dz
  
  def get_alpha(self):
    return self.__alpha

  def get_a(self):
    return self.__a

  def get_b(self):
    return self.__b

  def get_r(self):
    return self.__r

  def get_c(self):
    return self.__c

  def get_d(self):
    return self.__d

  def get_e(self):
    return self.__e

  def set_dz(self, dz):
    self.__dz = dz

  def get_parameters(self):
    r = [self.__alpha, self.__a, self.__b, self.__r, self.__c, self.__d, self.__e, self.__dz]
    return r

  def get_arms_lenght(self):
    return self.__a + self.__b + self.__d + self.__e

  def get_max_reach(self):
    return self.__a + self.__b + self.__d + self.__e + self.__r

  def t_points_half_sphere(self, radius, n_points):
    l= radius
    n_temp=0
    xt = []
    yt = []
    zt = []
    while(n_temp < n_points):
      x_temp = random.uniform(-l,l)
      y_temp = random.uniform(-l,l)
      if(x_temp * x_temp + y_temp * y_temp <= l):
        xt.append(x_temp)
        yt.append(y_temp)
        n_temp = n_temp + 1
      zt.append(random.uniform(0,l))
    return xt, yt, zt

  def valid_point(self,x,y,z):
    r_b = self.__r
    a = self.__a  
    b = self.__b
    c = self.__c
    d = self.__d
    e = self.__e  
    dz = self.__dz
  
    alpha_values=[0, radians(self.__alpha), radians(self.__alpha*2)]

    for i in range(3):
      alpha=alpha_values[i]
      pu = x * cos(alpha) + y * sin(alpha) - r_b
      pv = -x * sin(alpha) + y * cos(alpha)
      pw = z - dz

      if b == 0:
        return False
      aux = pv / b
      if aux < -1 or aux > 1:
        return False

      teta_3 = acos(aux)
      temp_1 = pw**2 + pu**2 + 2*c*pu -b**2* sin(teta_3)**2 - 2*b*e* sin(teta_3)-2*b*d* sin(teta_3)-2*d*e+a**2+c**2-d**2-e**2
                
      l0= temp_1-2*a*c-2*a*pu
      l1= -4*a*pw
      l2= temp_1+2*a*c+2*a*pu
      discriminant= l1**2-4*l2*l0

      if discriminant <= 0:
        return False

    return True

  def ik(self,x,y,z,return_mode=0):
        r_b = self.__r
        a = self.__a  
        b = self.__b
        c = self.__c
        d = self.__d
        e = self.__e  
        dz = self.__dz
        
        alpha_values=[0, radians(self.__alpha), radians(self.__alpha*2)]
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
