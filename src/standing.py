#Standing waves in Python
#TODO: reduce number of terms
import numpy as np
from math import *
import itertools as it
from scipy.integrate import solve_ivp
import collections


#Helper functions:

compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

def modes(sol):
	return np.float64([[sol['y'][i][n] for i in range(len(sol['y']))] for n in range(len(sol['y'][0]))])

def phase(sol):
	import matplotlib.pyplot as plt
	mode = modes(sol)
	size = int(len(mode)/2)
	fig, axes = plt.subplots(nrows=size,figsize=(10,5*size))

	for i in range(size):
		j=2*i
		axes[i].plot(mode[j],mode[j+1])
		axes[i].set_xlabel('a'+str(i+1))
		axes[i].set_ylabel('b'+str(i+1))

	return fig, axes

#Main Class:

class StandingSim:
	def __init__(self,max_order=2,mode=1,harmonics=2,
		params={'h':1000.0,'g0':1.0,'sigma':1.0,'k0':pi}):

		self.mode = mode
		self.max_n = max_order
		self.harmonics = harmonics
		self._js = [j for j in range(mode,mode*(harmonics+1),mode)]
		self.js = [-j for j in range(mode*harmonics,0,-mode)]+[j for j in range(mode,mode*(harmonics+1),mode)]
		self._h = params['h']
		self.g0 = params['g0']
		self.sigma = params['sigma']
		self.k0 = params['k0'] #unitless base wavenumber

		#------------------#
		# Differential eqs:#
		#------------------#

		#Initialize non-zero combinations using delta function. Initialize ODE coefficients

		non_zero = {}
		constant = {}
		for n in range(2,self.max_n+1):
			perms = [list(ele) for ele in it.product(self.js,repeat=n)] #Find all permutations of a given length and elements
			cs={}
			tp_dict={}
			for j in self._js:
				c = []
				tp_list=[]
				for perm in perms:
					if j-sum(perm)==0:  #delta filtering function
						c.append(perm)
				for onec in c:
					vector = list([j]+onec)
					tp_list.append(np.float64([self.alpha(n,vector),self.beta(n,vector),self.gamma(n,vector)])) #Find appropiate coeffs for a given perm

				if n>=4: #Hamiltonian simplification
					i=0
					while i<len(c):

					    l = i+1
					    while l<len(c):
					        if compare(list(map(abs,c[i][2:])), list(map(abs,c[l][2:]))) and (list(map(abs,c[i][:2]))==list(map(abs,c[l][:2]))):
					            tp_list[i] += tp_list[l]
					            del c[l]
					            del tp_list[l]
					        else:
					            l = l+1
					    i=i+1

				cs[j] = c
				tp_dict[j]=tp_list
			non_zero[n] = cs
			constant[n] = tp_dict			
		self.deltas = non_zero
		self.constants = constant


		#--------------#
		# Hamiltonian: #
		#--------------#

		#Initialize non-zero combinations using delta function

		non_zero_h = {}
		for n in range(2,self.max_n+2):
			perms = [list(ele) for ele in it.product(self.js,repeat=n)]
			c_h=[]
			for perm in perms:
				if sum(perm)==0:
					c_h.append(perm)
			non_zero_h[n] = c_h			
		self.deltas_h = non_zero_h

		#Initialize coefficients (h,s)

		constant_h = {}
		for n in range(3,self.max_n+2):
			tp_list_h=[]
			for delta_h in self.deltas_h[n]:
				vector = delta_h
				tp_list_h.append([self.h(n,vector),self.s(n,vector)])
			constant_h[n] = tp_list_h
		self.constants_h = constant_h


	#Define physical quantities: 
	def k(self,vec):
		#wavenumber
		return self.k0*vec
	def _k(self,vec):
		#absolute value of wave_number
		try: 
			len(vec)
			return np.absolute(sum([self.k(v) for v in vec]))
		except: return np.absolute(self.k(vec)) #TODO: absolute function not what we want
	def T(self,vec):
		#period
		return tanh(self._k(vec)*self._h)
	def g(self,vec):
		#gravity
		return self.g0*(1+self.sigma/self.g0*self._k(vec)**2)
	def G(self,vec,t):
		#gravity with forcing
		return self.g(vec)+sum([j[0]*cos(j[1]*t) for j in self.forcing]) 

	def omega(self,vec):
		#frequency
		return sqrt(self.g0*(1+self.sigma/self.g0*self._k(vec)**2)*self._k(vec)*self.T(vec))

	#Define coefficient functions of ODEs
	def h(self,n,vec):
		if n == 3:
			h=-(np.dot(self.k(vec[0]),self.k(vec[1]))+ \
				self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*self.T(vec[1]))
		elif n ==4:
			h=-(self._k(vec[0])*self.T(vec[0])*self._k(vec[1])**2-\
				self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*\
				self.T(vec[1])*self._k([vec[1],vec[2]])*self.T([vec[1],vec[2]]))
		elif n == 5:
			h=-((1/6.0*self._k(vec[1])**2-0.5*self._k([vec[1],vec[2]])**2+\
				self._k([vec[0],vec[2]])*self.T([vec[0],vec[2]])*self._k([vec[1],vec[3]])*self.T([vec[1],vec[3]]))*\
				self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*self.T(vec[1])-\
				self._k(vec[1])**2*self._k(vec[0])*self.T(vec[0])*self._k([vec[1],vec[2],vec[3]])*\
				self.T([vec[1],vec[2],vec[3]])+1/3.0*self._k(vec[0])**2*self._k(vec[1])**2)
		elif n == 6:
			h=((self._k(vec[0])*self.T(vec[0])*self._k([vec[0],vec[2]])*self.T([vec[0],vec[2]])-self._k(vec[0])**2)*\
				self._k(vec[1])*self.T(vec[1])*self._k([vec[1],vec[3]])*self.T([vec[1],vec[3]])*self._k([vec[1],vec[3],vec[4]])*self.T([vec[1],vec[3],vec[4]])+\
				(1/3.0*self._k(vec[1])**2-self._k([vec[1],vec[3]])**2)*self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*self.T(vec[1])*\
				self._k([vec[0],vec[2]])*self.T([vec[0],vec[2]])+0.25*self._k(vec[0])**2*self._k(vec[1])**2*self._k([vec[1],vec[4],vec[5]])*\
				self.T([vec[1],vec[4],vec[5]])+0.5*self._k(vec[0])**2*self._k([vec[0],vec[2],vec[3]])**2*self._k(vec[1])*self.T(vec[1])\
				-1/12.0*self._k(vec[0])**4*self._k(vec[1])*self.T(vec[1]))
		else: h=0.0
		return np.float64(h) 
	def s(self,n,vec):
		if n == 3:
			s=0.0
		elif n == 4:
			s=-1/4*self.sigma*np.dot(self.k(vec[0]),self.k(vec[1]))*np.dot(self.k(vec[2]),self.k(vec[3]))
		elif n == 5:
			s=0.0
		elif n == 6:
			s=-1/8.0*self.sigma*np.dot(self.k(vec[0]),self.k(vec[1]))*np.dot(self.k(vec[2]),self.k(vec[3]))*np.dot(self.k(vec[4]),self.k(vec[5]))
		else: s=0.0
		return np.float64(s)
	def alpha(self,n,vec):
		h1=list(vec)
		h2=list(vec)
		h1[0],h2[0],h2[1]=-vec[0],vec[1],-vec[0]
		return 1/2*(self.h(n+1,h1)+self.h(n+1,h2))
	def beta(self,n,vec):
		beta=0
		for i in range(2,len(vec)):
			vector=vec[1:i+1]+[-vec[0]]+vec[i+1:]
			beta=beta+self.h(n+1,vector)
		return -1/2*beta
	def gamma(self,n,vec):
		gamma=0
		for i in range(1,len(vec)+1):
			vector=vec[1:i]+[-vec[0]]+vec[i:]
			gamma=gamma+self.s(n+1,vector)
		return -1/2*gamma

	def prod(self,x,vec,ab):
		#Prod function allows to compute product of arbitrary number of solution elements
		prod = np.float64(1)
		for v,a in zip(np.absolute(vec),ab):
			prod = prod*x[int(v/self.mode-1)][a]
		return prod
	'''
	def prod(self,x,vec,ab):
		#Prod function allows to compute product of arbitrary number of solution elements
		prod = np.prod([x[int(v/self.mode-1)][a] for v,a in zip(vec,ab)])
		return prod
	'''
	def a_prime(self,x,j,t):
		linear = np.float64(self._k(j)*self.T(j)*x[int(j/self.mode-1)][1])
		non_linear = np.float64(0)
		for n in range(2,self.max_n+1):
			for i in range(len(self.deltas[n][j])):
				c = list(self.deltas[n][j][i])
				non_linear = non_linear + self.constants[n][j][i][0]*np.float64(self.prod(x,c,[1]+[0 for l in range(len(c)-1)]))
		return linear+non_linear	

	def b_prime(self,x,j,t):
		linear = - np.float64(self.G(j,t)*x[int(j/self.mode-1)][0])
		non_linear = np.float64(0)
		for n in range(2,self.max_n+1):
			for i in range(len(self.deltas[n][j])):
				c1,c2 = list(self.deltas[n][j][i]),list(self.deltas[n][j][i])
				non_linear = non_linear + self.constants[n][j][i][1]*np.float64(self.prod(x,c1,[1,1]+[0 for l in range(len(c1)-2)]))
				non_linear = non_linear + self.constants[n][j][i][2]*np.float64(self.prod(x,c2,[0 for l in range(len(c2))]))
		return linear+non_linear

	def f(self,x,t):
		#Main differential function
		#TODO: Can we split this into j processes??
		x = [[x[n],x[n+1]] for n in range(0,len(x),2)]
		f = [[0,0] for n in range(len(self._js))]
		for j in range(len(f)):
			f[j]=[self.a_prime(x,self.mode*(j+1),t),self.b_prime(x,self.mode*(j+1),t)]
		return np.array(f).flatten()

	def hamiltonian(self,x,t):
		x = [[x[n],x[n+1]] for n in range(0,len(x),2)]
		hamiltonian = np.float64(0.0)
		for n in range(2,self.max_n+2):
			if n==2:
				for vec in self.deltas_h[n]:
					hamiltonian = hamiltonian + (self.g0-self.sigma*np.dot(self.k(vec[0]),self.k(vec[1])))*self.prod(x,vec,[0,0]) +\
					 self._k(vec[0])*self.T(vec[0])*self.prod(x,vec,[1,1])
			else:
				for vec in self.deltas_h[n]:
					hamiltonian = hamiltonian + self.h(n,vec)*self.prod(x,vec,[1,1]+[0 for l in range(len(vec)-2)]) +\
					self.s(n,vec)*self.prod(x,vec,[0 for l in range(len(vec))])

		return 0.5*hamiltonian

	def solve(self,dt=0.01,N=1000,x0=np.array(None),
		forcing=[[0,0]],h_freq=100,method='homemade'):
		#Solved by employing RK4 method
		
		#Forcing parameter allows to input amplitude, frequency pairs for a linear combination of cosines that act as forcing function
		self.forcing = forcing

		ts = np.linspace(0.0,(N-1)*dt,num=N,dtype='float64')

		if method=='homemade':
			if x0.all()==None: 
				x0 = [0.001,0.0]+[0.0 for p in range(2*len(self._js)-2)]
				x0=np.float64(x0)
			hamiltonian = [self.hamiltonian(x0,ts[0])]
			sol = [x0]

			for q in range(1,N):

				ti = ts[q]
				k1 = dt*self.f(sol[q-1],ti)
				k2 = dt*self.f(sol[q-1]+0.5*k1,ti+dt*0.5)
				k3 = dt*self.f(sol[q-1]+0.5*k2,ti+dt*0.5)
				k4 = dt*self.f(sol[q-1]+k3,ti+dt)

				xi = sol[q-1] + 1/6.0*(k1+2.0*k2+2.0*k3+k4)
				sol.append(xi)
				if q%h_freq==0:
					hamiltonian.append(self.hamiltonian(xi,ti))
				else:
					hamiltonian.append(hamiltonian[q-1])

			return {'y':sol,'t':ts,'H':hamiltonian}


class Standing:
	def __init__(self,max_order=2,mode=1,harmonics=2,
		params={'h':1000.0,'g0':1.0,'sigma':1.0,'k0':pi}):

		self.mode = mode
		self.max_n = max_order
		self.harmonics = harmonics
		self._js = [j for j in range(mode,mode*(harmonics+1),mode)]
		self.js = [-j for j in range(mode*harmonics,0,-mode)]+[j for j in range(mode,mode*(harmonics+1),mode)]
		self._h = params['h']
		self.g0 = params['g0']
		self.sigma = params['sigma']
		self.k0 = params['k0'] #unitless base wavenumber

		#------------------#
		# Differential eqs:#
		#------------------#

		#Initialize non-zero combinations using delta function
		'''
		non_zero = {}
		for n in range(2,self.max_n+1):
			perms = [list(ele) for ele in it.product(self.js,repeat=n)]
			cs={}
			for j in self._js:
				c = []
				for perm in perms:
					if j-sum(perm)==0:
						c.append(perm)
				cs[j] = c
			non_zero[n] = cs			
		self.deltas = non_zero

		#Initialize coefficients (alphas,betas,gammas)

		constant = {}
		for n in range(2,self.max_n+1):
			tp_dict={}
			for j in self._js:
				tp_list=[]
				for delta in self.deltas[n][j]:
					vector = [j]+delta
					tp_list.append([self.alpha(n,vector),self.beta(n,vector),self.gamma(n,vector)])
				tp_dict[j]=tp_list
			constant[n] = tp_dict
		self.constants = constant
		'''

		non_zero = {}
		constant = {}
		for n in range(2,self.max_n+1):
			perms = [list(ele) for ele in it.product(self.js,repeat=n)] #Find all permutations of a given length and elements
			cs={}
			tp_dict={}
			for j in self._js:
				c = []
				tp_list=[]
				for perm in perms:
					if j-sum(perm)==0:  #delta filtering function
						c.append(perm)
				for onec in c:
					vector = list([j]+onec)
					tp_list.append(np.float64([self.alpha(n,vector),self.beta(n,vector),self.gamma(n,vector)])) #Find appropiate coeffs for a given perm

				cs[j] = c
				tp_dict[j]=tp_list
			non_zero[n] = cs
			constant[n] = tp_dict			
		self.deltas = non_zero
		self.constants = constant
		#--------------#
		# Hamiltonian: #
		#--------------#

		#Initialize non-zero combinations using delta function

		non_zero_h = {}
		for n in range(2,self.max_n+2):
			perms = [list(ele) for ele in it.product(self.js,repeat=n)]
			c_h=[]
			for perm in perms:
				if sum(perm)==0:
					c_h.append(perm)
			non_zero_h[n] = c_h			
		self.deltas_h = non_zero_h

		#Initialize coefficients (h,s)

		constant_h = {}
		for n in range(3,self.max_n+2):
			tp_list_h=[]
			for delta_h in self.deltas_h[n]:
				vector = delta_h
				tp_list_h.append([self.h(n,vector),self.s(n,vector)])
			constant_h[n] = tp_list_h
		self.constants_h = constant_h


	#Define physical quantities: 
	def k(self,vec):
		#wavenumber
		return self.k0*vec
	def _k(self,vec):
		#absolute value of wave_number
		try: 
			len(vec)
			return np.absolute(sum([self.k(v) for v in vec]))
		except: return np.absolute(self.k(vec)) #TODO: absolute function not what we want
	def T(self,vec):
		#period
		return tanh(self._k(vec)*self._h)
	def g(self,vec):
		#gravity, g_j
			return self.g0*(1+self.sigma/self.g0*self._k(vec)**2)

	def G(self,vec,t):
		#gravity
		try:
			if (self.force_eps and self.force_del)==0:
				return self.g(vec)
			else:
				return self.g(vec)+self.force_amp*cos(self.force_omega*t)
		except:
			return self.g(vec)

	def omega(self,vec):
		#frequency
		return sqrt(self.g(vec)*self._k(vec)*self.T(vec))

	#Define coefficient functions of ODEs
	def h(self,n,vec):
		if n == 3:
			h=-(np.dot(self.k(vec[0]),self.k(vec[1]))+ \
				self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*self.T(vec[1]))
		elif n ==4:
			h=-(self._k(vec[0])*self.T(vec[0])*self._k(vec[1])**2-\
				self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*\
				self.T(vec[1])*self._k([vec[1],vec[2]])*self.T([vec[1],vec[2]]))
		elif n == 5:
			h=-((1/6.0*self._k(vec[1])**2-0.5*self._k([vec[1],vec[2]])**2+\
				self._k([vec[0],vec[2]])*self.T([vec[0],vec[2]])*self._k([vec[1],vec[3]])*self.T([vec[1],vec[3]]))*\
				self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*self.T(vec[1])-\
				self._k(vec[1])**2*self._k(vec[0])*self.T(vec[0])*self._k([vec[1],vec[2],vec[3]])*\
				self.T([vec[1],vec[2],vec[3]])+1/3.0*self._k(vec[0])**2*self._k(vec[1])**2)
		elif n == 6:
			h=((self._k(vec[0])*self.T(vec[0])*self._k([vec[0],vec[2]])*self.T([vec[0],vec[2]])-self._k(vec[0])**2)*\
				self._k(vec[1])*self.T(vec[1])*self._k([vec[1],vec[3]])*self.T([vec[1],vec[3]])*self._k([vec[1],vec[3],vec[4]])*self.T([vec[1],vec[3],vec[4]])+\
				(1/3.0*self._k(vec[1])**2-self._k([vec[1],vec[3]])**2)*self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*self.T(vec[1])*\
				self._k([vec[0],vec[2]])*self.T([vec[0],vec[2]])+0.25*self._k(vec[0])**2*self._k(vec[1])**2*self._k([vec[1],vec[4],vec[5]])*\
				self.T([vec[1],vec[4],vec[5]])+0.5*self._k(vec[0])**2*self._k([vec[0],vec[2],vec[3]])**2*self._k(vec[1])*self.T(vec[1])\
				-1/12.0*self._k(vec[0])**4*self._k(vec[1])*self.T(vec[1]))
		else: h=0.0
		return np.float64(h) 
	def s(self,n,vec):
		if n == 3:
			s=0.0
		elif n == 4:
			s=-1/4*self.sigma*np.dot(self.k(vec[0]),self.k(vec[1]))*np.dot(self.k(vec[2]),self.k(vec[3]))
		elif n == 5:
			s=0.0
		elif n == 6:
			s=-1/8.0*self.sigma*np.dot(self.k(vec[0]),self.k(vec[1]))*np.dot(self.k(vec[2]),self.k(vec[3]))*np.dot(self.k(vec[4]),self.k(vec[5]))
		else: s=0.0
		return np.float64(s)
	def alpha(self,n,vec):
		h1=list(vec)
		h2=list(vec)
		h1[0],h2[0],h2[1]=-vec[0],vec[1],-vec[0]
		return 1/2*(self.h(n+1,h1)+self.h(n+1,h2))
	def beta(self,n,vec):
		beta=0
		for i in range(2,len(vec)):
			vector=vec[1:i+1]+[-vec[0]]+vec[i+1:]
			beta=beta+self.h(n+1,vector)
		return -1/2*beta
	def gamma(self,n,vec):
		gamma=0
		for i in range(1,len(vec)+1):
			vector=vec[1:i]+[-vec[0]]+vec[i:]
			gamma=gamma+self.s(n+1,vector)
		return -1/2*gamma

	def prod(self,x,vec,ab):
		#Prod function allows to compute product of arbitrary number of solution elements
		prod = np.float64(1)
		for v,a in zip(np.absolute(vec),ab):
			prod = prod*x[int(v/self.mode-1)][a]
		return prod
	'''
	def prod(self,x,vec,ab):
		#Prod function allows to compute product of arbitrary number of solution elements
		prod = np.prod([x[int(v/self.mode-1)][a] for v,a in zip(vec,ab)])
		return prod
	'''
	def a_prime(self,x,j,t):
		linear = np.float64(self._k(j)*self.T(j)*x[int(j/self.mode-1)][1])
		non_linear = np.float64(0)
		for n in range(2,self.max_n+1):
			for i in range(len(self.deltas[n][j])):
				c = list(self.deltas[n][j][i])
				non_linear = non_linear + self.constants[n][j][i][0]*np.float64(self.prod(x,c,[1]+[0 for l in range(len(c)-1)]))
		return linear+non_linear	

	def b_prime(self,x,j,t):
		linear = - np.float64(self.G(j,t)*x[int(j/self.mode-1)][0])
		non_linear = np.float64(0)
		for n in range(2,self.max_n+1):
			for i in range(len(self.deltas[n][j])):
				c1,c2 = list(self.deltas[n][j][i]),list(self.deltas[n][j][i])
				non_linear = non_linear + self.constants[n][j][i][1]*np.float64(self.prod(x,c1,[1,1]+[0 for l in range(len(c1)-2)]))
				non_linear = non_linear + self.constants[n][j][i][2]*np.float64(self.prod(x,c2,[0 for l in range(len(c2))]))
		return linear+non_linear

	def f(self,x,t):
		#Main differential function
		#TODO: Can we split this into j processes??
		x = [[x[n],x[n+1]] for n in range(0,len(x),2)]
		f = [[0,0] for n in range(len(self._js))]
		for j in range(len(f)):
			f[j]=[self.a_prime(x,self.mode*(j+1),t),self.b_prime(x,self.mode*(j+1),t)]
		return np.array(f).flatten()

	def hamiltonian(self,x,t):
		x = [[x[n],x[n+1]] for n in range(0,len(x),2)]
		hamiltonian = np.float64(0.0)
		for n in range(2,self.max_n+2):
			if n==2:
				for vec in self.deltas_h[n]:
					hamiltonian = hamiltonian + (self.g0-self.sigma*np.dot(self.k(vec[0]),self.k(vec[1])))*self.prod(x,vec,[0,0]) +\
					 self._k(vec[0])*self.T(vec[0])*self.prod(x,vec,[1,1])
			else:
				for vec in self.deltas_h[n]:
					hamiltonian = hamiltonian + self.h(n,vec)*self.prod(x,vec,[1,1]+[0 for l in range(len(vec)-2)]) +\
					self.s(n,vec)*self.prod(x,vec,[0 for l in range(len(vec))])

		return 0.5*hamiltonian

	def solve(self,dt=0.01,N=1000,x0=np.array(None),
		forcing={'epsilon':0,'delta':0,'j':1},h_freq=100,method='homemade'):
		#Solved by employing RK4 method
		self.force_eps = forcing['epsilon']
		self.force_del = forcing['delta']
		if (forcing['epsilon'] and forcing['delta'])!=0:
			self.force_j = forcing['j']
			self.force_omega = self.omega(self.force_j)/sqrt(self.force_del)
			self.force_amp = self.force_eps*self.force_omega**2/(self._k(self.force_j)*self.T(self.force_j))

		ts = np.linspace(0.0,(N-1)*dt,num=N,dtype='float64')

		if method=='homemade':
			if x0.all()==None: 
				x0 = [0.001,0.0]+[0.0 for p in range(2*len(self._js)-2)]
				x0=np.float64(x0)
			hamiltonian = [self.hamiltonian(x0,ts[0])]
			sol = [x0]

			for q in range(1,N):

				ti = ts[q]
				k1 = dt*self.f(sol[q-1],ti)
				k2 = dt*self.f(sol[q-1]+0.5*k1,ti+dt*0.5)
				k3 = dt*self.f(sol[q-1]+0.5*k2,ti+dt*0.5)
				k4 = dt*self.f(sol[q-1]+k3,ti+dt)

				xi = sol[q-1] + 1/6.0*(k1+2.0*k2+2.0*k3+k4)
				sol.append(xi)
				if q%h_freq==0:
					hamiltonian.append(self.hamiltonian(xi,ti))
				else:
					hamiltonian.append(hamiltonian[q-1])

			return {'y':sol,'t':ts,'h':hamiltonian}



















