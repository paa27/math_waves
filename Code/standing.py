#Standing waves in Python
import numpy as np
import pandas as pd
from math import *
import itertools as it

class Standing:
	def __init__(self,max_order=1,mode=1,harmonics=1,params={'h':1000.0,'g0':1.0,'sigma':1.0,'k0':pi}):
		self.max_n = max_order
		self._js = [j for j in range(mode,mode*(harmonics+1),mode)]
		self.js = [-j for j in range(mode*harmonics,0,-mode)]+[j for j in range(mode,mode*(harmonics+1),mode)]
		self._h = params['h']
		self.g0 = params['g0']
		self.sigma = params['sigma']
		self.k0 = params['k0']

		#Initialize non-zero terms in sums

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

		#Initialize constants (alphas,betas,gammas)

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


	def k(self,vec):
		return self.k0*vec
	def _k(self,vec):
		return np.absolute(self.k(vec)) #TODO: absolute function not what we want
	def T(self,vec):
		return tanh(self._k(vec)*self._h)
	def g(self,vec,t):
		return self.g0*(1*self.sigma/self.g0)*self._k(vec)**2 ## TODO: add time dependence
	def omega(self,vec):
		return sqrt(self.g(vec)*self._k(vec)*self.T(vec))

	def h(self,n,vec):
		if n == 3:
			h=-(np.dot(self.k(vec[0]),self.k(vec[1]))+ \
				self._k(vec[0])*self.T(vec[0])*self._k(vec[1])*self.T(vec[1]))
		else: h=0
		return np.float64(h) 
	def s(self,n,vec):
		if n == 3:
			s=np.float64(0)
		else: s=0
		return s
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
		for i in range(1,len(vec)):
			vector=vec[1:i]+[-vec[0]]+vec[i:]
			gamma=gamma+self.s(n+1,vector)
		return -1/2*gamma

	def prod(self,x,vec,ab):
		prod = np.float64(1)
		for v,a in zip(vec,ab):
			prod = prod*x[v-1][a]
		return prod

	def a_prime(self,x,j,t):
		linear = np.float64(self._k(j)*self.T(j)*x[j-1][1])
		non_linear = np.float64(0)
		for n in range(2,self.max_n+1):
			for i in range(len(self.deltas[n][j])):
				c = list(self.deltas[n][j][i])
				non_linear = non_linear + self.constants[n][j][i][0]*np.float64(self.prod(x,c,[1]+[0 for l in range(len(c)-1)]))
		return linear+non_linear

	def b_prime(self,x,j,t):
		linear = - np.float64(self.g(j,t)*x[j-1][0])
		non_linear = np.float64(0)
		for n in range(2,self.max_n+1):
			for i in range(len(self.deltas[n][j])):
				c1,c2 = list(self.deltas[n][j][i]),list(self.deltas[n][j][i])
				non_linear = non_linear + self.constants[n][j][i][1]*np.float64(self.prod(x,c1,[1,1]+[0 for l in range(len(c1)-2)]))
				non_linear = non_linear + self.constants[n][j][i][2]*np.float64(self.prod(x,c2,[0 for l in range(len(c2))]))
		return linear+non_linear

	def f(self,x,t):
		f = [[0,0] for n in range(len(self._js))]
		for j in range(len(f)):
			f[j]=[self.a_prime(x,j+1,t),self.b_prime(x,j+1,t)]
		return np.array(f)

	def solve(self,dt=0.01,N=200,x0=None):
		if x0==None: 
			x0=[[0.1,0.0]]+[[0,0] for p in range(len(self._js)-1)]
			x0=np.array(x0)
		sol = [x0]
		ts = np.linspace(0.0,(N-1)*dt,num=N)
		print(ts)

		for q in range(1,N):

			ti = ts[q]
			k1 = dt*self.f(sol[q-1],ti)
			k2 = dt*self.f(sol[q-1]+0.5*k1,ti+dt*0.5)
			k3 = dt*self.f(sol[q-1]+0.5*k2,ti+dt*0.5)
			k4 = dt*self.f(sol[q-1]+k3,ti+dt)

			xi = sol[q-1] + 1/6.0*(k1+2.0*k2+2.0*k3+k4)
			sol.append(xi)

		return sol
















