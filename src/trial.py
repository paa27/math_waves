#Trial of standing waves
#Completed tests. Code agrees with matlab hard-coded model (up to 2nd order)

from standing import Standing
import pdb
import matplotlib.pyplot as plt
from math import pi,sqrt,ceil,cos
import numpy as np
from scipy.signal import hilbert
from scipy.io import loadmat

wave = Standing(mode=1,max_order=3,harmonics=3,
	params={'h':3000.0,'g0':1.0,'sigma':1.0,'k0':pi})

sol,ts = wave.solve(dt=0.01,N=1000)

xs = np.linspace(0,pi/wave.k0,100)
plt.figure()

for time in range(0,len(sol),ceil(len(sol)/10)):
	surf = np.zeros(len(xs))
	for harmonic in range(wave.harmonics):
		coeff = sol[time][2*harmonic]
		surf = surf + coeff*np.cos(wave.k0*wave._js[harmonic]*xs)
	plt.plot(xs,surf)

f,(ax1,ax2) = plt.subplots(2,1)
for harmonic in range(wave.harmonics):
	ax1.plot(ts,[sol[i][2*harmonic] for i in range(len(sol))],label=str(harmonic+1))
ax1.set_xlabel('t')
ax1.legend()
#ax2.plot(ts,hamiltonian)
plt.show()
'''


operations = []
for order in range(1,6):

	wave = Standing(mode=1,max_order=order,harmonics=order,
		params={'h':3000.0,'g0':1.0,'sigma':1.0,'k0':pi})

	my_dict={}

	for n in wave.deltas:
		n_list=[]
		for m in wave.deltas[n]:
			n_list.append(len(wave.deltas[n][m]))
		my_dict[n]=n_list
	ops = 0
	for n in range(2,wave.max_n+1):
		ops = ops+n*sum(my_dict[n])
	ops = ops*2*wave.harmonics
	operations.append(ops)

plt.plot(range(1,6),operations)
plt.xlabel('Order')
plt.ylabel('# of operations')
plt.show()
'''


'''
for k in ks:
	for order in range(1,5):
		wave = Standing(mode=1,max_order=order,harmonics=order,
			params={'h':3000.0,'g0':1.0,'sigma':1.0,'k0':k})

		sol,ts,hamiltonian = wave.solve(dt=0.01,N=ceil(2*628/(0.01*wave.omega(1))),
			forcing={'epsilon':0.01,'delta':0.25,'j':1})



		env = np.abs(hilbert(np.array([sol[i][0][0] for i in range(len(sol))])))
		x = ts[1000:2000]
		y = np.log(env)[1000:2000]
		X = x - x.mean()
		Y = y - y.mean()

		slope = (X.dot(Y))/(X.dot(X))

		sigma1 = slope/(wave.force_eps*wave.force_omega)

		f,(ax1,ax2) = plt.subplots(2,1)
		for harmonic in range(wave.harmonics):
			ax1.plot(ts,[sol[i][harmonic][0] for i in range(len(sol))],label=str(harmonic+1))
		ax1.plot(ts[1000:2000],env[1000:2000],label='envelope')
		ax1.set_xlabel('t')
		ax1.legend()
		ax2.plot(ts,hamiltonian)

		f.suptitle('epsilon='+str(wave.force_eps)+' , sigma1='+str(sigma1) + ' , order='+str(wave.max_n))
		f.savefig('../fig/'+str(round(k,3))+'_'+str(order)+'.png')

'''




