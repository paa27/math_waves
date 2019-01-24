"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from standing import Standing
from math import sqrt

N=2000
dt = 0.01

wave = Standing(mode=1,max_order=2,harmonics=2,
	params={'h':30000.0,'g0':1.0,'sigma':1.0,'k0':np.pi})

sol,ts,h = wave.solve(dt=dt,N=N,x0=np.float64([0.01,0,0.01,0]))

sol = np.float64([[sol[i][n] for i in range(len(sol))]for n in range(2*wave.harmonics)])

a1 = sol[0]
a2 = sol[2]

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(20, 7))
fig.subplots_adjust(wspace=0.5)
ax1.plot(ts,a1)
ax1.set_xlabel('t',size=20)
ax1.set_ylabel('a',size=20)
ax1.set_title('j=1',size=20)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.tick_params(axis='both', which='minor', labelsize=10)

ax2.plot(ts,a2)
ax2.set_xlabel('t',size=20)
ax2.set_ylabel('a',size=20)
ax2.set_title('2j=2',size=20)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='minor', labelsize=10)

fig.suptitle(r'Fourier coeffs ($\zeta$)',size=20)

plt.show()


'''
fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(20, 7))
fig.subplots_adjust(wspace=0.5)
ax1.set_xlabel('x',size=20)
ax1.set_ylabel('y',size=20)
ax1.set_title('Total',size=20)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.tick_params(axis='both', which='minor', labelsize=10)

ax2.set_xlabel('x',size=20)
ax2.set_ylabel('y',size=20)
ax2.set_title('j=1',size=20)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='minor', labelsize=10)

ax3.set_xlabel('x',size=20)
ax3.set_ylabel('y',size=20)
ax3.set_title('2j=2',size=20)
ax3.tick_params(axis='both', which='major', labelsize=12)
ax3.tick_params(axis='both', which='minor', labelsize=10)

fig.suptitle(r'Free surface $\zeta$',size=20)

x = np.arange(0, np.pi/wave._k(1), 0.01)
line1, = ax1.plot(x,a1[0]*np.cos(wave._k(1)*x)+a2[0]*np.cos(wave._k(2)*x))
line2, = ax2.plot(x,a1[0]*np.cos(wave._k(1)*x))
line3, = ax3.plot(x,np.max(a2)*np.cos(wave._k(2)*x))


def animate(i):
	line1.set_ydata(a1[i]*np.cos(wave._k(1)*x)+a2[i]*np.cos(wave._k(2)*x))  # update the data
	line2.set_ydata(a1[i]*np.cos(wave._k(1)*x))  # update the data
	line3.set_ydata(a2[i]*np.cos(wave._k(2)*x))  # update the data
	return line1,line2,line3


# Init only required for blitting to give a clean slate.
def init():
	line1.set_ydata(np.ma.array(x, mask=True))
	line2.set_ydata(np.ma.array(x, mask=True))
	line3.set_ydata(np.ma.array(x, mask=True))
	return line1,line2,line3

ani = animation.FuncAnimation(fig, animate, np.arange(1, 20000, 20), init_func=init,
                              interval=20, blit=True,save_count=30000)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
ani.save('2modes21.mp4',writer=writer)
#plt.show()

'''



