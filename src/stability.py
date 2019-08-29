import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import colors
from math import cos,sin,sqrt,pi


def fptype(E,A,B,C):
	if E<0.0:
		return 0 #FP does not exist
	else:
		if B**2-4*A*C<0.0:
			return 1 #FP has complex roots
		else:
			if C/A<0.0:
				return 2 #FP has 2 real, 2 imag roots
			else:
				if -B/A<0.0:
					return 3 #FP has 4 imag roots
				elif -B/A>0.0:
					return 4 #FP has 4 real roots
				else:
					return 5 #Other




nphi,nbeta = (300,300)

phis = [(i/nphi)*pi/2 for i in range(nphi)]
#phis = [pi/2]
betas = [(1-j/nbeta)*2 for j in range(2*nbeta)]
#betas=[0]

pp, bb = np.meshgrid(phis,betas)

fp = np.zeros((2*nbeta,nphi))


for pn,p in enumerate(phis):
	for bn,b in enumerate(betas):

		# #Origin

		# E=1
		# A=1
		# B=5*b**2-1
		# C=(b**2-cos(p)**2)*(4*b**2-sin(p)**2)

		#\Delta=0, q_j=0

		E= (2*b-sin(p))*(b-cos(p))
		A= 1
		B= 12*b**2-4*b*sin(p)-8*b*cos(p)-sin(p)**2+4*sin(p)*cos(p)
		C= (-4*b*cos(p)-4*b*sin(p)
			+2*sin(p)*cos(p))*(4*b**2-2*b*sin(p)-
			4*b*cos(p)+2*sin(p)*cos(p))

		# #\Delta=pi, q_j=0

		# E= (2*b+sin(p))*(b-cos(p))
		# A= 1
		# B= 12*b**2+4*b*sin(p)-8*b*cos(p)-sin(p)**2-4*sin(p)*cos(p)
		# C= (-4*b*cos(p)+4*b*sin(p)
		# 	-2*sin(p)*cos(p))*(4*b**2+2*b*sin(p)-
		# 	4*b*cos(p)-2*sin(p)*cos(p))

		# #\Delta=0, p_j=0

		# E= (2*b-sin(p))*(b+cos(p))
		# A= 1
		# B= 12*b**2-4*b*sin(p)+8*b*cos(p)-sin(p)**2-4*sin(p)*cos(p)
		# C= (4*b*cos(p)-4*b*sin(p)
		# 	-2*sin(p)*cos(p))*(4*b**2-2*b*sin(p)+
		# 	4*b*cos(p)-2*sin(p)*cos(p))

		# # #\Delta=pi, p_j=0

		# E= (2*b+sin(p))*(b+cos(p))
		# A= 1
		# B= 12*b**2+4*b*sin(p)+8*b*cos(p)-sin(p)**2+4*sin(p)*cos(p)
		# C= (4*b*cos(p)+4*b*sin(p)
		# 	+2*sin(p)*cos(p))*(4*b**2+2*b*sin(p)+
		# 	4*b*cos(p)+2*sin(p)*cos(p))

		# #\Delta=0, \phi=pi/2

		# E=0.5*b*(2*b+1)
		# A=1
		# B=12*b**2+4*b-1
		# C=8*b**2*(1+2*b)

		# #\Delta=pi, \phi=pi/2

		# E=0.5*b*(2*b-1)
		# A=1
		# B=12*b**2+4*b-1
		# C=8*b**2*(1-2*b)

		# #\beta=0

		# E=1
		# A=1
		# B=4*sin(p)*cos(p)-sin(p)**2
		# C=4*sin(p)**2*cos(p)**2



		fp[bn,pn] = int(fptype(E,A,B,C))

for pn in range(len(phis)):
	for bn in range(len(betas)-1):
		if fp[bn,pn]!=fp[bn+1,pn]:
			fp[bn,pn]=5


# cmap = colors.ListedColormap(['lightgrey','ghostwhite','gold',
# 	'aquamarine','salmon','grey'])

# cmap = colors.ListedColormap(['black','dimgrey','darkgrey',
# 	'lightgrey','whitesmoke','white'])

cmap = colors.ListedColormap(['white','white','white',
	'white','white','black'])

extent = (min(phis),max(phis),min(betas),max(betas))
bounds = [0,1,2,3,4,5,6]
norm = colors.BoundaryNorm(bounds,cmap.N)

img = plt.imshow(fp, cmap=cmap, norm=norm, extent=extent,aspect=0.5)

# cbar = plt.colorbar(img,cmap=cmap,norm=norm,boundaries=bounds,aspect=5)

# #legend

# cbar.ax.get_yaxis().set_ticks([])
# for j, lab in enumerate(['$DNE$','$Comp.$','$Re+Im$','$Im$','$Re$','$Other$']):
#     cbar.ax.text(.5, (2 * j + 1) / 12.0, lab, ha='center', va='center')
# cbar.ax.get_yaxis().labelpad = 15
# cbar.ax.set_ylabel('Type of Polynomial Roots', rotation=270)

plt.xlabel(r'$\phi$')
plt.ylabel(r'$\beta$')
#plt.yticks([0])
#plt.title('Origin')
plt.title(r'$\Delta=0$, $q_j=0$')
#plt.title(r'$\beta=0$')

plt.show()


