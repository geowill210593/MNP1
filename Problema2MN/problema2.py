# Planck’s law of black body radiation
from math import *
from roots import *
import matplotlib.pyplot as plt

hP = 6.62606957e-34
kB = 1.3806488e-23
c = 2.99792458e8

def u(lam):
	global lam0, K

	if (lam == 0e0): return 0e0
	x = lam / lam0
	return K / (pow(x,5) * (exp(1e0/x) - 1e0))

def du(lam):
	global lam0, K

	if (lam == 0e0): return 0e0
	x = lam / lam0
	expx = exp(1e0/x)
	exp1 = expx - 1e0
	return K * lam0 * (5e0*x + (1e0-5e0*x)*expx) / (pow(x,7) * exp1*exp1)

lam_plot = 3e-6
h = 1e-8
n = int(lam_plot/h + 0.5) + 1
nmax = 3 * n

x = [0]*(nmax+1)
y = [0]*(nmax+1)
z = [0]*(nmax+1)

for j in range(0,3):
	T = 4000e0 + j * 1000e0
	lam0 = hP*c / (kB*T)
	K = 8e0 * pi * kB * T / pow(lam0,4)
	for i in range(1,n+1):
		lam = (i-1)*h
		x[i+j*n] = lam * 1e6
		y[i+j*n] = u(lam)
		z[i+j*n] = du(lam)

fmax = 0e0

for i in range(1,nmax+1):
	if (y[i] > fmax): fmax = y[i]
for i in range(1,nmax+1): y[i] /= fmax

fmax = 0e0

for i in range(1,nmax+1):
	if (z[i] > fmax): fmax = z[i]
for i in range(1,nmax+1): z[i] /= fmax

nn = [0, n, 2*n, 3*n]
col = ["", "blue", "green", "red"]
sty = [0, 1, 1, 1]

plt.plot(x,y, color = 'green')
plt.axis([0,3,0,1.1])
plt.show()

plt.plot(x,z, color = 'red')
plt.axis([0,3,-0.31,1.1])
plt.show()

T = 5778e0

lam0 = hP*c / (kB*T)

K = 8e0 * pi * kB * T / pow(lam0,4)

fmin = 1e10 
fmax = -1e10
for i in range(1,n+1):
	lam = (i-1)*h
	f = du(lam)
	if (f > fmax): fmax = f; lam_1 = lam
	if (f < fmin): fmin = f; lam_2 = lam
lam_max = 0.5e0 * (lam_1 + lam_2)
(lam_max,ierr) = NewtonNumDrv(du,lam_1,lam_2,lam_max)
print("Wien’s law:")
print("T = {0:4.0f} K".format(T))
print("lambda_max = {0:e} m".format(lam_max))
print("lambda_max * T = {0:e} m K".format(lam_max * T))
