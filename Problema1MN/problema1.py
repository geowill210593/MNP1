from math import *
from roots import *
import matplotlib.pyplot as plt

pi2 = 2e0 * pi

def Func(E):
	global e, M
	f = E - e * sin(E) - M
	df = 1e0 - e *cos(E)
	return(f, df)

n = 101
x = [0] * (n + 1)
y = [0] * (n + 1)

e = 0.9672671e0
T = 75.96e0
t0 = 1986.1113e0

t = 1986.2491e0
M = pi2 / T * (t - t0)
Emax = 1e0
h = Emax / (n - 1)

for i in range (1,n+1):
	E = (i - 1) * h
	x[i] = E
	(y[i],df) = Func(E)

plt.plot(x,y)

plt.show()

(E,ierr) = Newton(Func, 0e0, Emax, E)
print("E = ", E, " rad at t = ", t, "years")

h = T / (n - 1)
x[1] = 1986.1113e0
y[1] = E = 0e0
for i in range(2,n + 1):
	t = t0 + (i - 1) * h
	M = pi2 / T * (t - t0)
	(E,ierr) = Newton(Func,0e0,pi2,E)
	x[i] = t
	y[i] = E

plt.plot(x,y)
plt.axis([1900,2100,0,7])

plt.show()