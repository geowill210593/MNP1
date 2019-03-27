from math import *
from roots import *

def Func(f, x, n):

	f[1] = pow(x[1]-1,2) + x[2]*x[2] - 16e0
	f[2] = x[1]*x[1] - x[2] - 3e0

n = 2
f = [0]*(n+1)
x = [0]*(n+1)
dx = [0]*(n+1)
x[1] = -5e0
x[2] = 5e0

(dx,ierr) = NewtonSys(Func,x,n)
Func(f,x,n)

print("Solution 1:")
print("x dx f")
for i in range(1,n+1):
	print("{0:d} {1:15.7e} {2:7.0e} {3:7.0}".format(i,x[i],dx[i],f[i]))

x[1] = 5e0
x[2] = 5e0
(dx,ierr) = NewtonSys(Func,x,n)
Func(f,x,n)

print("\nSolution 2:")
print("x dx f")
for i in range(1,n+1):
	print("{0:d} {1:15.7e} {2:7.0e} {3:7.0}".format(i,x[i],dx[i],f[i]))