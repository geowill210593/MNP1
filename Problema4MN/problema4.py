from random import *
from linsys import *
from matutil import *

n = 100
a = [[0]*(n+1) for i in range(n+1)]
b = [[0]*(n+1) for i in range(n+1)]
c = [[0]*(n+1) for i in range(n+1)]
d = [[0]*(n+1) for i in range(n+1)]

for i in range(1,n+1):
	for j in range(1,n+1):
		a[i][j] = random()
		b[i][j] = random()

MatProd(a,b,c,n,n,n)
det = GaussJordan(c,d,n,0)

det = GaussJordan(a,d,n,0)
det = GaussJordan(b,d,n,0)
MatProd(b,a,d,n,n,n)

MatDiff(c,d,d,n,n)
print("(A B)^(-1) - B^(-1)A^(-1) (sample):\n")
MatPrint(d,5,5)

err = MatNorm(d,n,n)
print("\nMaximum error = ","{0:8.2e}".format(err))