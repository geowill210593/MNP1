

def NewtonNumDrv(Func, a, b, x):

eps = 1e-10
# relative precision of root
itmax = 100
# max. no. of iterations
for it in range(1,itmax+1):
	f = Func(x)
	dx = eps*fabs(x) if x else eps
# derivation step
	df = (Func(x+dx)-f)/dx
# numerical derivative
	dx = -f/df if fabs(df) > eps else -f
# root correction
	x += dx
# new approximation
	if ((x < a) or (x > b)): return (x,1)
# [a,b] does not contain a root
	if (fabs(dx) <= eps*fabs(x)): return (x,0)
# check convergence
print("NewtonNumDrv: max. no. of iterations exceeded !"); return (x,2)