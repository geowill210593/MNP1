import numpy 
from gauss_elim import gaussElimin


a = numpy.array([[8, -6,  2],
			    [-4, 11, -7],
			    [4,  -7,  6]
			   ])

b = numpy.array([28, -40, 33])

x = gaussElimin(a,b)

print(x)