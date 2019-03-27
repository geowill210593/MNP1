from math import *
from linsys import *


def Newton(Func, a, b, x):
   #----------------------------------------------------------------------------
   #  Determina una raíz real x de la función Func aislada en el intervalo [a, b] por
   #  el metodo de Newton-Raphson usando la derivada analítica. x es una aprox
   #  initial.
   #  Error code: 0 - ejecucion normal
   #              1 - intervalo no contiene una raiz
   #              2 - numero maximo de iteraciones excedido
   #----------------------------------------------------------------------------
   eps = 1e-10                                   # precision de la raiz
   itmax = 100                                       # max. no. de iteraciones

   for it in range(1,itmax+1):
      (f,df) = Func(x)                              # funcion y derivada
      dx = -f/df if fabs(df) > eps else -f                  # correccion de raiz
      x += dx                                             # new approximation
      if ((x < a) or (x > b)): return (x,1)   # [a,b] no tiene una raiz
      if (fabs(dx) <= eps*fabs(x)): return (x,0)          # verifica convergencia

   print("Newton: max. no. of iterations exceeded !"); return (x,2)



def Jacobian(x, jac, n, Func):

   #----------------------------------------------------------------------------
   #  Calcula el jacobiano de un sistema de n funciones reales con n variables.
   #  usando diferencias finitas centrales
   #  x[]     - Punto en el que se evalúa al jacobiano.
   #  jac[][] - Jacobiano
   #  n       - dimensión espacial
   #  Func    - Función de usuario que devuelve los valores de función en el punto x
   #               Func(f, x, n)
   #----------------------------------------------------------------------------
   eps = 1e-10

   fm = [0]*(n+1)
   fp = [0]*(n+1)

   for j in range(1,n+1):                             # bucle sobre coordenadas
      x0 = x[j]                                                  # store x[j]
      h = eps*fabs(x0) if x0 else eps                             # step-size
      x[j] = x0 - h; Func(fm,x,n)                            # decremento x[j]
      x[j] = x0 + h; Func(fp,x,n)                            # incremento x[j]
      h2 = 1e0/(2e0*h)
      for i in range(1,n+1): jac[i][j] = (fp[i] - fm[i]) * h2      # Jacobiano
      x[j] = x0                                                # restore x[j]


def NewtonSys(Func, x, n):
   #----------------------------------------------------------------------------
   #  Determina un cero real n-dimensional de un sistema de n funciones reales mediante
   #  el método Newton-Raphson.
   #  Func - user function returning the function values f[] for arguments x[]
   #            Func(f, x, n)
   #  x[]  -approximation inicial (input), solucion (output)
   #  dx[] - estimaciones de error de los componentes de la  solucion (output)
   #  n    - orden del sistema
   #  ierr - error code: 0 - ejecucion normal
   #                     1 - numero maximo de iteraciones excedido
   #  Calls: Jacobian - calcula Jacobiano
   #         MatInv   - invierte matriz de (n x n)
   #----------------------------------------------------------------------------
   eps = 1e-14                                # precisión para error acumulado
   itmax = 200                                       # max. no. de iteraciones

   f   = [0]*(n+1)
   dx  = [0]*(n+1)
   jac = [[0]*(n+1) for i in range(n+1)]

   for it in range(1,itmax+1):
      Func(f,x,n)                                                 # funciones
      Jacobian(x,jac,n,Func)                                       # Jacobiano
      det = MatInv(jac,n)                                  # Jacobiano inverso

      if (det):                                                 # correciones
         for i in range(1,n+1):                       # Jacobiano no-singular
            sum = 0e0
            for j in range(1,n+1): sum -= jac[i][j] * f[j]
            dx[i] = sum
      else:
         for i in range(1,n+1): dx[i] = -f[i]             # Jacobiano Singular

      err = 0e0
      for i in range(1,n+1):
         x[i] += dx[i]                                    # nueva aproximación
         err += fabs(f[i])                                   # error acumulado
      if (err <= eps): break                              # checa convergencia

   ierr = 0
   if (it >= itmax):
      ierr = 1; print("NewtonSys: max. no. of iterations exceeded !")
   return (dx,ierr)
