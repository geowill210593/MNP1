
from math import *

#============================================================================
def MatInv(a, n):
#----------------------------------------------------------------------------
#  Calcula la inversa de una matriz de (n x n) mediante la factorización LU.
#  a   - matriz (n x n), (input); a^(-1) (output)
#  det - determinante de la matriz de coeficiente (output).
#  Calls: LUFactor, LUSystem.
#----------------------------------------------------------------------------
   ainv = [[0]*(n+1) for i in range(n+1)]     # Almacenamiento temporal de la inversa
   b = [0]*(n+1)
   ipivot = [0]*(n+1)                                     # almacena filas de pivote

   det = LUFactor(a,ipivot,n)                         # factorizacion LU de
   if (det == 0e0):                                         # una matriz singular
      print("MatInv: singular matrix !"); return 0e0

   for j in range(1,n+1):                  # bucle sobre columnas de matriz unitaria
      for i in range(1,n+1): b[i] = 0e0                            # columna j
      b[j] = 1e0                          
      LUSystem(a,ipivot,b,n)                                   # resuelve el sistema
      for i in range(1,n+1): ainv[i][j] = b[i]          # columna j de inversa

   for j in range(1,n+1):                                 # copia la inversa en a
      for i in range(1,n+1): a[i][j] = ainv[i][j]

   return det


