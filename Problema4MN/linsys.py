from math import *

def LUFactor(a, ipivot, n):

   det = 1e0
   for j in range(1,n+1):
# loop over columns
      for i in range(1,j):
# elements of matrix U
         sum = a[i][j]
         for k in range(1,i): sum -= a[i][k]*a[k][j]
         a[i][j] = sum

      amax = 0e0
      for i in range(j,n+1):
         sum = a[i][j]
         for k in range(1,j): sum -= a[i][k]*a[k][j]
         a[i][j] = sum
# elements of matrix L
# undivided by pivot
# determine pivot
         if (amax < fabs(a[i][j])): amax = fabs(a[i][j]); imax = i

      if (amax == 0e0): print("LUFactor: singular matrix !"); return 0e0

      ipivot[j] = imax
# store pivot row index
# interchange rows imax and j
# to put pivot on diagonal
      if (imax != j):
         det = -det
         for k in range(1,n+1):
            t = a[imax][k]
            a[imax][k] = a[j][k]
            a[j][k] = t
      det *= a[j][j]
# multiply determinant with pivot
      t = 1e0/a[j][j]
      for i in range(j+1,n+1): a[i][j] *= t
# divide elements of L by pivot
   return det

#============================================================================
def LUSystem(a, ipivot, b, n):

   for i in range(1,n+1):
# solves Ly = b
      sum = b[ipivot[i]]
      b[ipivot[i]] = b[i]
      for j in range(1,i): sum -= a[i][j]*b[j]
      b[i] = sum
   
   for i in range(n,0,-1):
      sum = b[i]
      for j in range(i+1,n+1): sum -= a[i][j]*b[j]
      b[i] = sum/a[i][i]

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

#============================================================================
def GaussJordan(a, b, n, m):

   ipivot = [0]*(n+1); jpivot = [0]*(n+1)
# stores pivot rows and columns
   npivot = [0]*(n+1)
# marks used pivot columns
   det = 1e0
   for k in range(1,n+1):
# FORWARD ELIMINATION
      amax = 0e0
# determine pivot
      for i in range(1,n+1):
# loop over rows
         if (npivot[i] != 1):
            for j in range(1,n+1):
# loop over columns
               if (npivot[j] == 0):
# pivoting not yet done?
                  if (amax < fabs(a[i][j])):
                     amax = fabs(a[i][j])
                     imax = i
                     jmax = j
               if (npivot[j] > 1):
                  print("GaussJordan: singular matrix 1 !"); return 0e0
      if (amax == 0e0): print("GaussJordan: singular matrix 2 !"); return 0e0
      ipivot[k] = imax; jpivot[k] = jmax
      npivot[jmax] = npivot[jmax] + 1
# store pivot row and column
# mark used pivot column
      if (imax != jmax):
# interchange rows imax and jmax
         det = -det
# to put pivot on diagonal
         for j in range(1,n+1):
            (a[imax][j],a[jmax][j]) = (a[jmax][j],a[imax][j])
         for j in range(1,m+1):
            (b[imax][j],b[jmax][j]) = (b[jmax][j],b[imax][j])
      det *= a[jmax][jmax]
# multiply determinant with pivot
      t = 1e0/a[jmax][jmax]
# divide pivot row by pivot
      a[jmax][jmax] = 1e0
# diagonal element of unit matrix
      for j in range(1,n+1): a[jmax][j] *= t
      for j in range(1,m+1): b[jmax][j] *= t
      
      for i in range(1,n+1):
# reduce non-pivot rows
         if (i != jmax):
            t = a[i][jmax]
            a[i][jmax] = 0e0
# non-diagonal element of unit matrix
            for j in range(1,n+1): a[i][j] -= a[jmax][j]*t
            for j in range(1,m+1): b[i][j] -= b[jmax][j]*t
   for k in range(n,0,-1):
# rearrange columns of inverse
      imax = ipivot[k]; jmax = jpivot[k]
      if (imax != jmax):
         for i in range(1,n+1):
            (a[i][imax],a[i][jmax]) = (a[i][jmax],a[i][imax])

   return det
