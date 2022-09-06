import numpy as np
import math as mt

## Función auxiliar que regresa el
## índice del elemento máximo en valor absoluto.
# La función tiene dos argumentos, el primero será
# la lista donde se buscará, y el segundo será el
# índice desde donde se comienza a buscar
def ind_max(v,i): 
    v = np.array(v)
    v = np.absolute(v)
    v1 = np.zeros(len(v)-i)
    for j in range(i,len(v)):
        v1[j-i] = v[j] 
    M = max(v1)
    ind = np.where(v1==M)[0][0]
    ind += i
    return ind,M
# La función devuelve el índice en el que se encuentra
# el valor buscado y el máximo


## Función auxilar que permite intercambiar
## renglones en una matriz
# La función tiene 3 argumentos, el primero es la
# matriz a la que se le hará la permutación
# el segundo y tercero son los índices de los
# renglones
def permut(A,i,j):
    A = np.array(A)
    B = A.copy()
    
    A[i],A[j] = B[j],B[i]
    
    return A
# La función regresa la matriz con las permutaciones


#Mini código para verificar si la matriz es cuadrada
def comp_cuad(matrix):
    if len(matrix[0]) == len(matrix[:,0]):
        return True
    else:
        return False

#Código para verificar la simetría de una matriz
def comp_sim(matrix):
    mt = matrix.transpose()
    if np.array_equal(mt,matrix):
        return True
    else:
        return False
    
    
    
##Suma de diagonal en Cholesky
def sum_cho_1(L,ind):
    suma=0.0
    if ind==0:
        return suma
    else:
        for k in range(ind):
            suma += L[ind][k]
    return suma


##Suma de elementos fuera de diagonal en Cholesky
def sum_cho_2(L,ind1,ind2):
    suma=0.0
    if ind2==0:
        return suma
    else:
        for k in range(ind2):
            suma += L[ind1][k]*L[ind2][k]
    return suma

