from MCD.Fact import Factorizations as Fact
from MCD.Fact.AuxFunc import comp_cuad
from MCD.Sol import AuxFunc
from MCD.Sol.AuxFunc import solv_TrInf, solv_TrSup, diag_zero, euc_dist
import numpy as np


##Solucionador de Crout
def Crout(A,b):
    
    A = np.array(A)
    b = np.array(b)
    
    P = Fact.Crout(A)[0]
    L = Fact.Crout(A)[1] #Obtencion de la factorización
    U = Fact.Crout(A)[2]
    b = np.transpose(P)@b
    
    Y = np.zeros(len(b))#Vectores auxiliares para 
    X = np.zeros(len(b))#almacenar soluciones
    
    Y = solv_TrInf(L,b)
    X = solv_TrSup(U,Y)
    
    return X

##Solucionador Cholesky
def Cholesky(A,b):
    
    A = np.array(A)
    b = np.array(b)
    
    L = Fact.Cholesky(A)[0] #Obtencion de la factorización
    LT = Fact.Cholesky(A)[1]
    
    Y = np.zeros(len(b))#Vectores auxiliares para 
    X = np.zeros(len(b))#almacenar soluciones
    
    Y = solv_TrInf(L,b)
    X = solv_TrSup(LT,Y)
    
    return X

## Solucionador Jacobi
# Función del método de Jacobi
# La función tiene como argumentos a la matriz A
#  vector b, aproximación x_0, iteraciones n y
# tolerancia e
def Jacobi(A,b,x_0,n,e):
    
    A = np.array(A)
    b = np.array(b)
    m = len(A[0])
    i=0 #indice de iteracion
    d=0 #distancia entre vectores 
    
    x1 = np.array(x_0)#vectores auxiliares de aproximación
    x2 = np.zeros(m)
    
    if comp_cuad(A):
        if diag_zero(A):
    
            D = np.zeros((m,m))
            E = np.zeros((m,m)) #Matrices auxiliares
            T = np.zeros((m,m))
    
            for i in range(m): #Llenado de matrices E y D
                for j in range (m):
                    if i==j:
                        D[i][i] = (1/A[i][i])
                    else:
                        E[i][j] = A[i][j]
            
            T = np.matmul(D,E) # Tranformación T
            
            while i<n:
                
                x1=x2
                x2 = np.matmul(D,b) - np.matmul(T,x1)
                d = euc_dist(x2,x1)
                i += 1
                
                if d <= e:
                    print("Converge")
                    return x2,i
                elif i==n:
                    return print("No converge o se necesitan más iteraciones")
            
        else:
            print("A tiene al menos un cero en su diagonal.")
    else:
        return print("La matriz no es cuadrada.")
    
# Función del método de Gauss-Seidel
# La función tiene como argumentos a la matriz A
#  vector b, aproximación x_0, iteraciones n y
# tolerancia e
def GS(A,b,x_0,n,e):
    
    A = np.array(A)
    b = np.array(b)
    m = len(A[0])
    i=0 #indice de iteracion
    d=0 #distancia entre vectores 
    
    x1 = np.array(x_0)#vectores auxiliares de aproximación
    x2 = np.zeros(m)
    x3 = np.zeros(m)
    
    if comp_cuad(A):
        if diag_zero(A):
    
            D = np.zeros((m,m))
            E = np.zeros((m,m)) #Matrices auxiliares
            T = np.zeros((m,m))
    
            for i in range(m): #Llenado de matrices E y D
                for j in range (m):
                    if i==j:
                        D[i][i] = (1/A[i][i])
                    else:
                        E[i][j] = A[i][j]
            
            T = np.matmul(D,E) # Tranformación T
            
            while i<n:
                
                x1=x2
                
                for j in range(m):
                    x3 = x2
                    x2[j] = np.matmul(D[j],b) + np.matmul(T[j],x3)
                
                d = euc_dist(x2,x1)
                i += 1
                
                if d <= e:
                    print("Converge")
                    return x2,i-1
                elif i==n:
                    return print("No converge o se necesitan más iteraciones")
            
        else:
            print("A tiene al menos un cero en su diagonal.")
    else:
        return print("La matriz no es cuadrada.")