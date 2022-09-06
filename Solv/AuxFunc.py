import numpy as np

##Función que da la suma para solucionar
##Sistemas Triangulares Superiores
def sum_TrSup(A,X,k):
    A=np.array(A)
    X=np.array(X)
    #definimos una suma auxiliar
    if k==0:
        return 0
    else:
        suma=0
        for i in range(k):
            suma += A[k][i]*X[i]
        return suma
    

##Función que da la suma para solucionar
##Sistemas Triangulares Inferiores
def sum_TrInf(A,X,k):
    A=np.array(A)
    X=np.array(X)
    suma=0
    if k==0:
        return 0
    else:
        for i in range(k):
            suma+= A[len(X)-1-k][len(X)-1-i]*X[len(X)-1-i]
    return suma

##Función que sluciona sistemas
##Triangulares Inferiores
def solv_TrInf(A,b):
    
    A = np.array(A)
    b = np.array(b)
    X = np.zeros(len(b))
    
    for i in range(len(b)):
        X[i] = (1/A[i][i])*(b[i]-sum_TrInf(A,X,i))
    return X

##Función que sluciona sistemas
##Triangulares Superiores
def solv_TrSup(A,b):
    
    A = np.array(A)
    b = np.array(b)
    X = np.zeros(len(b))
    
    for i in range(len(b)):
        j=len(b)-i-1
        X[j] = (1/A[j][j])*(b[j]-sum_TrSup(A,X,i))
    return X

#código para verificar ceros en la diagonal
def diag_zero(A):
    for i in range(len(A[0])):
        if A[i][i]!=0:
            value=True
        else:
            return False
    return value

#código para calcular la dsitancia de dos vectores
def euc_dist(a,b):
    a=np.array(a)
    b=np.array(b)
    d = np.linalg.norm(a-b)
    
    return d