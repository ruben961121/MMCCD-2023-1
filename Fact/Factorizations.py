from MCD.Fact.AuxFunc import ind_max,permut,comp_cuad,comp_sim,sum_cho_1,sum_cho_2
import numpy as np

## Código principal para la factorización PLU
# La función tedrá como única entrada a la matriz
def Crout(A):
    
    A = np.array(A)
    if comp_cuad(A): #Verificación de ser cuadrada
        
        n = len(A[0])
        #Fijamos valores iniciales de PLU
        P = np.eye(n)
        L = np.zeros((n,n))
        U = A
        
        for j in range(n-1): #Ciclo para recorrer col
            
            index = ind_max(U[:,j],j)[0] #índice de máximo
            
            if ind_max(U[:,j],j)[1]==0:
                return print("La matriz es singular") #Cond. de singularidad
            
            else:
            
                P = permut(P,j,index)
                L = permut(L,j,index) #permutación de renglones
                U = permut(U,j,index)
            
                for i in range(j+1,n): # ciclo de renglones
                
                    L[i][j] = (U[i][j]/U[j][j]) # Llenado de L
                
                    U[i] = U[i] - L[i][j]*U[j] 
                    # Eliminación debajo de la diagona en col j
                
        L += np.eye(n) # Suma de Id a L
        P = np.transpose(P) # Transponer a P
        
        return P,L,U
        
    else:
        print("La matriz no es cuadrada")
        
        

def Cholesky(M):
    M=np.array(M)
    if comp_cuad(M):
        if comp_sim(M):
        
            n = len(M)
            L = np.zeros((n,n))
            LT = np.zeros((n,n))
        
            for j in range(n):
                for i in range(j,n):
                    
                    if i==j:
                        aux = M[i][i] - sum_cho_1(L,i)
                        if aux>0:
                            L[i][i] = np.sqrt(aux)
                        else:
                            return print("La matriz no es positivo-definida")
                    else:
                        L[i][j] = (1/L[j][j])*(M[i][j] - sum_cho_2(L,i,j))
            
            LT = L.transpose()
            
            return L,LT
                    
        else:
            print("La matriz no es simétrica")
    else:
        return print("La matriz no es cuadrada")
    
def QR(A):
    A = np.array(A)
    m,n=A.shape
    Q = np.zeros((m,n))
    R = np.zeros((n,n))
    
    for j in range(n):
        Q[:,j] = A[:,j] 
        for i in range(j):
            R[i][j] = Q[:,i]@A[:,j]
            Q[:,j] = Q[:,j] - R[i][j]*Q[:,i]
        
        R[j][j] = np.linalg.norm(Q[:,j])
        if R[j][j]==0:
            return print("La matriz no es de rango completo")
        else:
            Q[:,j] = Q[:,j]/R[j][j]
        
    return Q,R