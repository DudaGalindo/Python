import numpy as np

class Stock:
    def Gauss_Elimination(J,R):
        n = np.size(J[1])
        for k in range(0,n-1):
            for i in range(k+1,n):
                if J[i,k]!= 0:
                    R[i] = R[i]-(J[i,k]/J[k,k])*R[k]
                    J[i] = J[i]-(J[i,k]/J[k,k])*J[k]
        return J,R
