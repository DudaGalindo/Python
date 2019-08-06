import numpy as np
from math import sqrt
from Python_Codes.swap import swapRows

class Msolve:
    def Gauss_Elimination(J,R):
        n = np.size(J[1])
        for k in range(0,n-1):
            for i in range(k+1,n):
                if J[i,k]!= 0:
                    R[i] = R[i]-(J[i,k]/J[k,k])*R[k]
                    J[i] = J[i]-(J[i,k]/J[k,k])*J[k]
        return J,R

    def Backward_Elimination(J,R):
        n = len(J[0])
        m = np.size(R[0])
        if m > 1:
            x = np.zeros((n,m))
            for i in range (0,m):
                x[n-1,i] = R[n-1,i]/J[n-1,n-1]
                for k in range(n-2,-1,-1):
                    x[k,i]=(R[k,i]-np.dot(x[k+1:n,i],J[k,k+1:n]))/J[k,k]
        else:
            x = np.zeros(n)
            x[n-1] = R[n-1]/J[n-1,n-1]
            for k in range(n-2,-1,-1):
                x[k] = (R[k]-np.dot(x[k+1:n],J[k,k+1:n]))/J[k,k]

        return x

    def Doolittle(J):
        n = np.size(J[1])
        for k in range(0,n-1):
            for i in range(k+1,n):
             if J[i,k]!=0:
                lamb = (J[i,k]/J[k,k])
                J[i,k+1:n] = J[i,k+1:n]-lamb*J[k,k+1:n]
                J[i,k] = lamb
        return J


    def  Foward_Elimination(J,R):
        n = np.size(J[1])
        m = np.size(R[0])
        y = np.zeros((n,m))
        for i in range(0,m):
            y[0,i]=R[0,i]/J[0,0]
            for k in range(1,n):
                y[k,i] = (R[k,i]-np.dot(y[0:k,i],J[k,0:k]))
        return y

    def condition(a):
        det_a = np.linalg.det(a)
        if round(det_a,10) == 0:
            ans = "Singular Matrix"
        else:
            inva = np.linalg.inv(a)
            det_inv_a = np.linalg.det(inva)
            condition_number = np.dot(det_a,det_inv_a)
            if round(condition_number) == 1:
                ans = "The matrix is well conditioned"
            else:
                ans = "The matrix is ill conditioned"
        return ans

    def Choleski(J):
        n = np.size(J[1])
        for k in range(n):
            if J[k,k] - np.dot(J[k,0:k],J[k,0:k]) < 0:
                ValueError('Matrix is not positive definite')
            else:
                J[k,k] = sqrt(J[k,k] - np.dot(J[k,0:k],J[k,0:k]))
            for i in range(k+1,n):
                J[i,k] = (J[i,k] - np.dot(J[i,0:k],J[k,0:k]))/J[k,k]
        for k in range(1,n):
            J[0:k,k] = 0.0
        return J

    def solveCholeski(J,R):
        n = np.size(R)
        for k in range(n):
            R[k] = (R[k] - np.dot(J[k,0:k],R[0:k]))/J[k,k]
        for k in range(n-1,-1,-1):
            R[k] = (R[k] - np.dot(J[k+1:n,k],R[k+1:n]))/J[k,k]
        return R

    def identify_matrix(A):
        n = np.size(A[1])
        for k in range(0,n-1):
            for i in range(k+1,n):
                if A[i,k] == 0:
                    Matrix = A
                    a = 'Its a U matrix'
                elif A[k,i] == 0:
                    Matrix = A
                    a = 'Its a L matrix'
                return Matrix,a

    def finding_L(U):
        n = np.size(U[1])
        L = np.zeros((n,n))
        for k in range(0,n-1):
            L[k,k] = 1
            for i in range(k+1,n):
                L[i,k] = U[k,i]/U[k,k]
        L[n-1,n-1] = 1
        return L

    def finding_A(L,U):
        n = np.size(U[1])
        A = np.zeros((n,n))
        A = np.matmul(L,U)
        return A

    def A_to_LU(A): # to be used after Doolittle's
        n = np.size(A[1])
        L = np.zeros((n,n))
        U = np.zeros((n,n))
        for k in range(0,n-1):
            L[k,k] = 1
            U[k,k] = A[k,k]
            for i in range(k+1,n):
                L[i,k] = A[i,k]
                U[k,i] = A[k,i]
        L[n-1,n-1] = 1
        U[n-1,n-1] = A[n-1,n-1]
        return L,U

    def Gauss_pivoting(A,b):
        n = np.size(A[0 ])
        m = np.zeros(n)
        for i in range(0,n):
            m[i] = max(abs(A[i,0:n]))
        for k in range(0,n-1):
            j = np.argmax(abs(A[k:n,k])/m[k:n])+k
            if j!= k:
                swapRows(m,k,j)
                swapRows(A,k,j)
                swapRows(b,k,j)
            for i in range(k+1,n):
                if A[i,k]!= 0:
                    b[i] = b[i]-(A[i,k]/A[k,k])*b[k]
                    A[i] = A[i]-(A[i,k]/A[k,k])*A[k]
        #b[n-1] = b[n-1]/A[n-1,n-1]
        #for k in range(n-2,-1,-1):
        #    b[k] = (b[k] - np.dot(A[k,k+1:n],b[k+1:n]))/A[k,k]
        return A,b

    def mat_inv(A,I):
        n = len(A[0])
        aInv = np.identity(n)
        A,I = Msolve.Gauss_pivoting(A,I)
        L,U = Msolve.A_to_LU(A)
        y = Msolve.Foward_Elimination(L,I)
        Ainv = Msolve.Backward_Elimination(U,y)
        return Ainv

    def LUdecomp3(c,d,e):
        n = len(d)
        for k in range(1,n):
            lam = c[k-1]/d[k-1]
            d[k] = d[k] - lam*e[k-1]
            c[k-1] = lam
        return c,d,e

    def LUsolve3(c,d,e,b):
        n = len(d)
        for k in range(1,n):
            b[k] = b[k] - c[k-1]*b[k-1]
            b[n-1] = b[n-1]/d[n-1]
        for k in range(n-2,-1,-1):
            b[k] = (b[k] - e[k]*b[k+1])/d[k]
        return b

    def LUdecomp5(d,e,f):
        n = len(d)
        for k in range(n-2):
            lam = e[k]/d[k]
            d[k+1] = d[k+1] - lam*e[k]
            e[k+1] = e[k+1] - lam*f[k]
            e[k] = lam
            lam = f[k]/d[k]
            d[k+2] = d[k+2] - lam*f[k]
            f[k] = lam
            lam = e[n-2]/d[n-2]
            d[n-1] = d[n-1] - lam*e[n-2]
            e[n-2] = lam
        return d,e,f

    def LUsolve5(d,e,f,b):
        n = len(d)
        b[1] = b[1] - e[0]*b[0]
        for k in range(2,n):
            b[k] = b[k] - e[k-1]*b[k-1] - f[k-2]*b[k-2]
            b[n-1] = b[n-1]/d[n-1]
            b[n-2] = b[n-2]/d[n-2] - e[n-2]*b[n-1]
        for k in range(n-3,-1,-1):
            b[k] = b[k]/d[k] - e[k]*b[k+1] - f[k]*b[k+2]
        return b

    def Gauss_Seidel_for_normal_matrix(A,b,tol,iter_max,relax,n):  ##it can be optimized
        k = 10
        p = 1
        w = 1
        x_old = np.zeros(n)
        x = np.ones(n)
        z = 0
        prod1 = 0
        for z in range(iter_max):
            for i in range(0,n):
                prod = 0
                dx = 0
                for j in range(0,n):
                    prod1 = 0
                    if j!=i:
                        prod1 = np.dot(A[i,j],x[j])
                    prod = prod1 + prod
                x[i] = w*(b[i] - prod)/A[i,i] + (1-w)*x[i]
                dx = x[i] - x_old[i]
                if abs(dx)<tol: return z,x,wf
                x_old[i] = np.copy(x[i])
                wf = np.copy(w)
                if i == k: dx1 = dx
                if i == k + p:
                    dx2 = dx
                    w = (2.0/(1.0 + sqrt(1.0 - (dx2/dx1)**(1.0/p))))*(relax-1) + relax
            z = z + 1;
        print('Gauss-Seidel failed to converge')

# Method below was implemented for tridiagonal matrix
    def gaussSeidel_wNOT1_not_regular_matrix(iterEqs_choose,x,tol,n,iter_max):
        omega = 1.0
        k = 10
        p = 1
        for i in range(1,iter_max):
            xOld = np.copy(x)
            x = iterEqs_choose(x,n,omega)
            dx = sqrt(np.dot(x-xOld,x-xOld))
            if dx < tol:
                return x,i,omega
            # Compute relaxation factor after k+p iterations
            if i == k: dx1 = dx
            if i == k + p:
                dx2 = dx
                omega = 2.0/(1.0 + sqrt(1.0 - (dx2/dx1)**(1.0/p)))
        print('Gauss_Seidel failed to converge')

    def iterEqs3(x,n,w):
        x[0] = w*(x[1] - x[n-1])/2. + (1. - w)*x[0]
        for i in range(1,n-1):
            x[i] = w*(x[i-1] + x[i+1])/2. + (1. - w)*x[i]
        x[n-1] = w*(1. - x[0] + x[n-2])/2. + (1. - w)*x[n-1]
        return x
