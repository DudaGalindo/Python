import numpy as np
from numpy.linalg import inv, det
from ..msolve import Msolve
import unittest
#OPTMIZE CODE AND TRY TO USE ONE GAUSS SEIDEL FUNCTION FOR ALL TESTS
class Test_chapter2(unittest.TestCase):
    def test1_Gauss_Elimin(self):
        A = np.array([[6.,-4.,1.],[-4.,6.,-4.],[1.,-4.,6.]])
        b = np.array([[-14.,22.],[36.,-18.],[6.,7.]])
        sol_x = np.array([[10.,3.],[22.,-1.],[14.,0.]])
        #A_ = array([[4.,-2.,1.],[-2.,4.,-2.],[1.,-2.,4.]]);
        #b_ = array([[11.],[-16.],[17.]])
        #sol_x = np.array([[1.],[-2.],[3.]])
        # n = np.size(x)
        #for j in range(n):
            #self.assertEqual(x[i],sol_x[i],'ValueError: Failed') # ver se tem jeito mais simples de comparar vetor
        A,b = Msolve.Gauss_Elimination(A,b)
        x = np.round(Msolve.Backward_Elimination(A,b))
        n = np.size(x[:,0])
        m = np.size(x[0,:])
        for i in range(n):
            for j in range(m):
                self.assertEqual(x[i,j],sol_x[i,j],'ValueError: Failed') # ver se tem jeito mais simples de comparar vetor

    def test2_condition_number(self):
        #a = array([[-2,-1,0],[-1,2,-1],[0,-1,2]])
        a = np.array([[2.1,-0.6,1.1],[3.2,4.7,-0.8],[3.1,-6.5,4.1]])
        ans = Msolve.condition(a)
        sol_ans = "Singular Matrix"
        self.assertEqual(ans,sol_ans)

    def test3_Doolitle(self):
        LU = np.array([[1.,4.,1.],[1.,6.,-1.],[2.,-1.,2.]])
        b = np.array([[7.],[13.],[5.]])
        J = Msolve.Doolittle(LU)
        y = Msolve.Foward_Elimination(J,b)
        x = Msolve.Backward_Elimination(J,y)
        sol_x = np.array([5.,1.,-2.])
        n = np.size(x)
        for i in range(0,n):
            self.assertEqual(x[i],sol_x[i],'ValueError: Failed')

    def test4_Choleski(self):
        A = np.array([[1.44,-0.36,5.52,0.00],[-0.36,10.33,-7.78,0.00],[5.52,-7.78,28.40,9.00],[0.00,0.00,9.00,61.]])
        b = np.array([[0.04],[-2.15],[0.],[0.88]])
        #A = array([[4.,-2.,2.],[-2.,2.,-4.],[2.,-4.,11.]])
        L = Msolve.Choleski(A)
        n = np.size(L[0,:])
        x = np.round(Msolve.solveCholeski(L,b))
        sol_x = np.array([[3.09212567],[-0.73871706],[-0.8475723],[0.13947788]])
        for i in range(n):
            self.assertEqual(x[i],np.round(sol_x[i]),'ValueError: Failed')

    def test5_identify_matrix_and_found_A(self):  #incomplete test, meaning that works just for this specific entry matrix
        S = np.array([[4.,-2.,1.,0.],[0.,3.,-1.5,1.],[0.,0.,3.,-1.5],[0.,0.,0.,2.9166666666666667]])
        M,a = Msolve.identify_matrix(S)
        if a == 'Its a U matrix':  # a its a statement
          L = Msolve.finding_L(M)
          A = Msolve.finding_A(L,M)
        sol_A = np.array([[4.,-2.,1.,0.],[-2.,4.,-2.,1.],[1.,-2.,4.,-2],[0.,1.,-2.,4.]])
        n = np.size(A[0])
        for i in range(n):
            for j in range(n):
                self.assertEqual(np.round(A[i,j]),np.round(sol_A[i,j]),'ValueError: Failed')

    def test6_Doolittle_LDLt(self):
        A = np.array([[3.,-3.,3.],[-3.,5.,1.],[3.,1.,10.]])
        A = Msolve.Doolittle(A)
        L,U = Msolve.A_to_LU(A)
        Lt = np.transpose(L)
        D = np.matmul(U,np.linalg.inv(Lt))
        sol_D = np.array([[3.,0.,0.],[0.,2.,0.],[0.,0.,-1.]])
        n = np.size(D[0])
        for i in range(0,n):
            for j in range(0,n):
                self.assertEqual(D[i,j],sol_D[i,j],'ValueError: Failed')

    def test7_Gauss_pivoting(self):
        A = np.array([[2.,-2.,6.],[-2.,4.,3.],[-1.,8.,4.]])
        b = np.array([[16.],[0.],[-1.]])
        #A = array([[4.,-2.,1.],[-2.,1.,-1.],[-2.,3.,6.]])
        #b = array([[2.],[-1.],[0.]])
        J,R = Msolve.Gauss_pivoting(A,b)
        x = Msolve.Backward_Elimination(J,R)
        sol_x = np.array([1.,-1.,2.])
        n = np.size(x)
        for i in range(n):
            self.assertEqual(x[i],sol_x[i],'ValueError: Failed') # ver se tem jeito mais simples de comparar vetor

    def test8_inverse_matrix(self):
        A = np.array([[0.6,-0.4,1.],[-0.3,0.2,0.5],[0.6,-1.,0.5]])
        Aorigin = np.copy(A)
        n = np.size(A[1])
        I = np.identity(n)
        Iorigin = np.copy(I)
        Ainv = Msolve.mat_inv(A,I)
        I_sol = (np.matmul(Aorigin,Ainv))
        I_sol = np.round(I_sol)
        for i in range(n):
            for j in range(n):
                self.assertEqual(Iorigin[i,j],I_sol[i,j],'ValueError: Failed')

    def test9_Gauss_Seidel_w1(self): #n ta rodando
        A = np.array([[4.,-1.,1.],[-1.,4.,-2.],[1.,-2.,4.]])
        b = np.array([[12.],[-1.],[5.]])
        n = np.size(b)
        z,x,w = Msolve.Gauss_Seidel_for_normal_matrix(A,b,1.0e-9,500,1,n) #k is the number of iterations, in case the user want to know
        sol_x = np.array([3.,1.,1.])
        n = np.size(A[0])
        for i in range(n):
            self.assertEqual(round(x[i]),sol_x[i],'ValueError: Failed')

    def test10_Gauss_Seidel_wnot1_tridiagonal(self):
        n = 20
        d = np.ones((n))*2.
        c = np.ones((n-1))*(-1.)
        b = np.zeros((n))
        b[n-1] = 1.
        e = c.copy()
        x = np.zeros(n)
        x,k,w = Msolve.gaussSeidel_wNOT1_not_regular_matrix(Msolve.iterEqs3,x,1.e-9,n,500)
        sol_k = 259
        sol_w = 1.705
        sol_x = np.array([-4.5,-4.,-3.50,-3.,-2.50,-2.,-1.5,-9.99999997e-01,\
    -4.99999998e-01,2.14047151e-09,5.00000002e-01,1.00000000e+00,1.50000000e+00,\
    2.00000000e+00,2.50000000e+00,3.00000000e+00,3.50000000e+00,4.,4.5,5.])
        for i in range(0,n):
            self.assertAlmostEqual(round(x[i],5),round(sol_x[i],5),4,'ValueError: Failed')
        self.assertAlmostEqual(round(w,3),sol_w,4,'ValueError: Failed')
        self.assertEqual(k,sol_k,'ValueError: Failed')


        #k,x = Msolve.Gauss_Seidel(A,b,1e-9)




    #def test9_tridiagonal(self): NOT WORKING - DONT KNOW WHY
    #    d = np.ones((5))*2.0
    #    c = np.ones((4))*(-1.0)
    #    b = np.array([5.0, -5.0, 4.0, -5.0, 5.0])
    #    e = np.copy(c)
    #    c,d,e = Msolve.LUdecomp3(c,d,e)
    #    x = Msolve.LUsolve3(c,d,e,b)
    #    print(x)
    #    sol_x = np.array([2., -1., 1., -1., 2.])
    #    n = len(x)
    #    for i in range(n):
    #        self.assertEqual(x[i],sol_x[i], 'ValueError: Failed')
