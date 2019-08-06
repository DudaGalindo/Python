import numpy as np
import math
from numpy.linalg import inv, det
import unittest
from ..rsolve import rSolve

class Test_chapter4(unittest.TestCase):
    def test1_rootsearch(self):
        def f(x): return x**3 - 10.*x**2 + 5.
        x1 = 0.
        x2 = 1.
        for i in range(4): #4 stages, each stage having 10 search intervals
            dx = (x2 - x1)/10.
            x1,x2 = rSolve.root_search(f,x1,x2,dx)
        x = (x1+x2)/2.
        self.assertAlmostEqual(x,0.7346,4, 'ValueError: Failed')

    def test2_bisection(self):
        def f(x): return x**3 - 10.*x**2 + 5.
        x = rSolve.bisection(f,0.0,1.,1, 1.e-4)
        self.assertAlmostEqual(x,0.7346,4, 'ValueError: Failed')

    def test3_bisection(self):
        def g(x): return x - math.tan(x)
        a = 0.
        b = 20.
        dx = 0.01
        x1 = 0.
        i = 0
        while x1 != None:
            x1,x2 = rSolve.root_search(g,a,b,dx) #retorna o intervalo das respostas
            a = x2
            root = rSolve.bisection(g,x1,x2,1,1e-4)
            if root!= None:
                i = i + 1
                if i == 1: self.assertAlmostEqual(root,0.0,4, 'ValueError: Failed')
                if i == 2: self.assertAlmostEqual(root,4.493409458100745,4, 'ValueError: Failed')
                if i == 3: self.assertAlmostEqual(root,7.725251837074637,4, 'ValueError: Failed')
                if i == 4: self.assertAlmostEqual(root,10.904121659695917,4, 'ValueError: Failed')
                if i == 5: self.assertAlmostEqual(root,14.06619391292308,4, 'ValueError: Failed')
                if i == 6: self.assertAlmostEqual(root,17.220755272209537,4, 'ValueError: Failed')
            else: break

    def test4_ridder(self):
        def f(x): return x**3 - 10*x**2 + 5
        a = 0.6
        b = 0.8
        x = rSolve.ridder(f,a,b,1e-9)
        self.assertAlmostEqual(x,0.7346,4, 'ValueError: Failed')

    def test5_ridder(self):
        def f(x): return 1/((x - 0.3)**2 + 0.01) - 1/((x - 0.8)**2 + 0.04)
        a = 0.
        b = 1.
        x = rSolve.ridder(f,a,b,1e-9)
        self.assertAlmostEqual(x,0.5800,4, 'ValueError: Failed')

    def test6_newtonRaphson(self):
        def f(x): return x**2 - 2
        def df(x): return 2*x
        a = 1.0
        b = 1.5
        x,i = rSolve.newton_raphson(f,df,a,b,1e-9,1)
        self.assertAlmostEqual(x,1.4142135,5, 'ValueError: Failed')

    def test7_newtonRaphson_for_a_double_root(self):
        def f(x): return x**4 - 6.4*x**3 + 6.45*x**2 + 20.538*x - 31.752
        def df(x): return 4.*x**3 - 19.2*x**2 + 12.9*x + 20.538
        a = 2.
        b = 0 #don't need it is this case
        x,i = rSolve.newton_raphson(f,df,a,b,1e-9,2)
        self.assertAlmostEqual(x,2.0999999786199406,5,'ValueError: Failed' )

    def test8_newtonRaphson_system_two_equations(self):
        def f(x):
            f = np.zeros(len(x))
            f[0] = x[0]**2 + x[1]**2 - 3
            f[1] = x[0]*x[1] - 1
            return f
        x = np.array([0.5,1.5])
        sol = rSolve.newton_raphson_system(f,x,1e-9)
        x_sol = np.array([0.61803, 1.61803])
        for i in range (len(x)):
            self.assertAlmostEqual(sol[i],x_sol[i],5,'ValueError: Failed')

    def test9_newtonRaphson_system_three_equations(self):
        def f(x):
            f = np.zeros(len(x))
            f[0] = math.sin(x[0]) + x[1]**2 + math.log(x[2])-7
            f[1] = 3*x[0] + 2**(x[1]) - x[2]**3 + 1
            f[2] = x[0] + x[1] + x[2] - 5
            return f
        x = np.array([1.,1.,1.])
        sol = rSolve.newton_raphson_system(f,x,1e-9)
        x_sol = np.array([0.59905376, 2.3959314, 2.00501484])
        for i in range (len(x)):
            self.assertAlmostEqual(sol[i],x_sol[i],5,'ValueError: Failed')

    def test10_polynomial_roots(self):
        a = np.array([-250.0,155.0,-9.0,-5.0,1.0]) #starting value
        roots = rSolve.polyRoots(a,1e-9)
        sol = np.array([2.+0.j, 4.-3.j, 4.+3.j, -5.+0.j])
        #p,dp,ddp = rSolve.evalPoly(a,roots)
        for i in range(len(roots)):
            self.assertAlmostEqual(roots[i],sol[i],9, 'ValueError: Failed')
