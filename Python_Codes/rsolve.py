import numpy as np
from Python_Codes.Error import error
import math
from Python_Codes.msolve import Msolve
from random import random
import cmath

class rSolve:
    def root_search(f,a,b,dx):
        x1 = a
        x2 = a + dx
        f1 = f(a)
        f2 = f(x2)
        while np.sign(f1) == np.sign(f2):
            if x1 >= b:
                return None, None
            x1 = x2
            f1 = f2
            x2 = x1 + dx
            f2 = f(x2)
        return x1,x2

    def bisection(f,x1,x2,switch,tol):
        f1 = f(x1)
        f2 = f(x2)
        if f1 == 0.: return x1
        if f2 == 0.: return x2
        if np.sign(f1) == np.sign(f2):
            error.err('Root is not bracketed')
        n = int(math.ceil(math.log(abs(x2-x1)/tol)/math.log(2.))) #the ceiling of n is the smallest integer of greater than n
        #n is the number of bissections required to reach determined tolerance
        for i in range(n):
            x3 = 0.5*(x1 + x2)
            f3 = f(x3)
            if (switch) and (abs(f3) > abs(f1)) and (abs(f3) > abs(f2)):
                return None
            if f3 == 0.: return x3
            if np.sign(f2) != np.sign(f3):
                x1 = x3; f1 = f3
            else: f2 = f3; x2 = x3
        return (x1 + x2)/2

    def ridder(f,x1,x2,tol):
        f1 = f(x1)
        f2 = f(x2)
        if f1 == 0.: return x1
        if f2 == 0.: return x2
        if np.sign(f1) == np.sign(f2):
            error.err('Root is not bracketed')
        for i in range(80): #critério colocado para convergir com tolerância de 10-9 em um exemplo (avacalhei talvez um pouco)
            x3 = (x1+x2)/2.; f3 = f(x3)
            if np.sign(f2) != np.sign(f3): x1 = x3; f1 = f3
            s = math.sqrt(f3**2 - f1*f2)
            if s == 0.0: return None
            if f1 > f2:
                x4 = x3 - (x3 - x1)*f3/s
            else:
                x4 = x3 + (x3 - x1)*f3/s
            f4 = f(x4)
            #testing for convergence:
            if i > 0:
                if abs(x4 - x4old) < tol*max(abs(x4),1.): #porque eu n sei -  acho que é pra tornar o erro proporcional à magnitude de x4
                    return x4
            x4old = x4
            if np.sign(f3) == np.sign(f4):
                if np.sign(f1) != np.sign(f4): x2 = x4; f2 = f4
                else: x1 = x4; f1 = f4
            else: x1 = x3; f1 = f3; x2 = x4; f2 = f4
        return None

    def newton_raphson(f,df,x1,x2,tol,m): #m is the multiplicity of the root(if there is a double root for example) - it speeds the convergence
        f1 = f(x1)
        f2 = f(x2)
        if f1 == 0.: return x1,0
        if f2 == 0.: return x2,0
        if m == 1 and np.sign(f1) == np.sign(f2):
            error.err('Root is not bracketed')
            print('ERROR')
        x = 0.5*(x1 + x2)
        for i in range(30):
            fx = f(x)
            if fx == 0.: return x,i
            if np.sign(f1) != np.sign(fx): x2 = x
            else: x1 = x
            dfx = df(x)
            try: dx = -fx/dfx
            except ZerodivisionError: dx = x2 - x1
            x = x + m*dx
            # if x is out of the brackets, wich means that x<a else x>b:
            if m == 1 and (x2 - x)*(x - x1) < 0.0:
                dx = 0.5*(x2 - x1)
                x = (x1 + x2)/2 # or just: x = (x1 + x2)/2.
            if abs(dx)<tol*max(abs(x),1.): return x,i
        return None,i

    def newton_raphson_system(f,x,tol):
        def Jacobian(f,x): # from the finite difference approximation
            h = 1.0e-4
            n = len(x)
            jac = np.zeros((n,n))
            f0 = f(x)
            for i in range(n):
                temp = x[i]
                x[i] = temp + h
                f1 = f(x)
                x[i] = temp
                jac[:,i] = (f1 - f0)/h
            return jac,f0

        for i in range(30):
            jac,f0 = Jacobian(f,x)
            if math.sqrt(np.dot(f0,f0)/len(x)) < tol: return x
            J,f0 = Msolve.Gauss_pivoting(jac,-f0)
            dx = Msolve.Backward_Elimination(J,f0)
            x = x + dx#resolver o que merda ta acontecendo
            if np.sqrt(np.dot(dx[0],dx[0])) < tol*max(np.max(np.absolute(x)),1.0):
                return x
        print('too many iterations') #too many iterations

    def evalPoly(a,x): #where a are the polynomial coefficients
        n = len(a)-1
        p = a[n]
        dp = 0. + 0.0j
        ddp = 0. + 0.0j
        for i in range(1,n+1):
            ddp = ddp*x + 2.0*dp
            dp = dp*x + p
            p = p*x + a[n-i]
        return p,dp,ddp

    def polyRoots(a,tol):
        def laguerre(a,tol):
            x = random()
            n = len(a) - 1
            for i in range(30):
                p,dp,ddp = rSolve.evalPoly(a,x)
                if abs(p) < tol: return x
                g = dp/p
                h = g*g - ddp/p
                f = cmath.sqrt((n - 1)*(n*h - g*g))
                if abs(g + f) > abs(g - f): dx = n/(g + f)
                else: dx = n/(g - f)
                x = x - dx
                if abs(dx) < tol: return x
            #print(’Too many iterations’)

        def deflPoly(a,root):
            # Deflates a polynomial b0 = a1...b(i) = a(i+1)
            n = len(a)-1
            b = [(0.0 + 0.0j)]*n
            b[n-1] = a[n]
            for i in range(n-2,-1,-1):
                b[i] = a[i+1] + root*b[i+1]
            return b

        n = len(a) - 1
        roots = np.zeros((n),dtype=complex)
        for i in range(n):
            x = laguerre(a,tol)
            if abs(x.imag) < tol: x = x.real
            roots[i] = x
            a = deflPoly(a,x)
        return roots
