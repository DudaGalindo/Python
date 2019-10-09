import numpy as np

class Kelem:
    def Kel(n,A,E,he):
        if n == 3:
            k = 1/3*np.array([[7, -8, 1],[-8, 16, -8],[1, -8, 7])
        if n == 2:
            k = np.array([[1 -1],[-1, 1]])
        Kel = A*E/he*k
        return Kel
class FatorForma:
    def N1(xe,xe1,xe2):return 1 - 3*((xe-xe1)/(xe2-xe1))**2 + 2*((xe-xe1)/(xe2-xe1))**3
    def N2(xe,xe1,xe2):return -3*(xe-xe1)*(1-(xe-xe1)/(xe2-xe1))**2
    def N3(xe,xe1,xe2):return 3*((xe-xe1)/(xe2-xe1))**2 - 2*((xe-xe1)/(xe2-xe1))**3
    def N4(xe,xe1,xe2):return -(xe-xe1)*(((xe-xe1)/(xe2-xe1))**2 - (xe-xe1)/(xe2-xe1))
