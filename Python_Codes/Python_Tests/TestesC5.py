import numpy as np
from Python_Codes.difsolve import differenciate
import unittest

class Test_chapter5(unittest.TestCase):
    def test1_derivatives_O2(self):
        x = [0.,0.1,0.2,0.3,0.4]
        f = [0.0000,0.0819,0.1341,0.1646,0.1797]
        x1 = 0.0
        x2 = 0.2
        n = 2
        xd = [x1,x2]
        h = abs(x[1] - x[0]) #the intervals h should be equal for any consecutives values of x
        df1 = np.zeros(n)
        df2 = np.zeros(n)
        df1 = differenciate.fowardDiff_secOrder(f,x,x1,n,h);
        df2 = differenciate.centeredDiff_secOrder(f,x,x2,n,h);
        sol_df1 = np.array([0.9675,-3.77])
        sol_df2 = np.array([0.4135,-2.17])
        for i in range (n):
            self.assertAlmostEqual(df1[i],sol_df1[i],4,'ValueError: Failed')
            self.assertAlmostEqual(df2[i],sol_df2[i],4,'ValueError: Failed')

    def test2_Richardson(self):
        h1 = 0.1;
        h2 = 0.2;
        x = [0,0.1,0.2,0.3,0.4]
        f = [0.0000,0.0819,0.1341,0.1646,0.1797]
        x1 = 0.0
        p = 2.
        n = 1
        G = differenciate.Richard_extrapolation(f,x,h1,h2,n,x1,p);
        self.assertEqual(G,0.99275,'ValueError: Failed')
