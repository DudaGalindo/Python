import numpy as np


class differenciate:
    def fowardDiff_secOrder(f,x,xr,n,h): #n is the degree of the last derivative you want (ex. if n=2, you want the first and second derivatives)
        df = np.zeros(n)
        pos_xr = x.index(xr)
        pos_xrh = x.index(round(xr+h,2))
        pos_xr2h = x.index(round(xr+2.*h,2))
        if n >= 1:
            df[0] = (-3.*f[pos_xr] +  4.*f[pos_xrh] -  1.*f[pos_xr2h])/(2.*h)
            if n >= 2:
                pos_xr3h = x.index(round(xr+3.*h,2))
                df[1] = ( 2.*f[pos_xr] + -5.*f[pos_xrh] +  4.*f[pos_xr2h] -  1.*f[pos_xr3h])/(h**2)
                if n >= 3:
                    pos_xr4h = x.index(x+4.*h)
                    df[2] = (-5.*f[pos_xr] + 18.*f[pos_xrh] - 24.*f[pos_xr2h] + 14.*f[pos_xr3h] - 3.*f[pos_xr4h])/(2*h**3)
                    if n == 4:
                        pos_xr5h = x.index(x+5.*h)
                        df[3] = ( 3.*f[pos_xr] - 14.*f[pos_xrh] + 26.*f[pos_xr2h] - 24.*f[pos_xr3h] + 11.*f[pos_xr4h] - 2.*f[pos_xr5h])/(h**4)
        return df

    def centeredDiff_secOrder(f,x,xr,n,h):
        df = np.zeros(n)
        pos_xr = x.index(xr)
        pos_xr_h = x.index(round(xr-h,2))
        pos_xrh = x.index(round(xr+h,2))
        pos_xr2h = x.index(round(xr+2.*h,2))
        df[0] = (-1.*f[pos_xr_h] +  0.*f[pos_xr] +  1.*f[pos_xrh])/(2.*h)
        df[1] = ( 1.*f[pos_xr_h] - 2.*f[pos_xr] +  1.*f[pos_xrh])/(h**2)
        if n >= 3:
            pos_xr_2h = x.index(x-2.*h)
            df[2] = (-1.*f[pos_xr_2h] + 2.*f[pos_xr_h] - 0.*f[pos_xr] - 2.*f[pos_xrh] + 1.*f[pos_xr2h])/(2*h**3)
        elif n == 4:
            df[3] = (1.*f[pos_xr_2h] - 4.*f[pos_xr_h] + 6.*f[pos_xr] - 4.*f[pos_xrh] + 1.*f[pos_xr2h])/(h**4)
        return df

    def Richard_extrapolation(f,x,h1,h2,n,x1,p):
        g1 = np.zeros(n)
        g2 = np.zeros(n)
        pos = x.index(x1)
        g1 = differenciate.fowardDiff_secOrder(f,x,x1,n,h1)
        g2 = differenciate.fowardDiff_secOrder(f,x,x1,n,h2)
        G = ((2.**p)*g1 - g2)/(2.**p - 1.)
        return G
