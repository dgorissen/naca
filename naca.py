"""
Python code to generate 4 and 5 digit NACA profiles

Pots of the Matlab code available here:
    http://www.mathworks.com/matlabcentral/fileexchange/19915-naca-4-digit-airfoil-generator
    http://www.mathworks.com/matlabcentral/fileexchange/23241-naca-5-digit-airfoil-generator
    
Copyright (C) 2011 by Dirk Gorissen <dgorissen@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import math

"""
Emulate Matlabs linspace
"""
def linspace(start,stop,np):

    delta = (stop - start) / (np - 1.0)
    
    p = [start]
    for i in range(np-1):
        p.append(p[-1] + delta)

    return p

"""
A cubic spline interpolation on a given set of points (x,y)
Recalculates everything on every call which is far from efficient but does the job for now
should eventually be replaced by an external helper class
"""
def Interpolate(xa,ya,queryPoints):
    
    # PreCompute() from Paint Mono which in turn adapted:
    # NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
    # ISBN 0-521-43108-5, page 113, section 3.3.
    # http://paint-mono.googlecode.com/svn/trunk/src/PdnLib/SplineInterpolator.cs
    
    #number of points
    n = len(xa)
    u = [0]*n
    y2 = [0]*n
    
    u[0] = 0;
    y2[0] = 0;
    
    for i in range(1,n-1):
    
        # This is the decomposition loop of the tridiagonal algorithm.
        # y2 and u are used for temporary storage of the decomposed factors.
    
        wx = xa[i + 1] - xa[i - 1]
        sig = (xa[i] - xa[i - 1]) / wx
        p = sig * y2[i - 1] + 2.0
        
        y2[i] = (sig - 1.0) / p
    
        ddydx = (ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]) - (ya[i] - ya[i - 1]) / (xa[i] - xa[i - 1])
    
        u[i] = (6.0 * ddydx / wx - sig * u[i - 1]) / p;
     
    
    y2[n - 1] = 0;
    
    # This is the backsubstitution loop of the tridiagonal algorithm
    #((int i = n - 2; i >= 0; --i):
    for i in range(n-2,-1,-1):
        y2[i] = y2[i] * y2[i + 1] + u[i];

    # Interpolate() adapted from Paint Mono which in turn adapted:
    # NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
    # ISBN 0-521-43108-5, page 113, section 3.3.
    # http://paint-mono.googlecode.com/svn/trunk/src/PdnLib/SplineInterpolator.cs

    results = [0]*n;

    #loop over all query points
    for i in range(len(queryPoints)):
        # bisection. This is optimal if sequential calls to this
        # routine are at random values of x. If sequential calls
        # are in order, and closely spaced, one would do better
        # to store previous values of klo and khi and test if

        klo = 0;     
        khi = n - 1; 

        while (khi - klo > 1):
            k = (khi + klo) >> 1;
            if (xa[k] > queryPoints[i]):
                khi = k; 
            else:
                klo = k;

        h = xa[khi] - xa[klo];
        a = (xa[khi] - queryPoints[i]) / h;
        b = (queryPoints[i] - xa[klo]) / h;

        # Cubic spline polynomial is now evaluated.
        results[i] = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2[klo] + (b * b * b - b) * y2[khi]) * (h * h) / 6.0;

    return results;

"""
Returns n points (for EACH HALF) in [0 1] for the given 4 digit NACA number string
"""
def naca4(number, n):
    
    finite_TE = False
    half_cosine_spacing = False
    
    m = float(number[0])/100.0;
    p = float(number[1])/10.0;
    t = float(number[2:])/100.0;

    a0 = 0.2969;
    a1 = -0.1260;
    a2 = -0.3516;
    a3 = 0.2843;

    if finite_TE:
        a4 = -0.1015 # For finite thick TE
    else:
        a4 = -0.1036  # For zero thick TE

    if half_cosine_spacing:
        beta = linspace(0.0,math.pi,n+1)
        x = map(lambda xx : (0.5*(1.0-math.cos(xx))),beta)  # Half cosine based spacing
    else:
        x = linspace(0.0,1.0,n+1)
    
    yt = map(lambda xx: (t/0.2)*(a0*math.sqrt(xx)+a1*xx+a2*math.pow(xx,2)+a3*math.pow(xx,3)+a4*math.pow(xx,4)),x);

    xc1 = filter(lambda xx: xx <= p, x);
    xc2 = filter(lambda xx: xx > p, x);
    xc = xc1 + xc2;

    if p == 0:
        xu = x;
        yu = yt;

        xl = x;
        yl = [-xx for xx in yt];
    
        zc = [0]*len(xc)
    else:
        yc1 = map(lambda xx: (m/math.pow(p,2))*(2*p*xx-math.pow(xx,2)),xc1);
        yc2 = map(lambda xx: (m/math.pow((1-p),2))*((1-2*p)+2*p*xx-math.pow(xx,2)),xc2);
        zc = yc1 + yc2

        dyc1_dx = map(lambda xx: (m/math.pow(p,2))*(2*p-2*xx),xc1);
        dyc2_dx = map(lambda xx: (m/math.pow((1-p),2))*(2*p-2*xx),xc2);
        dyc_dx = dyc1_dx + dyc2_dx;
        
        theta = [math.atan(xx) for xx in dyc_dx];
    
        xu = map(lambda xx,yy,zz: xx-yy * math.sin(zz),x,yt,theta);
        yu = map(lambda xx,yy,zz: xx+yy * math.cos(zz),zc,yt,theta);
    
        xl = map(lambda xx,yy,zz: xx + yy * math.sin(zz),x,yt,theta);
        yl = map(lambda xx,yy,zz: xx - yy * math.cos(zz),zc,yt,theta);
        
    X = xu[::-1] + xl[1:]
    Z = yu[::-1] + yl[1:]
    
    pts = zip(X,Z)

    return pts
    
"""
Returns n points (for EACH HALF) in [0 1] for the given 5 digit NACA number string
"""
def naca5(number,n):

    finite_TE = False
    half_cosine_spacing = False
    
    naca1 = int(number[0]);
    naca23 = int(number[1:3]);
    naca45 = int(number[3:]);

    cld = naca1*(3.0/2.0)/10.0;
    p = 0.5*naca23/100.0;
    t = naca45/100.0;

    a0 =  0.2969;
    a1 = -0.1260;
    a2 = -0.3516;
    a3 = 0.2843;

    if finite_TE:
        a4 = -0.1015; # For finite thickness trailing edge
    else:
        a4 = -0.1036;  # For zero thickness trailing edge
    
    if half_cosine_spacing:
        beta = linspace(0.0,math.pi,n+1)
        x = map(lambda x : (0.5*(1.0-math.cos(x))),beta)  # Half cosine based spacing
    else:
        x = linspace(0.0,1.0,n+1)

    yt = map(lambda xx : (t/0.2)*(a0*math.sqrt(xx)+a1*xx+a2*math.pow(xx,2)+a3*math.pow(xx,3)+a4*math.pow(xx,4)),x);

    P = [0.05,0.1,0.15,0.2,0.25]
    M = [0.0580,0.1260,0.2025,0.2900,0.3910]
    K = [361.4,51.64,15.957,6.643,3.230]

    m = Interpolate(P,M,[p])[0]
    k1 = Interpolate(M,K,[m])[0]
    
    xc1 = filter(lambda xx: xx <= m, x);
    xc2 = filter(lambda xx: xx > m, x);
    xc = xc1 + xc2;
    
    if p == 0:
        xu = x
        yu = yt
        
        xl = x
        yl = [-x for x in yt];
    
        zc = [0]*len(xc)
    else:
        yc1 = map(lambda xx : (1.0/6.0)*k1*( math.pow(xx,3)-3*m*math.pow(xx,2)+ math.pow(m,2)*(3-m)*xx),xc1)
        yc2 = map(lambda xx : (1.0/6.0)*k1*math.pow(m,3)*(1-xx),xc2)
        zc = map(lambda xx : (cld/0.3) * xx, yc1 + yc2)

        dyc1_dx = map(lambda xx : cld/0.3*(1.0/6.0)*k1*( 3*math.pow(xx,2)-6*m*xx+math.pow(m,2)*(3-m) ), xc1);
        dyc2_dx = [cld/0.3*(1.0/6.0)*k1*math.pow(m,3)]*len(xc2);
        
        dyc_dx = dyc1_dx + dyc2_dx;
        theta = [math.atan(xx) for xx in dyc_dx];

        xu = map(lambda xx,yy,zz: xx - yy * math.sin(zz),x,yt,theta);
        yu = map(lambda xx,yy,zz: xx + yy * math.cos(zz),zc,yt,theta);

        xl = map(lambda xx,yy,zz: xx + yy * math.sin(zz),x,yt,theta)
        yl = map(lambda xx,yy,zz: xx - yy * math.cos(zz),zc,yt,theta)

    X = xu[::-1] + xl[1:];
    Z = yu[::-1] + yl[1:];
  
    pts = zip(X,Z)

    return pts

if __name__ == "__main__":
    
    # Examples
    pts = naca4("0012",60)
    #pts = naca5("23015",60)
    
    for p in pts:
        print p[0],p[1]
