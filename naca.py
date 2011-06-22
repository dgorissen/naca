"""
Python code to generate 4 and 5 digit NACA profiles

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
Returns n points in [0 1] for the given 4 digit NACA number string
"""
def naca4(number, n):

    naca1 = int(number[0]);
    naca2 = int(number[1]);
    naca34 = int(number[2:]);
    
    ml = naca1/100.0;
    pl = naca2/10.0;
    dl = naca34/100.0;

    m = n/2+1;
    
    # Airfoil-Drop
    a0 = 1.4845;
    a1 = -0.63;
    a2 = -1.7580;
    a3 = 1.4215;
    a4 = -0.5075;

    # Create the point objects
    ptList = [[0,0] for i in range(n)];
    
    #Calculate the points ...

    #TrailingEdge
    ptList[0] = [1.0,0.0]
               
    # Between
    for i in range(1,m-1):
        x=1.0-0.5*(1.0-math.cos(math.pi*(1.0*i)/(1.0*m-1.0)))
        yt=dl*(a0*math.sqrt(x)+x*(a1+x*(a2+x*(a3+x*a4))))
        
        if (x<=pl):
            y = (ml/math.pow(pl,2))*(2.0*pl-x)*x;
            dy = (ml/math.pow(pl,2))*(2.0*pl-2.0*x)
        else:
            y = (ml/math.pow((1.0-pl),2))*((1.0-2.0*pl)+2.0*pl*x-math.pow(x,2))
            dy = (ml/math.pow((1.0-pl),2))*(2.0*pl-2.0*x)

        theta = math.atan(dy);
        xu = x-yt*math.sin(theta);
        yu = y+yt*math.cos(theta);
        xl = x+yt*math.sin(theta);
        yl = y-yt*math.cos(theta);

        # upper point
        ptList[i] = [xu,yu];
        
        #lower point
        ptList[n-i] = [xl,yl];

    # LeadingEdge point
    ptList[m-1] = [0.0,0.0];

    return ptList;

"""
Returns n points in [0 1] for the given 5 digit NACA number string
"""
def naca5(number,n):

    naca1 = int(number[0]);
    naca23 = int(number[1:3]);
    naca45 = int(number[3:]);

    cld = naca1*0.75/10.0;
    pl = 0.5*naca23/100.0;
    dl = naca45/100.0;

    m = n/2+1;

    # Airfoil-Drop
    a0 = 1.4845;
    a1 = -0.63;
    a2 = -1.7580;
    a3 = 1.4215;
    a4 = -0.5075;

    # Create the point objects
    ptList = [[0.0,0.0] for i in range(n)];

    # Calculate the points ...

    # TrailingEdge
    ptList[0] = [1.0,0.0]

    P = [0.05,0.1,0.15,0.2,0.25]
    M = [0.0580,0.1260,0.2025,0.2900,0.3910]
    K = [361.4,51.64,15.957,6.643,3.230]
    mll = Interpolate(P,M,[pl])

    ml = mll[0];

    k1l = Interpolate(M,K,[ml]);
    k1 = k1l[0];
    
    for i in range(1,m-1):
        x=1.0-0.5*(1.0-math.cos(math.pi*(1.0*i)/(1.0*m-1.0)));
        yt=dl*(a0*math.sqrt(x)+x*(a1+x*(a2+x*(a3+x*a4))));
        
        if (x<=pl):
            y = cld/0.3 * (1.0/6.0)*k1*( math.pow(x,3)-3.0*ml*math.pow(x,2)+ml*ml*(3.0-ml)*x );
            dy = (1.0/6.0)*k1*( 3.0*math.pow(x,2)-6.0*ml*x+ml*ml*(3.0-ml) );
        else:
            y = cld/0.3 * (1.0/6.0)*k1*math.pow(ml,3)*(1.0-x);
            dy = -(1.0/6.0)*k1*math.pow(ml,3);

        theta = math.atan(dy);
        xu = x-yt*math.sin(theta);
        yu = y+yt*math.cos(theta);
        xl = x+yt*math.sin(theta);
        yl = y-yt*math.cos(theta);
        
        # upper point
        ptList[i] = [xu,yu];

        # lower point
        ptList[n-i] = [xl,yl];

    # LeadingEdge point
    ptList[m-1] = [0.0,0.0]

    return ptList;

if __name__ == "__main__":
    
    # Examples
    pts = naca4("2441",100)
    #pts = naca5("23015",100)
    
    for p in pts:
        print p[0],p[1]