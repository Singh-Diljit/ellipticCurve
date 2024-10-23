"""Implement EllipticCurve class representing elliptic curves over Q."""

from Coordinate import Coordinate, inf
from helperFuncs import makeFrac, makeStr_EC, sqRt_Q, isSmooth

class EllipticCurve:
    """Implements elliptic curves over Q."""

    def __init__(self, a, b, tup=None):
        """Initialize an elliptic curve.

        Parameters
        ----------
        a, b : *     : Data for coefficients in E : (y**2 = x**3 + ax + b).
        tup  : tuple : For backend init, a and b are confirmed type 'Fraction'.

        Initializes
        -----------
        self.a   : Fraction : Coefficient of linear term in associated cubic.
        self.b   : Fraction : Constant term in associated cubic.
        self.tup : tuple    : Tuple of self.x and self.y

        Raises
        ------
        Exception : If the point does not satisify the EC.

        Notes
        -----
        * See documentation for 'makeFrac' for accepted inputs and examples.

        """
        if tup != None:
            self.a, self.b = tup
        else:
            self.a, self.b = makeFrac(a), makeFrac(b)

        self.tup = (self.a, self.b)

        if not(isSmooth(self.a, self.b)):
            raise Exception('The curve is not smooth.')

    # Discriminant related properties
    
    @property
    def discriminant(self):
        """Return the discriminant.

        Notes
        -----
        Not only is the discriminant dependant on the equation, but different
        authors use different expressions. Common expressions include:
            (1) -16(4a**3 + 27b**2)
            (2) -(4a**3 + 27b**2)
            (3) 4a**3 + 27b**2 [Cassels]
        This class uses option (1).

        See Also
        --------
        See the 'jInvariant' property for an invariant related to the
        discriminant.

        """
        return -64*self.a**3 - 432*self.b**2

    @property
    def cubicDisc(self):
        """Return the discriminant of the associated cubic.

        Notes
        -----
        The result is -(4a**3 + 27b**2).

        See Also
        --------
        See the 'jInvariant' and 'discriminant' properties.
        
        """
        return -4*self.a**3 - 27*self.b**2

    @property
    def jInvariant(self):
        """Return the j-invariant."""
        fourACubed = 4 * self.a**3
        return 1728 * fourACubed / (fourACubed + 27*self.b**2)
    
    @property
    def numberComponents(self):
        """Return the number of components the graph has.

        Notes
        -----
        The number of components is determined by the sign of the discriminant.
        If positive the curve has two components, otherwise the graph is
        connected and has one component.
        
        """
        return 2 if (27*self.b**2 < -4*self.a**3) else 1

    # Basic Classification, Representation, and Accessing

    def isomorphic(self, other):
        """Return if self and other are isomorphic."""
        return self.jInvariant == other.jInvariant

    @property
    def integral(self):
        """Return if (self.a, self.b) is in Z^2."""
        return (self.a.denominator == 1) and (self.b.denominator == 1)
            
    def __eq__(self, other):
        """Return if self and other represent the same curve."""
        return self.tup == other.tup

    def __ne__(self, other):
        """Return if self and other do not represent the same curve."""
        return not(self.__eq__(other))
    
    def __repr__(self):
        """Return repr(self)."""
        return f'EllipticCurve(a={repr(self.a)}, b={repr(self.b)})'
    
    def __str__(self):
        """Return str(self)."""
        return makeStr_EC(self.a, self.b)

    # Features As A Cubic

    @property
    def cubic(self):
        """Return the associated cubic as a callable function."""
        return lambda x: (x**3 + self.a*x + self.b)

    @property
    def derivative(self):
        """Return dy/dx as a callable function.

        Returns
        -------
        res : func : Function accepting a point in Q^2 as an input. 

        Notes
        -----
        This is an application of implicit differentiation.
        
        """
        return lambda P: (3*P.x**2 + self.a) / (2*P.y)

   # Point Verification
   
    def onCurve(self, P):
        """Return if a point in Q^2 is on the curve."""
        return P.inf or (P.y**2 == self.cubic(P.x))

    def findY(self, x):
        """Return the value(s) of y so that P = (x, y) are on the curve.

        Returns
        -------
        y : set : Solution(s) to sqrt(x^3 + ax + b).

        Raises
        ------
        Exception: 'Input is not in the domain.'

        See Also
        --------
        'inDomain' method can be useful.

        """
        ySq = self.cubic(x)
        y = sqRt_Q(ySq)
        if y != -1: #y == -1 iff ySq not in Q^2.
            res = {y, -y}
        else:
            raise Exception('Input is not in the domain.')

        return res

    def inDomain(self, x):
        """Return if a value is in the domain of the curve."""
        if x == float('inf'):
            res = True
        else:
            res = sqRt_Q(self.cubic(x))

        return res != -1

    # Group Arithmetic

    def tangent(self, P):
        """Return slope of the line tangent to self and a given point.

        Parameters
        ----------
        P : Coordinate : A point on the curve.
        
        """
        if P.y == float('inf'):
            res = float('inf')
        else:
            res = (3*P.x**2 + self.a) / (2*P.y)
            
        return res

    def __add_sumExists(self, P, Q, slope):
        """Given the slope of the secant/tangent line, return P + Q.

        Parameters
        ----------
        P, Q  : Coordinate : Points in Q^2.
        slope : Fraction   : Slope the secant (or tangent) line.

        Returns
        -------
        - : Point : The sum of P and Q.
        
        """
        resX = slope**2 - P.x - Q.x
        resY = slope*(P.x - resX) - P.y
        
        return Coordinate(resX, resY)

    def double(self, P):
        """Return [2]P.

        Parameters
        ----------
        P : Coordinate : A point on the curve.

        Returns
        -------
        res : Coordinate : [2]P.
        
        """      
        if P == P.reflect():
            res = inf

        else:
            tangentSlope = self.tangent(P)
            res = self.__add_sumExists(P, P, tangentSlope)

        return res

    def add(self, P, Q):
        """Return the sum of two points.

        Parameters
        ----------
        P, Q  : Coordinate : Point on the curve.

        Returns
        -------
        res : Coordinate: P+Q.
        
        """    
        if P == Q:
            res = self.double(P)

        elif (P == inf) or (Q == inf):
            res = P if (Q == inf) else Q

        elif P == Q.reflect():
            res = inf
            
        else:
            slope_ = P.slope(Q)
            res = self.__add_sumExists(P, Q, slope_)

        return res

    def sub(self, P, Q):
        """Return the difference of two points.

        Parameters
        ----------
        P, Q  : Coordinate : Point on the curve.

        Returns
        -------
        res : Coordinate: P-Q.
        
        """
        return self.add(P, Q.reflect())
                    
    def mult(self, P, k):
        """Return k[P].

        Parameters
        ----------
        P : Coordinate : Point on the curve.
        k : int        : Scaling factor.

        Returns
        -------
        res : Coordinate : [k]P.
        
        """
        if k < 0:
            P, k = P.reflect(), -k

        if k in {0, 1}:
            res = inf if (k == 0) else P

        elif P.inf: res = P

        else:
            res = inf
            for bit in bin(k)[:1:-1]:
                if bit == '1':
                    res = self.add(P, res)
                
                P = self.double(P)
                
        return res

    # Torsion, Nagell-Lutz, Mazur's Theorem
    
    def nagellLutz(self, P):
        """Return if a point is a candidate for having finite order."""
        if P == inf:
            return True
        torsionCand = False
        if self.integral and P.integral:
            if P.y == 0:
                torsionCand = True
            else:
                D = self.cubicDisc
                torsionCand = (D % P.y**2 == 0)
                
        elif not self.integral:
            if P.integral:
                torsionCand = True
            else:
                torsionCand = (P.x.denominator == 4 and P.y.denominator == 8)
            
        return torsionCand
        
    def order(self, P):
        """Return the order of a point."""
        orderP = float('inf')
        if P == inf:
            orderP = 1
        elif self.nagellLutz(P):
            Q = Coordinate(inf=True)
            for k in range(1, 13):
                Q = self.add(Q, P)
                if Q.inf:
                    orderP = k
                    break

        return orderP

    def isTorsion(self, P):
        """Return is a point has finite order."""
        return self.order(P) < 13

    def isOrder(self, P, k):
        """Return the if a value is the order of a point."""
        res = False
        if P.inf:
            res = (k==1)
            
        elif self.mult(P, k) != inf:
            res = False
            
        elif k in {2, 3, 5, 7, 11}:
            res = True

        elif k == 9:
            res = not(self.isOrder(P, 3))

        elif k in {4, 6, 10}:
            res = not(self.isOrder(P, 2))

        elif k == 8:
            res = not(self.isOrder(P, 2) or self.mult(P, 4)==inf)

        elif k == 12:
            if self.mult(P, 6) != inf:
                res = not(self.mult(P, 4) == inf)
            else:
                res = False

        return res

            
            
