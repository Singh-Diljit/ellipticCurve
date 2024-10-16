"""Implement 'Coordinate' class representing points in Q^2."""

from helperFuncs import makeFrac, primID

class Coordinate:
    """Implements points of the form: (x, y) for x, y in Q."""
                  
    def __init__(self, x=0, y=0, inf=False, tup=None):
        """Initialize a point.

        Parameters
        ----------
        x, y : *     : Data for rational numbers.
        inf  : bool  : If representing the point at infinity.
        tup  : tuple : For backend init, x and y are confirmed type 'Fraction'.

        Initializes
        -----------
        self.inf : bool             : If point represents 'point at infinity'. 
        self.x   : Fraction, float  : x-coordinate of point.
        self.y   : Fraction, float  : y-coordinate of point.
        self.tup : tuple            : Tuple of self.x and self.y

        Notes
        -----
        * See documentation for 'makeFrac' for accepted inputs and examples.
        
        """
        self.inf = inf
        if tup != None:
            self.x, self.y = tup
                       
        elif self.inf:
            self.x, self.y = float('inf'), float('inf')
            
        else:
            self.x = makeFrac(x)
            self.y = makeFrac(y)
            
        self.tup = (self.x, self.y)
            
    def reflect(self, axis=0):
        """Return the reflection of self in the x- or y- axis.

        Parameters
        ----------
        axis : int : Either 0 (for reflection in the x-axis) or 1 (y-axis).

        Returns
        -------
        P : Coordinate : Either point at infinity or the reflected point.

        """
        if self.inf:
            res = Coordinate(inf=True)
        else:
            newX = self.x if (axis==0) else -self.x 
            newY = self.y if (axis==1) else -self.y
            
            res = Coordinate(tup=(newX, newY))

        return res

    @property
    def integral(self):
        """Return if self is in Z^2.

        Notes
        -----
        The point at infinity is considered in Z^2.

        """
        if self.inf:
            res = True
        else:
            res = (self.x.denominator == 1) and (self.y.denominator == 1)
            
        return res

    def __repr__(self):
        """Return repr(self)."""
        return f'Coordinate(x={repr(self.x)}, y={repr(self.y)}, inf={self.inf})'
    
    def __str__(self):
        """Return str(self)."""
        return f'({str(self.x)}, {str(self.y)})'

    def __getitem__(self, index):
        """Allow access to x (index=0) and y (index=1) components."""
        return self.tup[index]

    def __eq__(self, other):
        """Return if self and other represent the same point in space."""
        
        return self.tup == other[:2]
    
    def __ne__(self, other):
        """Return if self and other do not represent the same point in space."""
        return not(self.__eq__(other))

    def slope(self, other):
        """Return the slope of the line connecting self and other."""
        if self.inf or other.inf or (self.x == other.x):
            res = float('inf')
        else:
            res = (self.y - other.y) / (self.x - other.x)
            
        return res

    def isCurveInd(self, Q=primID, k=1):
        """Return if [k]self + Q is curve-independent under the EC group law.

        Parameters
        ----------
        Q : Coordinate : Point in Q^2.
        k : int        : Scaling factor.

        Returns
        -------
        - : bool : If [k]self + Q can be computed without curve data.

        """
        if Q == primID:
            Q = Coordinate(inf=True)
            
        trivAddition = (self.inf or Q.inf or self.reflect() == Q)
        trivMult = (k in {-1, 0, 1})

        return trivAddition and trivMult
    
    def resCurveInd(self, Q=primID, k=1):
        """Return [k]self + Q given sum is EC-independent.

        Parameters
        ----------
        Q : Coordinate : Point in Q^2.
        k : int        : Scaling factor.

        Returns
        -------
        res : Coordinate : Is equal to [k]P + Q.

        Notes
        -----
        The logic in this function is based off of backend use cases, namely
        it is never the case: (k != 1) and (Q.inf == False).

        """
        if Q == primID:
            Q = Coordinate(inf=True)
        
        if self.inf and k == 1:
            res = Q
        elif (self.reflect() == Q) or (k == 0):
            res = Coordinate(inf=True)
        elif k == -1:
            res = self.reflect()
        if self.inf:
            res = Q

        return res
        
inf = Coordinate(inf=True)
