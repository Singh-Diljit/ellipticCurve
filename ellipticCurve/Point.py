"""Implement 'Point' class representing points in Q^2 on an Elliptic Curve."""

from Coordinate import Coordinate, inf
from EllipticCurve import EllipticCurve
from helperFuncs import *

class Point(Coordinate):
    """Implements points of the form: (x, y) for x, y in Q and on curve E."""

    def __init__(self, curve, x=0, y=0, inf=False, tup=None, coord=False):
        """Initialize a point.

        Parameters
        ----------
        x, y : *     : Data for rational numbers.
        inf  : bool  : If representing the point at infinity.
        tup  : tuple : For backend init, x and y are confirmed type 'Fraction'.

        Initializes
        -----------
        self.curve : EllipticCurve    : Elliptic Curve point lives on.
        self.inf   : bool             : If point represents 'point at infinity'. 
        self.x     : Fraction, float  : x-coordinate of point.
        self.y     : Fraction, float  : y-coordinate of point.
        self.tup   : tuple            : Tuple of self.x and self.y
        self.coord : Coordinate       : Point as seen without ambient curve.

        Raises
        ------
        Exception : If the point does not satisfy the EC.

        Notes
        -----
        * See documentation for 'makeFrac' for accepted inputs and examples.
        
        """
        self.curve = curve
        if coord:
            super().__init__(tup=coord.tup)
            self.coord = coord
        else:
            super().__init__(x, y, inf, tup)
            self.coord = Coordinate(tup=self.tup)

        if not(self.curve.onCurve(self.coord) or self.inf):
            raise Exception('The point is not on the curve.')

    @property
    def ID(self):
        """If point represents the additive ID."""
        return self.inf

    @property
    def order(self):
        """Return the order of a point."""
        return self.curve.order(self.coord)

    def isOrder(self, k):
        """Return the if the order of a point is equal to 'k'."""
        return self.curve.isOrder(self.coord, k)

    def isFinite(self):
        """Return if order is finite."""
        return self.curve.isTorsion(self.pos)

    def __repr__(self):
        """Return repr(self)."""
        
        curve_ = f'curve={repr(self.curve)}'
        x_ = f'x={repr(self.x)}'
        y_ = f'y={repr(self.y)}'
        
        return f'Point({x_}, {y_}, {curve_})'

    @property
    def inverse(self):
        """Return -E (under the group law)."""
        return Point(curve=self.curve, coord=self.reflect())

    def __neg__(self):
        return self.inverse    

    def __add__(self, other):
        """Return the sum of self and other."""
        res = self.curve.add(self.coord, other)
        return Point(curve=self.curve, coord=res)

    def __radd__(self, other):
        """Return the sum of self and other."""
        return self.__add__(other)

    def __iadd__(self, other):
        """Return the sum of self and other."""
        return self.__add__(other)

    def __sub__(self, other):
        """Return the self - other."""
        return self.__add__(-other)

    def __rsub__(self, other):
        """Return the other - self."""
        return other.__add__(self.__neg__)

    def __isub__(self, other):
        """Return the self - other."""
        return self.__neg__(other)

    def __mul__(self, k):
        """Return the self*other."""
        res = self.curve.mult(self.coord, k)
        return Point(curve=self.curve, coord=res)
    
    def __rmul__(self, k):
        """Return the other*self."""
        return self.__mul__(k)
    
    def __imul__(self, k):
        """Return the self*other."""
        return self.__mul__(k)
