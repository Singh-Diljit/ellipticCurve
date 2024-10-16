"""Functions in support of the 'Coordinate' and 'EC' classes."""

from fractions import Fraction

primID = (float('inf'), float('inf'))

def sqRt_Z(N):
    """Return if input is in Z^2, if so return the principal square root.

    Parameters
    ----------
    N : int : Integer whose root is to be inspected.

    Returns
    -------
    res : int : Either the principal square root or -1.
    
    """
    res = -1
    if N in {0, 1}:
        res = N

    #Binary search to find square root if possible.
    elif (N > 1) and (N%4 in {0, 1}):
        lo, hi = 0, N
        while lo <= hi:
            mid = (lo + hi) // 2
            guess = mid ** 2     
            if guess > N: 
                hi = mid - 1
            elif guess < N: 
                lo = mid + 1
            else:
                res = mid
                break
            
    return res

def sqRt_Q(X):
    """Return if input is in Q^2, if so return the principal square root.

    Parameters
    ----------
    N : Fraction : Rational whose root is to be inspected.

    Returns
    -------
    res : Fraction : Either the principal square root or -1.
    
    """
    res = Fraction(-1)
    if X in {0, 1}:
        res = X

    elif X > 0:
        num = sqRt_Z(X.numerator)
        den = -1
        if num != -1: #flag for num not in Z^2
            den = sqRt_Z(X.denominator)
        if -1 not in {num, den}:
            res = Fraction(num, den)

    return res

def gcd(n, k):
    """Return the GCD of two integers via the Euclidean algorithm."""
    while k:
        n, k = k, n%k

    return abs(n)

def lcm(n, k):
    """Return the LCM of two integers."""
    return (n*k) // gcd(n, k)

def isSmooth(a, b):
    """Return if (y**2 = x**3 + ax + b) for a, b in Q is smooth.

    Parameters
    ----------
    a, b : int : Coefficients used to define: (y**2 = x**3 + ax + b).

    Returns
    -------
    - : bool : If the defined curve is smooth.
    
    """
    negDisc = 64 * a**3 + 432 * b**2
    return (negDisc != 0)

def makeStr_EC(a, b):
    """Represent E : (y**2 = x**3 + ax + b) as a string.

    Parameters
    ----------
    a, b : Fraction : Coefficient of EC's associated cubic.

    Returns
    -------
    res : str : Represents E : (y**2 = x**3 + ax + b) as a string.
    
    Example(s)
    ----------
    >>> makeStr_EC(Fraction(-1), Fraction(0))
    >>> E : (y**2 = x**3 - x)

    >>> makeStr_EC(Fraction(0), Fraction(1, 5))
    >>> E : (y**2 = x**3 + 1/5)

    >>> makeStr_EC(Fraction(-2, 3), Fraction(1, 5))
    >>> E : (y**2 = x**3 - 2/3x + 1/5)
    
    """
    prefix = 'y**2 = x**3'
    linear, constant = '', ''
    if a != 0:
        sgn = ' + ' if (a > 0) else ' - '
        a = abs(a)
        linear = (sgn+'x') if (a == 1) else (sgn+str(a)+'x')
    if b != 0:
        sgn = ' + ' if (b > 0) else ' - '
        b = abs(b)
        constant = sgn + str(b)

    return f'E : ({prefix + linear + constant})'

def makeFrac(X):
    """Convert input data to an instance of 'Fraction'.

    Parameters
    ----------
    X : * : Data for Fraction.__new__.

    Returns
    -------
    res : Fraction : Represents an element in Q.

    Notes
    -----
    * Accepted objects are: int, float, or have a __getitem__ attribute for
    index equals 0 or 1. In the last case, X[0] and X[1] should subclass to
    instances of 'fractions.Rational'.

    See Also
    --------
    Documentation for 'fractions.Fraction.__new__'.
    
    Example(s)
    ----------
    >>> res = makeFrac(-2.5)
    >>> repr(res)
    >>> Fraction(-5, 2)

    >>> res = makeFrac((3, 5))
    >>> repr(res)
    >>> Fraction(3, 5)

    >>> res = makeFrac(4, 9)
    >>> repr(res)
    >>> Fraction(4, 9)

    >>> res = makeFrac(3.0)
    >>> repr(res)
    >>> Fraction(3, 1)

    >>> makeFrac((3.0, 1))
    >>> TypeError: both arguments should be Rational instances
    
    """
    try:
        res = Fraction(X)
    except:
        res = Fraction(X[0], X[1])
    return res
