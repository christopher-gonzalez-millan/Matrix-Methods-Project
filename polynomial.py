import fieldmath


# Polynomial operations without using class
# This although more painful to implement into code should speed things up a huge amount
# Using the Polynomial class seems to be slightly slower than PolynomialOperations
class PolynomialOperations(object):
    def __init__(self, field):
        self.f = field

    def poly_call(self, poly, x):
        res = self.f.zero()
        for c in range(len(poly)-1, -1, -1):
            res = self.f.add(self.f.multiply(res, x), poly[c])
        return res

    def poly_degree(self, poly):
        """
        Return the degree of the polynomial. ex 3 + 3x + x^4 has degree 4
        """
        self.poly_trim(poly)
        if poly:
            return len(poly) - 1
        else:
            return 0

    def poly_trim(self, poly):
        """
        removes all leading 0 terms from the coefficients array
        """
        while poly and poly[-1] == self.f.zero():
            poly.pop()
        if len(poly) == 0:
            poly = [self.f.zero()]

    def poly_lead(self, poly):
        """
        Returns the leading term. ex x^4 + 3x + 3 will return 1. returns 0 for the polynomial 0.
        """
        self.poly_trim(poly)
        for i in range(len(poly) - 1, -1, -1):
            if poly[i] != 0:
                return poly[i]
        return self.f.zero()
    
    def set(self, poly, degree_term, coeff_val):
        """
        Sets coefficient of the degree term to the coeff value. ex set(2, 2) for x+2 will result in the polynomial being
        2x^2 + x + 2
        """
        if len(poly) <= degree_term:
            poly += [self.f.zero()] * (degree_term - len(poly) + 1)
        poly[degree_term] = coeff_val

    def poly_divmod(self, poly1, poly2):
        """
        Heavily modified version of a polynomial long division originally found here:
        https://rosettacode.org/wiki/Polynomial_long_division
        This might need to be changed to no affect poly1.
        """
        # top = self.__class__(*self.coefficients, f=self.f, flipped=True)
        # bottom = self.__class__(*other.coefficients, f=other.f, flipped=True)
        deg_top = self.poly_degree(poly1)
        deg_bot = self.poly_degree(poly2)

        if deg_bot < 0:
            raise ZeroDivisionError

        if deg_top >= deg_bot:
            q = [self.f.zero()] * deg_top
            while deg_top >= deg_bot and self.poly_lead(poly1) != self.f.zero():  # and top.lead() != self.f.zero():
                d = [self.f.zero()] * (deg_top - deg_bot) + poly2
                inv = self.f.reciprocal(self.poly_lead(d))
                multiplier = self.f.multiply(self.poly_lead(poly1), inv)
                self.set(q, deg_top - deg_bot, multiplier)
                d = [self.f.multiply(multiplier, c) for c in d]
                poly1 = self.poly_subtract(poly1, d)
                deg_top = self.poly_degree(poly1)
        else:
            q = [self.f.zero()]
        return q, poly1

    def poly_add(self, poly1, poly2):
        return [self.f.add(t1, t2) for t1, t2 in zip_longest(poly1, poly2, fillvalue=self.f.zero())]

    def poly_subtract(self, poly1, poly2):
        return [self.f.subtract(t1, t2) for t1, t2 in zip_longest(poly1, poly2, fillvalue=self.f.zero())]

    def poly_equal(self, poly1, poly2):
        self.poly_trim(poly1)
        self.poly_trim(poly2)
        if poly1 == poly2:
            return True
        else:
            return False

    def poly_not_equal(self, poly1, poly2):
        return not self.poly_equal(poly1, poly2)

    def poly_gcd(self, poly1, poly2):
        if self.poly_degree(poly1) > self.poly_degree(poly2):
            return self.poly_gcd(poly2, poly1)

        if self.poly_equal(poly1, [self.f.zero()]):
            return poly2

        q, r = self.poly_divmod(poly2, poly1)
        self.poly_trim(r)

        return self.poly_gcd(r, poly1)
# Poly class


def poly_gcd(a, b):
    if not isinstance(a, Polynomial) or not isinstance(b, Polynomial):
        raise Exception("poly_gcd requires two Polynomial objects.")
    if a.f != b.f:
        raise Exception("Polynomial objects do not have the same fields.")

    if a.degree() > b.degree():
        return poly_gcd(b, a)

    if a == Polynomial(0, f=a.f):
        return b

    q, r = divmod(b, a)
    r.trim()

    return poly_gcd(r, a)


class Polynomial(object):
    """
    This is very loosely built off the implementation here, much has been changed and added:
    https://www.python-course.eu/polynomial_class_in_python.php

    This entire class now works entirely in the space F[x] where F is the provided field.
    To use this function for R[x] a new field should be created representing the real numbers.
    I have added the following functions to the Polynomial class:
            degree, lead, set, __divmod__, __mod__, __mul__, __rmul__, __neg__, __eq__, __ne__, trim
    And I have modified:
            __init__, __call__, __add__, __sub__
    I have also added the following additional function to this class:
            poly_gcd

    """
    def __init__(self, *coefficients, f, flipped=False):
        """
        example implementation Polynomial(1, 2, 3, f=field) -> 1 + 2x + 2x^2.
        Note: all polynomial operations with 2 different polynomials must be done in the same field.
        :param coefficients: default in the order a_0, a_1, a_2, ..., a_n
        :param f: Specifies finite field. All operations will be down in that field. refer to fieldmath.py
        :param flipped: True if coeff order a_n, a_n-1, a_n-2, ..., a_0
        """

        # Notice elements are stores in the order a_n, a_n-1, a_n-2, ..., a_0. this is to make some operations easier
        if flipped:
            self.coefficients = list(coefficients)
        else:
            self.coefficients = list(coefficients)[::-1]

        if isinstance(f, fieldmath.Field):
            self.f = f
        else:
            raise Exception("Please provide a valid Field")

    def __repr__(self):
        """
        method to return the canonical string representation
        of a polynomial.
        """
        return "Polynomial" + str(self.coefficients)

    def __call__(self, x):
        res = self.f.zero()
        for coeff in self.coefficients:
            res = self.f.add(self.f.multiply(res, x), coeff)
        return res

    def degree(self):
        """
        Return the degree of the polynomial. ex x^4 + 3x + 3 has degree 4
        """
        self.trim()
        for i in range(len(self.coefficients)):
            if self.coefficients[i] != 0:
                return len(self.coefficients)-1-i
        return 0

    def lead(self):
        """
        Returns the leading term. ex x^4 + 3x + 3 will return 1. returns 0 for the polynomial 0.
        """
        self.trim()
        for i in range(len(self.coefficients)):
            if self.coefficients[i] != 0:
                return self.coefficients[i]
        return self.f.zero()

    def set(self, degree_term, coeff_val):
        """
        Sets coeffiecent of the degree term to the coeff value. ex set(2, 2) for x+2 will result in the polynomial being
        2x^2 + x + 2
        """
        if len(self.coefficients) <= degree_term:
            self.coefficients = [self.f.zero()] * (degree_term - len(self.coefficients) + 1) + self.coefficients
        self.coefficients[-degree_term-1] = coeff_val

    def __divmod__(self, other):
        """
        This function overrides the divmod(a, b) function.
        Heavily modified version of a polynomial long division originally found here:
        https://rosettacode.org/wiki/Polynomial_long_division
        """
        if self.f != other.f:
            raise Exception("Polynomial objects do not have the same fields.")

        top = self.__class__(*self.coefficients, f=self.f, flipped=True)
        bottom = self.__class__(*other.coefficients, f=other.f, flipped=True)
        deg_top = top.degree()
        deg_bot = bottom.degree()

        if deg_bot < 0:
            raise ZeroDivisionError

        if deg_top >= deg_bot:
            q = self.__class__(self.f.zero(), f=self.f)
            while deg_top >= deg_bot and top.lead() != self.f.zero():
                d = self.__class__(*(other.coefficients + [self.f.zero()] * (deg_top - deg_bot)), f=self.f,
                                   flipped=True)
                inv = self.f.reciprocal(d.lead())
                multiplier = self.f.multiply(top.lead(), inv)
                q.set(deg_top - deg_bot, multiplier)
                d = multiplier * d
                top = top - d
                deg_top = top.degree()
        else:
            q = self.__class__(0, f=self.f)
        return q, top
    
    def __mod__(self, other):
        _, r = self.__divmod__(other)
        return r

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            if self.f != other.f:
                raise Exception("Polynomial objects do not have the same fields.")
            c1 = self.coefficients
            c2 = other.coefficients
            res = [self.f.zero()]*(len(c1)+len(c2)-1)
            for i, c1_coeff in enumerate(c1):
                for j, c2_coeff in enumerate(c2):
                    res[i + j] = self.f.add(res[i + j], self.f.multiply(c1_coeff, c2_coeff))
            return self.__class__(*res, f=self.f, flipped=True)
        else:
            # self.f.multiply will automatically check that everything is being down in the field f.
            cs = self.coefficients
            return self.__class__(*[self.f.multiply(other, c) for c in cs], f=self.f, flipped=True)

    __rmul__ = __mul__

    def __add__(self, other):
        if not isinstance(other, Polynomial):
            raise Exception("Can only add Polynomial to another Polynomial")
        if self.f != other.f:
            raise Exception("Polynomial objects do not have the same fields.")
        c1 = self.coefficients[::-1]
        c2 = other.coefficients[::-1]

        return self.__class__(*[self.f.add(t1, t2) for t1, t2 in zip_longest(c1, c2, fillvalue=self.f.zero())],
                              f=self.f)

    def __sub__(self, other):
        if not isinstance(other, Polynomial):
            raise Exception("Can only subtract a Polynomial from another Polynomial")
        if self.f != other.f:
            raise Exception("Polynomial objects do not have the same fields.")
        c1 = self.coefficients[::-1]
        c2 = other.coefficients[::-1]
        return self.__class__(*[self.f.subtract(t1, t2) for t1, t2 in zip_longest(c1, c2, fillvalue=0)], f=self.f)

    def __neg__(self):
        return self.__class__(*[self.f.negate(coeff) for coeff in self.coefficients], f=self.f, flipped=True)

    def __str__(self):
        res = ""
        degree = len(self.coefficients) - 1
        if degree > 0:
            res += str(self.coefficients[0]) + "x^" + str(degree)
            for i in range(1, len(self.coefficients) - 1):
                coeff = self.coefficients[i]
                if coeff < 0:
                    res += " - " + str(-coeff) + "x^" + str(degree - i)
                else:
                    res += " + " + str(coeff) + "x^" + str(degree - i)

            if self.coefficients[-1] < 0:
                res += " - " + str(-self.coefficients[-1])
            else:
                res += " + " + str(self.coefficients[-1])
        else:
            res = str(self.coefficients[-1])

        return res

    def __eq__(self, other):
        if isinstance(other, Polynomial):
            self.trim()
            other.trim()
            if self.coefficients == other.coefficients and self.f == other.f:
                return True
            else:
                return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def trim(self):
        """
        removes all leading 0 terms from the coefficients array
        """
        while self.coefficients and self.coefficients[0] == 0:
            self.coefficients.pop(0)
        if len(self.coefficients) == 0:
            self.coefficients = [0]


def zip_longest(iter1, iter2, fillvalue=None):
    # I did not create this. Found with original framework for Polynomial class. used for __add__ and __sub__
    for i in range(max(len(iter1), len(iter2))):
        if i >= len(iter1):
            yield (fillvalue, iter2[i])
        elif i >= len(iter2):
            yield (iter1[i], fillvalue)
        else:
            yield (iter1[i], iter2[i])
        i += 1


if __name__ == '__main__':
    f = fieldmath.BinaryField(0x11d)
    p = PolynomialOperations(f)
    polyn1 = [0, 1, 2, 3, 0, 0]
    polyn2 = [3, 2, 3, 2, 0, 0, 1]
    poly1_class = Polynomial(*polyn1, f=f)
    poly2_class = Polynomial(*polyn2, f=f)

    print(poly_gcd(poly2_class, poly1_class))
    print(p.poly_gcd(polyn2, polyn1))

