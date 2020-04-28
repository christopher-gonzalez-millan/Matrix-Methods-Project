
# Modification of fieldmath.py, modified by Ian Jorquera
# Author: Nayuki
# Original Project Title: Reed-Solomon error-correcting code decoder (Python)
# Copyright (c) 2020 Project Nayuki
# All rights reserved. Contact Nayuki for licensing.
# https://www.nayuki.io/page/reed-solomon-error-correcting-code-decoder

# Modified by Ian Jorquera
# I have personally added to following function to the Field class:
#       __eq__ and __ne__
# and the following to the Matrix class:
#       get_sub_matrix, to_list, copy, __mul__, transpose, any, and kernel_space
# I have also added these addition function:
#       pow_over_field, identity_n, create_matrix, augmented_a_b_matrix, solve_ax_b and solve_lstsq
# Slight modification may have also been made to other existing functions not listed above


# ---- Field abstract class ----
class Field:
    """An abstract base class representing a field in abstract algebra. Every field must
    satisfy all these axioms, where x, y, z are arbitrary elements of the field:
    - 0 is an element of the field, and 0 + x = x. (Existence of additive identity)
    - 1 is an element of the field, and 1 * x = x. (Existence of multiplicative identity)
    - 0 != 1. (Distinctness of additive and multiplicative identities)
    - x + y = y + x. (Commutativity of addition)
    - x * y = y * x. (Commutativity of multiplication)
    - (x + y) + z = x + (y + z). (Associativity of addition)
    - (x * y) * z = x * (y * z). (Associativity of multiplication)
    - x * (y + z) = (x * y) + (x * z). (Distributivity of multiplication over addition)
    - -x is an element of the field, such that x + (-x) = 0. (Existence of additive inverse)
    - x^-1 is an element of the field, such that x * (x^-1) = 1. (Existence of multiplicative inverse)
    Each Field object should be stateless and immutable. The field element objects should be immutable too."""

    # -- Constant values --

    def zero(self):
        """Returns the additive identity constant of this field."""
        raise AssertionError("Not implemented")

    def one(self):
        """Returns the multiplicative identity constant of this field."""
        raise AssertionError("Not implemented")

    # -- Comparison --

    def equals(self, x, y):
        """Tests whether the two given elements are equal.
        Note that the elements are not required to implement their own __eq__() correctly.
        This means x == y is allowed to mismatch f.equals(x, y)."""
        raise AssertionError("Not implemented")

    # -- Addition/subtraction --

    def negate(self, x):
        """Returns the additive inverse of the given element."""
        raise AssertionError("Not implemented")

    def add(self, x, y):
        """Returns the sum of the two given elements."""
        raise AssertionError("Not implemented")

    def subtract(self, x, y):
        """Returns the difference of the two given elements.
        A correct default implementation is provided."""
        return self.add(x, self.negate(y))

    # -- Multiplication/division --

    def reciprocal(self, x):
        """Returns the multiplicative inverse of the given non-zero element."""
        raise AssertionError("Not implemented")

    def multiply(self, x, y):
        """Returns the product of the two given elements."""
        raise AssertionError("Not implemented")

    def divide(self, x, y):
        """Returns the quotient of the given elements.
        A correct default implementation is provided."""
        return self.multiply(x, self.reciprocal(y))

    def __eq__(self, other):
        raise AssertionError("Not implemented")

    def __ne__(self, other):
        return not self.__eq__(other)


# ---- PrimeField class ----
class PrimeField(Field):
    """A finite field of the form Z_p, where p is a prime number.
    Each element of this kind of field is an integer in the range [0, p).
    Both the field and the elements are immutable and thread-safe."""

    def __init__(self, mod):
        """Constructs a prime field with the given modulus. The modulus must be a
        prime number, but this crucial property is not checked by the constructor."""
        if mod < 2:
            raise ValueError("Modulus must be prime")
        # The modulus of this field, which is also the number of elements in this finite field. Must be prime.
        self.modulus = mod

    def zero(self):
        return 0

    def one(self):
        return 1

    def equals(self, x, y):
        return self._check(x) == self._check(y)

    def negate(self, x):
        return -self._check(x) % self.modulus

    def add(self, x, y):
        return (self._check(x) + self._check(y)) % self.modulus

    def subtract(self, x, y):
        return (self._check(x) - self._check(y)) % self.modulus

    def multiply(self, x, y):
        return (self._check(x) * self._check(y)) % self.modulus

    def reciprocal(self, w):
        # Extended Euclidean GCD algorithm
        x = self.modulus
        y = self._check(w)
        if y == 0:
            raise ValueError("Division by zero")
        a = 0
        b = 1
        while y != 0:
            q, r = x // y, x % y
            x, y = y, r
            a, b = b, (a - q * b)
        if x == 1:
            return a % self.modulus
        else:  # All non-zero values must have a reciprocal
            raise AssertionError("Field modulus is not prime")

    # Checks if the given object is the correct type and within
    # the range of valid values, and returns the value itself.
    def _check(self, x):
        if not isinstance(x, int):
            raise TypeError()
        if not (0 <= x < self.modulus):
            raise ValueError("Not an element of this field: " + str(x))
        return x


# ---- BinaryField class ----
class BinaryField(Field):
    """A Galois field of the form GF(2^n/mod). Each element of this kind of field is a
    polynomial of degree less than n where each monomial coefficient is either 0 or 1.
    Both the field and the elements are immutable and thread-safe."""

    def __init__(self, mod):
        """Constructs a binary field with the given modulus. The modulus must have
        degree at least 1. Also the modulus must be irreducible (not factorable) in Z_2,
        but this critical property is not checked by the constructor."""
        if mod <= 1:
            raise ValueError("Invalid modulus")

        # The modulus of this field represented as a string of bits in natural order.
        # For example, the modulus x^5 + x^1 + x^0 is represented by the integer value 0b100011 (binary) or 35 (decimal)
        self.modulus = mod

        # The number of (unique) elements in this field. It is a positive power of 2, e.g. 2, 4, 8, 16, etc.
        # The size of the field is equal to 2 to the power of the degree of the modulus.
        self.size = 1 << (mod.bit_length() - 1)

        self.mults = [[None] * self.size for _ in range(self.size)]
        self.recips = [None] * self.size
        self.pows = [[None] * self.size for _ in range(self.size)]
        self.init_tables()

    def init_tables(self):
        for i in range(0, self.size):
            for j in range(0, self.size):
                self.mults[i][j] = self.multiply(i, j)

        for i in range(1, self.size):
            self.recips[i] = self.reciprocal(i)

        for i in range(self.size):
            for j in range(self.size):
                self.pows[i][j] = pow_over_field(i, j, self)

    def zero(self):
        return 0

    def one(self):
        return 1

    def equals(self, x, y):
        return self._check(x) == self._check(y)

    def negate(self, x):
        return self._check(x)

    def add(self, x, y):
        return self._check(x) ^ self._check(y)

    def subtract(self, x, y):
        return self.add(x, y)

    def multiply(self, x, y):
        self._check(x)
        self._check(y)
        if self.mults[x][y] is not None:
            return self.mults[x][y]
        result = 0
        while y != 0:
            if y & 1 != 0:
                result ^= x
            x <<= 1
            if x >= self.size:
                x ^= self.modulus
            y >>= 1
        return result

    def reciprocal(self, w):
        # Extended Euclidean GCD algorithm
        x = self.modulus
        y = self._check(w)
        if self.recips[w] is not None:
            return self.recips[w]

        if y == 0:
            raise ValueError("Division by zero")
        a = 0
        b = 1
        while y != 0:
            q, r = self._divide_and_remainder(x, y)
            if q == self.modulus:
                q = 0
            x, y = y, r
            a, b = b, (a ^ self.multiply(q, b))
        if x == 1:
            return a
        else:  # All non-zero values must have a reciprocal
            raise AssertionError("Field modulus is not irreducible")

    # Returns a new tuple containing the pair of values (x div y, x mod y).
    def _divide_and_remainder(self, x, y):
        quotient = 0
        ylen = y.bit_length()
        for i in reversed(range(x.bit_length() - ylen + 1)):
            if x.bit_length() == ylen + i:
                x ^= y << i
                quotient |= 1 << i
        return (quotient, x)

    # Checks if the given object is the correct type and within the
    # range of valid values, and returns the same value.
    def _check(self, x):
        if not isinstance(x, int):
            raise TypeError()
        if not (0 <= x < self.size):
            raise ValueError("Not an element of this field: " + str(x))
        return x

    def __eq__(self, other):
        if isinstance(other, BinaryField):
            if self.modulus == other.modulus:
                return True
        return False


# ---- Additional field operations ----
# Created by Ian Jorquera
def pow_over_field(base, exp, field):
    if exp < 0:
        raise ValueError("Exp can not be negative.")
    if not isinstance(field, Field):
        raise Exception("Must provide a valid Field.")
    if isinstance(field, BinaryField) and field.pows[base][exp] is not None:
        return field.pows[base][exp]
    if isinstance(field, BinaryField) and field.pows[base][exp-1] is not None:
        return field.multiply(base, field.pows[base][exp-1])

    res = field.one()
    for _ in range(exp):
        res = field.multiply(base, res)
    return res


# ---- Matrix class ----
# I have added a decent amount of additional function using the exiting ones
class Matrix:
    """Represents a mutable matrix of field elements, supporting linear algebra operations.
    Note that the dimensions of a matrix cannot be changed after construction. Not thread-safe."""

    def __init__(self, rows, cols, field, zeros=False):
        """Constructs a blank matrix with the given number of rows and columns,
        with operations from the given field. All the elements are initially None."""
        if rows <= 0 or cols <= 0:
            raise ValueError("Invalid number of rows or columns")
        if not isinstance(field, Field):
            raise TypeError()

        # The field used to operate on the values in the matrix.
        self.f = field
        # The values of the matrix stored in row-major order, with each element initially None.
        if zeros:
            self.values = [[self.f.zero()] * cols for _ in range(rows)]
        else:
            self.values = [[None] * cols for _ in range(rows)]

    # -- Basic matrix methods --

    def row_count(self):
        """Returns the number of rows in this matrix, which is a positive integer."""
        return len(self.values)

    def column_count(self):
        """Returns the number of columns in this matrix, which is a positive integer."""
        return len(self.values[0])

    def get(self, row, col):
        """Returns the element at the given location in this matrix. The result may be None."""
        if not (0 <= row < len(self.values) and 0 <= col < len(self.values[row])):
            raise IndexError("Row or column index out of bounds")
        return self.values[row][col]

    def get_sub_matrix(self, row_i, row_t, col_i, col_t):
        if row_i is None:
            row_i = 0
        if row_t is None:
            row_t = self.row_count()
        if col_i is None:
            col_i = 0
        if col_t is None:
            col_t = self.column_count()

        if row_t <= row_i or col_t <= col_i:
            raise Exception("Invalid parameters. Terminator can not be leq init.")
        result = Matrix(row_t - row_i, col_t - col_i, self.f)
        for r in range(row_i, row_t):
            for c in range(col_i, col_t):
                result.set(r - row_i, c - col_i, self.get(r, c))
        return result

    def set(self, row, col, val):
        """Stores the given element at the given location in this matrix. The value to store can be None."""
        if not (0 <= row < len(self.values) and 0 <= col < len(self.values[row])):
            raise IndexError("Row or column index out of bounds")
        self.values[row][col] = val

    def to_list(self, single=False):
        lst = []
        for r in range(self.row_count()):
            row = []
            for c in range(self.column_count()):
                if single:
                    lst.append(self.get(r, c))
                else:
                    row.append(self.get(r, c))
            if not single:
                lst.append(row)
        return lst

    def __str__(self):
        """Returns a string representation of this matrix. The format is subject to change."""
        result = "["
        for (i, row) in enumerate(self.values):
            if i > 0:
                result += ",\n "
            result += "[" + ", ".join(str(val) for val in row) + "]"
        return result + "]"

    def __mul__(self, other):
        if not isinstance(other, Matrix):
            raise TypeError()
        if self.f != other.f:
            raise Exception("Fields do align.")

        if self.column_count() != other.row_count():
            raise Exception("Can not multiple matrices, inner dimensions do no align.")

        result = Matrix(self.row_count(), other.column_count(), self.f)
        for r in range(result.row_count()):
            for c in range(result.column_count()):
                val = self.f.zero()
                for i in range(self.column_count()):
                    val = self.f.add(val, self.f.multiply(self.get(r, i), other.get(i, c)))
                result.set(r, c, val)
        return result

    def transpose(self):
        result = Matrix(self.column_count(), self.row_count(), self.f)
        for r in range(result.row_count()):
            for c in range(result.column_count()):
                result.set(r, c, self.get(c, r))
        return result

    def any(self):
        for r in range(self.row_count()):
            for c in range(self.column_count()):
                if not self.f.equals(self.get(r, c), self.f.zero()):
                    return True
        return False

    def copy(self):
        result = Matrix(self.row_count(), self.column_count(), self.f)
        result.values = self.values
        return result

    # -- Simple matrix row operations --

    def swap_rows(self, row0, row1):
        """Swaps the two given rows of this matrix. If the two row indices are the same, the swap is a no-op.
        Any matrix element can be None when performing this operation."""
        if not (0 <= row0 < len(self.values) and 0 <= row1 < len(self.values)):
            raise IndexError("Row index out of bounds")
        self.values[row0], self.values[row1] = self.values[row1], self.values[row0]

    def multiply_row(self, row, factor):
        """Multiplies the given row in this matrix by the given factor. In other words, row *= factor.
        The elements of the given row should all be non-None when performing this operation."""
        if not (0 <= row < len(self.values)):
            raise IndexError("Row index out of bounds")
        self.values[row] = [self.f.multiply(val, factor) for val in self.values[row]]

    def add_rows(self, srcrow, destrow, factor):
        """Adds the first given row in this matrix multiplied by the given factor to the second given row.
        In other words, destdow += srcrow * factor. The elements of the given two rows
        should all be non-None when performing this operation."""
        if not (0 <= srcrow < len(self.values) and 0 <= destrow < len(self.values)):
            raise IndexError("Row index out of bounds")
        self.values[destrow] = [self.f.add(destval, self.f.multiply(srcval, factor))
                                for (srcval, destval) in zip(self.values[srcrow], self.values[destrow])]

    # -- Advanced matrix operations --
    def row_echelon_form(self):
        rows = self.row_count()
        cols = self.column_count()

        # Compute row echelon form (REF)
        numpivots = 0
        for j in range(cols):  # For each column
            if numpivots >= rows:
                break
            pivotrow = numpivots
            while pivotrow < rows and self.f.equals(self.get(pivotrow, j), self.f.zero()):
                pivotrow += 1
            if pivotrow == rows:
                continue  # Cannot eliminate on this column
            self.swap_rows(numpivots, pivotrow)
            pivotrow = numpivots
            numpivots += 1

            # Simplify the pivot row
            # Eliminate rows below
            for i in range(pivotrow + 1, rows):
                self.add_rows(pivotrow, i, self.f.negate(self.f.multiply(self.f.reciprocal(self.get(pivotrow, j)), self.get(i, j))))

    def reduced_row_echelon_form(self, only_n_pivots=None):
        """Converts this matrix to reduced row echelon form (RREF) using Gauss-Jordan elimination.
        All elements of this matrix should be non-None when performing this operation.
        Always succeeds, as long as the field follows the mathematical rules and does not raise an exception.
        The time complexity of this operation is O(rows * cols * min(rows, cols))."""
        rows = self.row_count()
        if only_n_pivots:
            cols = only_n_pivots
        else:
            cols = self.column_count()

        # Compute row echelon form (REF)
        numpivots = 0
        for j in range(cols):  # For each column
            if numpivots >= rows:
                break
            pivotrow = numpivots
            while pivotrow < rows and self.f.equals(self.get(pivotrow, j), self.f.zero()):
                pivotrow += 1
            if pivotrow == rows:
                continue  # Cannot eliminate on this column
            self.swap_rows(numpivots, pivotrow)
            pivotrow = numpivots
            numpivots += 1

            # Simplify the pivot row
            self.multiply_row(pivotrow, self.f.reciprocal(self.get(pivotrow, j)))

            # Eliminate rows below
            for i in range(pivotrow + 1, rows):
                self.add_rows(pivotrow, i, self.f.negate(self.get(i, j)))

        # Compute reduced row echelon form (RREF)
        for i in reversed(range(numpivots)):
            # Find pivot
            pivotcol = 0
            while pivotcol < cols and self.f.equals(self.get(i, pivotcol), self.f.zero()):
                pivotcol += 1
            if pivotcol == cols:
                continue  # Skip this all-zero row

            # Eliminate rows above
            for j in range(i):
                self.add_rows(i, j, self.f.negate(self.get(j, pivotcol)))

    def kernel_space(self):
        """
        I used a rather simple method for doing this as explained in the Computation by Gaussian elimination section
        found here https://en.wikipedia.org/wiki/Kernel_(linear_algebra).
        """
        ai_matrix = augmented_a_b_matrix(self.transpose(), identity_n(self.column_count(), self.f))
        ai_matrix.reduced_row_echelon_form(only_n_pivots=self.row_count())
        res = []
        for r in range(ai_matrix.row_count()):
            valid = True
            for c in range(self.row_count()):
                if not self.f.equals(ai_matrix.get(r, c), self.f.zero()):
                    valid = False
                    break
            if valid:
                res.append([ai_matrix.get(r, c) for c in range(self.row_count(), ai_matrix.column_count())])

        if len(res) == 0 or len(res[0]) == 0:
            return 0
        result = Matrix(len(res), len(res[0]), self.f)
        for r in range(len(res)):
            for c in range(len(res[0])):
                result.set(r, c, res[r][c])
        return result.transpose()


# ---- Additional Functions for Matrix class ----
# Created by Ian Jorquera
def identity_n(n, field):
    """
    returns a nxn identity matrix
    """
    if not isinstance(field, Field):
        raise TypeError()

    result = Matrix(n, n, field)
    for r in range(result.row_count()):
        for c in range(result.column_count()):
            if r == c:
                result.set(r, c, field.one())
            else:
                result.set(r, c, field.zero())
    return result


def create_matrix(lst, field):
    """
    Helper function to more easily initialize a matrix from a list
    """
    rows = len(lst)
    if rows == 0:
        raise Exception("Invalid Input")
    if isinstance(lst[0], list):
        columns = len(lst[0])
    else:
        columns = 1

    result = Matrix(rows, columns, field)
    for r in range(rows):
        for c in range(columns):
            result.set(r, c, lst[r][c])
    return result


def augmented_a_b_matrix(a, b):
    """
    Simply combines to matrices as a augmented matrix. Returns (A | B)
    """
    if not isinstance(a, Matrix) or not isinstance(b, Matrix):
        raise TypeError()
    if a.f != b.f:
        raise Exception("Fields do no align.")

    if a.row_count() != b.row_count():
        raise Exception(
            f"Matrices do no align. a: {a.row_count()}x{a.column_count()} and b: {b.row_count()}x{b.column_count()}")

    axb = Matrix(a.row_count(), a.column_count() + b.column_count(), a.f)
    for r in range(axb.row_count()):
        for c in range(axb.column_count()):
            if c >= a.column_count():
                val = b.get(r, c - a.column_count())
            else:
                val = a.get(r, c)
            axb.set(r, c, val)
    return axb


def solve_ax_b(a, b):
    """
    Creates an augmented matrix and find the RREF. After that parses solution. Note this does not check for singular
    matrices and tries to provide a particular solution, may have unexpected behavior.
    """

    if not isinstance(a, Matrix) or not isinstance(b, Matrix):
        raise TypeError
    if b.column_count() != 1 or b.row_count() != a.row_count():
        raise Exception("Matrix b must be nx1.")

    axb = augmented_a_b_matrix(a, b)
    axb.reduced_row_echelon_form()

    # now we got to find the solution
    c = 0
    result = Matrix(a.column_count(), 1, a.f, zeros=True)
    for r in range(axb.row_count()):
        c = 0
        while c < axb.column_count():
            if axb.f.equals(axb.get(r, c), axb.f.zero()):
                c += 1
            elif c == axb.column_count() - 1:
                raise Exception("Inconsistent linear system, or singular matrix A")
            else:
                break
        if c < axb.column_count() - 1:
            result.set(c, 0, axb.get(r, axb.column_count() - 1))
    return result


def solve_ax_b_fast(a, b):
    """
    Creates an augmented matrix and find the RREF. After that parses solution. Note this does not check for singular
    matrices and tries to provide a particular solution, may have unexpected behavior.
    """

    if not isinstance(a, Matrix) or not isinstance(b, Matrix):
        raise TypeError
    if b.column_count() != 1 or b.row_count() != a.row_count():
        raise Exception("Matrix b must be nx1.")

    axb = augmented_a_b_matrix(a, b)
    axb.row_echelon_form()

    # now we got to find the solution
    c = 0
    result = Matrix(a.column_count(), 1, a.f, zeros=True)
    if axb.row_count() < axb.column_count() - 1:
        raise Exception("Can not solve, too few rows")
    for r in range(axb.row_count()-1, -1, -1):
        c = (min(axb.column_count() - 2, r))
        if c != r:
            if axb.get(r, c) != 0 or axb.get(r, c + 1) != 0:
                raise Exception("Inconsistent linear system, or singular matrix A")
        else:
            to_sub = axb.f.zero()
            for c2 in range(c + 1, axb.column_count() - 1):
                to_sub = axb.f.add(to_sub, axb.f.multiply(axb.get(r, c2), result.get(c2, 0)))
            result.set(c, 0, axb.f.multiply(axb.f.reciprocal(axb.get(r, c)), axb.f.subtract(axb.get(r, axb.column_count() - 1), to_sub)))
    return result


def solve_lstsq(a, b):
    """
    This function solves the normal equations A^T*A*x = A^T*b. Not this function considers only the 'euclidean dot
    product'
    This really doesnt work with the Binary field. In fact this doesnt really work at all for fields
    """
    if not isinstance(a, Matrix) or not isinstance(a, Matrix):
        raise TypeError

    ata = a.transpose() * a
    atb = a.transpose() * b

    x = solve_ax_b(ata, atb)
    return x


"""Here we will constrct a program to correct error in a (k,n,3) code"""
from scipy.linalg import null_space
import numpy as np


bf= BinaryField(2)
def hamming(a):
    b = str(a)
    cw=list(b)
    k,m=len(cw),0
    while k+m+1>2**m:
        m+=1
    Ik= np.identity(k)
    PP=np.ones((k,m))
    li = False
    while (li == False):
        PP = np.ones((k, m))
        for x in range(0, m):
            z = 0
            for z in range(0, k - m):
                PP[np.random.randint(0, k - 1), x] = 0
                z += 1
        b = null_space(PP)
        li= np.array_equal(b,np.empty((m,0)))
    G=np.hstack((Ik,PP))
    a_matrix = create_matrix(G.tolist(), bf)
    print(a_matrix)
    """Here we create the Generator Matrix that will out put the original k message bits with additional m parity bits
       appended. We construct the functions for each parity bit by generating m linearly independent vectors that will contain
       all possible combinations of bit pairs as its rows.The first four columns will return the original 4 bits and thus
       are the I(k) matrix, the resulting columns will define the bits and will be the aforementioned linearly independent
       basis for all of the parity check bits"""
    x_matrix=create_matrix(cw,bf)
    x_matrix=x_matrix.transpose()

    print(x_matrix)
    bitsfortransmission = x_matrix*a_matrix
    print(bitsfortransmission)

def decode(d):
    r = str(d)
    e = list(r)

    received= np.asarray(cw, dtype=int)
    Pcheck=np.hstack((np.matrix.transpose(PP),np.identity(m)))
    """We create the Parity Check Matrix, the orthogonal complement to the space defined by G. This is expedited using 
    the definition provided by _________"""
    parity=np.matmul(Pcheck,np.transpose(received))
    if parity == np.empty((m,0)):
        print(received[0:3]) ####print first four elements of the input hamming code word
    else:
        w=Pcheck.index(np.transpose(parity))
        if received[w]==1:
            received[w]=0
            print(received[0:3])
        elif received[w]==0:
            received[w]=1
            print(received[0:3])




hamming(1101)
"""Idk how to fix leading decimal issue either for integers when we are putting them into the function,
tried to put other stuff into it but it hasn't been working"""
