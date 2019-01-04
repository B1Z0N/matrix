from itertools import permutations
from copy import deepcopy
from random import random
from math import sqrt
from collections.abc import Sequence
from functools import singledispatch
import numbers
import reprlib

#TODO: Comment and test
#TODO: write rank function

"""Exceptions"""
#Blocks of exception classes
class MatrixError(Exception):
    """Base class for exceptions in this module"""
#1
class RectangularError(MatrixError):
    """Thrown when the matrix is not rectangular"""
    def __init__(self, mat):
        print("Not rectangular matrix {0}".format(mat))
#2
class SquareError(MatrixError):
    """Thrown when the matrix is not square, but function is 
            defined on square matrices only"""
    def __init__(self, m, reason):
        self.m = m
        self.reason = reason

    def __str__(self):

        return("Not able to take {0} of non-square matrix {1}".format( \
            self.reason, self.m))
#4
class ArithmeticMatrixError(MatrixError):
    """General class for mul and add dimension error"""
    def __init__(self, m1, m2, operation):
        self.m1 = m1
        self.m2 = m2
        self.operation = operation
    
    def __str__(self):
        return ("It is impossible to {0} {1}x{2} and {3}x{4} matrices".format( \
            self.operation, self.m1.rows(), self.m1.cols(), \
                self.m2.rows(), self.m2.cols()))
#5
class SOLEError(MatrixError): #A * x = b
    """General class for SOLE errors"""
    def __init__(self, A, b, err):
        self.A = A
        self.b = b
        self.operation = err
    
    def __str__(self):
        return ("The error is '{0}':\n {1} * \n(X vector)\n = {3}".format( \
            self.err, self.A, self.x, self.b))
#6    
class NotConvergableError(MatrixError):
    def __init__(self, m, helpinf, operation):
        self.m = m
        self.help = helpinf
        self.operation = operation

    def __str__(self):

        return("Not able to perform {0} with {1} help".format( \
            self.operation, self.m, self.help))

"""Generation functions"""
def seq_to_type(seq, totype = float):
    if not issubclass(totype, numbers.Number):
        raise TypeError

    return [[totype(j) for j in i] for i in seq]

def gen_id_mat(y, x = None, scalartype = float):
    if not x:
        x = y
    if not issubclass(scalartype, numbers.Number):
        raise TypeError

    check_dim(y, x)

    mat = [[scalartype(i == j) for i in range(x)] for j in range(y)]

    return(mat)

def gen_zero_mat(y, x = None, scalartype = float):
    if not x:
        x = y
    if not issubclass(scalartype, numbers.Number):
        raise TypeError

    check_dim(y, x)

    mat = [[scalartype(0) for i in range(x)] for j in range(y)]

    return(mat)

def gen_rand_mat(y, x = None, start = -100, end = 100, scalartype = float):
    if not x:
        x = y
    if not issubclass(scalartype, numbers.Number):
        raise TypeError

    check_dim(y, x)

    randrange = lambda: random() * (end - start) + start
    if scalartype == complex:
        matrix = [[scalartype(randrange(), randrange()) for j in range(x)] for i in range(y)]
    else:
        matrix = [[scalartype(randrange()) for j in range(x)] for i in range(y)]
                      
    return(matrix)
"""Tests"""
def is_rectangular(lst):
    len(lst) #checking if it has a len method
    width = len(lst[0])
    for i in lst:
        if len(i) != width:
            raise RectangularError(lst)

def check_dim(row, col):
    if type(row) != int:
        raise TypeError("Rows can't be of type %s" % type(row))
    if type(col) != int:
        raise TypeError("Columns can't be of type %s" % type(col))

    if row <= 0:
        raise ValueError("Invalid number of rows %d" % row)
    if col <= 0:
        raise ValueError("Invalid number of columns %d" % col)

"""Classes"""
####################################################################
class MatrixMathMixin:
    """LU solution"""
    def pivotize(self):#create P matrix of PLU decomposition
        if self.cols() != self.rows():
            raise SquareError(self, "pivotize")
        
        n = self.cols()
        P = matr(gen_id_mat(n))
        swaps = 0
        for j in range(n):
            row = max(range(j, n), key = lambda i: abs(self[i][j]))
            if j != row:
                P[j], P[row] = P[row], P[j]
                swaps += 1

        return (P, swaps)
        
    def PLU(self):#PA = LU
        P, _  = self.pivotize()
        A = P @ self
        L, U = A.LU()
        
        return(P, L, U)
        
    def LU(self):
        if self.rows() != self.cols():
            raise SquareError(self, "LU")
        
        n = self.rows()
        L = matr(n)
        U = matr(n)
        
        for i in range(n):
            for k in range(i, n):
                s = sum([L[i][j] * U[j][k] for j in range(i)])
                U[i][k] = self[i][k] - s
            
            for k in range(i, n):
                if i == k:
                    L[i][i] = 1.0
                else:
                    if U[i][i] == 0.0:
                        raise ZeroDivisionError(self)
                    s = sum([L[k][j] * U[j][i] for j in range(i)])
                    L[k][i] = (self[k][i] - s) / U[i][i]

        return (L, U)
    
    def LU_solve(self, right): #A * X = B
        if right.rows() != self.rows():
            raise SOLEError(self, right, "wrong dimensions of the righthand vector")
        
        y = right.rows()
        x = right.cols()#dims
        P, L, U = self.PLU()
        
        solution = matr(y, x)
        res_col = 0
        
        for col in right.get_cols():
        #(P^-1)LUx = b, (P^-1)x1 = b, x1 = P * b
            x1 = P @ matr(col)

            #LUx = x1, Lx2 = x1
            x2 = matr(y, 1)
            for i in range(y):
                s = sum([L[i][k] * x2[k][0] for k in range(i)])
                x2[i][0] = x1[i][0] - s

            #Ux = x2
            for i in range(y - 1, -1, -1):
                s = sum([U[i][k] * solution[k][res_col] for k in range(y - 1, i, -1)])
                solution[i][res_col] = (x2[i][0] - s) / U[i][i]
                
            res_col += 1

        return(solution)
    """Iterative methods"""
    def residual(self, X, B):
        
        return(self @ X - B)
    
    def Gauss_Seidel(self, right, eps = 10e-5): #self * X = right
        if right.rows() != self.rows():
            raise SOLEError(self, right, "wrong dimensions of the righthand vector")

        totype = type(self[0][0])
        
        n = right.rows()
        x = matr([[totype(0.0)] for i in range(n)])

        norms = {}
        norm = self.spectral_norm()
        if  norm >= totype(1):
            norms['spec'] = norm
            norm = self.vec_norm1()
            if norm >= totype(1):
                norms['vec1'] = norm        #collecting norms info
                norm = self.vec_norm2()
                if norm >= totype(1):
                    norms['vec2'] = norm
                    
        conv = False
        try:
            while not conv:
                nx = deepcopy(x)
                for i in range(n):
                    s1 = sum([self[i][j] * nx[j][0] for j in range(i)])
                    s2 = sum([self[i][j] * x[j][0] for j in range(i + 1, n)])
                    nx[i][0] = (right[i][0] - s1 - s2) / self[i][i]

                conv = self.residual(nx, right).spectral_norm() <= eps
                x = nx
        except OverflowError:
            raise NotConvergableError(self, str(norms), "Gauss Seidel")

        return(x)
    
    def spectral_norm(self):
        y = self.rows()
        x = self.cols()

        return(sqrt(sum(self[i][j] ** 2 for i in range(y) \
                    for j in range(x))))
    
    def vec_norm1(self):
        y = self.rows()
        x = self.cols()
        
        return( \
            max( \
                sum( \
                    abs(self[i][j]) for j in range(x) \
                    ) \
                for i in range(y)))
            
    
    def vec_norm2(self):
        y = self.rows()
        x = self.cols()
        
        return( \
            max( \
                sum( \
                    abs(self[i][j]) for i in range(y) \
                    ) \
                for j in range(x)))
    """Eigenvalues"""
    def power_iteration(self, accuracy = 0.001):
        def eigen_value(matr, acc):
            y1 = matr @ gen_id_vec(matr.rows())
            max_ev = 0.0
            try:
                while(True):
                    y2 = matr @ y1
                    eps = abs(max_ev - y2[0][0]/y1[0][0])
                    max_ev = y2[0][0]/y1[0][0]
                    if eps < acc:
                        return(max_ev)
                    y1 = y2
                      
                return max_ev
            except OverflowError as OE:
                norm_msg = "{'norm' : {} }"
                raise NotConvergableError(self, norm_msg.format(str(norm)), "Power Iteration")

        max_ev = eigen_value(self, accuracy)
        min_ev = 1 / eigen_value(~self, accuracy * 0.1)       

        return (max_ev, min_ev)

    def Frobenius_form(self):
        A = self
        dim = self.rows()

        for n in reversed(range(self.rows() - 1)):
            similar = matr(gen_id_mat(dim))
            for i in range(dim):
                similar[n][i] = -A[n + 1][i] / A[n + 1][n]
            similar[n][n] = 1 / A[n + 1][n]

            B = matr(gen_id_mat(dim))
            for i in range(dim):
                for j in range(dim):
                    B[i][j] = A[i][j] + A[i][n] * similar[n][j]
            for i in range(dim):
                B[i][n] = A[i][n] * similar[n][n]

            reversed_similar = matr(gen_id_mat(dim))
            for i in range(dim):
                reversed_similar[n][i] = A[n + 1][i]
    
            A = reversed_similar * B
            print(A)

        return(A)

    def Frob_charach(self, decimals = 3):#eigenvalues
        A = self.Frobenius_form()
        dim = A.rows()
        sign = (-1) ** dim

        cursign = lambda x: ('-', '+')[x > 0]
        delone = lambda x: (cursign(x) + ' ' + cursign(round(x, decimals)), cursign(round(x, decimals)))[abs(x) == 1]
        
        charach = (delone(sign) + 'x^'  + str(dim))[1:]


        i = dim
        for element in A.get_row(0)[0]:
            charach += ' ' + delone(sign * element) + 'x^'  + str(i)
            i -= 1

        return(charach, A.get_row(0)[0])

    def find_evecs(self, eigen):
        dim = self.rows()
        A = self - eigen * matr(gen_id_mat(dim))
        evec = A.LU_solve(matr([0.0] * dim))
            
        return evec
##################################################################################
class CustomConstructorsMixin:
    @classmethod
    def from_string(cls, string, delimiter = ';'):
        return cls([m.split() for m in string.split(delimiter)[:-1]])
    
    @classmethod
    def from_file(cls, fpath, delimiter = '\n'):
        s = ''
        with open(fpath) as f:
            for line in f:
                s += line

        return cls.from_string(s, delimiter = '\n')

    def to_string(self, delimiter = ';'):
        s = ''
        for i in self:
            for j in i:
                s += str(j) + ' '
            s = s[:-1] + delimiter
        return s
    
    def to_file(self, fpath):
        s = self.to_string(delimiter = '\n')
        with open(fpath, 'w') as f:
            f.write(s)
###############################################################################################
class matr(MatrixMathMixin, CustomConstructorsMixin):
    """Different constructors"""
    def __init__(self, x, y = None, scalartype = float):
        @singledispatch
        def matr_init(arg):
            return NotImplemented
        @matr_init.register(numbers.Number)
        def _(arg):
            try:
                return gen_id_mat(arg, scalartype = scalartype)
            except TypeError:
                raise TypeError("argument must be of positive integer type")
        @matr_init.register(Sequence)
        def _(arg):
            try:
                is_rectangular(arg)
            except TypeError:
                raise TypeError("argument must be two-dimensional iterable")
            return [list(i) for i in arg]

        if y == None:
                self.matrix = matr_init(x)
                if self.matrix == NotImplemented:
                    raise TypeError("wrong arguments of type {!r} and {!r}".format(type(x), type(y)))
        else:
            try:
                self.matrix = gen_id_mat(x, y, scalartype = scalartype)
            except TypeError:
                raise TypeError("arguments must be of positive integer type")
       
    """Get elements"""
    def mat_to_file(self, fpath):
        y = self.rows()
        x = self.cols()
        try:
            with open(fpath, "w") as f:
                for i in range(y):
                    for j in range(x):
                        if j < x - 1:
                            print(self[i][j], end = ' ', file = f)
                        elif i < y - 1:
                            print(self[i][j], end ='\n', file = f)
                        else:
                            print(self[i][j], end = '', file = f)
        except IOError as er:
            print(er)

        return
                    
    def cols(self):
        return(len(self.matrix[0]))

    def rows(self):
        return(len(self.matrix))

    def minor(self, num_lst):  
    #num_lst is list of points(rows and columns) which we need to slice
        row = (num_lst[j][0] for j in range(len(num_lst)))#to get our minor
        col = (num_lst[j][1] for j in range(len(num_lst)))

        minor = [[self.matrix[i][j] \
                for j in range(self.cols()) if j not in col] \
                for i in range(self.rows()) if i not in row]

        return(matr(minor))

    def mat_2_lst(self):
        
        return(self.matrix)

    def get_col(self, num):
        col = [[self[i][num]] for i in range(self.rows())]

        return col

    def get_row(self, num):
    
        return ([self[num]])
    def get_cols(self):
        for i in range(self.cols()):
            yield self.get_col(i)
            
    def get_rows(self):
        for i in range(self.rows()):
            yield self.get_row(i)
    """Operations"""
    
    def rank():
        pass
        
    def trace (self):
        if self.rows() != self.cols():
            raise SquareError(self, "trace")

        return( sum( [ int(i == j) * self[j][i] 
                for i in range(self.cols()) 
                for j in range(self.rows()) ] ) )
    
    def det(self):
        if self.rows() != self.cols():
            raise SquareError(self, "det")

        P , swaps  = self.pivotize()
        
        try:
            _ , U = (P @ self).LU()
        except ZeroDivisionError:
            return(0)
        
        det = 1
        for i in range(self.rows()):
            det *= U[i][i]
        
        return(det*((-1)**(swaps)))

    def transpose(self):
        is_rectangular(self)

        trans = [[self[j][i] for j in \
                  range(self.rows())] for i in \
                  range(self.cols())]
        
        return(matr(trans))

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        if self.rows() != other.rows() \
            or self.cols() != other.cols():
                raise ArithmeticMatrixError(self, other, "equality checking")
        
        rows = self.rows()
        cols = self.cols()

        return all(self[i][j] == other[i][j] for j in range(cols) for i in range(rows))

    def __repr__(self):
        if self.rows() >= 20:
            end = 9
            start = self.rows() - 9
        else:
            end = self.rows()
            start = end
            
        s = '\n'
        for i in range(end):
            s1 = reprlib.repr(self[i])
            s1 = 7 * ' ' + s1
            s += s1 + '\n'
        
        if start != end:
            mid = int((len(s1) - 7) / 2 + 7)
            
            s += mid * ' ' + '...\n'
            
        for i in range(start, self.rows()):
            s1 = reprlib.repr(self[i])
            s1 = 7 * ' ' + s1
            s += s1 + '\n'
            
        return "Matrix ({}".format(s) + 7 * ' ' + ')'
    
    def __mul__(self, other):
        if not isinstance(other, numbers.Number):
            return NotImplemented
        m = matr([[self[i][j] * other for j in range(self.cols())] for i in range(self.rows())])
        return matr(m)

    def __rmul__(self, other):
        
        return self * other

    def __matmul__(self, other):
        if not isinstance(other, matr) or isinstance(other, Sequence):
            return NotImplemented
        else:
            if self.cols() != other.rows():
                raise ArithmeticMatrixError(self, other, "matrix multiplication")

            m = [
                    [ \
                        sum([self[i][k] * other[k][j] for k in range(self.cols())]) \
                        for j in range(other.cols()) \
                    ] \
                        for i in range(self.rows()) \
                ]
                        
        return(matr(m))

    def __rmatmul__(self, other):
        return matr(other) @ self
    
    def __add__(self, other):
        if self.rows() != other.rows() \
            or self.cols() != other.cols():
                raise ArithmeticMatrixError(self, other, "matrix addition")

        m = matr([[self[i][j] + other[i][j] for j in range(self.cols())] for i in range(self.rows())])

        return(m)

    def __sub__(self, other):
        return self + -other
    
    def __neg__(self):
        return(-1 * self)   
    
    def __getitem__(self, index):
        cls = type(self)
        if isinstance(index, slice):
            return matr(self.matrix[index])
        elif isinstance(index, numbers.Integral):
            return self.matrix[index]
        else:
            print("{cls.__name__} indices must be integers")
            raise TypeError(cls = cls)

    def __setitem__(self, row, col, val):
        if not isinstance(val, numbers.Number):
            return NotImplemented
        self.matrix[row][col] = val

    def __setitem__(self, row, val):
        if not isinstance(val, numbers.Number):
            return NotImplemented
        self.matrix[row] = val
        
    def __str__(self):
        s = '\n'
        for i in self.matrix:
            s += '['
            s1 = ''
            for j in i:
                s1 += ', ' + str(j) 
            s1 = s1[2:]
            s += s1 + ']\n'
            
        return(s)

    def __invert__(self):
        if self.rows() != self.cols():
            raise SquareError(self, "inversion")
        if self.det() == 0:
            raise ZeroDivisionError
        return(self.LU_solve(matr(gen_id_mat(self.rows()))))

    def __pow__(self, num):
        if not isinstance(num, numbers.Integral):
            raise TypeError("can't raise {} to nonintegral type: {}".format(type(self).__name__, type(num)))

        if num == 0:
            temp = matr(gen_id_mat(self.rows(), self.cols()), scalartype = type(self[0][0]))
        elif num > 0:
            temp = deepcopy(self)
            for i in range(num - 1):
                temp = temp * self
        else:
            temp = (~self) ** (-num) ##temp = (A^-1) ^ (-num) when num < 0
        
        return(temp)
