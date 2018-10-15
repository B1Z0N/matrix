from itertools import permutations
from copy import deepcopy
from random import random

"""Exceptions"""
class Matrix_Error(Exception):
    """Base class for exceptions in this module"""

class Rectangular_Error(Matrix_Error):
    """Thrown when the matrix is not rectangular"""
    def __init__(self, mat):
        print("Not rectangular matrix %s" % mat)
        
class Square_Error(Matrix_Error):
    """Thrown when the matrix is not square, but function is 
            defined on square matrices only"""
    def __init__(self, m, reason):
        self.m = m
        self.reason = reason

    def __str__(self):

        return("Not able to take %s of non-square matrix %s" % \
        self.reason, self.m)

class Det_Error(Square_Error):
    """Thrown when we are not able to take det"""
    def __init__(self, m):
        Square_Error.__init__(m, "determinant")

class Trace_Error(Square_Error):
    """Thrown when we are not able to take trace"""
    def __init__(self, m):
        Square_Error.__init__(m, "trace")

class Arithmetic_Matrix_Error(Matrix_Error):
    """General class for mul and add dimension error"""
    def __init__(self, m1, m2, operation):
        self.m1 = m1
        self.m2 = m2
        self.operation = operation
    
    def __str__(self):
        return ("It is impossible to %s %dx%d and %dx%d matrices" % \
        self.operation, self.m1.rows(), self.m1.cols(), \
        self.m2.rows(), self.m2.cols())

class Addition_Matrix_Error(Arithmetic_Matrix_Error):
    def __init__(self, m1, m2):
        Arithmetic_Matrix_Error.__init__(self, m1, m2, "add")

class Multiplication_Matrix_Error(Arithmetic_Matrix_Error):
    def __init__(self, m1, m2):
        Arithmetic_Matrix_Error.__init__(self, m1, m2, "muliply")


"""Tests"""
def is_rectangular(mat):
    width = mat.cols()
    for i in mat.matrix:
        if len(i) != width:
            raise Rectangular_Error(mat)

def check_dim(row, col):
    if type(row) != int:
        raise TypeError("Rows value can't be of type %s" % type(row))
    if type(col) != int:
        raise TypeError("Cols value can't be of type %s" % type(col))
    
    if row <= 0:
        raise ValueError("Invalid number of rows %d" % row)

    if col <= 0:
        raise ValueError("Invalid number of columns %d" % col)

"""Class"""
class matr:                             #SLAY in Russian
    """Init, operators, iterators"""
    def __init__(self, *args):
        if len(args) == 1:#one integer
            if isinstance(args[0], list):
                self.matrix = args[0]
                is_rectangular(self)
            else:
                check_dim(args[0], args[0])
                self.matrix = gen_zero_mat(args[0])
        elif len(args) == 2:#two int dims
            rows, cols = args
            check_dim(rows, cols)

            self.matrix = gen_zero_mat(rows, cols)
        else:
            raise ValueError("Matrix takes 1 or 2 args, not %d" % len(args))

        self.index = 0 ##index for iterators in matrix

    """def __init__(self, row, col = -1):## is this ok -1?
        self.matrix = gen_id_mat(row, col)
        self.index = 0"""#make inside of main init
    
    def __iter__(self):

        return(self)
    
    def __next__(self):
        width = self.cols()
        height = self.rows()
        if self.index == width * height:
            self.index = 0
            raise StopIteration
        y = self.index // width
        x = self.index % width
        self.index += 1
        
        return (self[y][x])
    
    def __mul__(self, other):
        m = []

        if is_scalar(other):
            if type(other) != type(self[0][0]):
                raise TypeError("Matrix and scalar must be the same type, \
                not %s and %s" % type(self[0][0]), type(other))
            
            m = [ 
                  [ 
                    self[i][j] * other \
                    for j in range(self.cols()) \
                  ] \
                  for i in range(self.rows()) \
                ]
        else:
            if type(other[0][0]) != type(self[0][0]):
                raise TypeError("Matrices must be the same type, \
                not %s and %s" % type(self[0][0]), type(other[0][0]))
            
            if self.cols() != other.rows():
                raise Multiplication_Matrix_Error(self, other)

            m = [
                    [ \
                        sum([self[i][k] * other[k][j] for k in range(self.cols())]) \
                        for j in range(other.cols()) \
                    ] \
                        for i in range(self.rows()) \
                ]
                        
        return(matr(m))
    
    def __rmul__(self, other):
        
        return(self * other)

    def __add__(self, other):
        if self.rows() != other.rows() \
            or self.cols() != other.cols():
            raise Addition_Matrix_Error
        if not isinstance(other, matr):
            raise TypeError("Not able to add matrix to %s" % type(other))

        m = [ 
              [ 
                self[i][j] + other[i][j] \
                for j in range(self.cols()) \
              ] \
                for i in range(self.rows()) \
            ]

        return(matr(m))

    def __sub__(self, other):
        return(self + -other)
    
    def __neg__(self):
        return(-1 * self)    
    
    def __getitem__(self, index):

        return(self.matrix[index])

    def __setitem__(self, row, col, val):
        self.matrix[row][col] = val
        
    def __str__(self):
        s = ''
        for i in self.matrix:
            s += str(i) + '\n'

        return(s)

    """def __invert__(self):

        for i in self.matrix:
            for j in self.matrix:"""

    def __pow__(self, num):
        if isinstance(num, int) == False:
            raise TypeError("Can't raise matrix to type %s" % type(num))

        if num == 0:
            temp = matr(gen_id_mat(self.rows(), self.cols()))
        elif num > 0:
            temp = deepcopy(self)
            for i in range(num - 1):
                temp = temp * self
        else:
            temp = (~self) ** (-num) ##temp = (A^-1)^(-num) when num < 0
        
        return(temp)

    """Get elements"""
    def cols(self):
        return(len(self.matrix[0]))

    def rows(self):
        return(len(self.matrix))

    def minor(self, num_lst):  
    #num_lst is list of points, rows and columns of which we need to slice
        row = [num_lst[j][0] for j in range(len(num_lst))]#to get our minor
        col = [num_lst[j][1] for j in range(len(num_lst))]

        minor = [[self.matrix[i][j] \
                for j in range(self.cols()) if j not in col] \
                for i in range(self.rows()) if i not in row]

        return(matr(minor))    
    """Operations"""
    def trace (self):
        if self.rows() != self.cols():
            raise Trace_Error(self)

        return( sum( [ int(i == j) * self[j][i] 
                for i in range(self.cols()) 
                for j in range(self.rows()) ] ) )
    
    def det(self):
        if self.rows() != self.cols():
            raise Det_Error(self)
        
        det = 1
        LU = self.LU()
        for i in range(self.rows()):
            det *= (LU[1])[i][i]

        return(det)
        
    def LU(self):
        pass
    
    def transpose(self):
        is_rectangular(self)

        trans = [[self[j][i] for j in \
                  range(self.rows())] for i in \
                  range(self.cols())]
        
        return(matr(trans))
    ##others

"""Support functions""" 
def is_scalar(num):
    
    return(isinstance(num, int)  \
        or isinstance(num, float))

def gen_id_mat(y, x = -1):
    if x == -1:
        x = y
    
    check_dim(y, x)

    mat = [[int(i == j) for i in range(x)] for j in range(y)]

    return(mat)

def gen_zero_mat(y, x = -1):
    if x == -1:
        x = y
        
    check_dim(y, x)

    mat = [[0 for i in range(x)] for j in range(y)]

    return(mat)

def gen_rand_mat(y, x = -1, start = -100, end = 100, isint = 0):
    if x == -1:
        x = y

    check_dim(y, x)

    matrix = []
    for i in range(y):
        matrix.append([(random() * (end - start)) \
                         + start for j in range(x)])
    if isint:
        matrix = [[int(matrix[j][i]) for i in \
        range(len(matrix[0]))] for j in range(len(matrix))]
        
    return(matrix)

class SOLE(matr):
    def __init__(self, matrix):
        matr.__init__(self, matrix)
        self.LU = self.LU()
        self.det = self.det()
        
m1 = matr(gen_rand_mat(y=3,start = 0, end = 5, isint=1))
m2 = matr(gen_rand_mat(y=3,start = 0, end = 5, isint=1))
print(m1)
print(m1.minor([[0,2]]))
print(m2)
print(m2.minor([[2,2],[0,0]]))
print(m1 + -m2)
print((m1 + -m2).transpose())
print(m1 * m2)
print((m1 * m2).transpose())
