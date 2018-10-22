from itertools import permutations
from copy import deepcopy
from random import random
from math import sqrt
CONST_TEST_MAT = [[7.14, 1.02, 1.05, 1.12, -0.96], \
           [1.18, 3.28, 1.3, -1.63, -0.93], \
           [0.89, -2.46, 6.32, 2.1, -0.467],\
           [1.36, 0.16, 2.1, 5.22, 8], \
           [1.44, 0.27, 0.733, -4, 12]]
CONST_TEST_VEC = [[2.1], [0.96], [-1.72], [3.68], [-1.92]]
"""Exceptions"""
#Blocks of exception classes
class Matrix_Error(Exception):
    """Base class for exceptions in this module"""
#1
class Rectangular_Error(Matrix_Error):
    """Thrown when the matrix is not rectangular"""
    def __init__(self, mat):
        print("Not rectangular matrix {0}".format(mat))
#1
#2
class Square_Error(Matrix_Error):
    """Thrown when the matrix is not square, but function is 
            defined on square matrices only"""
    def __init__(self, m, reason):
        self.m = m
        self.reason = reason

    def __str__(self):

        return("Not able to take {0} of non-square matrix {1}".format( \
            self.reason, self.m))

class Det_Error(Square_Error):
    """Thrown when we are not able to take det"""
    def __init__(self, m):
        Square_Error.__init__(self, m, "determinant")

class Inv_Error(Square_Error):
    """Thrown when we are not able to take inverse matrix"""
    def __init__(self, m):
        Square_Error.__init__(self, m, "inverse matrix")

class Trace_Error(Square_Error):
    """Thrown when we are not able to take trace"""
    def __init__(self, m):
        Square_Error.__init__(self, m, "trace")
        
class LU_Square_Error(Square_Error):
    """Thrown when we are not able to take LU"""
    def __init__(self, m):
        Square_Error.__init__(self, m, "LU")

class Pivotize_Error(Square_Error):
    """Thrown when we are not able to pivotize"""
    def __init__(self, m):
        Square_Error.__init__(self, m, "pivotizing")

#2
#3
class LU_Error(Matrix_Error):
    """General class for all exceptions related to LU decomposition"""
    def __init__(self, m, reason, operation):
        self.m = m
        self.reason = reason
        self.operation = operation
        
    def __str__(self):

        return("Not able to perform " + \
                   "{0} with {1} because {2}".format( \
                       self.operation, self.m, self.reason))
        

class LU_Zero_Division_Error(LU_Error):
    def __init__(self, m):
        LU_Error.__init__(self, m, "can't div on zero", \
                          "LU decompose")
#3
#4
class Arithmetic_Matrix_Error(Matrix_Error):
    """General class for mul and add dimension error"""
    def __init__(self, m1, m2, operation):
        self.m1 = m1
        self.m2 = m2
        self.operation = operation
    
    def __str__(self):
        return ("It is impossible to {0} {1}x{2} and {3}x{4} matrices".format( \
            self.operation, self.m1.rows(), self.m1.cols(), \
                self.m2.rows(), self.m2.cols()))

class Addition_Matrix_Error(Arithmetic_Matrix_Error):
    def __init__(self, m1, m2):
        Arithmetic_Matrix_Error.__init__(self, m1, m2, "add")

class Multiplication_Matrix_Error(Arithmetic_Matrix_Error):
    def __init__(self, m1, m2):
        Arithmetic_Matrix_Error.__init__(self, m1, m2, "muliply")
#4
#5
class SOLE_Error(Matrix_Error): #A * x = b
    """General class for SOLE errors"""
    def __init__(self, A, x, b, err):
        self.A = A
        self.x = x
        self.b = b
        self.operation = err
    
    def __str__(self):
        return ("The error is '{0}':\n {1} X {2} = {3}".format( \
            self.err, self.A, self.x, self.b))

class Right_Vector_Error(SOLE_Error):#A * x = b
    """If dims of right vector incompatible with dims of matrix"""
    def __init__(self, A, x, b):
        SOLE_Error.__init__(self, A, x, b, "wrong right vector dims")
#5
#6
class Gauss_Seidel_Error(Matrix_Error):
    """General class for all exceptions related to LU decomposition"""
    def __init__(self, m, reason, operation):
        self.m = m
        self.reason = reason
        self.operation = operation
        
    def __str__(self):

        return("Not able to perform {0} with {1} because {2}".format( \
            self.operation, self.m, self.reason))
               
class Not_Convergable_Error(Gauss_Seidel_Error):
    def __init__(self, m, norm):
        Gauss_Seidel_Error.__init__(self, m, \
        "norms is " + str(norm) + " ( >= 1 )", \
        "Gauss-Seidel itertive method")
#6
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
class matr:
    """Init, operators, iterators"""
    def __init__(self, *args, rows_or_cols = 1):# <= 0 for rows, > 0 for cols
        #trying to guess what arguments were passed
        if isinstance(args[0], list) and is_int_float(args[0][0]):#if there are lists in args
            for i in args:
                if isinstance(i, list):
                    if is_scalar(i) == False:
                        raise TypeError("Neither float nor int type of list")
                else:
                    raise TypeError("Can't join list with %s" % type(i))
                
            if rows_or_cols < 0:    #means rows
                self.matrix = []
                for i in args:
                    self.matrix.append(i)
                    
            else: #means cols
                self.matrix = [[ i[j] for i in args] for j in range(len(args[0]))]

            is_rectangular(self)
        else:    
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
                raise ValueError("Matrix takes 1, 2 args or many lists to join into Matrix, not %d args" % len(args))

            #self.P, self.L, self.U = self.PLU()
        
        self.index = 0 ##index for iterators in matrix
    
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

    def __truediv__(self, other):
        return(self * ~other)
        
    def __mul__(self, other):
        m = []

        if is_int_float(other):
            print(type(other))
            print(type(self[0][0]))
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
        return(-1.0 * self)    
    
    def __getitem__(self, index):

        return(self.matrix[index])

    def __setitem__(self, row, col, val):
        self.matrix[row][col] = val

    def __setitem__(self, row, val):
        self.matrix[row] = val
        
    def __str__(self):
        s = '\n'
        for i in self.matrix:
            s += str(i) + '\n'#[round(i[j], 2) for j in range(len(i))]) + '\n'

        return(s)

    def to_float(self):
        return( matr([[float(self[i][j]) \
                for j in range(self.cols())] \
                for i in range(self.rows())]))
    
    def to_int(self):
        return( matr([[float(int(self[i][j])) \
                for j in range(self.cols())] \
                for i in range(self.rows())]))

    def __invert__(self):
        if self.rows() != self.cols():
            raise Inv_Error(self.matrix)
        if self.det() == 0:
            raise Det_Error(self.matrix)
        return(self.LU_solve(matr(gen_id_mat(self.rows()))))

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
        row = [num_lst[j][0] for j in range(len(num_lst))]#to get our minor
        col = [num_lst[j][1] for j in range(len(num_lst))]

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
            raise Trace_Error(self)

        return( sum( [ int(i == j) * self[j][i] 
                for i in range(self.cols()) 
                for j in range(self.rows()) ] ) )
    
    def det(self):
        if self.rows() != self.cols():
            raise Det_Error(self)

        P , swaps  = self.pivotize()
        
        try:
            _ , U = (P * self).LU()
        except LU_Zero_Division_Error:
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

    """Solving"""
    """LU solution"""
    def pivotize(self):#create P matrix of PLU decomposition
        if self.cols() != self.rows():
            raise Pivotize_Error(self)
        
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
        A = P * self
        L, U = A.LU()
        
        return(P, L, U)
        
    def LU(self):
        if self.rows() != self.cols():
            raise LU_Square_Error(self)
        
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
                        raise LU_Zero_Division_Error(self)
                    s = sum([L[k][j] * U[j][i] for j in range(i)])
                    L[k][i] = (self[k][i] - s) / U[i][i]

        return (L, U)
    
    def LU_solve(self, right): #A * X = B
        if right.rows() != self.rows():
            raise Right_Vector_Error
        
        y = right.rows()
        x = right.cols()#dims
        P, L, U = self.PLU()
        
        solution = matr(y, x)
        res_col = 0
        
        for col in right.get_cols():
        #(P^-1)LUx = b, (P^-1)x1 = b, x1 = P * b
            x1 = P * matr(col)

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
        
        return(self * X - B)
    
    def Gauss_Seidel(self, right, eps = 10e-5): #self * X = right
        if right.rows() != self.rows():
            raise Right_Vector_Error

        n = right.rows()
        x = matr([[0.0] for i in range(n)])

        norms = {}
        norm = self.spectral_norm()
        if  norm >= 1.0:
            norms['spec'] = norm
            norm = self.vec_norm1()
            if norm >= 1.0:
                norms['vec1'] = norm        #collecting norms info
                norm = self.vec_norm2()
                if norm >= 1.0:
                    norms['vec2'] = norm
                    
        conv = False
        try:
            while not conv:
                nx = deepcopy(x)
                for i in range(n):
                    s1 = sum([self[i][j] * nx[j][0] for j in range(i)])
                    s2 = sum([self[i][j] * x[j][0] for j in range(i + 1, n)])
                    nx[i][0] = (right[i][0] - s1 - s2) / A[i][i]

                conv = self.residual(nx, right).spectral_norm() <= eps
                x = nx
        except OverflowError:
            raise Not_Convergable_Error(self, norms)

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
            
"""Support functions""" 
def is_scalar(val):
    if isinstance(val, list):
        for i in val:
            if is_int_float(i) == False:
                return(False)
    else:
        return(is_int_float(val))

def is_int_float(num):
    return(isinstance(num, int)  \
        or isinstance(num, float))

def gen_id_mat(y, x = -1):
    if x == -1:
        x = y
    
    check_dim(y, x)

    mat = [[float(i == j) for i in range(x)] for j in range(y)]

    return(mat)

def gen_zero_mat(y, x = -1):
    if x == -1:
        x = y
        
    check_dim(y, x)

    mat = [[0.0 for i in range(x)] for j in range(y)]

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
        matrix = [[float(int(matrix[j][i])) for i in \
        range(len(matrix[0]))] for j in range(len(matrix))]
        
    return(matrix)

def file_to_mat(fpath):
    try:
        with open(fpath) as f:
            l = [[float(num) for num in line.split(' ')] \
                 for line in f]
    except IOError as er:
        print(er)
        return
        
    return(matr(l))

class SOLE(matr):
    def __init__(self, matrix):
        matr.__init__(self, matrix)
        try:
            self.LU = self.LU()
        except ZeroDivisionError:
            print("Can't find LU decomposition")
        self.det = self.det()

A = matr(CONST_TEST_MAT)
A.mat_to_file('matr2.txt')
file_to_mat('mat2r.txt')
