from array import array
import math
import reprlib
import numbers
import functools
import operator
import itertools

class Vector:
    typecode = 'd'
    
    def __init__(self, components):
        self._components = array(self.typecode, components)
        
    @classmethod
    def frombytes(cls, octets):
        typecode = chr(octets[0])
        memv = memoryview(octets[1:]).cast(typecode)
        return cls(memv)

    def angle(self, n):
        r = math.sqrt(sum(x * x for x in self[n:]))
        a = math.atan2(r, self[n - 1])
        if (n == len(self) - 1) and (self[-1] < 0):
            return math.pi * 2 - a
        else:
            return a

    def angles(self):
        return(self.angle(n) for n in range(1, len(self)))

    def __iter__(self):
        return iter(self._components)

    def __repr__(self):
        components = reprlib.repr(self._components)
        components = components[components.find('['):-1]
        return 'Vector({})'.format(components)

    def __str__(self):
        return str(tuple(self))

    def __bytes__(self):
        return (bytes([ord(self.typecode)]) +
        bytes(self._components))    

    def __bool__(self):
        return bool(abs(self))
    
    def __format__(self, fmt_spec=''):
        if fmt_spec.endswith('h'):
            fmt_spec = fmt_spec[:-1]
            coords = itertools.chain([abs(self)], self.angles())
            outer_fmt = '<{}>'
        else:
            coords = self
            outer_fmt = '({})'
        components = (format(c, fmt_spec) for c in coords)
        return outer_fmt.format(', '.join(components))
    def __len__(self):
        return len(self._components)

    def __getitem__(self, index):
        cls = type(self)
        if isinstance(index, slice):
            return cls(self._components[index])
        elif isinstance(index, numbers.Integral):
            return self._components[index]
        else:
            print('{cls.__name__} indices must be integers')
            raise TypeError(cls = cls)

    shortcuts = 'xyzt'
     
    def __getattr__(self, name):
        cls = type(self)
        if len(name) == 1:
            pos = cls.shortcuts.find(name)
            if 0 <= pos < len(self.shortcuts):
                return self._components[pos]
        msg = '{.__name__!r} objects has no attribute {!r}'
        raise AttributeError(msg.format(cls, name))

    def __setattr__(self, name, value):
        cls = type(self)
        if len(name) == 1:
            if name in self.shortcuts:
                error = 'readonly attribute {attr_name!r}'
            elif name.islower():
                error = 'can`t set attributes from "a" to "z" in {cls_name!r}'
            else:
                error = ''
            if error:
                msg = error.format(cls_name = cls.__name__, attr_name = name)
                raise AttributeError(msg)
        super().__setattr__(name, value)

    def __hash__(self):
        hashes = map(hash, self._components)
        return functools.reduce(operator.xor, hashes, 0)
    
    def __eq__(self, other):
        if isinstance(other, Vector):
            return len(self) == len(other) and all(x == y for x,y in zip(self, other))
        else:
            return NotImplemented
        
    def __abs__(self):
        return math.sqrt(sum(x * x for x in self))

    def __neg__(self):
        return Vector(-x for x in self)

    def __pos__(self):
        return Vector(self)

    def __add__(self, other):
        try:
            pairs = itertools.zip_longest(self, other, fillvalue = 0.0)
            return Vector(a + b for a, b in pairs)
        except TypeError:
            return NotImplemented
        
    def __radd__(self, other):
        return self + other

    def __mul__(self, scalar):
        if isinstance(scalar, numbers.Real):
            return Vector(x * scalar for x in self)
        else:
            return NotImplemented

    def __rmul__(self, scalar):
        return self * scalar
        
    def __matmul__(self, other):
        try:
            return sum(a * b for a, b in zip(self, other))
        except TypeError:
            return NotImplemented
        
    def __rmatmul__(self, other):
        return self @ other

    def __sub__(self, other):
        return -(other - self)

    def __rsub__(self, other):
        return (other + -self)
    
