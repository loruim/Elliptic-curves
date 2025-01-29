import copy
import random
from sympy import factorint

def Exponentiation(value, degree, n):
    bin_degree = bin(degree)[2:][::-1]
    t = len(bin_degree)

    b = 1
    if degree == 0:
        return b
    
    A = copy.copy(value)
    if bin_degree[0] == 1:
        b = copy.copy(value)

    for i in range(1, t):
        A = A**2 % n
        if (int(bin_degree[i]) == 1):
            b = (A * b) % n

    return b

def Fermat_test(value : int, iterations : int) -> bool:
    if (value == 0):
        raise ValueError(f"You put zero")

    for i in range(1, iterations):
        a = random.randint(2, value - 1)
        r = Exponentiation(a, value-1, value)
        if (r != 1):
            return False
    
    return True

def Euclidean_algorithm(value1, value2):
    y2 = 0
    y1 = 1
    while (value2 > 0):
        q = value1 // value2
        r = value1 - q * value2
        y = y2 - q * y1

        value1 = copy.copy(value2)
        value2 = copy.copy(r)
        y2 = copy.copy(y1)
        y1 = copy.copy(y)

    nod = value1
    y = y2

    if nod == 1:
        return y2
    
    raise ValueError(f"NOD > 1") 

def factorize_number(value):
    factors = factorint(value)
    factors_set = list(dict.keys(factors))
    factors_set.insert(0, 1)
    factors_set.append(value)
    return factors_set