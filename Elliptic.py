from general_functions import Fermat_test, Euclidean_algorithm, factorize_number
from sympy.ntheory.residue_ntheory import is_quad_residue
from sympy.ntheory import sqrt_mod
import random
import math
import copy


class Eleptic_curves:
    def __init__(self, a : int, b : int, p : int):
        if Fermat_test(p, 20) != True:
            raise ValueError(f"The number {p} is not prime")
        if (-4 * (a) ** 3 - 27 * (b) ** 2) % p == 0:
            raise ValueError(f"The equation y**2 = x**3 + {a} * x + {b} is not a curve equation")
        self.a, self.b, self.p = a, b, p # y**2 = x^3 + a * x + b

        self.Lower_limit = math.ceil(self.p + 1 - 2 * (self.p) ** 0.5)
        self.Upper_limit = math.floor(self.p + 1 + 2 * (self.p) ** 0.5)

    def Sum_of_two_points(self, one_point : list, two_point : list) -> list:  
        if (one_point == [0, 0]): return two_point
        elif (two_point == [0, 0]): return one_point

        self.__check_point(one_point), self.__check_point(two_point)

        if (one_point[0] == two_point[0] and one_point[1] == -two_point[1]):
            return [0, 0]
        elif (one_point == two_point):
            return self.__identical(one_point)
        elif (one_point[0] != two_point[0]):
            return self.__different(one_point, two_point)
    
    def Exponentiation(self, degree : int, value : list):
        bin_degree = bin(degree)[2:][::-1]
        t = len(bin_degree)
        value[0] %= self.p
        value[1] %= self.p

        b = [0, 0]
        if degree == 0:
            return b

        A = copy.copy(value)
        if int(bin_degree[0]) == 1:
            b = copy.copy(value)

        for i in range(1, t):
            A = self.Sum_of_two_points(A, A)
            if (int(bin_degree[i]) == 1):
                b = self.Sum_of_two_points(A, b)

        if b == None:
            b =[0, 0]
        return b
    
    def find_the_generator(self):
        order_one = self.__find_point_order(self.__point_generation())
        order_two = self.__find_point_order(self.__point_generation())
        nok = math.lcm(order_one, order_two)
        loop_flag = False
        while not loop_flag:
            count = 0
            for i in range(self.Lower_limit, self.Upper_limit + 1):
                if (i % nok == 0):
                    count += 1
                    if (count > 1):
                        break
                    N = i

            if (count == 1):
                self.N = N
                print(f"E({self.a}, {self.b})(F_{self.p}), orfer: {N}")
                return N
            nok = math.lcm(nok, self.__find_point_order(self.__point_generation()))

    def finding_points(self):
        order_divisors = factorize_number(self.N)
        order_divisors.remove(1)
        loop_flag = False
        while not loop_flag:
            point = self.__point_generation()
            Q = self.Exponentiation((self.N // order_divisors[0]), point)
            if(Q == [0, 0]):
                    continue
            if self.Exponentiation(order_divisors[0], Q) == [0, 0]:
                print(f"Point {Q}, order: {order_divisors[0]}")
                order_divisors.pop(0)
            
            if (len(order_divisors) == 0):
                loop_flag = True

    def __identical(self, point : list):
        x_one = point[0]
        y_one = (point[1]) % self.p

        first_part = (3 * (x_one**2) + self.a) % self.p
        second_part = Euclidean_algorithm(self.p, 2 * y_one)

        x_three = (((first_part * second_part) ** 2) - (2 * x_one)) % self.p
        y_three = ((first_part * second_part) * (x_one - x_three) - y_one) % self.p
        return [x_three, y_three]

    def __different(self, point_one : list, point_two : list):
        x_one, x_two = point_one[0], point_two[0]
        y_one, y_two = point_one[1], point_two[1]

        first_part = y_two - y_one
        second_part = (x_two - x_one) % self.p

        x_three = ((first_part * Euclidean_algorithm(self.p, second_part)) ** 2 - x_one - x_two) % self.p
        y_three = ((first_part * Euclidean_algorithm(self.p, second_part)) * (x_one - x_three) - y_one) % self.p

        return [x_three, y_three]

    def __check_point(self, point : list):
        if (type(point) == int or len(point) == 0):
            raise ValueError(f"Complete the coordinates with zeros according to x and y")
        elif (len(point) > 2):
            raise ValueError(f"Point have only x, y coordinates")
        elif (point[0] % self.p or point[1] % self.p) not in range(self.p):
            raise ValueError(f"The point value does not lie within the field F({self.p})")
        elif (point[1]**2 % self.p != (point[0]**3 + self.a * point[0] + self.b) % self.p):
            raise ValueError(f"The point {point} does not belong to the eleptic line")

    def __find_point_order(self, point : list):
        point[0] %= self.p
        point[1] %= self.p

        mP_Set = []
        m = math.ceil(self.p ** 0.25)
        for j in range(m+1):
            mP_Set.append(self.Exponentiation(j, point))
        
        giant_step = self.Exponentiation(self.Lower_limit, point)
        new_order = self.Lower_limit
        if giant_step in mP_Set:
            l = next(j for j in range(len(mP_Set)) if giant_step == mP_Set[j])
            order = new_order - l
            order = self.__order_check(order, point)
            return order
        
        for i in range(self.Lower_limit, self.Upper_limit):
            giant_step = self.Sum_of_two_points(giant_step, mP_Set[len(mP_Set) - 1])
            new_order += (len(mP_Set) - 1)
            if giant_step in mP_Set:
                l = next(j for j in range(len(mP_Set)) if giant_step == mP_Set[j])
                Order = new_order - l
                Order = self.__order_check(Order, point)
                return Order
            elif new_order > self.Upper_limit:
                raise RuntimeError(f"Not work")

    def __order_check(self, order, point):
        factorize_order = factorize_number(order)
        for i in factorize_order:
            check_point = self.Exponentiation(i, point)
            if (check_point == [0, 0]):
                return i
            
    def __point_generation(self):
        loop_flag = False
        while not loop_flag:
            x_coord = random.randint(0, self.p-1)
            square_y_coord = ((x_coord) ** 3 + self.a * x_coord + self.b) % self.p
            if (is_quad_residue(square_y_coord, self.p)):
                point, loop_flag = [x_coord, sqrt_mod(square_y_coord, self.p)], True
        return point