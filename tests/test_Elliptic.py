import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from Elliptic import Eleptic_curves

@pytest.fixture
def ppass(): # For some reason, the first function is not considered a test, because of this I added this useless function
    pass

def test_generator():
    ellip_one = Eleptic_curves(-2, 1, 7)
    ellip_two = Eleptic_curves(-2, 1, 50873)
    ellip_three = Eleptic_curves(2, 3, 1125899906843)
    assert ellip_one.find_the_generator() == 12
    assert ellip_two.find_the_generator() == 50560
    assert ellip_three.find_the_generator() == 1125900059286

def test_point():
    ellip_one = Eleptic_curves(-2, 1, 7)
    ellip_two = Eleptic_curves(-2, 1, 51803)
    ellip_three = Eleptic_curves(2, 3, 1125899906843)

    assert ellip_one.Sum_of_two_points([-3, -1], [-3, -1]) == [3, 6]
    with pytest.raises(ValueError):
        ellip_one.Sum_of_two_points([-3, 0], [-3, 1])

    assert ellip_two.Exponentiation(13037, [291, 1035]) == [0, 0]
    assert ellip_two.Exponentiation(52148, [409, 1522]) == [0, 0]

    assert ellip_three._Eleptic_curves__find_point_order([939250958713, 348791869945]) == 2143
    assert ellip_three._Eleptic_curves__find_point_order([442629937818, 1061215920018]) == 46901