import unittest
import math as m
import StationaryMeasurments2


class TestStationaryMeasurements2(unittest.TestCase):
    def test_degtorad(self):
        result1 = StationaryMeasurments2.degtorad(180)
        self.assertEqual(result1, m.pi)
        result2=StationaryMeasurments2.degtorad(30)
        self.assertEqual(result2,m.pi/6)
        result3 = StationaryMeasurments2.degtorad(270)
        self.assertEqual(result3,3*m.pi/2)

    def test_lbstokg(self):
        result1 = StationaryMeasurments2.lbstokg(15)
        self.assertEqual(result1, 6.8038856)

    def test_total_weight(self):
        result1=StationaryMeasurments2.total_weight(0)
        self.assertEqual(result1,)
