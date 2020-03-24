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
        self.assertEqual(result1,
                         (StationaryMeasurments2.OEW+StationaryMeasurments2.weight_payload+
                          StationaryMeasurments2.weight_fuel)*StationaryMeasurments2.g0)
        result2=StationaryMeasurments2.total_weight(StationaryMeasurments2.weight_fuel)
        self.assertEqual(round(result2,6),
                         round((StationaryMeasurments2.OEW+StationaryMeasurments2.weight_payload)*
                               StationaryMeasurments2.g0,6))

    def test_fttom(self):
        result2=StationaryMeasurments2.fttom(20)
        self.assertEqual(result2,6.096)

    def test_density(self):
        result1=StationaryMeasurments2.density(1000,2,281.65-273.15)
        self.assertEqual(round(result1,3),1.112) #comparing to ISA values
        result2 = StationaryMeasurments2.density(3000, 2, 268.65 - 273.15)
        self.assertEqual(round(result2,3), 0.909) #comparing to ISA values