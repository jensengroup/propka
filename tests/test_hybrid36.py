import unittest

import propka.hybrid36 as hybrid36

class Hybrid36Test(unittest.TestCase):
    def testDecode(self):
        test_values = {
            "99999": 99999,
            "A0000": 100000,
            "0": 0,
            "9": 9,
            "A": 10,
            "  ZZZZY": 43770014,
            "ZZZZZ": 43770015, # ZZZZZ - A0000 + 100000
            "a0000": 43770016,
            "zzzzz": 87440031,
            "zzzzy": 87440030,
            "99": 99,
            "A0": 100,
            "ZZ": 1035,
            "zz": 1971,
            "-99999": -99999,
            "-A0000": -100000,
            "-0": 0,
            "-9": -9,
            "-A": -10,
            "-ZZZZY": -43770014,
            "-ZZZZZ": -43770015, # ZZZZZ - A0000 + 100000
            "-a0000": -43770016,
            "-zzzzz": -87440031,
            "-zzzzy": -87440030,
            "-99": -99,
            "-A0": -100,
            "-ZZ": -1035,
            "-zz": -1971,
            "PROPKA": 954495146,
            "A001Z": 100071,
            "B0000": 1779616,
        }

        for k, v in test_values.items():
            self.assertEqual(hybrid36.decode(k), v)

    def testErrors(self):
        test_values = [
            "99X99",
            "X9-99",
            "XYZa",
            "",
            "-",
            "!NotOk",
        ]

        for v in test_values:
            with self.assertRaises(ValueError) as e:
                hybrid36.decode(v)
            self.assertTrue(v in str(e.exception))