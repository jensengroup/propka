"""Test the hybrid36 module."""
import unittest
import propka.hybrid36 as hybrid36


class Hybrid36Test(unittest.TestCase):
    """Test class for hybrid36."""

    def test_decode(self):
        """Test decoding functions."""
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
        for key, value in test_values.items():
            self.assertEqual(hybrid36.decode(key), value)

    def test_errors(self):
        """Test values that should raise errors."""
        test_values = [
            "99X99",
            "X9-99",
            "XYZa",
            "",
            "-",
            "!NotOk",
        ]
        for value in test_values:
            with self.assertRaises(ValueError) as err:
                hybrid36.decode(value)
            self.assertTrue(value in str(err.exception))
