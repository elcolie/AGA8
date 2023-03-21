import unittest
import typing as typ
import detail

class TestStringMethods(unittest.TestCase):

    def test_simple(self):
        """Test naive case."""
        detail.SetupDetail()
        x: typ.List[float] = [0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005,
        0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009, 0.004, 0.005,
        0.002, 0.0001, 0.0025, 0.007, 0.001]
        x.insert(0, 0.0)

        mm: float = 0.0
        mm = detail.MolarMassDetail(x)

        T: float = 400
        P: float = 50000
        D: float = 1e10
        detail.DensityDetail(T, P, x, D)

        self.assertEqual('foo'.upper(), 'FOO')


if __name__ == '__main__':
    unittest.main()
