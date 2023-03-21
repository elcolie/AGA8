import typing as typ
import unittest

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
        print(f"Molar Mass: {mm}")

        T: float = 400
        P: float = 50000
        D: float = 1e10
        D_of_T_and_P, ierr_no, herr = detail.DensityDetail(T, P, x, D)

        _P, _Z, _dPdD, _d2PdD2, _d2PdTD, _dPdT, _U, _H, _S, _Cv, _Cp, _W, _G, _JT, _Kappa = detail.PropertiesDetail(T,
                                                                                                                    D_of_T_and_P,
                                                                                                                    x)
        print(f"D_of_T_and_P: {D_of_T_and_P}")
        print(f"Speed of sound [m/s] : {_W}")
        self.assertEqual(_W, 712.63936840579)  # [m/s]


if __name__ == '__main__':
    unittest.main()
