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

        _P, _Z, _dPdD, _d2PdD2, _d2PdTD, _dPdT, _U, _H, _S, _Cv, _Cp, _W, _G, _JT, _Kappa = detail.PropertiesDetail(
            T, D_of_T_and_P, x)

        mm_reference: float = 20.54333051
        D_reference: float = 12.80792403648801
        P_reference: float = 50000.0
        Z_reference: float = 1.173801364147326
        dPdD_reference: float = 6971.387690924090
        d2PdD2_reference: float = 1118.803636639520
        dPdT_reference: float = 235.6641493068212
        U_reference: float = -2739.134175817231
        H_reference: float = 1164.699096269404
        s_reference: float = -38.54882684677111
        Cv_reference: float = 39.12076154430332
        Cp_reference: float = 58.54617672380667
        W_reference: float = 712.6393684057903
        G_reference: float = 16584.22983497785
        JT_reference: float = 7.432969304794577E-05
        Kappa_reference: float = 2.672509225184606

        print("Inputs------")
        print(f"Temperature [K]: 400.00 != {T}")
        print(f"Pressure [kPa]: 50000.00 != {P}")
        print("Outputs------")
        if abs(mm_reference - mm) > 1.0e-8:
            print(f"Molar mass [g/mol]:                 {mm_reference} != {mm}")
            raise Exception("Molar mass error")
        if abs(D_reference - D_of_T_and_P) > 1.0e-8:
            print(f"Molar density [mol/l]:              {D_reference} != {D_of_T_and_P}")
            raise Exception("Molar density error")
        if abs(P_reference - _P) > 1.0e-8:
            print(f"Pressure [kPa]:                     {P_reference} != {_P}")
            raise Exception("Pressure error")
        if abs(Z_reference - _Z) > 1.0e-8:
            print(f"Compressibility factor:             {Z_reference} != {_Z}")
            raise Exception("Compressibility error")
        if abs(dPdD_reference - _dPdD) > 1.0e-8:
            print(f"d(P)/d(rho) [kPa/(mol/l)]:          {dPdD_reference} != {_dPdD}")
            raise Exception("d(P)/d(rho) error")
        if abs(d2PdD2_reference - _d2PdD2) > 1.0e-8:
            print(f"d^2(P)/d(rho)^2 [kPa/(mol/l)^2]:    {d2PdD2_reference} != {_d2PdD2}")
            raise Exception("d^2(P)/d(rho)^2 error")
        if abs(dPdT_reference - _dPdT) > 1.0e-8:
            print(f"d(P)/d(T) [kPa/K]:                  {dPdT_reference} != {_dPdT}")
            raise Exception("d(P)/d(T) error")
        if abs(s_reference - _S > 1.0e-8):
            print(f"Energy [J/mol]:                     {U_reference} != {_U}")
            raise Exception("Energy error")
        if abs(H_reference - _H) > 1.0e-8:
            print(f"Enthalpy [J/mol]:                   {H_reference} != {_H}")
            raise Exception("Enthalpy error")
        if abs(s_reference - _S) > 1.0e-8:
            print(f"Entropy [J/mol-K]:                  {s_reference} != {_S}")
            raise Exception("Entropy error")
        if abs(Cv_reference - _Cv) > 1.0e-8:
            print(f"Isochoric heat capacity [J/mol-K]:  {Cv_reference} != {_Cv}")
            raise Exception("Isochoric heat capacity error")
        if abs(Cp_reference - _Cp) > 1.0e-8:
            print(f"Isobaric heat capacity [J/mol-K]:   {Cp_reference} != {_Cp}")
            raise Exception("Isobaric heat capacity error")
        if abs(W_reference - _W) > 1.0e-8:
            print(f"Speed of sound [m/s]:               {W_reference} != {_W}")
            raise Exception("Speed of sound error")
        if abs(G_reference - _G) > 1.0e-8:
            print(f"Gibbs energy [J/mol]:               {G_reference} != {_G}")
            raise Exception("Gibbs energy error")
        if abs(JT_reference - _JT) > 1.0e-8:
            print(f"Joule-Thomson coefficient [K/kPa]:  {JT_reference} != {_JT}")
            raise Exception("Joule-Thomson coefficient error")
        if abs(Kappa_reference - _Kappa) > 1.0e-8:
            print(f"Isentropic exponent:                {Kappa_reference} != {_Kappa}")
            raise Exception("Isentropic exponent error")


if __name__ == '__main__':
    unittest.main()
