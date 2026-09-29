"""Regression checks for scientific interpretation and lossless preparation."""
import tempfile
from pathlib import Path
import unittest
from bound_states import prepare, read_states
from coupled import fields
from barrier_spectrum import spectrum_ceiling


class BoundTests(unittest.TestCase):
    def test_below_limit_barrier_does_not_truncate_bound_spectrum(self):
        self.assertGreater(spectrum_ceiling((2.1,49750),55564),55564)
        self.assertGreater(spectrum_ceiling((2.55,61500),55564),61500)

    def test_native_columns_do_not_confuse_v_and_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'sample.states'
            p.write_text('142 55617.402696 36 1 + f A1Pi 41 1 0 1 b 1.89829 7.2e-27\n'
                         '118 55151.383139 36 1 + f A1Pi 17 1 0 1 u 5.20976 0.00054\n')
            a,b=read_states(p)
            self.assertEqual(a['v'],41)
            self.assertEqual(a['flag'],'b')
            self.assertGreater(a['energy'],55163.877082)
            self.assertEqual(b['flag'],'u')
            self.assertLess(b['energy'],55163.877082)
            self.assertGreater(b['mean_r'],4)

    def test_preparation_preserves_curves_and_removes_fitting(self):
        source='''jrot 0 - 98
grid
npoints 601
range .7 5
end
contraction
vib
vmax 60 60
end
poten 1
type EMO
values
V0 0
RE 1.65 fit
end
FITTING
itmax 5
energies (J parity NN energy)
0 + 1 0 1 0 0 0 0 1
end
dipole 1 1
type grid
values
1.6 1.0
end
'''
        result=prepare(source)
        self.assertNotIn('FITTING',result)
        self.assertEqual(result.count('INTENSITY'),1)
        self.assertIn('states_only',result)
        self.assertIn('dipole 1 1',result)
        self.assertEqual(fields(source)['POTEN',1,1].values,fields(result)['POTEN',1,1].values)

    def test_unknown_states_layout_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'sample.states';p.write_text('1 0 8 0 0\n')
            with self.assertRaises(ValueError):read_states(p)


if __name__=='__main__':unittest.main()
