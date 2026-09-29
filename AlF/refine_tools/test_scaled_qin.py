"""Energy-zero and bound-character regressions for the scaled Qin experiment."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scale_qin_reference import scale
from bound_states import summarize
from fit_scaled_qin import remove_A_reference


class ScaledQinTests(unittest.TestCase):
    def test_old_reference_removed_without_changing_model(self):
        old='poten 2\ntype EMO\nvalues\nV0 44000\nend\n'
        ref='(old reference)\nabinitio poten 2\ntype grid\nvalues\n1.5 45000\nend\n'
        cleaned=remove_A_reference(old+ref+'FITTING\nenergies\nend\n')
        self.assertTrue(cleaned.startswith(old))
        self.assertNotIn('abinitio',cleaned)
        self.assertNotIn('old reference',cleaned)
        self.assertIn('FITTING',cleaned)

    def test_scale_preserves_radii_and_every_relative_energy_ratio(self):
        p=np.array([[1.5,10.],[1.7,-10.],[2.5,2.],[9.,-1.]])
        q,s=scale(p,100.,120.,0.)
        np.testing.assert_array_equal(q[:,0],p[:,0])
        np.testing.assert_allclose(q[:,1]-100,s*(p[:,1]+10))
        self.assertEqual(s,2.)
        self.assertEqual(q[-1,1],118.) # finite endpoint must not become infinity silently

    def test_endpoint_convention_maps_endpoint_to_requested_limit(self):
        p=np.array([[1.5,10.],[1.7,-10.],[2.5,2.],[9.,-1.]])
        q,s=scale(p,100.,120.,-1.)
        self.assertAlmostEqual(q[-1,1],120.)
        self.assertAlmostEqual(s,20/9)

    def test_invalid_depth_or_radii_rejected(self):
        for p,te,limit,src in [([[1,0],[2,-1],[3,0]],2,1,0),
                              ([[1,0],[1,-1],[3,0]],0,2,0),
                              ([[1,0],[2,-1],[3,0]],0,2,-2)]:
            with self.assertRaises(ValueError):scale(p,te,limit,src)

    def test_grid_summary_uses_absolute_limit_minus_zpe(self):
        source='''poten 2
name "A1Pi"
type grid
values
1.5 45000
2.5 57000
9.0 55530
end
FITTING
energies
0 + 1 0 1 0 0 0 0 1
end
'''
        base=dict(J=1,parity='+',ef='f',state='A1Pi',mean_r=2.,tail_probability=0.,flag='b')
        rows=[dict(base,id=1,v=0,energy=55150),dict(base,id=2,v=1,energy=55200)]
        with tempfile.TemporaryDirectory() as tmp:
            s=summarize(source,rows,400,Path(tmp),dissociation=55564.)
            self.assertEqual(s['dissociation_above_X_ground_cm'],55164.)
            self.assertEqual(s['categories']['below_limit_localized'],1)
            self.assertEqual(s['categories']['above_limit_localized_candidate'],1)
            with self.assertRaises(ValueError):summarize(source,rows,400,Path(tmp))


if __name__=='__main__':unittest.main()
