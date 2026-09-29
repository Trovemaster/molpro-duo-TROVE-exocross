"""Scientific data-integrity checks for reference PEC preparation."""
import unittest
import numpy as np
from add_abinitio_reference import convert, block_text
from prepare_qin_reference import extract_pair


class ReferenceTests(unittest.TestCase):
    def test_each_state_uses_its_own_radius(self):
        text = ('R/Angstrom X1Sigma+ R/Angstrom A1Pi\n'
                '1.72\t-100\t1.74\t-20\n'
                '1.74\t-99\t1.76\t-19\n'
                '1.76\t-98\t\t\n')
        rows, cols = extract_pair(text)
        self.assertEqual(cols, [3, 4])
        self.assertEqual(rows, [('1.74', '-20', 2), ('1.76', '-19', 3)])

    def test_missing_internal_pair_does_not_shift_other_columns(self):
        text = ('R/Angstrom X1Sigma+ R/Angstrom A1Pi R/Angstrom 21Pi\n'
                '1.7\t-100\t\t\t1.9\t200\n'
                '1.8\t-99\t2.0\t-19\t2.1\t210\n')
        self.assertEqual(extract_pair(text)[0], [('2.0', '-19', 3)])

    def test_half_missing_pair_is_an_error(self):
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            extract_pair('R/Angstrom A1Pi\n1.7\t\n')

    def test_vertical_alignment_preserves_shape_in_physical_units(self):
        raw = np.array([[2, -.1], [3, -.2], [4, -.15]])
        r, e, minimum, shift = convert(raw, 43000, 'bohr', 'hartree')
        np.testing.assert_allclose(r, [1.058354421088, 1.587531631632, 2.116708842176], atol=1e-12, rtol=0)
        np.testing.assert_allclose(e, [64947.46313632, 43000, 53973.73156816], atol=1e-9, rtol=0)
        self.assertAlmostEqual(minimum, -43894.92627264)
        self.assertAlmostEqual(shift, 86894.92627264)

    def test_invalid_numeric_data_are_not_silently_repaired(self):
        for points in [[[1, 0], [1, 1], [2, 3]], [[1, 0], [2, float('nan')], [3, 4]]]:
            with self.assertRaises(ValueError):
                convert(points, 43000, 'angstrom', 'cm-1')

    def test_explicit_weights_do_not_enable_live_ps1997(self):
        block = block_text(2, 'A1Pi', 1, 1, [1, 2, 3], [5, 0, 3], 100, .001, 18000, [1, .5, .1])
        self.assertNotIn('\nWeighting ', block)
        self.assertEqual(block.count('\nvalues\n'), 1)
        self.assertIn('abinitio poten 2\n', block)
        self.assertIn('fit_factor 100\n', block)


if __name__ == '__main__':
    unittest.main()
