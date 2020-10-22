import unittest
from exposure_calculator import converters as convert


class TestConversions(unittest.TestCase):

    def test_magab_to_photons(self):
        photons = convert.magab_to_photons(25, 500)
        self.assertAlmostEqual(0.010959, photons, 4,
                               msg="Mag25 should be ~0.011 photons at 500nm")

    def test_jansky_to_photons(self):
        photons = convert.jansky_to_photons(3.6307805477010036e-07, 500)
        self.assertAlmostEqual(0.010959, photons, 4,
                               msg=" should be ~0.011 photons at 500nm")

    def test_photons_to_magab(self):
        magab = convert.photons_to_magab(0.01095907669405222, 500)
        self.assertAlmostEqual(25, magab, 1,
                               msg="Mag25 should be ~0.011 photons at 500nm")

    def test_photons_to_jansky(self):
        jansky = convert.photons_to_jansky(0.01095907669405222, 500)
        self.assertAlmostEqual(3.6307805477010036e-07, jansky, 1,
                               msg=" should be ~0.011 photons at 500nm")


if __name__ == '__main__':
    unittest.main()
