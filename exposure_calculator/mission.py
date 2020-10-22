# Imports
import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt
from pprint import pprint
from . import common_tools as ct


class Mission:
    """Struct like class to handle mission data neatly"""

    def __init__(self, detector=None):

        # Physical propertoes
        self.D = None  # aperture diamater
        self.obscuration = None  # percentage of the aperture obsecured by the obscuration
        self.pixel_scale = None  # milliarcseconds
        self.jitter = 0  # 1 sigma pointing jitter in milliarcseconds

        # Mission bandpass
        self.bands = np.array([])  # low end of filter
        self.llo = np.array([])  # low end of filter
        self.lhi = np.array([])  # high end of filter
        self.lc = np.array([])  # filter central wavelength

        # a range that spans the sensitivty of the mission
        self.wavelengths = np.array([])
        self.throughputs = np.array([])  # throughput at every wavelength

        # pretty plot params
        # jitter convolved psf at the center wavelength of band
        self.psf = np.array([])
        self.mag = np.array([])  # 5 sigma snr 1 hour integration per pixel

    # Missons collecting area

    @property
    def collecting_area(self):
        if self.D is None or self.obscuration is None:
            print("Mission does not have a defined diameter or obscuration")
            return None
        else:
            return np.pi*(self.D/2)**2 * (1 - self.obscuration)

    def show_vars(self):
        pprint(vars(self))

    def plot_instrument(self, color='b', psf=True):
        for i in range(len(self.llo)):
            plt.plot([self.llo[i], self.lhi[i]], [
                     self.mag[i], self.mag[i]], color=color)
            if psf:
                plt.text(self.lc[i], self.mag[i], str(self.psf[i]))


# Define some standard missions

### superbit ###
################
superbit = Mission()

# Physical attributes
superbit.D = 0.5  # meters
superbit.obscuration = 0.384  # fraction of aperture that is obscured
superbit.jitter = 2.355*50  # milliarcs from perfromance paper approx
superbit.pixel_scale = 206  # milliarcs

# bands
superbit.bands = np.array(['u', 'b', 'g', 'r', 'i', 'lum', 'ha'])
superbit.llo = np.array([300, 350, 475, 550, 700, 375])
superbit.lhi = np.array([400, 500, 550, 800, 1100, 700])
superbit.lc = (superbit.llo+superbit.lhi)/2.0
superbit.wavelenghts = np.arange(310, 1101)
superbit.throughputs = np.zeros((7, 791))

for i in range(len(superbit.bands)):
    superbit.throughputs[i, :] = np.loadtxt(
        "../data/superbit_throughput/" + superbit.bands[i] + ".csv", delimiter=',')


# Predicted PSF performance
superbit.psf = np.array(np.sqrt(ct.get_psf(
    superbit.lc*u.nm, superbit.D*u.m)**2+superbit.jitter**2), dtype=np.int32)

superbit.psf_nojitter = np.array(
    ct.get_psf(superbit.lc*u.nm, superbit.D*u.m), dtype=np.int32)

# special atribute includes the background mag per pixel as measured by ajay
superbit.mag_b = np.array([24.489, 25.642, 24.721, 24.624, 22.273, 25.022])


### gigabit ###
################
gigabit = Mission()

# physical attributes
gigabit.D = 1.3
gigabit.jitter = 2.355*20  # milliarcs from perfromance paper approx
gigabit.obscuration = 0.3
gigabit.pixel_scale = 100  # milliarcs

# bands
gigabit.llo = np.array([300, 350, 475, 550, 700, 375])
gigabit.lhi = np.array([400, 500, 550, 800, 1100, 700])
gigabit.lc = (gigabit.llo+gigabit.lhi)/2.0
gigabit.bands = np.array(['u', 'b', 'g', 'r', 'i', 'lum', 'ha'])
gigabit.wavelenghts = np.arange(310, 1101)
gigabit.throughputs = np.zeros((7, 791))

for i in range(len(gigabit.bands)):
    gigabit.throughputs[i, :] = np.loadtxt(
        "../data/superbit_throughput/" + gigabit.bands[i] + ".csv", delimiter=',')

gigabit.psf = np.array(np.sqrt(
    ct.get_psf(gigabit.lc*u.nm, gigabit.D*u.m)**2+gigabit.jitter**2), dtype=np.int32)
gigabit.psf_nojitter = np.array(
    ct.get_psf(gigabit.lc*u.nm, gigabit.D*u.m), dtype=np.int32)
gigabit.aperture = np.pi*(1.3/2)**2
gigabit.aperture_effective = gigabit.aperture*(1-gigabit.obscuration)

gigabit.mag_b = np.array([24.489, 25.642, 24.721, 24.624, 22.273, 25.022])


### NGRST ###
################
ngrst = Mission()
ngrst.D = 1.3
ngrst.jitter = 2.355*14
ngrst.llo = np.array([480, 760, 927, 1131, 1380, 1683, 927])
ngrst.lhi = np.array([760, 977, 1192, 1454, 1774, 2000, 2000])
ngrst.lc = (ngrst.lhi+ngrst.llo)/2.0

ngrst.psf = np.array(
    np.sqrt(ct.get_psf(ngrst.lc*u.nm, ngrst.D*u.m)**2+ngrst.jitter**2), dtype=np.int32)
ngrst.psf_nojitter = np.array(
    ct.get_psf(ngrst.lc*u.nm, ngrst.D*u.m), dtype=np.int32)
ngrst.mag = [28.5, 28.02, 27.95, 27.87, 27.81,
             27.32, 28.33]  # 5 sigma snr 1 hour integration
