# Imports
import numpy as np
from astropy import units as u
from astropy import constants as const
from matplotlib import pyplot as plt
from pprint import pprint
import matplotlib.patches as mpatches


class Detector:
    def __init__(self, name, res, depth, n_dark, n_read, x_dim=4000, y_dim=4000, pix_size=5, sensitivity=35):

        self.name = name  # sensor name
        self.res = 2**res*u.adu  # adc res 2^res
        self.depth = depth * u.electron  # well depth in e-
        self.n_dark = n_dark * u.electron/u.second  # dark current e-/s
        self.n_read = n_read * u.electron  # read noise e-
        self.gain = self.depth/(self.res)  # e-/ADU
        self.x_dim = x_dim * u.pix  # pixels
        self.y_dim = y_dim * u.pix  # pixels
        self.pix_size = pix_size * u.micrometer  # microns
        self.sensitivity = sensitivity * u.microVolt/u.electron

    def add_qe_from_csv(self, path):
        from scipy.interpolate import interp1d
        temp = np.loadtxt(path, delimiter=',')
        wave_sens = temp[:, 0]
        qe_sens = temp[:, 1]
        # snippet taken form ajays code to interpolate between points and IR
        # used same code as ajay to ensure consistency
        f_sens = interp1d(wave_sens, qe_sens,
                          kind='quadratic')  # interpolated QE
        wave_sens_intp = np.arange(310, 950, 1)
        qe_sens_intp = f_sens(wave_sens_intp)
        i_overall = qe_sens_intp
        y2, x2 = i_overall[len(i_overall)-1], 949
        y1, x1 = 0, 1100

        slope = (y2-y1) / (x2-x1)
        b = y2 - (slope * x2)

        wave_extra = np.arange(950, 1101, 1)
        ir_extra = np.zeros(len(wave_extra))
        rest_extra = np.zeros(len(wave_extra))

        for i in range(len(wave_extra)):
            ir_extra[i] = (slope * wave_extra[i]) + b

        wave_tot_ir = np.append(wave_sens_intp, wave_extra)
        i_overall_intp = np.append(i_overall, ir_extra)
        self.wave = wave_tot_ir
        self.qe = i_overall_intp

    def show_vars(self):
        pprint(vars(self))

    def plot_qe(self):
        plt.figure(figsize=(14, 8))
        plt.plot(self.wave, self.qe)
        plt.ylabel("QE")
        plt.xlabel("Wavelength nm")
        plt.title(self.name + " QE Curve")
