# Imports
import numpy as np
from astropy import units as u
from astropy import constants as const
from . import converters as convert
from . import common_tools as ct


def get_mag_from_SNR_slow(
    snr,
    mission,
    band,
    n_dark,
    n_read,
    exptime,
    mag_background,
    mags=np.arange(23, 28),
    n=3,
    per_pixel=False,
    psf_fwhm=None
):

    pts = []
    for mag in mags:
        pts.append(get_SNR_from_mag(mission, band, mag,
                                    mag_background, n_dark, n_read, exptime, per_pixel, psf_fwhm))

    # fit mag vs snr relation ship and find the mag at which snr = given
    z = np.polyfit(np.log(pts), mags, n)
    p = np.poly1d(z)
    return p(np.log(snr))


def get_SNR_from_mag(
    mission, band, mag_source, mag_background, n_dark, n_read, exptime, per_pixel=False, psf_fwhm=None
):

    # the index that represnets the band name
    band = np.argwhere(mission.bands == band)[0][0]

    # throughput and wavelngths of the mission at the specified band
    e_nu = mission.throughputs[band]
    lam = mission.wavelenghts  # nm

    # Intialize accumulators for background and source electron counts
    n_background = 0
    n_source = 0

    # Handle provided magnitude SEDs
    if isinstance(mag_source, list) or isinstance(mag_source, np.ndarray):
        print("You have provided a source SED")
    else:
        print("Source SED will be assumed as flat with valuev = ", mag_source)
        mag_source = np.ones(len(lam))*mag_source

    if isinstance(mag_background, list) or isinstance(mag_background, np.ndarray):
        print("You have provided a backgroudn SED")
    else:
        print("Background SED will be assumed as flat with valuev = ", mag_background)
        mag_background = np.ones(len(lam))*mag_background

    for i in range(len(lam)):
        if e_nu[i] > 1e-2:  # Cut off at throughput of 0.1%

            # TODO: Package the bakground and source number of electrons calc into a function so it can be used in conjuction with t from SNR

            # if psf is not provided assume diffraction limited
            if psf_fwhm is None:
                psf_fwhm = ct.get_psf(lam[i] * u.nm, mission.D * u.m)

            # handle calculating SNR per pixel vs for the entire source
            if per_pixel:
                pix_illum = ct.pixel_illumnation_fraction(
                    psf_fwhm, mission.pixel_scale)
                approx_src_size = 1
            else:
                pix_illum = 1
                approx_src_size = ((3*(psf_fwhm/2.35482)) /
                                   mission.pixel_scale)**2  # 3*sigma pixels wide

            # Calculate number of electrons per wavelenght contribution from background
            background_flux_density = convert.magab_to_photons(
                mag_background[i], lam[i]
            )  # assumes flat background SED

            background_electrons = (
                background_flux_density *
                e_nu[i] * mission.collecting_area*approx_src_size
            )

            n_background = background_electrons + n_background

            # Calculate number of electrons per wavelength contribution from sources
            source_flux_density = convert.magab_to_photons(
                mag_source[i], lam[i]
            )

            source_electrons = (
                source_flux_density
                * e_nu[i]
                * mission.collecting_area
                * pix_illum
            )
            n_source = source_electrons + n_source

    return SNR(exptime, n_source, n_background, n_dark*approx_src_size, n_read)


# CCD equation in all three arangments for snr, t, source_count


def SNR(t, source, background, dark, read):
    """use ccd equation to return snr"""
    return (
        source
        * t
        / np.sqrt(source * t + background * t + dark * t + read ** 2)
    )


def t_from_SNR(snr, source, background, dark, read):
    """use ccd equation to return exposure time for given snr"""
    root = (snr ** 4) * (
        -background - source - dark
    ) ** 2 + 4 * source ** 2 * read ** 2 * snr ** 2
    num = (
        np.sqrt(root)
        + background * snr ** 2
        + source * snr ** 2
        + dark * snr ** 2
    )
    denum = 2 * source ** 2
    return num / denum


def source_count_from_snr(snr, t, background, dark, read):
    """use ccd equation to return source photon count for a given snr"""
    root = 4 * background * t + 4 * dark * t + 4 * read ** 2 + snr ** 2
    num = snr ** 2 * t + snr * t * np.sqrt(root)
    denum = 2 * t ** 2
    return num / denum
