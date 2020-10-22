# This file houses useful flux density related converters
from astropy import units as u


def MAGAB_to_photons(magab, central_wavelength):
    """Converts magAB flux density to photons/nm/s/m^2 at a given wavelength

    Args:
        magab (float): value of flux density in AB magnitude.
        central_wavelength (float): value of photon wavelength in nanometers

    Returns:
        float: value of flux density in photons/nm/s/m^2
    """
    f_nu = magab * u.ABmag
    g_nu = f_nu.to(
        u.photon / u.m ** 2 / u.s / u.nm,
        equivalencies=u.spectral_density(central_wavelength * u.nm),
    )
    return g_nu.value  # photons/nm/s/m^2


def Jansky_to_photons(jy, central_wavelength):
    """Converts Janskys to photons/nm/s/m^2 at a given wavelength

    Args:
        jy (float): value of flux density in jansky
        central_wavelength (float): value photon wavelength in nanometers

    Returns:
        float: value of flux density in photons/nm/s/m^2
    """
    f_nu = jy * u.Jansky
    g_nu = f_nu.to(
        u.photon / u.m ** 2 / u.s / u.nm,
        equivalencies=u.spectral_density(central_wavelength * u.nm),
    )
    return g_nu.value  # photons/nm/s/m^2


def photons_to_MAGAB(photons_per_second, central_wavelength):
    """Converts photon flux density from photons/nm/s/m^2 to AB magnitude at a given wavelength

    Args:
        photons_per_second (float): value of photon flux density in photons/nm/s/m^2
        central_wavelength (float): value of photon wavelength in nanometers

    Returns:
        float: value of flux density in AB magnitude
    """
    g_nu = photons_per_second * u.photon / u.m ** 2 / u.s / u.nm
    f_nu = g_nu.to(
        u.ABmag, equivalencies=u.spectral_density(central_wavelength * u.nm)
    )
    return f_nu.value  # mag AB


def photons_to_Jansky(photons_per_second, central_wavelength):
    """[summary]

    Args:
        photons_per_second (float): value of photon flux density in photons/nm/s/m^2
        central_wavelength (float): value of photon wavelength in nanometers

    Returns:
        float: value of flux density in Jansky
    """
    """Converts photons/nm/s/m^2 values to magAB at a given wavelength"""
    g_nu = photons_per_second * u.photon / u.m ** 2 / u.s / u.nm
    f_nu = g_nu.to(
        u.Jansky, equivalencies=u.spectral_density(central_wavelength * u.nm)
    )
    return f_nu.value  # mag AB
