# This file houses useful flux density related converters
from astropy import units as u


def convert_MAGAB_to_flux_photons(magab, central_wavelength):
    """Converts magAB values to photons/nm/s/m^2 at a given wavelength"""
    f_nu = magab * u.ABmag
    g_nu = f_nu.to(
        u.photon / u.m ** 2 / u.s / u.nm,
        equivalencies=u.spectral_density(central_wavelength * u.nm),
    )
    return g_nu.value  # photons/nm/s/m^2


def convert_MAGAB_to_Jansky(magab, central_wavelength):
    """Converts magAB values to janskay at a given wavelength"""
    f_nu = magab * u.ABmag
    f_nu = f_nu.to(
        u.Jansky, equivalencies=u.spectral_density(central_wavelength * u.nm)
    )
    return f_nu.value  # photons/nm/s/m^2


def convert_Jansky_to_flux_photons(jy, central_wavelength):
    """Converts Jansky values to photons/nm/s/m^2 at a given wavelength"""
    f_nu = jy * u.Jansky
    g_nu = f_nu.to(
        u.photon / u.m ** 2 / u.s / u.nm,
        equivalencies=u.spectral_density(central_wavelength * u.nm),
    )
    return g_nu.value  # photons/nm/s/m^2


def convert_Jansky_to_MAGAB(jy, central_wavelength):
    """Converts Jansky values to magab at a given wavelength"""
    f_nu = jy * u.Jansky
    g_nu = f_nu.to(
        u.ABmag, equivalencies=u.spectral_density(central_wavelength * u.nm)
    )
    return g_nu.value  # photons/nm/s/m^2


def convert_flux_photons_to_MAGAB(photons_per_second, central_wavelength):
    """Converts photons/nm/s/m^2 values to magAB at a given wavelength"""
    g_nu = photons_per_second * u.photon / u.m ** 2 / u.s / u.nm
    f_nu = g_nu.to(
        u.ABmag, equivalencies=u.spectral_density(central_wavelength * u.nm)
    )
    return f_nu.value  # mag AB


def convert_flux_photons_to_Jansky(photons_per_second, central_wavelength):
    """Converts photons/nm/s/m^2 values to magAB at a given wavelength"""
    g_nu = photons_per_second * u.photon / u.m ** 2 / u.s / u.nm
    f_nu = g_nu.to(
        u.Jansky, equivalencies=u.spectral_density(central_wavelength * u.nm)
    )
    return f_nu.value  # mag AB
