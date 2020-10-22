import random
from IPython.display import HTML
from matplotlib import pyplot as plt
import numpy as np
from astropy import units as u
from astropy import constants as const
plt.rcParams.update({'font.size': 14})


def hide_toggle(for_next=False, ttext="Code"):
    """Useful tool for jupyter notebooks in order to hide code"""
    this_cell = """$('div.cell.code_cell.rendered.selected')"""
    next_cell = this_cell + '.next()'

    toggle_text = ttext  # text shown on toggle link
    target_cell = this_cell  # target cell to control with toggle
    # bit of JS to permanently hide code in current cell (only when toggling next cell)
    js_hide_current = ''

    if for_next:
        target_cell = next_cell
        toggle_text += ' next cell'
        js_hide_current = this_cell + '.find("div.input").hide();'

    js_f_name = 'code_toggle_{}'.format(str(random.randint(1, 2**64)))

    html = """
        <script>
            function {f_name}() {{
                {cell_selector}.find('div.input').toggle();
            }}

            {js_hide_current}
        </script>

        <a href="javascript:{f_name}()">{toggle_text}</a>
    """.format(
        f_name=js_f_name,
        cell_selector=target_cell,
        js_hide_current=js_hide_current,
        toggle_text=toggle_text
    )

    return HTML(html)


def makeGaussian(size, sig=3, center=None):
    """ Make a square gaussian kernel."""
    x = np.arange(0, size, 1, float)
    y = x[:, np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    g = np.exp(-.5 * ((x-x0)**2 + (y-y0)**2) / sig**2)
    return g/np.sum(g)


def get_psf(lam, D):
    """return psf fwhm in milli arcseconds given diameter and wavelength with astrop units"""
    ans = (1.22 * lam / D).decompose() * u.rad
    return ans.to(u.milliarcsecond).value


def pixel_illumnation_fraction(psf_fwhm, pixel_scale):
    """given a psf fwhm and a plate scale assume gaussian psf and return illumnation fraction of the pixel"""
    pixel = np.zeros((pixel_scale * 5, pixel_scale * 5))
    pixel[
        2 * pixel_scale: 3 * pixel_scale + 1,
        2 * pixel_scale: 3 * pixel_scale + 1,
    ] = 1

    return np.sum(
        makeGaussian(pixel_scale * 5, sig=psf_fwhm / 2.35482) * pixel
    )


def is_fits(name):
    """Checks if the image is a fits image"""
    return ".fits" in name


def indices_array_generic(m, n):
    """ Generates an m by n numpy array of index positions """
    r0 = np.arange(m)  # Or r0,r1 = np.ogrid[:m,:n], out[:,:,0] = r0
    r1 = np.arange(n)
    out = np.empty((m, n, 2), dtype=int)
    out[:, :, 0] = r0[:, None]
    out[:, :, 1] = r1
    return out


def plotimg(img, sig=3, title="", cb=True, report=True):
    mini = np.mean(img)-sig*np.std(img)
    maxi = np.mean(img)+sig*np.std(img)
    plt.figure(figsize=(14, 8))
    plt.imshow(img, vmin=mini, vmax=maxi, cmap="viridis")
    if cb:
        plt.colorbar()
    plt.title(title)
    if report:
        print(np.mean(img), np.median(img), np.std(img))


def plothist(img, title="", bins=10000, sig=1, nfig=True, label="", report=True):
    m = np.mean(img)
    s = np.std(img)
    hist, edge = np.histogram(img, bins=bins)
    if nfig:
        plt.figure(figsize=(14, 8))
    plt.plot(edge[1:], hist, label=label)
    plt.xlim([m-(sig*s), m+(sig*s)])
    plt.title(title)
    plt.legend()
    temp = sorted(hist, key=lambda x: abs(x-np.max(hist)/2))
    print("Peak =", np.max(hist), "FWHM = ", np.abs(
        edge[np.where(hist == temp[0])] - edge[np.where(hist == temp[1])])[0])
    if report:
        print(m, np.median(img), s)
