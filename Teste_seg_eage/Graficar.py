import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rc('font', size=16)
mpl.rc('figure', figsize=(8, 6))

def plot_shotrecord(rec, origin, domain_size, t0, tn, factor= 10., cmap = "gray", colorbar=True):
    """
    Plot a shot record (receiver values over time).
    Parameters
    ----------
    rec :
        Receiver data with shape (time, points).
    model : Model
        object that holds the velocity model.
    t0 : int
        Start of time dimension to plot.
    tn : int
        End of time dimension to plot.
    """
    scale = np.max(rec) / factor
    extent = [1e-3*origin[0], 1e-3*origin[0] + 1e-3*domain_size[0],
              1e-3*tn, t0]

    plt.figure()
    plot = plt.imshow(rec, vmin=-scale, vmax=scale, cmap=cmap, extent=extent, aspect = 'auto')
    plt.xlabel('X position (km)')
    plt.ylabel('Time (s)')
    

    # Create aligned colorbar on the right
    if colorbar:
        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(plot, cax=cax)
    
    plt.tight_layout()

def plot_image(data, extent, vmin=None, vmax=None, colorbar=True, cmap="gray"):
    """
    Plot image data, such as RTM images or FWI gradients.
    Parameters
    ----------
    data : ndarray
        Image data to plot.
    cmap : str
        Choice of colormap. Defaults to gray scale for images as a
        seismic convention.
    """
    plt.figure()
    plot = plt.imshow(data,
                      vmin=vmin or 0.9 * np.min(data),
                      vmax=vmax or 1.1 * np.max(data),
                      cmap=cmap, extent= extent)

    # Create aligned colorbar on the right
    if colorbar:
        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(plot, cax=cax)

    plt.tight_layout()


def plot_seismic_traces(traces, t0, tn):
    """
    Plot seismic traces
    Parameters
    ----------
    traces : list
        List of traces to plot
    t0 :
        initial time
    tn :
        final time
    """
    plt.figure()
    for trace in traces:
        x = np.linspace(t0, tn, trace.shape[0])
        plt.plot(trace, x)
        plt.ylim(tn, t0)

    plt.tight_layout()


# #p_devito_600 = np.load("p_devito_ref_t600_Test5.npy")
# p_devito_1200 = np.load("p_devito_ref_t1200_Test5.npy")
# rec_devito = np.load("rec_devito_ref_t1200_Test5.npy")

# #plot_shotrecord(rec_devito, (3000, 0), (7000, 4200), 0, 1200, factor=100)
# #plot_image(p_devito_600.T)
# plot_image(p_devito_1200.T, [3000, 10000, 4200, 0])
# traces = [p_devito_1200[5760, :]]
# plot_seismic_traces(traces, 0, 1200)

# plt.show()