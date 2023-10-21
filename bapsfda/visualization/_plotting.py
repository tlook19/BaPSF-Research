__all__ = ["plot_3d_surf"]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator


def set_ax_labels(ax, **kwargs):
    ax.set_xlabel(kwargs.get("xlabel", "x"))
    ax.set_ylabel(kwargs.get("ylabel", "y"))
    try:
        ax.set_zlabel(kwargs.get("zlabel", "z"))
    except AttributeError:
        pass


def plot_3d_surf(X, Y, Z, **kwargs):
    """
    Plots a 3D surface colormap plot.

    Parameters
    ----------
    X, Y, Z : Array-like
        Input data. See matplotlib.Axes.contour for supported data shapes.

    file_name : str, optional
        Filename to save plot to. If not provided, plot is not saved.

    title : str, optional
        Title for plot.

    contour : bool, optional
        Whether to plot contours on the x-y, x-z, and y-z planes.

    xlabel, ylabel, zlabel : str, optional
        Axes labels. Defaults to "x", "y", "z".

    """
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    X, Y = np.meshgrid(X, Y)

    # Plot the surface.
    surf = ax.plot_surface(
        X, Y, Z, cmap=cm.inferno, linewidth=0, antialiased=False, alpha=0.7
    )
    if kwargs.get("contour"):
        ax.contour(X, Y, Z, zdir="z", offset=0, cmap="coolwarm")
        ax.contour(X, Y, Z, zdir="x", offset=X[0, 0], cmap="coolwarm")
        ax.contour(X, Y, Z, zdir="y", offset=Y[-1, 0], cmap="coolwarm")
    # Customize the z axis.
    # ax.set_zlim(0, 5)
    # ax.zaxis.set_major_locator(LinearLocator(11))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter("{x:.01f}")
    set_ax_labels(ax, **kwargs)
    if kwargs.get("title") is not None:
        ax.set_title(kwargs.get("title"))
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=20, pad=0.1)
    if kwargs.get("file_name") is not None:
        plt.savefig(kwargs.get("file_name"))
    plt.show()
    plt.close()


# import matplotlib
# import scaleogram as scg


# def debug_plotter(xarr, yarr, debug_func=None, debug_vline=None, debug_hline=None):
#     plt.plot(xarr, yarr)
#     if debug_func is not None:
#         if hasattr(debug_func, "__iter__"):
#             for f in debug_func:
#                 if callable(f):
#                     plt.plot(xarr, f(xarr))
#         elif callable(debug_func):
#             plt.plot(xarr, debug_func(xarr))
#     if debug_vline is not None:
#         if hasattr(debug_vline, "__iter__"):
#             for v in debug_vline:
#                 plt.axvline(v)
#         else:
#             plt.axvline(debug_vline)
#     if debug_hline is not None:
#         if hasattr(debug_hline, "__iter__"):
#             for h in debug_hline:
#                 plt.axhline(h)
#         else:
#             plt.axhline(debug_hline)
#     plt.show()


# font = {"family": "Helvetica", "weight": "normal", "size": 18}

# matplotlib.rc("font", **font)


# def plot_stack(
#     xcoords,
#     ydata,
#     labels,
#     xlabel,
#     ylabel,
#     save_plot=False,
#     save_name="plot.pdf",
#     save_dir="/data/tlook/analysis/plots/",
#     dark=False,
#     plt_style="default",
# ):
#     """
#     Plots a stack of lines on the same plot.

#     Args:
#         xcoords (_type_): _description_
#         ydata (_type_): _description_
#         labels (_type_): _description_
#         xlabel (_type_): _description_
#         ylabel (_type_): _description_
#         save_plot (bool, optional): _description_. Defaults to False.
#         save_name (str, optional): _description_. Defaults to "plot.pdf".
#         save_dir (str, optional): _description_. Defaults to "/data/tlook/analysis/plots/".
#         dark (bool, optional): _description_. Defaults to False.
#         plt_style (str, optional): _description_. Defaults to "default".
#     """
#     if dark:
#         plt.style.use("dark_background")
#     else:
#         plt.style.use("default")

#     color_bwr = plt.cm.bwr(np.linspace(0, 1, 8))
#     for i, y in enumerate(ydata):
#         plt_args = (xcoords, y)
#         plt_kwargs = {"label": f"{labels[i]}", "color": color_bwr[i], "marker": "."}
#         if plt_style == "loglog":
#             plt.loglog(*plt_args, **plt_kwargs)
#         elif plt_style == "semilogx":
#             plt.semilogx(*plt_args, **plt_kwargs)
#         elif plt_style == "semilogy":
#             plt.semilogy(*plt_args, **plt_kwargs)
#         else:
#             if not plt_style == "default":
#                 print(f"Warning: unknown plt_style {plt_style}, using default")
#             plt.plot(*plt_args, **plt_kwargs)
#     plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
#     plt.xlabel(xlabel=xlabel)
#     plt.ylabel(ylabel=ylabel)
#     if save_plot:
#         plt.savefig(save_dir + save_name, bbox_inches="tight")
#     plt.show()
