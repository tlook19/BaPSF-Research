import numpy as np
import matplotlib.pyplot as plt

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
