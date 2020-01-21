from math import pi, sin, cos, tan, atan2, log, copysign
from traceback import format_exc

import tkinter as tk
import tkinter.ttk as ttk
import tkinter.messagebox as msgbox
from tkinter.filedialog import askopenfilename

import numpy as np

from matplotlib import use
use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as Navigation


def gz_poly(xs, zs, xv, zv, density, con=13.3464e0):
    """
    Compute the vertical component of gravitational acceleration due to a
    polygon in a two-dimensional space (X, Z) at multiple stations.

    Note: X-axis is positive-down, Z axis is positive right.

    Parameters
    ----------
    xs : list
        X-coordinates of the stations.
    zs : list
        Z-coordinates of the stations.
    xv : list
        X-coordinates of the polygon vertexes (clockwise).
    zv : list
        Z-coordinates of the polygon vertexes (clockwise).
    density : float
        Polygon density.
    con : float
        Conversion factor (dependent on units). Example below:

    con = 13.3464e+0  # g/cm**3, mgals, km
    con = 13.3464e-3  # g/cm**3, mgals, m
    con = 4.0680e+0   # g/cm**3, mgals, kft
    con = 13.3464e-8  # kg/m**3, mm/s**2, m
    con = 13.3464e-11 # kg/m**3, m/s**2, m

    If the chosen units for gravity, length and density differ from the
    associated SI units by factors C1, C2 and C3 respectively, then the
    parameter `con` should be set to (C2*C3/C1) * 13.3464e-11.

    Returns
    -------
    grav_z : list
        Gravitational acceleration at each station due to the polygon

    Original code in Fortran by Won and Bevis (1987)
    """

    assert len(xs) == len(zs)
    assert len(xv) == len(zv)

    nstn = len(xs)
    nvert = len(xv)

    gravz = []

    for i in range(nstn):
        grav = 0.0
        xst = xs[i]
        zst = zs[i]
        for ic in range(nvert):
            gz = 0.0
            x1 = xv[ic] - xst
            z1 = zv[ic] - zst

            if ic == nvert-1:  # last loop
                x2 = xv[0] - xst
                z2 = zv[0] - zst
            else:
                x2 = xv[ic+1] - xst
                z2 = zv[ic+1] - zst

            if x1 == 0 and z1 == 0:
                continue
            else:
                th1 = atan2(z1, x1)

            if x2 == 0 and z2 == 0:
                continue
            else:
                th2 = atan2(z2, x2)

            if copysign(1, z1) != copysign(1, z2):
                test = x1*z2 - x2*z1
                if test > 0:
                    if z1 >= 0:
                        th2 += 2*pi
                elif test < 0:
                    if z2 >= 0:
                        th1 += 2*pi
                else:  # station is on polygon side
                    continue

            t12 = th1 - th2
            z21 = z2 - z1
            x21 = x2 - x1
            xz1 = x1*x1 + z1*z1
            xz2 = x2*x2 + z2*z2

            if x21 == 0:
                gz = 0.5 * x1 * log(xz2/xz1)
            else:
                gz = x21 * (x1*z2 - x2*z1) * (t12 + 0.5 * (z21/x21) * log(xz2/xz1)) / (x21*x21 + z21*z21)

            grav += gz
        gravz.append(con*density*grav)
    return gravz


def anomaly(polygons, stations, con=13.3464e0):
    """
    Return the anomaly from `polygons` at `stations`.
    """
    out = [0] * len(stations)

    xs, ys = zip(*stations)

    for points, density in polygons:
        xv, yv = zip(*points)
        for i, a in enumerate(gz_poly(xs, ys, xv, yv, density, con)):
            out[i] += a
    return out


def anomalyfromcircle(r, x, z, x1, z1, density):
    """Return the anomaly from a circle at (x1, z1) with radius r and density at a station at (x, z)"""
    f = 6.67408e-11  # grav constant
    return 2 * pi * f * density * r*r * (z1-z) / ( (x-x1)**2 + (z-z1)**2 )


def loadpoly(file):
    """Return a list [([(x1, y1), ...], density), ...] from file."""
    out = []
    with open(file) as f:
        while True:
            try:
                density = float(f.readline())
            except ValueError:
                break
            l = []
            try:
                while True:
                    line = f.readline().replace(",", " ").split()
                    l.append((float(line[0]), float(line[1])))
            except (ValueError, IndexError):
                out.append((l, density))
    if not out:
        raise ValueError(f"Could not read '{file}'")
    return out


def loadstations(file):
    """Return a list [(x1, y1), ...] from file."""
    out = []
    with open(file) as f:
        for line in f:
            try:
                line = line.replace(",", " ").split()
                out.append((float(line[0]), float(line[1])))
            except (ValueError, IndexError):
                continue
    return out


fig = Figure()
ax = fig.add_subplot(111)


class _MplCanvas(ttk.Frame):
    """
    A matplotlib canvas for tkinter.

    Can be updated, but remember to call .draw()
    """

    def __init__(self, master):
        ttk.Frame.__init__(self, master)

        self.canvas = FigureCanvasTkAgg(fig, master=self)

        self.canvas.get_tk_widget().pack(fill="both", expand=1)

        self.draw = self.canvas.draw
        self.draw()

        Navigation(self.canvas, self)


class Window(tk.Tk):
    window_title = "Gravity Modelling"
    units_con = {
        "g/cm**3, mgals, km": 13.3464e+0,
        "g/cm**3, mgals, m": 13.3464e-3,
        "g/cm**3, mgals, kft": 4.0680e+0,
        "kg/m**3, mm/s**2, m": 13.3464e-8,
        "kg/m**3, m/s**2, m": 13.3464e-11}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title(self.window_title)

        self.polygons = []
        self.stations = []

        self.canvas = _MplCanvas(self)
        self.canvas.grid(row=0, column=0, sticky="nsew")

        options = ttk.Frame()
        options.grid(row=1, column=0, sticky="ew", padx=3, pady=3)

        self.unit = tk.StringVar()

        ttk.Label(options, text="Units:").grid(row=0, column=0, sticky="nw")
        ttk.OptionMenu(options, self.unit, next(iter(self.units_con)), *self.units_con.keys()
            ).grid(row=0, column=1, sticky="nw")

        ttk.Button(options, text="Load Polygon File", command=self.loadpoly).grid(row=0, column=3, sticky="ne")
        ttk.Button(options, text="Load Station File", command=self.loadstations).grid(row=0, column=4, sticky="ne")

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        options.grid_columnconfigure(2, weight=1)

        ax.set_aspect("equal")

        for widget in options.winfo_children():
            widget.grid_configure(pady=2, padx=2)


    def loadpoly(self):
        try:
            self.polygons = loadpoly(askopenfilename())
        except Exception:
            msgbox.showerror(self.window_title, "An error occurred.\n\n" + format_exc())
            raise

        self.update()

    def loadstations(self):
        try:
            self.stations = loadstations(askopenfilename())
        except Exception:
            msgbox.showerror(self.window_title, "An error occurred.\n\n" + format_exc())
            raise

        self.update()

    def update(self):
        ax.clear()
        patches = []
        densities = []

        ax.plot(*zip(*self.stations), color="black", marker=".", linestyle="")

        if self.polygons:
            for polygon in self.polygons:
                patches.append(Polygon(polygon[0], False))
                densities.append(polygon[1])
            p = PatchCollection(patches)
            p.set_array(np.array(densities))
            ax.add_collection(p)
            ax.autoscale_view()
            if self.stations:
                an = anomaly(self.polygons, self.stations, con=self.units_con[self.unit.get()])
                for station, a in zip(self.stations, an):
                    ax.text(*station, "{:.2f}".format(a),
                            verticalalignment="bottom", horizontalalignment="center")
                    a *= 5 / max(an)
                    ax.annotate("", xy=(station[0], station[1]+a), xytext=station,
                        arrowprops=dict(arrowstyle="->"))



        ax.invert_yaxis()
        self.canvas.draw()


if __name__ == '__main__':
    #print(anomaly(loadpoly("circle.csv"), loadstations("originstation.csv"), con=13.3464e-11))
    #print(anomalyfromcircle(30, 0, 0, 0, 20, 10))

    root = Window()
    root.mainloop()
