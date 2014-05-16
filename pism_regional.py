#!/usr/bin/env python
from netCDF4 import Dataset as NC

import numpy as np
import sys
import dbg

def permute(variable, output_order = ('time', 'z', 'zb', 'y', 'x')):
    '''
    Permute dimensions of a NetCDF variable to match the output
    storage order.

    Parameters
    ----------
    variable : a netcdf variable
               e.g. thk = nc.variables['thk']
    output_order: dimension tuple (optional)
                  default ordering is ('time', 'z', 'zb', 'y', 'x')

    Returns
    -------
    var_perm : array_like
    '''
    input_dimensions = variable.dimensions

    # filter out irrelevant dimensions
    dimensions = filter(lambda(x): x in input_dimensions,
                        output_order)

    # create the mapping
    mapping = map(lambda(x): dimensions.index(x),
                  input_dimensions)

    if mapping:
        return np.transpose(variable[:], mapping)
    else:
        return variable[:]              # so that it does not break processing "mapping"

def find_coordinate_variables(input_file):
    "Find names of coordinate variables in input_file."
    # set defaults:
    x_name = "x"
    y_name = "y"
    for name in input_file.variables:
        variable = input_file.variables[name]
        if getattr(variable, "standard_name", "") == "projection_x_coordinate":
            x_name = name
        if getattr(variable, "standard_name", "") == "projection_y_coordinate":
            y_name = name

    return x_name, y_name

def load_data(input_file):
    """Loads data from an input file.

    An input file has to contain variables 'x', 'y', 'usurf', 'thk'.
    """

    nc = NC(input_file)

    xdim, ydim = find_coordinate_variables(nc)

    dimension_order = ('time', 'z', 'zb', ydim, xdim)

    x = np.array(nc.variables[xdim][:], dtype=np.double)
    y = np.array(nc.variables[ydim][:], dtype=np.double)
    try:
        z = np.array(np.squeeze(permute(nc.variables['usurf'], dimension_order)), dtype=np.double, order='C')
    except:
        z = np.array(np.squeeze(permute(nc.variables['usrf'], dimension_order)), dtype=np.double, order='C')
    thk = np.array(np.squeeze(permute(nc.variables['thk'], dimension_order)), dtype=np.double, order='C')

    nc.close()

    return (x, y, z, thk)

def save_mask(input_file, output_file, result, cutout_command, history):
    """ Saves the computed drainage basin mask to a file.
    """

    print "Saving the mask to %s..." % output_file,

    nc_in = NC(input_file)

    xdim, ydim = find_coordinate_variables(nc_in)

    x_orig = nc_in.variables[xdim]
    y_orig = nc_in.variables[ydim]

    nc_out = NC(output_file, 'w', format='NETCDF3_64BIT')

    nc_out.createDimension('x', x_orig.size)
    nc_out.createDimension('y', y_orig.size)

    x = nc_out.createVariable("x", 'f8', ('x',))
    y = nc_out.createVariable("y", 'f8', ('y',))
    mask = nc_out.createVariable("ftt_mask", 'i4', ('y', 'x'))

    mask.long_name = "Drainage basin area for regional modeling"

    # copy attributes
    for var, old_var in zip([x,y], [x_orig, y_orig]):
        for attr in old_var.ncattrs():
            value = old_var.getncattr(attr)
            if isinstance(value, (str, unicode)):
                value = str(value.encode('ASCII', 'ignore'))
            var.setncattr(attr, value)

    # copy coordinate data
    x[:] = x_orig[:]
    y[:] = y_orig[:]

    mask[:] = result != 2

    nc_out.cutout_command = cutout_command
    nc_out.history = history

    nc_out.close()
    nc_in.close()

    print "done."

def initialize_mask(thk, x, y, terminus):

    mask = dbg.initialize_mask(thk)

    if terminus is not None:
        x_min, x_max, y_min, y_max = terminus
        for j in range(y.size):
            for i in range(x.size):
                inside = (x[i] >= x_min and
                          x[i] <= x_max and
                          y[j] >= y_min and
                          y[j] <= y_max)

                if inside:
                    mask[j,i] = 2
                else:
                    if mask[j,i] > 0:
                        mask[j,i] = 1

    return mask

def compute_bbox(input_file, mask, x, y, border):
    """Compute the bounding box around a drainage basin and return the NCO
    command what would cut it out of the bigger dataset.
    """
    x0 = x.size - 1
    x1 = 0
    y0 = y.size - 1
    y1 = 0

    nc = NC(input_file)
    xdim, ydim = find_coordinate_variables(nc)
    nc.close()

    for j in range(y.size):
        for i in range(x.size):
            if mask[j,i] == 2:
                if x[i] < x[x0]:
                    x0 = i

                if x[i] > x[x1]:
                    x1 = i

                if y[j] < y[y0]:
                    y0 = j

                if y[j] > y[y1]:
                    y1 = j


    x0 = np.maximum(x0 - border, 0)
    x1 = np.minimum(x1 + border, x.size - 1)

    y0 = np.maximum(y0 - border, 0)
    y1 = np.minimum(y1 + border, y.size - 1)

    return "ncks -d %s,%d,%d -d %s,%d,%d %s output.nc" % (xdim, x0, x1, ydim, y0, y1, input_file)


class App:
    """An application class containing methods of the drainage basin tool.
    """
    def __init__(self, master):
        self.input_file = None
        self.master = master
        self.fill = None
        self.terminus = None
        self.nc = None
        self.mask_computed = False
        self.cutout_command = ""
        self.Ncontours = 30

        self.create_widgets(master)

        self.load_data()


    def load_data(self):
        self.input_file = tkFileDialog.askopenfilename(parent=root,
                                                       filetypes = ["NetCDF .nc"],
                                                       title='Choose an input file')

        if len(self.input_file) == 0:
            print "No input file selected. Exiting..."
            sys.exit(0)

        self.x, self.y, self.z, self.thk = load_data(self.input_file)

        self.mask = initialize_mask(self.thk, self.x, self.y, None)
        print "Mask initialization: done"

        plt.figure(1)
        self.plot_mask(0, cmaps.binary)

        try:
            # block=False is not supported in some earlier versions of matplotlib
            plt.show(block=False)

            f = plt.get_current_fig_manager().window
            w, x0, y0 = f.winfo_width(), f.winfo_x(), f.winfo_y()

            self.master.geometry("+%d+%d" % (x0 + w + 10, y0))
            self.master.lift()
        except:
            plt.show()

    def save_results(self):

        self.output_file = self.get_output()

        if self.output_file is None:
            print "No output file selected; cannot proceed."
            return

        save_mask(self.input_file, self.output_file, self.mask, self.cutout_command,
                  self.compute_command())

    def plot_mask(self, threshold, colormap):
        """Plots mask > threshold using the given colormap.
        Only 2 colors in the colormap matter, though...
        """
        plt.pcolormesh(self.x, self.y, self.mask > threshold, cmap=colormap)
        plt.contour(self.x, self.y, self.z, self.Ncontours, colors='black')
        plt.axis('tight')
        plt.axes().set_aspect('equal')
        plt.xticks([])
        plt.yticks([])

    def get_output(self):
        """Asks the user for the name of the output file.
        """
        output = tkFileDialog.asksaveasfilename(parent=root,
                                                filetypes = ["NetCDF .nc"],
                                                title="Save the mask in...")
        if len(output) > 0:
            return output
        else:
            return None

    def create_widgets(self, master):
        """Creates all the widgets."""
        frame = Frame(master)
        frame.grid()

        # 1
        label = Label(master, text="1.")
        label.grid(padx=2, pady=2, row=1, column=1, sticky=E+W)

        button = Button(master, text="Select terminus rectangle", command=self.get_terminus)
        button.grid(padx=2, pady=2, row=1, column=2, columnspan=2, sticky=E+W)

        # 2
        label = Label(master, text="2.")
        label.grid(padx=2, pady=2, row=2, column=1, sticky=E+W)

        label = Label(master, text="Set border width (cells):")
        label.grid(padx=2, pady=2, row=2, column=2, sticky=W)

        self.entry = Entry(master, width=10)
        self.entry.grid(padx=2, pady=2, row=2, column=3)
        self.entry.insert(0, "5")

        # 3
        label = Label(master, text="3.")
        label.grid(padx=2, pady=2, row=3, column=1, sticky=E+W)

        button = Button(master, text="Compute the drainage basin mask", command=self.compute_mask)
        button.grid(padx=2, pady=2, row=3, column=2, columnspan=2, sticky=E+W)

        # 4
        label = Label(master, text="4.")
        label.grid(padx=2, pady=2, row=4, column=1, sticky=E+W)

        button = Button(master, text="Save the drainage basin mask", command=self.save_results)
        button.grid(padx=2, pady=2, row=4, column=2, columnspan=2, sticky=E+W)

        # 5
        label = Label(master, text="5.")
        label.grid(padx=2, pady=2, row=5, column=1, sticky=E+W)

        button = Button(master, text="Quit", command=master.quit)
        button.grid(padx=2, pady=5, row=5, column=2, columnspan=2, sticky=E+W)

        master.update()

        w, h = master.winfo_width(), master.winfo_height()
        sw, sh = master.winfo_screenwidth(), master.winfo_screenheight()

        master.geometry("+%d+%d" % ((sw - w) / 2, (sh - h) / 2))
        master.wm_resizable(False, False)

    def get_terminus(self):
        """Gets (and plots) the terminus rectangle.
        """
        from matplotlib.widgets import Cursor

        if self.mask_computed == True:
            self.mask = dbg.initialize_mask(self.thk)

            plt.clf()
            self.plot_mask(0, cmaps.binary)
            plt.draw()

        plt.setp(plt.gca(),autoscale_on=False)

        cursor = Cursor(plt.axes(), useblit=True, color='black', linewidth=2 )

        # remove the rectangle
        if self.fill is not None and self.mask_computed == False:
            for p in self.fill:
                p.remove()
            self.fill = None

        x_min, x_max, y_min, y_max = plt.axis()

        x0, y0 = plt.ginput(timeout=-1)[0]

        l1 = plt.plot([x0, x0], [y_min, y_max], color='black', lw=2)
        l2 = plt.plot([x_min, x_max], [y0, y0], color='black', lw=2)

        x1, y1 = plt.ginput(timeout=-1)[0]

        l3 = plt.plot([x1, x1], [y_min, y_max], color='black', lw=2)
        l4 = plt.plot([x_min, x_max], [y1, y1], color='black', lw=2)

        dx = x1 - x0
        dy = y1 - y0

        xs = [x0, x1, x1, x0]
        ys = [y0, y0, y1, y1]
        self.fill = plt.fill(xs, ys, fill=False, edgecolor='black', lw = 2, hatch='/')

        # remove guides
        for line in [l1, l2, l3, l4]:
            line[0].remove()

        x_min = np.minimum(x0, x1)
        x_max = np.maximum(x0, x1)

        y_min = np.minimum(y0, y1)
        y_max = np.maximum(y0, y1)

        self.terminus = (x_min, x_max, y_min, y_max)

        plt.draw()

        self.mask_computed = False

    def compute_mask(self):
        """Calls gbd.upslope_area() to compute the drainage basin mask (in place).
        """

        self.mask = initialize_mask(self.thk, self.x, self.y, self.terminus)

        dbg.upslope_area(self.x, self.y, self.z, self.mask)
        print "Drainage basin computation: done"
        self.mask_computed = True

        self.compute_bbox()

        self.plot_mask(1, cmaps.Blues)
        plt.show()

    def compute_bbox(self):
        """Computes the bounding box of the drainage basin and prints the NCO
        command that would cut it out of the whole-icesheet dataset.
        """

        try:
            self.border = int(self.entry.get())
        except:
            print "Invalid border width value: %s, using the default (5)." % self.entry.get()
            self.border = 5

        self.cutout_command = compute_bbox(self.input_file, self.mask, self.x, self.y, self.border)

        print "To cut out the drainage basin from the original dataset, run:"
        print self.cutout_command

    def compute_command(self):
        x_min, x_max, y_min, y_max = self.terminus

        ii = np.r_[0:self.x.size][(self.x >= x_min) & (self.x <= x_max)]
        jj = np.r_[0:self.y.size][(self.y >= y_min) & (self.y <= y_max)]

        i_min, i_max = ii[0], ii[-1]
        j_min, j_max = jj[0], jj[-1]

        cmd = "pism_regional.py -i %s -o %s -x %d,%d -y %d,%d -b %d" % (
            self.input_file, self.output_file, i_min, i_max, j_min, j_max, self.border)

        return cmd

def batch_process():
    """
    Process a file using command-line options (without a GUI).
    """

    from optparse import OptionParser

    parser = OptionParser()
    parser.usage = "usage: %prog -i foo.nc -o bar.nc --x_range ... --y_range ..."
    parser.description = "Computes the drainage basin mask given a DEM and a terminus location."

    parser.add_option("-x", "--x_range", dest="x_range", help="x_min,x_max (in grid indices)")
    parser.add_option("-y", "--y_range", dest="y_range", help="y_min,y_max (in grid indices)")
    parser.add_option("-b", "--border",  dest="border",  help="y_min,y_max (in grid indices)")
    parser.add_option("-i", dest="input", help="input file name")
    parser.add_option("-o", dest="output", help="output file name")

    (opts, args) = parser.parse_args()

    if (opts.input is None or opts.output is None or
        opts.x_range is None or opts.y_range is None):
        return

    sys.stderr.write("Loading data from %s..." % opts.input)
    x, y, z, thk = load_data(opts.input)
    sys.stderr.write("done.\n")

    i_min, i_max = map(lambda(x): int(x), opts.x_range.split(','))
    j_min, j_max = map(lambda(x): int(x), opts.y_range.split(','))

    sys.stderr.write("Initializing the mask...")
    mask = initialize_mask(thk, x, y, (x[i_min], x[i_max], y[j_min], y[j_max]))
    sys.stderr.write("done.\n")

    sys.stderr.write("Computing the drainage basin mask...")
    dbg.upslope_area(x, y, z, mask)
    sys.stderr.write("done.\n")

    sys.stderr.write("Computing the cutout command...")
    cutout_command = compute_bbox(opts.input, mask, x, y, int(opts.border))
    sys.stderr.write("done.\n")

    print "Command: %s\n" % cutout_command

    save_mask(opts.input, opts.output, mask, cutout_command,
              ' '.join(sys.argv))

    sys.exit(0)

if __name__ == "__main__":

    batch_process()             # calls sys.exit(0) if batch processing succeeded

    import matplotlib
    import matplotlib.cm as cmaps
    matplotlib.use('TkAgg')
    import pylab as plt

    from Tkinter import Tk, Frame, Label, Button, Entry, E, W
    import tkFileDialog

    root = Tk()
    root.wm_title("PISM drainage basin mask creator")

    a = App(root)

    root.mainloop()

