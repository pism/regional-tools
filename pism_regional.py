#!/usr/bin/env python
try:
    from netCDF3 import Dataset as NC
except:
    from netCDF4 import Dataset as NC

import matplotlib
matplotlib.use('TkAgg')

from Tkinter import Tk, Frame, Label, Button, Entry, E, W
import tkFileDialog

import numpy as np
import pylab as plt

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

class App:
    """An application class containing methods of the drainage basin tool.
    """
    def __init__(self, master):
        self.master = master
        self.fill = None
        self.terminus = None
        self.nc = None
        self.mask_computed = False
        self.Ncontours = 30

        self.create_widgets(master)

        self.load_data()

    def save_results(self):
        """ Saves the computed drainage basin mask to a file.
        """
        if self.nc is None:
            return

        output_file = self.get_output()

        if output_file is None:
            print "No output file selected; cannot proceed."
            return

        print "Saving the mask to %s" % output_file

        nc = NC(output_file, 'w')

        nc.createDimension('x', self.x.size)
        nc.createDimension('y', self.y.size)

        x = nc.createVariable("x", 'f8', ('x',))
        y = nc.createVariable("y", 'f8', ('y',))
        mask = nc.createVariable("ftt_mask", 'i4', ('y', 'x'))

        mask.long_name = "Drainage basin area for regional modeling"

        x_orig = self.nc.variables['x']
        y_orig = self.nc.variables['y']

        for var, old_var in zip([x,y], [x_orig, y_orig]):
            for attr in old_var.ncattrs():
                value = old_var.getncattr(attr)
                if isinstance(value, (str, unicode)):
                    value = str(value.encode('ASCII', 'ignore'))
                var.setncattr(attr, value)

        x[:] = self.x
        y[:] = self.y

        nc.variables['ftt_mask'][:] = (self.mask == 2)
        nc.close()

        print "Done."

    def load_data(self):
        """Loads data from an input file.

        An input file has to contain variables 'x', 'y', 'usurf', 'thk'.
        """
        self.input_file = tkFileDialog.askopenfilename(parent=root,
                                                       filetypes = ["NetCDF .nc"],
                                                       title='Choose an input file')

        if len(self.input_file) == 0:
            print "No input file selected. Exiting..."
            exit(0)

        self.nc = NC(self.input_file)
        nc = self.nc

        self.x = np.array(nc.variables['x'][:], dtype=np.double)
        self.y = np.array(nc.variables['y'][:], dtype=np.double)
        self.z = np.array(np.squeeze(permute(nc.variables['usurf'])), dtype=np.double, order='C')
        self.thk = np.array(np.squeeze(permute(nc.variables['thk'])), dtype=np.double, order='C')

        self.mask = dbg.initialize_mask(self.thk)
        print "Mask initialization: done"

        plt.figure(1)
        plt.pcolormesh(self.x, self.y, self.mask)
        plt.contour(self.x, self.y, self.z, self.Ncontours, colors='black')
        plt.axis('tight')
        plt.axes().set_aspect('equal')
        plt.xticks([])
        plt.yticks([])

        try:
            # block=False is not supported in some earlier versions of matplotlib
            plt.show(block=False)

            f = plt.get_current_fig_manager().window
            w, x0, y0 = f.winfo_width(), f.winfo_x(), f.winfo_y()

            self.master.geometry("+%d+%d" % (x0 + w + 10, y0))
            self.master.lift()
        except:
            plt.show()

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
            plt.pcolormesh(self.x, self.y, self.mask)
            plt.contour(self.x, self.y, self.z, self.Ncontours, colors='black')
            plt.axis('tight')
            plt.axes().set_aspect('equal')
            plt.draw()

        plt.setp(plt.gca(),autoscale_on=False)

        cursor = Cursor(plt.axes(), useblit=True, color='white', linewidth=1 )

        # remove the rectangle
        if self.fill is not None and self.mask_computed == False:
            for p in self.fill:
                p.remove()
            self.fill = None

        x_min, x_max, y_min, y_max = plt.axis()

        x0, y0 = plt.ginput(timeout=-1)[0]

        l1 = plt.plot([x0, x0], [y_min, y_max], color='white')
        l2 = plt.plot([x_min, x_max], [y0, y0], color='white')

        x1, y1 = plt.ginput(timeout=-1)[0]

        l3 = plt.plot([x1, x1], [y_min, y_max], color='white')
        l4 = plt.plot([x_min, x_max], [y1, y1], color='white')

        dx = x1 - x0
        dy = y1 - y0

        xs = [x0, x1, x1, x0]
        ys = [y0, y0, y1, y1]
        self.fill = plt.fill(xs, ys, 'white', lw = 2, alpha=0.5)

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
        x_min, x_max, y_min, y_max = self.terminus

        if self.terminus is not None:
            def correct_mask(mask, x, y):
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

            correct_mask(self.mask, self.x, self.y)

        dbg.upslope_area(self.x, self.y, self.z, self.mask)
        print "Drainage basin computation: done"
        self.mask_computed = True

        self.compute_bbox()

        plt.figure(1)
        plt.pcolormesh(self.x, self.y, self.mask)
        plt.contour(self.x, self.y, self.z, self.Ncontours, colors='black')
        plt.axis('tight')
        plt.axes().set_aspect('equal')
        plt.show()

    def compute_bbox(self):
        """Computes the bounding box of the drainage basin and prints the NCO
        command that would cut it out of the whole-icesheet dataset.
        """
        x = self.x
        y = self.y

        x0 = x.size - 1
        x1 = 0
        y0 = y.size - 1
        y1 = 0

        for j in range(y.size):
            for i in range(x.size):
                if self.mask[j,i] == 2:
                    if x[i] < x[x0]:
                        x0 = i

                    if x[i] > x[x1]:
                        x1 = i

                    if y[j] < y[y0]:
                        y0 = j

                    if y[j] > y[y1]:
                        y1 = j

        try:
            border = int(self.entry.get())
        except:
            print "Invalid border width value: %s, using the default (5)." % self.entry.get()
            border = 5

        x0 = np.maximum(x0 - border, 0)
        x1 = np.minimum(x1 + border, x.size - 1)

        y0 = np.maximum(y0 - border, 0)
        y1 = np.minimum(y1 + border, y.size - 1)

        print "To cut out the drainage basin from the original dataset, run:"
        print "ncks -d x,%d,%d -d y,%d,%d %s output.nc" % (x0, x1, y0, y1, self.input_file)


if __name__ == "__main__":
    root = Tk()
    root.wm_title("PISM drainage basin mask creator")

    a = App(root)

    root.mainloop()

