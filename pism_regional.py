#!/usr/bin/env python
try:
    from netCDF3 import Dataset as NC
except:
    from netCDF4 import Dataset as NC

import matplotlib
matplotlib.use('TkAgg')

from Tkinter import *
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
    def __init__(self, master):
        self.master = master
        self.ph = None
        self.pts = None
        self.nc = None
        self.mask_computed = False
        self.Ncontours = 30

        self.create_widgets(master)
        master.geometry('+%d+%d' % (400, 200))

        self.load_data()

    def save_results(self):
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
                var.setncattr(attr, old_var.getncattr(attr))

        x[:] = self.x
        y[:] = self.y

        nc.variables['ftt_mask'][:] = (self.mask == 2)
        nc.close()

        print "Done."

    def load_data(self):
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

        plt.show(block=False)

        f = plt.get_current_fig_manager().window
        w, x0, y0 = f.winfo_width(), f.winfo_x(), f.winfo_y()

        self.master.geometry("+%d+%d" % (x0 + w, y0))
        self.master.lift()

    def get_output(self):
        output = tkFileDialog.asksaveasfilename(parent=root,
                                                filetypes = ["NetCDF .nc"],
                                                title="Save the mask in...")
        if len(output) > 0:
            return output
        else:
            return None

    def create_widgets(self, master):
        # The frame

        self.frame = Frame(master)
        self.frame.grid()

        button = Button(master, text="Select terminus location", command=self.get_terminus)
        button.grid(pady=2, row=1, column=1, sticky=E+W)

        button = Button(master, text="Compute the drainage basin mask", command=self.compute_mask)
        button.grid(pady=2, row=2, column=1, sticky=E+W)

        button = Button(master, text="Save the drainage basin mask", command=self.save_results)
        button.grid(pady=2, row=3, column=1, sticky=E+W)

        master.wm_resizable(False, False)

    def get_terminus(self):
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

        if self.ph is not None and self.mask_computed == False:
            for p in self.ph:
                p.remove()
            self.ph = None

        pts = []
        while len(pts) < 4:
            pts = np.asarray( plt.ginput(4, timeout=-1) )

        self.ph = plt.fill(pts[:,0], pts[:,1], 'white', lw = 2, alpha=0.5)
        plt.draw()

        self.pts = pts

        self.mask_computed = False

    def compute_mask(self):
        import matplotlib.nxutils as nx

        if self.pts is not None:
            def correct_mask(mask, x, y, pts):
                    for j in range(y.size):
                        for i in range(x.size):
                            if mask[j,i] > 0:
                                if nx.pnpoly(x[i], y[j], pts):
                                    mask[j,i] = 2
                                else:
                                    mask[j,i] = 1

            correct_mask(self.mask, self.x, self.y, self.pts)

        dbg.upslope_area(self.x, self.y, self.z, self.mask)
        print "Drainage basin computation: done"
        self.mask_computed = True

        plt.figure(1)
        plt.pcolormesh(self.x, self.y, self.mask)
        plt.contour(self.x, self.y, self.z, self.Ncontours, colors='black')
        plt.axis('tight')
        plt.axes().set_aspect('equal')
        plt.show()

if __name__ == "__main__":
    root = Tk()
    root.wm_title("PISM drainage basin mask creator")

    a = App(root)

    root.mainloop()

