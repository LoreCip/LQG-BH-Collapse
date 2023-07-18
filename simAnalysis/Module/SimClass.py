import os
import re
from pathlib import Path

import warnings
import errno

import h5py

import numpy as np
import scipy.signal as sg
import scipy.optimize as op

import matplotlib.pyplot as plt

class _Sims():

    def __init__(self, _path):

        self.check(_path)

        self.root = _path
        self.sims = {}

        for path in Path(self.root).rglob('*.h5'):

            match = re.search(r"id[\d]_m[\d.]+_dx[\d.]+_xMax[\d.]+_tf[\d.]+_r[\d.]+_a0[\d.]+(?:_corrected(?:_finer)?)?", str(path))
            if match:
                key = match.group()

            self.sims[key] = _Sim(str(path)[:-18], _name=key)

        self.keylist = list(self.sims.keys())


    def check(self, _path):
        if not os.path.isdir(_path):
            raise NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), _path)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.sims[key]
        elif isinstance(key, int):
            return self.sims[self.keylist[key]]
        elif isinstance(key, float):
            key = int(key)
            return self.sims[self.keylist[key]]
        elif isinstance(key, list) or isinstance(key, np.ndarray):
            try:
                key = key.tolist()
            except AttributeError:
                pass
            return [self.__getitem__(k) for k in key]
        elif isinstance(key, slice):
            start = key.start or 0
            stop = key.stop or len(self.keylist)
            step = key.step or 1
            return [self.__getitem__(k) for k in range(start, stop, step)]

class _Sim():

    def __init__(self, _path, _name=None):

        self.check(_path)

        self.path = _path
        self.outpath = os.path.join(_path, 'outputs/output.h5')

        if _name == None:
            match = re.search(r"id[\d]_m[\d.]+_dx[\d.]+_xMax[\d.]+_tf[\d.]+_r[\d.]+_a0[\d.]+(?:_corrected(?:_finer)?)?", str(_path))
            if match:
                self.name = match.group()
        else:
            self.name = _name

        data = self.readParameters()
        self.id      = int(data[0])
        self.simtime = data[1]
        self.r0      = data[2]
        self.a0      = data[3]
        self.mass    = data[4]
        self.order   = int(data[5])
        self.xmax    = data[6]
        self.dx      = data[7]
        self.dSave   = int(data[8])
        self.dOut    = int(data[9])

        self.hor_loc = 2 * self.mass
        self.xb = self.xbounce()

        self.outfile = h5py.File(self.outpath, 'r')
        self.iterations = self.sort_groups()
        self.niter = len(self.iterations)

        self.xgrid = self.__getitem__('Xgrid')

        self.valid_keys = ['t', 'rho', 'e^b', 'B', 'beta', 'dt']
        

    def check(self, _path):
        if not os.path.isdir(_path):
            raise NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), _path)

        if not os.path.isfile(os.path.join(_path, 'ParameterFile.par')):
            raise FileNotFoundError(os.path.join(_path, 'ParameterFile.par'))

        if not os.path.isfile(os.path.join(_path, 'outputs/output.h5')):
            raise FileNotFoundError(os.path.join(_path, 'outputs/output.h5'))


    def __getitem__(self, key):
        item = self.outfile[key]
        if isinstance(item, h5py.Group):
            return GroupWrapper(item, self.xgrid)
        else:
            return item[()]

    def get(self, iteration, item):
        if isinstance(iteration, int):
            indices = [iteration]
        elif isinstance(iteration, list) or isinstance(iteration, np.ndarray):
            indices = iteration
        elif isinstance(iteration, slice):
            start = iteration.start or 0
            stop = iteration.stop or self.niter
            step = iteration.step or 1
            indices = range(start, stop, step)
        else:
            try:
                iteration = int(iteration)
            except ValueError:
                raise ValueError("Invalid iteration value: {} - It should be an integer or convertible to one.")
            indices = [iteration]

        result = [self.__getitem__(self.iterations[i])[item] for i in indices]

        if len(result) == 1:
            return result[0]
        return result

    def get_at_time(self, time, item, find_closest = True):
        
        if not isinstance(time, float):
            time = float(time)
        if time > self.simtime:
            raise ValueError("Requested time: {} - Maximum time available {}".format(time, self.simtime))
        if item not in self.valid_keys:
            raise IndexError("Invalid key: {} - Accepted keys: {}".format(item, self.valid_keys))

        # Max iteration possible if all timesteps where done with the largest timestep possible dt = dx/100
        max_iter = str(int(100 * time / self['dx']))
        # Find its index and add some iterations for safety
        try:
            idx = self.iterations.index(max_iter) + 20
        except ValueError as e:
            if not find_closest:
                raise e
            else:
                closest_idx = str(min(self.iterations, key=lambda x: abs(int(x) - int(max_iter))))
                idx = self.iterations.index(str(int(closest_idx))) + 20

        # Check distance from wanted time
        dt = np.abs(time - self.get(idx, 't')) 

        # Going backwards
        for i in range(idx-1, 0, -1):
            # Find new time
            t = self.get(i, 't')
            # Check its disntace from wanted time
            dt_tmp = np.abs(time - t)
            # If it is getting closer save it
            if dt_tmp < dt:
                dt = dt_tmp
            # If it is going further stop the cycle
            elif dt_tmp > dt:
                i = i + 1
                break

        return self.__getitem__(self.iterations[i])[item]

    def readParameters(self):
        comment_char = "/*"
        # Load the file using genfromtxt, ignoring comments
        return np.genfromtxt(
                                os.path.join(self.path, 'ParameterFile.par'),
                                comments=comment_char,
                                dtype = float
                            )
    
    def xbounce(self):

        def fxb(x, alpha, rs):
            return x**3 + alpha * x - rs
        
        try:
            return op.brentq(fxb, 
                        1e-6, 3*(self.hor_loc)*(1/3), 
                        args = (self.r0**2 / self.a0**2, self.hor_loc))
        except ValueError:
            psamp = np.linspace(self.xgrid[0], self.xgrid[-1], 1000)
            sampl = fxb(psamp, self.r0**2 / self.a0**2, self.hor_loc)
            for i in range(1000):
                if sampl[i] * samp[i+1] <= 0:
                    break
            return op.brentq(fxb, 
                        self.xgrid[i], self.xgrid[i+1], 
                        args = (self.r0**2 / self.a0**2, self.hor_loc))

    def sort_groups(self):
        group_names = list(self.outfile.keys())
        rl = []
        for i in range(len(group_names)):
            try:
                int(group_names[i])  # Try converting the element to an integer
            except ValueError:
                rl.append(i)  # Remove the element if it cannot be converted
        for r in rl[::-1]:
            group_names.pop(r)

        return sorted(group_names, key=int)

    def find_timeout(self):

        x = self.xgrid
        lx = len(x)

        for iter in range(self.niter-1, 0, -1):
            t = self.get(iter, 't')
            if t < 5:
                continue

            rho = self.get(iter, 'rho')

            cond1 = x >= self.xb
            skipped = len(cond1) - sum(cond1)
            cond2 = x <= self.hor_loc*2
            cond = cond1 & cond2

            rho = rho[cond]
            locmaxrho = self.find_peak(rho)
            locmaxrho += skipped

            if x[locmaxrho] <= self.hor_loc:
                return t

        return np.NaN

    def find_peak(self, rho, height = [1e-2]):
        idx_MAX = sg.find_peaks(rho, height=height, distance=2)[0][::-1]
        return idx_MAX[0]

    def plot(self, Y, X=None, iteration=None, time=None, color=None, linestyle=None, xrange=None, yrange=None, xlabel=None, ylabel=None, title=None, printTime=True, showHor=False, return_handles=False, savefig=False, name=None, savepath='.'):
        
        if iteration == None and time != None:
            qtt = self.get_at_time(time, Y)
            t = time
        elif time == None and iteration != None:    
            qtt = self.get(iteration, Y)
            t = self.get(iteration, 't')
        elif iteration == None and time == None:
            raise ValueError("Both iteration and time are missing. Please provide either iteration or time.")
        else:
            warnings.warn("Both iteration and time given. Using iteration.")
            qtt = self.get(iteration, Y)
            t = self.get(iteration, 't')

        if X == None:
            xgrid = self.xgrid
        else:
            xgrid = X

            if len(xgrid) != len(self.__getitem__('Xgrid')):
                raise ValueError("Length of provided grid is not equal to the length of the {} array.".format(item))

        if color != None:
            c = color
        else:
            c = 'k'
        if linestyle != None:
            ls = linestyle
        else:
            ls = '-'

        fig, ax = plt.subplots()
        plt.plot(xgrid, qtt, color = c, linestyle = ls)

        if showHor:
            plt.axvline(self.hor_loc, color = 'k', linestyle='--', alpha = 0.5)

        if xrange != None:
            if len(xrange) != 2:
                raise ValueError("Parameter xrange must be a tuple or list of lenght 2.")
            plt.xlim(xrange[0], xrange[1])
            
            if yrange == None:
                mask = (xgrid >= xrange[0]) & (xgrid <= xrange[1])
                bot = np.min(qtt[mask])
                top = np.max(qtt[mask])

                if np.abs(bot) > np.abs(top):
                    plt.ylim(bot*1.1, 0)
                else:
                    plt.ylim(0, top*1.1)

        if yrange != None:
            if len(yrange) != 2:
                raise ValueError("Parameter yrange must be a tuple or list of lenght 2.")
            plt.ylim(yrange[0], yrange[1])

        if xlabel != None:
            xl = str(xlabel)
        else:
            xl = 'Position'
        plt.xlabel(xl)
        
        if ylabel != None:
            yl = str(ylabel)
        else:
            yl = Y
        plt.ylabel(yl)

        if title != None:
            tl = str(title)
        else:
            tl = str(Y)
        plt.title(tl)

        if printTime:
            textstr = "Time: " + str(np.round(t, 3))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            plt.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

        if savefig:
            if name == None:
                name = f'{str(Y)}_time{np.round(t,3)}.png'
            if not name[:-4] in ['.png', '.jpg']:
                name = name + '.png'
            savepath = os.path.join(savepath, name)

            plt.savefig(savepath, format='png')

        if return_handles:
            return fig, ax

class GroupWrapper:
    def __init__(self, group, xgrid):
        self.group = group
        self.conv = False
        self.xgrid = xgrid

    def __getitem__(self, key):

        if key == 'beta':
            key = 'B'
            self.conv = True

        if not self.conv:
            return self.group[key][()]
        else:
            return self.group[key][()] / self.xgrid**2