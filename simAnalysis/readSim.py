import os
import re
from pathlib import Path
import errno

import h5py

import numpy as np
import scipy.signal as sg

class Sims():

    def __init__(self, _path):

        self.check(_path)

        self.root = _path
        self.sims = {}

        for path in Path(self.root).rglob('*.h5'):

            match = re.search(r"id0_m[\d.]+_dx[\d.]+_xMax[\d.]+_tf[\d.]+_r[\d.]+_a0[\d.]+(?:_corrected(?:_finer)?)?", str(path))
            if match:
                key = match.group()

            self.sims[key] = Sim(str(path)[:-18])

        self.simslist = list(self.sims.keys())


    def check(self, _path):
        if not os.path.isdir(_path):
            raise NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), _path)

    def __getitem__(self, key):
        return self.sims[key]



class Sim():

    def __init__(self, _path):

        self.check(_path)

        self.path = _path
        self.outpath = os.path.join(_path, 'outputs/output.h5')

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

        self.outfile = h5py.File(self.outpath, 'r')
        self.iterations = self.sort_groups()
        self.niter = len(self.iterations)

        self.valid_keys = ['t', 'rho', 'e^b', 'B', 'dt']
        

    def __getitem__(self, key):
        item = self.outfile[key]
        if isinstance(item, h5py.Group):
            return GroupWrapper(item)  # Wrap the group in a custom class
        else:
            return item[()]

    def get(self, iteration, item):
        if isinstance(iteration, slice):
            start = iteration.start or 0
            stop = iteration.stop or self.niter
            step = iteration.step or 1
            indices = range(start, stop, step)
        elif isinstance(iteration, int):
            indices = [iteration]
        else:
            try:
                iteration = int(iteration)
                indices = [iteration]
            except ValueError:
                raise ValueError("Invalid iteration value: {} - It should be an integer or convertible to one.")

        result = [self.__getitem__(self.iterations[i])[item] for i in indices]

        if len(result) == 1:
            return result[0]
        return result


    # def get(self, iteration, item):

    #     if not isinstance(iteration, int):
    #         try:
    #             iteration = int(iteration)
    #         except ValueError:
    #             raise ValueError("Invalid iteration value: {} - It should be an integer or convertible to one.".format(iteration))

    #     if iteration >= self.niter: 
    #         raise IndexError("Invalid index: {} - Array length: {}".format(iteration, self.niter))
    #     if item not in self.valid_keys:
    #         raise IndexError("Invalid key: {} - Accepted keys: {}".format(item, self.valid_keys))
        
    #     return self.__getitem__(self.iterations[iteration])[item]

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


    def check(self, _path):
        if not os.path.isdir(_path):
            raise NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), _path)

        if not os.path.isfile(os.path.join(_path, 'ParameterFile.par')):
            raise FileNotFoundError(os.path.join(_path, 'ParameterFile.par'))

        if not os.path.isfile(os.path.join(_path, 'outputs/output.h5')):
            raise FileNotFoundError(os.path.join(_path, 'outputs/output.h5'))


    def readParameters(self):
        comment_char = "/*"
        # Load the file using genfromtxt, ignoring comments
        return np.genfromtxt(
                                os.path.join(self.path, 'ParameterFile.par'),
                                comments=comment_char,
                                dtype = float
                            )
        
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

        hor_loc = 2 * self.mass
        x_min = hor_loc**(1/3)

        x = self['Xgrid']
        lx = len(x)
        x = x[:lx//2]

        for iter in range(self.niter-1, 0, -1):
            t = self.get(iter, 't')
            if t < 5:
                continue

            rho = self.get(iter, 'rho')[:lx//2]

            cond1 = x > x_min
            skipped = len(cond1) - sum(cond1)
            cond2 = x < hor_loc*1.3
            cond = cond1 & cond2

            rho = rho[cond]
            locmaxrho = self.find_peak(rho)
            locmaxrho += skipped

            if x[locmaxrho] <= hor_loc:
                return t

        return np.NaN

    def find_peak(self, rho):
        idx_MAX = sg.find_peaks(rho, height=[1e-2], distance=2)[0][::-1]
        return idx_MAX[0]

class GroupWrapper:
    def __init__(self, group):
        self.group = group

    def __getitem__(self, key):
        return self.group[key][()]