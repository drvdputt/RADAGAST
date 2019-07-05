"""Collection of functions that grab some commonly needed quantities
when plotting the benchmark results."""

from pathlib import Path
import numpy as np
import pandas as pd


class BenchmarkResult:
    def __init__(self, bench_directory):
        self.d = Path(bench_directory)
        self.parameters = np.loadtxt(self.d / 'parameters.dat')
        self.cloudy = CloudyResult(self.d)
        self.gasmodule = GasModuleResult(self.d / 'MRNdust')

    def get_nh_tc_lum(self):
        """get the parameters for which the benchmark was run: nH, Tc and luminosity"""
        return self.parameters


class CloudyResult:
    def __init__(self, directory):
        self.d = Path(directory)
        self.ovr = pd.read_csv(self.d / 'hsphere.ovr', sep='\t')

    def get_densities(self, depth_index=0, numpy=False):
        hden = self.ovr['hden'][depth_index]
        vals = [self.ovr['eden'][depth_index],
                self.ovr['HII'][depth_index] * hden,
                self.ovr['HI'][depth_index] * hden,
                self.ovr['2H_2/H'][depth_index] * hden / 2]
        if numpy:
            return np.array(vals)
        else:
            keys = ['e', 'H+', 'H', 'H2']
            return {key: value for (key, value) in zip(keys, vals)}


class GasModuleResult:
    def __init__(self, directory):
        self.d = Path(directory)
        self.ovr = np.loadtxt(self.d / 'overview.dat')

    def get_densities(self, numpy=False):
        if numpy:
            return self.ovr[1:5]
        else:
            return {'e': self.ovr[1],
                    'H+': self.ovr[2],
                    'H': self.ovr[3],
                    'H2': self.ovr[4]}
