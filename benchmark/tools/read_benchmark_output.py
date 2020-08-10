"""Collection of functions that grab some commonly needed quantities
when plotting the benchmark results."""

from pathlib import Path
import numpy as np
import pandas as pd


class BenchmarkResult:
    def __init__(self, bench_directory):
        self.d = Path(bench_directory)
        self.parameters = np.loadtxt(self.d / "parameters.dat")
        self.cloudy = CloudyResult(self.d)
        self.gasmodule = GasModuleResult(self.d / "MRNDust")

    def get_nh_tc_lum(self):
        """get the parameters for which the benchmark was run: nH, Tc and luminosity"""
        return self.parameters


class CloudyResult:
    def __init__(self, directory):
        self.d = Path(directory)
        self.ovr = pd.read_csv(self.d / "overview.dat", sep="\t")
        self.ok = len(self.ovr) >= 1 and self.ovr.to_numpy().dtype == np.float64
        if not self.ok:
            print("Warning: problem with cloudy for {}".format(directory))

    def get_temp(self, depth_index=0):
        if self.ok:
            t = self.ovr["Te"][depth_index]
        else:
            t = 0
        return t

    def get_densities(self, depth_index=0, numpy=False):
        if self.ok:
            hden = self.ovr["hden"][depth_index]
            vals = [
                self.ovr["eden"][depth_index],
                self.ovr["HII"][depth_index] * hden,
                self.ovr["HI"][depth_index] * hden,
                self.ovr["2H_2/H"][depth_index] * hden / 2,
            ]
        else:
            vals = [0.0] * 4

        if numpy:
            return np.array(vals)
        else:
            keys = ["e", "H+", "H", "H2"]
            return {key: value for (key, value) in zip(keys, vals)}

    def get_emissivity(self, depth_index=0):
        """returns wavelength, 4pi nu j_nu"""
        if self.ok:
            cloudy_diffuse = np.loadtxt(
                self.d / "all_4pi_nu_jnu.erg_cm-3_s-1.out", max_rows=2
            )
            return cloudy_diffuse[0], cloudy_diffuse[depth_index + 1]
        else:
            return [0.0], [0.0]

    def get_h_populations(self):
        """indexed arbitrarily, but should be the same for cloudy and gas module"""
        if self.ok:
            cloudy_pops = np.loadtxt(self.d / "hpopulations.cm-3.out", ndmin=2)
            y = cloudy_pops[0][3:]
            return y / sum(y)
        else:
            return np.zeros(1)

    def get_h2_populations(self):
        """indexed arbitrarily, but should be the same for cloudy and gas module"""
        if self.ok:
            cloudy_h2pops = pd.read_csv(self.d / "h2populations.frac.out", sep="\t")
            # these should already be population fractions
            pops_h2 = cloudy_h2pops["pops/H2"].to_numpy()[3:]
        else:
            pops_h2 = [0.0]

        return pops_h2

    def get_h2_rates(self):
        """get formation, dissocation rates"""
        cloudy_grainh2rate = pd.read_csv(self.d / "grainH2rate.out", sep="\t")
        h2form_g = sum(cloudy_grainh2rate.to_numpy()[0, 1:])

        cloudy_h2destruction = pd.read_csv(self.d / "h2destruction.out", sep="\t")
        h2dissoc = cloudy_h2destruction["PHOTON,H2=>H,H"][0]

        return h2form_g, h2dissoc

    def get_rates(self):
        """ in this order: h2 form, h2 dissoc, h ion, h rec """

        h2rates = self.get_h2_rates()
        hionization = pd.read_csv(
            self.d / "hionization.s-1.cm-3_s-1.out", sep="\t", nrows=1
        ).dropna()
        ne, n_p = self.get_densities(numpy=True)[[0, 1]]
        ion = hionization["gam1"][0]
        rec = hionization["RecTot"][0] / ne / n_p
        return np.array([h2rates[0], h2rates[1], ion, rec])

    def get_grain_temps(self):
        """get list of grain temperatures"""
        cloudy_graintemps = pd.read_csv(self.d / "graintemperature.out", sep="\t")
        return cloudy_graintemps.to_numpy()[0, 1:]

    def get_heat_or_cool(self, is_heat):
        if is_heat:
            fname = "heating.out"
            i_total = 2
            labels = ["H  1", "H2dH", "GrnP", "H2vH"]
        else:
            fname = "cooling.out"
            i_total = 3
            labels = [
                "H2ln 0.0",
                "dust 0.0",
                "FF c 0.0",
                "ISrcolH H  0.0",
                "ISclinH H  0.0",
                "H2cX 0.0",
            ]

        with open(self.d / fname, "r") as f:
            f.readline()  # skip legend
            data_line = f.readline()

        contributions = {}
        fragments = data_line.split("\t")
        total = float(fragments[i_total])
        for l in labels:
            contributions[l] = total * find_fraction(fragments, l)

        return contributions

    def get_heat(self):
        """ get heating in this order: h photo, H2 diss, grain photo, h2ln """
        dic = self.get_heat_or_cool(True)
        return np.array([dic[k] for k in ["H  1", "H2dH", "GrnP", "H2vH"]])

    def get_cool(self):
        """ get cooling in this order: H rec, H line, H ff, H2 line + H2cx, dust col """
        dic = self.get_heat_or_cool(False)
        result = [
            dic[k]
            for k in [
                "ISrcolH H  0.0",
                "ISclinH H  0.0",
                "FF c 0.0",
                "H2ln 0.0",
                "dust 0.0",
            ]
        ]
        result[3] += dic["H2cX 0.0"]
        return np.array(result)


class GasModuleResult:
    def __init__(self, directory):
        self.d = Path(directory)
        self.ovr = np.loadtxt(self.d / "overview.dat")

    def get_temp(self):
        return self.ovr[0]

    def get_densities(self, numpy=False):
        if numpy:
            return self.ovr[1:5]
        else:
            return {
                "e": self.ovr[1],
                "H+": self.ovr[2],
                "H": self.ovr[3],
                "H2": self.ovr[4],
            }

    def get_emissivity(self):
        """returns wavelength and 4pi nu j_nu"""
        optical = np.loadtxt(self.d / "opticalProperties.dat")
        return optical[:, 1], optical[:, 2]

    def get_opacity(self):
        optical = np.loadtxt(self.d / "opticalProperties.dat")
        return optical[:, 1], optical[:, 3]

    def get_h_populations(self):
        hpop = np.loadtxt(self.d / "hpopulations.dat")
        y = hpop[:, 2]
        return y / sum(y)

    def get_h2_populations(self):
        h2pop = np.loadtxt(self.d / "h2populations.dat")
        y = h2pop[:, 2]
        return y / sum(y)

    def get_h2_rates(self):
        rates = np.loadtxt(self.d / "rates.dat")
        return rates[:2]

    def get_rates(self):
        rates = np.loadtxt(self.d / "rates.dat")
        return rates[[0, 1, 2, 4]]

    def get_grain_temps(self):
        """get list of grain temperatures"""
        grainprops = np.loadtxt(self.d / "grainpop_0.dat")
        return grainprops[:, 3]

    def get_heat(self):
        """ get heating in this order: h photo, H2 diss, grain photo, h2ln """
        heats = np.loadtxt(self.d / "heat.dat")
        return heats[[1, 3, 4, 2]]

    def get_cool(self):
        """ get cooling in this order: H rec, H line, H ff, H2 line + H2cx, dust col """
        cools = np.loadtxt(self.d / "cool.dat")
        return cools[[2, 0, 3, 1, 4]]


def find_fraction(line_fragments, label):
    """Helper function for heating and cooling files, which extracts the
    heating/cooling fraction from a line of the file.

    line_fragments: list containing parts of the lines, after splitting
    on tab

    label: name of the heating contribution we want to find in the line

    """
    try:
        index = line_fragments.index(label)
    except ValueError:
        return 0

    return float(line_fragments[index + 1])
