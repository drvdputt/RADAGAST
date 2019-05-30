import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd

MRN_DIR = Path('../../../../run/MRNDust')
OUTPUT_DIR = Path('plots')


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    compare_grain_size_distribution()
    compare_incident_radiation_field()
    compare_equilibrium()

def compare_grain_size_distribution():
    plt.figure()

    # cloudy
    cloudy_abun = np.loadtxt('grainabundance.g_cm-3.out')
    # first column = depth, last column = total
    x = range(cloudy_abun.shape[1])[1:-1]
    y = cloudy_abun[0, 1:-1]
    plt.semilogy(x, y, marker='+')

    # mrn test
    mrn_abun = np.loadtxt(MRN_DIR / 'grainpop_0.dat')
    x = range(1, mrn_abun.shape[0] + 1)
    y = mrn_abun[:, 2]
    plt.semilogy(x, y, marker='x')

    plt.xlabel('grain pop number')
    plt.ylabel('abundance (g_cm-3)')
    plt.title('grain abundance')
    plt.savefig(OUTPUT_DIR / 'grain_size.pdf')


def compare_incident_radiation_field():
    plt.figure()
    cloudy_cont = pd.read_csv('continuum.out', sep='\t')
    x = cloudy_cont['#Cont  nu']
    y = cloudy_cont['incident']
    plt.loglog(x, y)

    mrn_cont = np.loadtxt(MRN_DIR / 'nu_jnu.dat')
    plt.loglog(mrn_cont[:, 0], mrn_cont[:, 1])

    plt.xlabel('$\\lambda$ (micron)')
    plt.ylabel('indicent radiation field $\\nu L_{\\nu}$ (erg / s)')
    plt.savefig(OUTPUT_DIR / 'radiation_field.pdf')


def compare_equilibrium():
    cloudy_ovr = pd.read_csv('hsphere.ovr', sep='\t')
    t = cloudy_ovr['Te'][0]
    heat = cloudy_ovr['Htot'][0]
    hden = cloudy_ovr['hden'][0]
    ne = cloudy_ovr['eden'][0]
    np = cloudy_ovr['HII'][0] * hden
    nh = cloudy_ovr['HI'][0] * hden
    nh2 = cloudy_ovr['2H_2/H'][0] * hden / 2
    print("Te = ", t)
    print("Htot = ", heat)
    print("eden = ", ne)
    print("H+ = ", np)
    print("HI = ", nh)
    print("H2 = ", nh2)


if __name__ == "__main__":
    main()
