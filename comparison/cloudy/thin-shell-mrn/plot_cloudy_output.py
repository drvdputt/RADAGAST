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
    compare_emission()
    compare_equilibrium()


def compare_grain_size_distribution():
    plt.figure()

    # cloudy
    cloudy_abun = np.loadtxt('grainabundance.g_cm-3.out', ndmin=2)
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


def compare_emission():
    plt.figure()
    cloudy_cont = pd.read_csv('continuum.out', sep='\t')
    x = cloudy_cont['#Cont  nu']
    y = cloudy_cont['DiffOut']
    plt.loglog(x, y, label='cloudy')

    mrn_opt = np.loadtxt(MRN_DIR / 'opticalProperties.dat')
    plt.loglog(mrn_opt[:, 1], mrn_opt[:, 2], label='gasmodule')
    plt.xlabel('$\\lambda$ (micron)')
    plt.ylabel('outward emission')
    plt.legend()
    plt.savefig(OUTPUT_DIR / 'emission.pdf')
    # plt.show()


def compare_equilibrium():
    cloudy_ovr = pd.read_csv('hsphere.ovr', sep='\t')
    t = cloudy_ovr['Te'][0]
    heat = cloudy_ovr['Htot'][0]
    hden = cloudy_ovr['hden'][0]
    ne = cloudy_ovr['eden'][0]
    np = cloudy_ovr['HII'][0] * hden
    nh = cloudy_ovr['HI'][0] * hden
    nh2 = cloudy_ovr['2H_2/H'][0] * hden / 2

    cloudy_gheat = pd.read_csv('grainheating.erg_cm-3_s-1.out', sep='\t')
    gheat = sum(cloudy_gheat.to_numpy()[0, 1:])

    cloudy_h2creation = pd.read_csv('h2creation.out', sep='\t')
    h2form = cloudy_h2creation['grn,H,H=>grn,H2'][0] / nh
    # I want to the total formation per H atom. Therefore, we do not use
    # the 'coef' options here. Using give rather large numbers (order of
    # 1). Probably has to do with grain number density.

    cloudy_h2destruction = pd.read_csv('h2destruction.out', sep='\t')
    h2dissoc = cloudy_h2destruction['PHOTON,H2=>H,H'][0]
    # h2destruction was saved with the 'coef' option. According to hazy
    # 1, a 2-body reaction would have a coefficient of unit cm3 s-1 when
    # this option is used. Multiplying with the densities of the 2
    # bodies yields cm-3 s-1. The destruction by photons is a 1-body
    # reaction, so this gives us something in s-1 for this specific
    # reaction. So using 'coef' for the destruction gives us the same
    # coefficient as my gas code

    cloudy_hionization = pd.read_csv(
        'hionization.s-1.cm-3_s-1.out', sep='\t', nrows=1)
    hphotoion = cloudy_hionization['gam1'][0]
    hcolion = cloudy_hionization['coll ion1'][0]
    hrec = cloudy_hionization['RecTot'][0]

    print("Htot = ", heat)
    print("grainHeat = ", gheat)
    print("eden = ", ne)
    print("H+ = ", np)
    print("HI = ", nh)
    print("H2 = ", nh2)

    print("h2form = ", h2form)
    print("h2dissoc = ", h2dissoc)
    print("hphotoion = ", hphotoion)
    print("hcolion = ", hcolion)
    print("hrec = ", hrec)

    print("Te = ", t)


if __name__ == "__main__":
    main()
