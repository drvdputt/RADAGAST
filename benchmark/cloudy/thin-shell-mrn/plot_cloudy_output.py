import pandas as pd
from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.rc('lines', linewidth=1)

MRN_DIR = Path('../../../../run/MRNDust')
OUTPUT_DIR = Path('plots')


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    compare_equilibrium()
    compare_grain_size_distribution()
    compare_incident_radiation_field()
    compare_emission()
    compare_opacity()
    compare_populations()


def compare_grain_size_distribution():
    plt.figure()

    # cloudy
    cloudy_abun = np.loadtxt('grainabundance.g_cm-3.out', ndmin=2)
    cloudy_size = np.loadtxt('grainsize.out', ndmin=2)
    cloudy_num = np.loadtxt('grainnumberdensity.out', ndmin=2)
    # first column = depth, last column = total
    x = cloudy_size[0, 1:]
    # y = cloudy_abun[0, 1:-1]
    y = cloudy_num[0, 1:]
    plt.loglog(x, y, marker='+', label='cloudy')

    # mrn test
    mrn_abun = np.loadtxt(MRN_DIR / 'grainpop_0.dat')
    x = mrn_abun[:, 0]
    y = mrn_abun[:, 1]
    plt.loglog(x, y, marker='x', label='MRNtest')

    plt.xlabel('average radius')
    plt.ylabel('number density')
    plt.title('grain abundance')
    plt.legend()
    plt.savefig(OUTPUT_DIR / 'grain_size.pdf')
    plt.close()


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
    plt.close()


def compare_emission():
    plt.figure()

    # cloudy_cont = pd.read_csv('continuum.out', sep='\t')
    # x = cloudy_cont['#Cont  nu']
    # y = cloudy_cont['DiffOut']

    # cloudy_diffuse = pd.read_csv('4pi_nu_jnu.erg_cm-3_s-1.out', sep='\t')
    # x = cloudy_diffuse['#energy/um']
    # y = cloudy_diffuse['Total']

    cloudy_diffuse = np.loadtxt('all_4pi_nu_jnu.erg_cm-3_s-1.out', max_rows=2)
    x = cloudy_diffuse[0]
    y = cloudy_diffuse[1]

    max_y = max(y)
    lower_limit = 1e-20 * max_y
    plt.loglog(x, y, label='cloudy')

    mrn_opt = np.loadtxt(MRN_DIR / 'opticalProperties.dat')
    plt.loglog(mrn_opt[:, 1], mrn_opt[:, 2], label='gasmodule', ls='dashed')
    plt.xlabel('$\\lambda$ (micron)')
    plt.ylabel('outward emission')
    plt.legend()
    plt.ylim([lower_limit, max_y])
    plt.xlim([5.e-2, 1.e6])
    plt.savefig(OUTPUT_DIR / 'emission.pdf')
    # plt.close()


def compare_opacity():
    plt.figure()
    cloudy_op = pd.read_csv('opacity.out', sep='\t')
    x = cloudy_op['#nu/um']
    y = cloudy_op['Tot opac']
    plt.loglog(x, y, label='cloudy')

    mrn_opt = np.loadtxt(MRN_DIR / 'opticalProperties.dat')
    plt.loglog(mrn_opt[:, 1], mrn_opt[:, 3], label='gasmodule')
    plt.xlabel('$\\lambda$ (micron)')
    plt.ylabel('opacity')
    plt.legend()
    plt.savefig(OUTPUT_DIR / 'opacity.pdf')
    plt.close()


def compare_populations():
    fig, axs = plt.subplots(1, 2)

    cloudy_pops = np.loadtxt('hpopulations.cm-3.out', ndmin=2)
    y = cloudy_pops[0, 3:]
    x = range(len(y))
    axs[0].semilogy(x, y, label='cloudy')

    hpop = np.loadtxt(MRN_DIR / 'hpopulations.dat')
    axs[0].semilogy(range(hpop.shape[0]), hpop[:, 2], label='gasmodule')

    axs[0].set_xlabel('index')
    axs[0].set_ylabel('population density (cm-3)')
    axs[0].set_title('H')

    cloudy_h2pops = pd.read_csv('h2populations.frac.out', sep='\t')
    pops_h2 = cloudy_h2pops['pops/H2'].to_numpy()[3:]
    energy = cloudy_h2pops['energy(wn)'].to_numpy()[3:]
    y = pops_h2
    x = range(len(y))
    axs[1].semilogy(x, y, label='cloudy')

    h2pop = np.loadtxt(MRN_DIR / 'h2populations.dat')
    y = h2pop[:, 2]
    y /= sum(y)
    x = range(len(y))

    axs[1].semilogy(x, y, label='gasmodule')

    axs[1].set_xlabel('index')
    axs[1].set_title('H2')

    axs[0].legend()
    fig.savefig(OUTPUT_DIR / 'hpopulations.pdf')
    plt.show()


def compare_equilibrium():
    cloudy_ovr = pd.read_csv('overview.dat', sep='\t')
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

    cloudy_grainh2rate = pd.read_csv('grainH2rate.out', sep='\t')
    h2form_g = sum(cloudy_grainh2rate.to_numpy()[0, 1:])
    # I want to the total formation per H atom. Therefore, we do not use
    # the 'coef' options here. Using give rather large numbers (order of
    # 1). Probably has to do with grain number density.

    # Alternatively, I could sum the formation rates from the grain
    # H2rate file. They are in units s-1 and the source code states that
    # multiplying them with hden will yield the total formation rate
    # cm-3 s-1.

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

    print("h2form = ", h2form_g)
    print("h2dissoc = ", h2dissoc)
    print("hphotoion = ", hphotoion)
    print("hcolion = ", hcolion)
    print("hrec = ", hrec)

    print("Te = ", t)


if __name__ == "__main__":
    main()
