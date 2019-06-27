import argparse
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('-x', default='nh')
    args = ap.parse_args()

    if args.x == 'nh':
        par_index = 0
        par_xlabel = 'nh (cm-3)'
    elif args.x == 'tc':
        par_index = 1
        par_xlabel = 'Tc (K)'
    elif args.x == 'lum':
        par_index = 2
        par_xlabel = 'lum (solar)'
    else:
        raise 'invalid -x option'

    dirs = [Path(d) for d in args.dirs]
    x = np.zeros(len(dirs))
    lab_curves = ['Te', 'e', 'p', 'H', 'H2']
    num_curves = len(lab_curves)
    y_cloudy = np.zeros((num_curves, len(dirs)))
    y_gasmod = np.zeros((num_curves, len(dirs)))

    # one point per dir
    for i, d in enumerate(dirs):
        par = np.loadtxt(d / 'parameters.dat')
        nh = par[par_index]
        x[i] = nh

        cloudy_ovr = pd.read_csv(d / 'hsphere.ovr', sep='\t')
        t = cloudy_ovr['Te'][0]
        heat = cloudy_ovr['Htot'][0]
        hden = cloudy_ovr['hden'][0]
        ne = cloudy_ovr['eden'][0]
        nhp = cloudy_ovr['HII'][0] * hden
        nh = cloudy_ovr['HI'][0] * hden
        nh2 = cloudy_ovr['2H_2/H'][0] * hden / 2

        y_cloudy[:, i] = (t, ne, nhp, nh, nh2)

        gasmod_ovr = np.loadtxt(d / 'MRNdust/overview.dat')
        y_gasmod[:, i] = gasmod_ovr[:]

    plt.figure()
    for i in range(num_curves):
        lc = plt.plot(x, y_cloudy[i], label='cloudy ' +
                      lab_curves[i], marker='x')
        plt.plot(x, y_gasmod[i], label='gasmod ' +
                 lab_curves[i], marker='+', color=lc[0].get_color())

    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel(par_xlabel)
    plt.legend()

    # one curve per dir
    # _________________

    # spectrum
    plt.figure()
    for d in dirs:
        cloudy_diffuse = np.loadtxt(
            d / 'all_4pi_nu_jnu.erg_cm-3_s-1.out', max_rows=2)
        x = cloudy_diffuse[0]
        y = cloudy_diffuse[1]

        par = np.loadtxt(d / 'parameters.dat')
        nh = par[par_index]
        x[i] = nh
        lc = plt.loglog(x, y, label='cloudy {}'.format(nh))

        mrn_opt = np.loadtxt(d / 'MRNdust' / 'opticalProperties.dat')
        plt.loglog(mrn_opt[:, 1], mrn_opt[:, 2], label='gasmodule',
                   ls='dashed', color=lc[0].get_color())

    plt.xlabel('$\\lambda$ (micron)')
    plt.ylabel('outward emission')
    plt.gcf().legend()

    # populations
#    plt.figure()
    fig, axs = plt.subplots(1, 2)
    for i, d in enumerate(dirs):
        MRN_DIR = d / 'MRNdust'

        par = np.loadtxt(d / 'parameters.dat')
        nh = par[par_index]

        cloudy_pops = np.loadtxt(d / 'hpopulations.cm-3.out', ndmin=2)
        y = cloudy_pops[0, 3:]
        x = range(len(y))
        # lc = axs[0].semilogy(x, y, label='cloudy {}'.format(nh), marker='x')

        color = lc[0].get_color()

        hpop = np.loadtxt(MRN_DIR / 'hpopulations.dat')
        # axs[0].semilogy(range(hpop.shape[0]), hpop[:, 2], label='gasmodule', marker='+', color=color)

        maxindex = min(len(y), len(hpop[:, 2]))
        axs[0].semilogy(range(maxindex), hpop[:maxindex, 2] /
                        y[:maxindex], label='ratio at {}'.format(nh), marker='+')

        cloudy_h2pops = pd.read_csv(d / 'h2populations.frac.out', sep='\t')
        pops_h2 = cloudy_h2pops['pops/H2'].to_numpy()[3:]
        energy = cloudy_h2pops['energy(wn)'].to_numpy()[3:]
        yc = pops_h2
        x = range(len(yc))
        # axs[1].semilogy(x, y, label='cloudy', marker='x', color=color)

        h2pop = np.loadtxt(MRN_DIR / 'h2populations.dat')
        yg = h2pop[:, 2]
        yg /= sum(yg)
        x = range(len(yg))
        # axs[1].semilogy(x, y, label='gasmodule', marker='+', color=color)

        maxindex = min(len(yg), len(yc))
        axs[1].semilogy(x[:maxindex], yg[:maxindex] /
                        yc[:maxindex], label='gasmodule', marker='+')

    axs[0].legend()
    axs[0].set_xlabel('index')
    axs[0].set_ylabel('population density (cm-3)')
    axs[0].set_title('H')

    axs[1].set_xlabel('index')
    axs[1].set_title('H2')

    plt.show()


if __name__ == '__main__':
    main()
