import argparse
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import pandas as pd

from read_benchmark_output import BenchmarkResult

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
    lab_curves = ['e', 'H+', 'H', 'H2']
    num_curves = len(lab_curves)
    y_cloudy = np.zeros((num_curves, len(dirs)))
    y_gasmod = np.zeros((num_curves, len(dirs)))

    # one point per dir
    # -----------------
    for i, d in enumerate(dirs):
        br = BenchmarkResult(d)
        x[i] = br.get_nh_tc_lum()[par_index]
        y_cloudy[:, i] = br.cloudy.get_densities(numpy=True)
        y_gasmod[:, i] = br.gasmodule.get_densities(numpy=True)

    fig, axs = plt.subplots(2, 1, sharex=True)
    for i in range(num_curves):
        lc = axs[0].plot(x, y_cloudy[i], label='cloudy ' +
                      lab_curves[i], marker='x')
        axs[0].plot(x, y_gasmod[i], label='gasmod ' +
                 lab_curves[i], marker='+', color=lc[0].get_color(), ls='dashed')

        axs[1].plot(x, y_gasmod[i] / y_cloudy[i], color=lc[0].get_color())

    for ax in axs:
        ax.set_xscale('log')
        ax.set_yscale('log')

    axs[1].set_xlabel(par_xlabel)
    axs[0].set_ylabel('density (cm-3)')
    axs[1].set_ylabel('ratio (gasmodule / cloudy)')
    fig.legend(*axs[0].get_legend_handles_labels())
    fig.subplots_adjust(hspace=0.03)
    fig.savefig('curves.pdf', bbox_inches='tight')

    # one curve per dir
    # _________________

    # spectrum
    plt.figure()
    for d in dirs:
        br = BenchmarkResult(d)
        par = br.get_nh_tc_lum()[par_index]

        x, y = br.cloudy.get_emissivity()
        lc = plt.loglog(x, y, label='cloudy {}'.format(par))

        x, y = br.gasmodule.get_emissivity()
        plt.loglog(x, y, label='gasmodule',
                   ls='dashed', color=lc[0].get_color())

    plt.xlabel('$\\lambda$ (micron)')
    plt.ylabel('outward emission')
    plt.gcf().legend()
    plt.savefig('emission.pdf')
    

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

    plt.savefig('populations.pdf')


if __name__ == '__main__':
    main()
