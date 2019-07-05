import argparse
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path

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
    plt.ylim(bottom=1e-30, top=max(y) * 1.1)
    plt.xlim(left=1e-2, right=3e5)
    plt.gcf().legend()
    plt.savefig('emission.pdf')
    

    # populations
#    plt.figure()
    fig, axs = plt.subplots(2, 2, sharex='col')
    for i, d in enumerate(dirs):
        br = BenchmarkResult(d)
        par = br.get_nh_tc_lum()[par_index]

        yc = br.cloudy.get_h_populations()
        lc = axs[0][0].semilogy(range(len(yc)), yc, label='cloudy {}'.format(par), marker='x')
        color = lc[0].get_color()

        yg = br.gasmodule.get_h_populations()
        axs[0][0].semilogy(range(len(yg)), yg, label='gasmodule', marker='+', color=color, ls='dashed')

        max_h_level = min(len(yc), len(yg))
        axs[1][0].semilogy(range(max_h_level), yg[:max_h_level] / yc[:max_h_level], color=color)

        y2c = br.cloudy.get_h2_populations()
        lc = axs[0][1].semilogy(range(len(y2c)), y2c, color=color)

        y2g = br.gasmodule.get_h2_populations()
        axs[0][1].semilogy(range(len(y2g)), y2g, color=color, ls='dashed')

        max_h2_level = min(len(y2c), len(y2g))
        axs[1][1].semilogy(range(max_h2_level), y2g[:max_h2_level] / y2c[:max_h2_level], color=color)

    fig.legend(*axs[0][0].get_legend_handles_labels())
    axs[1][0].set_xlabel('h level index')
    axs[1][1].set_xlabel('h2 level index')
    axs[0][0].set_ylabel('population fraction')
    axs[0][0].set_title('H')
    axs[0][1].set_title('H2')

    plt.savefig('populations.pdf')


if __name__ == '__main__':
    main()
