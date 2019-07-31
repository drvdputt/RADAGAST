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

    def prepare_curve_arrays(num_curves):
        curves_array = np.zeros((num_curves, len(dirs)))
        return curves_array, curves_array.copy()

    def plot_curves(labels, data_cloudy, data_gasmod, ylabel='value', fn='curves.pdf'):
        fig, axs = plt.subplots(2, 1, sharex=True)
        for i in range(len(labels)):
            lc = axs[0].plot(x, data_cloudy[i],
                             label='cloudy ' + labels[i], marker='d', ls='none')
            axs[0].plot(x, data_gasmod[i], label='gasmod ' +
                        labels[i], color=lc[0].get_color())
            axs[1].plot(x, data_gasmod[i] / data_cloudy[i],
                        color=lc[0].get_color())

        for ax in axs:
            ax.set_xscale('log')
            ax.set_yscale('log')
        axs[0].set_ylabel(ylabel)
        axs[1].set_xlabel(par_xlabel)
        axs[1].set_ylabel('ratio (gasmod / cloudy)')

        fig.legend(*axs[0].get_legend_handles_labels())
        fig.subplots_adjust(hspace=0.03)
        fig.savefig(fn, bbox_inches='tight')

    dens_curves = ['e', 'H+', 'H', 'H2']
    dens_cloudy, dens_gasmod = prepare_curve_arrays(len(dens_curves))

    rate_curves = ['h2form', 'h2dissoc']
    rate_cloudy, rate_gasmod = prepare_curve_arrays(len(rate_curves))

    # one point per dir
    # -----------------
    dens_fig, dens_axs = plt.subplots(2, 1, sharex=True)
    rate_fig, rate_axs = plt.subplots(2, 1, sharex=True)
    for i, d in enumerate(dirs):
        br = BenchmarkResult(d)
        x[i] = br.get_nh_tc_lum()[par_index]
        dens_cloudy[:, i] = br.cloudy.get_densities(numpy=True)
        dens_gasmod[:, i] = br.gasmodule.get_densities(numpy=True)
        rate_cloudy[:, i] = br.cloudy.get_h2_rates()
        rate_gasmod[:, i] = br.gasmodule.get_h2_rates()

    plot_curves(dens_curves, dens_cloudy, dens_gasmod,
                'density (cm-3)', 'density.pdf')
    plot_curves(rate_curves, rate_cloudy,
                rate_gasmod, 'rate (s-1)', 'rate.pdf')

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
        lc = axs[0][0].semilogy(
            range(len(yc)), yc, label='cloudy {}'.format(par), ls='none', marker='d')
        color = lc[0].get_color()

        yg = br.gasmodule.get_h_populations()
        axs[0][0].semilogy(range(len(yg)), yg, label='gasmodule',
                           color=color)

        max_h_level = min(len(yc), len(yg))
        axs[1][0].semilogy(range(max_h_level),
                           yg[:max_h_level] / yc[:max_h_level], color=color)

        y2c = br.cloudy.get_h2_populations()
        lc = axs[0][1].semilogy(range(len(y2c)), y2c, color=color)

        y2g = br.gasmodule.get_h2_populations()
        axs[0][1].semilogy(range(len(y2g)), y2g, color=color, ls='dashed')

        max_h2_level = min(len(y2c), len(y2g))
        axs[1][1].semilogy(range(max_h2_level),
                           y2g[:max_h2_level] / y2c[:max_h2_level], color=color)

    fig.legend(*axs[0][0].get_legend_handles_labels())
    axs[1][0].set_xlabel('h level index')
    axs[1][1].set_xlabel('h2 level index')
    axs[0][0].set_ylabel('population fraction')
    axs[0][0].set_title('H')
    axs[0][1].set_title('H2')

    plt.savefig('populations.pdf')


if __name__ == '__main__':
    main()
