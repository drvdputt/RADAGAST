import argparse
from pathlib import Path
import numpy as np

from read_benchmark_output import BenchmarkResult

COLUMNS = ['nH', 'Tc', 'lum', 'ne_reldiff', 'nH+_reldiff', 'nH_reldiff', 'nH2_reldiff', 'T_reldiff']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='+')
    args = ap.parse_args()

    dirs = [Path(d) for d in args.dirs]

    with open('benchmark_overview.txt', 'w') as f:
        write_header(f)
        for d in dirs:
            write_bench_overview(f, d)


def write_header(file_handle):
    line = '\t'.join(COLUMNS)
    file_handle.write('#' + line + '\n')


def write_bench_overview(file_handle, bench_dir):
    br = BenchmarkResult(bench_dir)

    line_data = np.zeros(len(COLUMNS))
    pos = 0

    # benchmark input parameters
    par = br.get_nh_tc_lum()
    line_data[pos:pos+len(par)] = par
    pos += len(par)

    # density reldiffs
    dens_cloudy = br.cloudy.get_densities(numpy=True)
    dens_gasmod = br.gasmodule.get_densities(numpy=True)
    dens_reldiffs = dens_gasmod / dens_cloudy - 1
    line_data[pos:pos+len(dens_reldiffs)] = dens_reldiffs
    pos += len(dens_reldiffs)

    T_reldiff = br.gasmodule.get_temp() / br.cloudy.get_temp() - 1
    line_data[pos] = T_reldiff
    pos += 1

    file_handle.write('\t'.join(('{:.2e}'.format(d) for d in line_data)) + '\n')


if __name__ == '__main__':
    main()
