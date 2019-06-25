import argparse
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='+')
    args = ap.parse_args()

    # one point per dir
    # x = range(len(dirs))
    # y_cloudy = np.zeros(len(dirs))
    # y_gasmod = np.zeros(len(dirs))
    # for i, d in enumerate(args.dirs):
        
    # one curve per dir
    plt.figure()
    for d in args.dirs:
        pd = Path(d)
        cloudy_diffuse = np.loadtxt(pd / 'all_4pi_nu_jnu.erg_cm-3_s-1.out', max_rows=2)
        x = cloudy_diffuse[0]
        y = cloudy_diffuse[1]

        plt.loglog(x, y, label='cloudy' + d)

        mrn_opt = np.loadtxt(pd / 'MRNdust' / 'opticalProperties.dat')
        plt.loglog(mrn_opt[:, 1], mrn_opt[:, 2], label='gasmodule', ls='dashed')

    plt.xlabel('$\\lambda$ (micron)')
    plt.ylabel('outward emission')

    plt.gcf().legend()
    plt.show()

if __name__ == '__main__':
    main()
