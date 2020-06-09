#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

UPPER_LIMIT = 6000


def plot_data(datafile):
    """Plot data file with the 1st column x, other columns y.
    """
    df = pd.read_csv(datafile, sep=',', comment='#', header=None)
    if df.shape[0] == 0:
        with open(datafile, 'r', encoding='utf-8') as f:
            for line in f:
                print(line)
        return

    print('columns:')
    for col in df.columns:
        print(str(col), end=' ')
    print()
    print('number of points: {}'.format(df.shape[0]))

    x = df.iloc[:, 0]
    ys = df.iloc[:, 1:]
    ncols = ys.shape[1]
    fig, axes = plt.subplots(ncols, 1, sharex=True)
    if x.shape[0] > UPPER_LIMIT:
        step_size = x.shape[0]//UPPER_LIMIT
    else:
        step_size = 1
    x = x[::step_size]
    ys = ys.iloc[::step_size, :]

    if ncols == 1:
        axes.plot(x, ys[df.columns[1]], 'b.-', label=str(df.columns[1]).strip())
    else:
        for i, colname in enumerate(df.columns[1:]):
            y = ys[colname]
            axes[i].plot(x, y, '.-', label=str(colname).strip())
            ymin = np.min(y)
            ymax = np.max(y)
            if np.abs(ymax - ymin) > 1.0e-14:
                y_range = ymax - ymin
                y_range_margin = 0.01 * y_range
                axes[i].set_ylim(ymin - y_range_margin, ymax + y_range_margin)

    plt.subplots_adjust(hspace=0)
    plt.legend(loc='best')
    plt.show()


def main():
    description = 'plot data generated from dnpsoup'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('datafile', help='name of the datafile')
    args = parser.parse_args()
    plot_data(args.datafile)


if __name__ == '__main__':
    main()

