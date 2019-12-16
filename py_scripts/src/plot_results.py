import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse


def plot_data(datafile):
    """Plot data file with the 1st column x, other columns y.
    """
    df = pd.read_csv(datafile, sep=',', comment='#')
    print('columns: {}'.format(df.columns))
    print('number of points: {}'.format(df.shape[0]))

    x = df.iloc[:, 0]
    y = df.iloc[:, 1:]
    fig = plt.figure()
    ncols = y.shape[1]
    for colname in df.columns[1:]:
        plt.plot(x, y[colname], '.-', label=str(colname).rstrip().lstrip())
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

