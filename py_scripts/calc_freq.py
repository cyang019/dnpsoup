#!/usr/bin/env python
import argparse


def calc_freq(field, spin_type, g):
    mu_b = 9.274009994e-24
    hbar = 1.0545718e-34
    gyro = {
        'H': 42.57747872e6,
        'e': mu_b/hbar * g / (2.0 * 3.1415927)
    }[spin_type]
    return field * gyro


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('field')
    parser.add_argument('spin')
    parser.add_argument('g')
    args = parser.parse_args()
    f = float(args.field)
    t = str(args.spin)
    g = float(args.g)
    freq = calc_freq(f, t, g)
    if freq > 0.1e9:
        print('{} GHz'.format(freq/1.0e9))
    else:
        print('{} MHz'.format(freq/1.0e6))



if __name__ == '__main__':
    main()

