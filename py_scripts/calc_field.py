#!/usr/bin/env python
import argparse


def calc_field(freq, spin_type, g):
    mu_b = 9.274009994e-24
    hbar = 1.0545718e-34
    gyro = {
        'H': 42.57747872e6,
        'e': mu_b/hbar * g / (2.0 * 3.1415927)
    }[spin_type]
    return freq / gyro


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('freq')
    parser.add_argument('spin')
    parser.add_argument('g')
    args = parser.parse_args()
    f = float(args.freq)
    t = str(args.spin)
    g = float(args.g)
    print(calc_field(f, t, g))


if __name__ == '__main__':
    main()
