#!/usr/bin/env python
import argparse


def calc_field(freq, spin_type):
    gyro = {
        'H': 42.57747872e6,
        'e': 28024951108.281105
    }[spin_type]
    return freq / gyro


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('freq')
    parser.add_argument('spin')
    args = parser.parse_args()
    f = float(args.freq)
    t = str(args.spin)
    print(calc_field(f, t))


if __name__ == '__main__':
    main()
