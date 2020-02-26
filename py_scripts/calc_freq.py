#!/usr/bin/env python
import argparse


def calc_freq(field, spin_type):
    gyro = {
        'H': 42.57747872e6,
        'e': 28024951108.281105
    }[spin_type]
    return field * gyro


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('field')
    parser.add_argument('spin')
    args = parser.parse_args()
    f = float(args.field)
    t = str(args.spin)
    freq = calc_freq(f, t)
    if freq > 0.1e9:
        print('{} GHz'.format(freq/1.0e9))
    else:
        print('{} MHz'.format(freq/1.0e6))



if __name__ == '__main__':
    main()

