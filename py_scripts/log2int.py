import argparse

def log2int(val):
    target = 0
    while val > 0:
            val //= 2
            target += 1
    return target

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('value')
    args = parser.parse_args()
    value = int(args.value)
    print(log2int(value))


if __name__ == '__main__':
    main()
