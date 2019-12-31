from src import run_dnpsoup, plot_data
import argparse


def main():
    parser = argparse.ArgumentParser(description="Run dnpsoup and plot results.")
    parser.add_argument('spinsys', help='js file for spin system')
    parser.add_argument('pulseseq', help='js file for pulse sequence')
    parser.add_argument('hardware', help='js file for hardware')
    parser.add_argument('output', help='output filename')
    args = parser.parse_args()
    run_dnpsoup(args.spinsys, args.pulseseq, args.hardware, args.output)
    plot_data(args.output)

if __name__ == '__main__':
    main()

