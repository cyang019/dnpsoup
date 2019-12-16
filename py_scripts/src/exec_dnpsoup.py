import subprocess
import argparse

dnpsoup_path = './build/dnpsoup_cli/dnpsoup_exec'


def run_dnpsoup(spinsys_js_file, pulseseq_js_file, hardware_js_file, result_file):
    subprocess.run([dnpsoup_path,
                    spinsys_js_file,
                    pulseseq_js_file,
                    hardware_js_file,
                    result_file], capture_output=True)

def main():
    descrption = (
        "run dnpsoup with a spinsys jason file, "
        "a pulse sequence jason file, a hardwaer jason file,"
        "and output filename")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('spinsys', help='js file for spin system')
    parser.add_argument('pulseseq', help='js file for pulse sequence')
    parser.add_argument('hardware', help='js file for hardware')
    parser.add_argument('output', help='output filename')
    args = parser.parse_args()
    run_dnpsoup(args.spinsys, args.pulseseq, args.hardware, args.output)

if __name__ == '__main__':
    main()

