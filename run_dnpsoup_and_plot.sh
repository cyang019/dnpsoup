#!/usr/bin/env bash

#python py_scripts/exec_and_run_dnpsoup.py examples/BDPA_spinsys.json examples/cw_pulse.json examples/hardware.json results.log
./build/dnpsoup_cli/dnpsoup_exec examples/spinsys/BDPA_spinsys.json examples/pulseseq/cw_pulse.json examples/hardware/hardware.json results.log
python py_scripts/src/plot_results.py results.log

