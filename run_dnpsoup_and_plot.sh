#!/usr/bin/env bash

SPINFILE=examples/spinsys/DPNO_BzP.json
#SPINFILE=examples/spinsys/BDPA_spinsys.json
#SPINFILE=examples/spinsys/eH_spinsys.json
#SPINFILE=examples/spinsys/TOTAPOL.json
#SPINFILE=examples/spinsys/trityl.json

#PULSEFILE=examples/pulseseq/cw_pulse_short.json
#PULSEFILE=examples/pulseseq/cw_pulse.json
#PULSEFILE=examples/pulseseq/topdnp.json
PULSEFILE=examples/pulseseq/NOVEL.json

#HARDWAREFILE=examples/hardware/400MHz263GHz.json
HARDWAREFILE=examples/hardware/14pt7MHz9pt7GHz0p35T.json

OUTPUT=results.log
if [ -f $OUTPUT ]; then
  rm $OUTPUT
fi

#python py_scripts/exec_and_run_dnpsoup.py examples/BDPA_spinsys.json examples/cw_pulse.json examples/hardware.json results.log
#./build/dnpsoup_cli/dnpsoup_exec examples/spinsys/TOTAPOL.json examples/pulseseq/cw_pulse.json examples/hardware/hardware.json results.log
./build/dnpsoup_cli/dnpsoup_exec $SPINFILE $PULSEFILE $HARDWAREFILE $OUTPUT

if [ ! -f $OUTPUT ]; then
  echo "$OUTPUT not found!"
else
  python py_scripts/src/plot_results.py $OUTPUT
fi

