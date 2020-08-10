#!/usr/bin/env bash

FILENAME=$1
CATEGORYNAME=$2
VERSION=$3

INPUTPATH=examples/$CATEGORYNAME/$FILENAME.json
APPENDIX=""

#OUTPUTPATH=results/v2/$CATEGORYNAME/$FILENAME.result$VERSION$APPENDIX

OUTPUTDIR=py_scripts/case_studies/dnpsoup_analytics/results
OUTPUTPATH=$OUTPUTDIR/$CATEGORYNAME/$FILENAME.result$VERSION$APPENDIX

date
echo input: $INPUTPATH
./fast_exec_dnpsoup $INPUTPATH $OUTPUTPATH
echo input: $INPUTPATH
echo output: $OUTPUTPATH
echo "finished."
date
echo ""

