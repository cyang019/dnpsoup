#!/usr/bin/env bash

FILENAME=$1
CATEGORYNAME=$2
VERSION=$3

INPUTPATH=dnpsoup_analytics/inputs/$CATEGORYNAME/$FILENAME.json
APPENDIX=""

#OUTPUTPATH=results/v2/$CATEGORYNAME/$FILENAME.result$VERSION$APPENDIX

OUTPUTDIR=dnpsoup_analytics/outputs
OUTPUTPATH=$OUTPUTDIR/$CATEGORYNAME/$FILENAME.result$VERSION$APPENDIX

date
echo input: $INPUTPATH
./build/dnpsoup_cli/dnpsoup_exec $INPUTPATH $OUTPUTPATH
echo input: $INPUTPATH
echo output: $OUTPUTPATH
echo "finished."
date
echo ""

