#!/usr/bin/env bash

export OMP_NUM_THREADS=1

# =====================================
#   FILENAME | CATECORYNAME | VERSION
# =====================================

#CATEGORYNAME=Top-DNP
#CATEGORYNAME=ISE_SSE
#CATEGORYNAME=CW_CrossEffect/v3
#CATEGORYNAME=NOVEL
CATEGORYNAME=NOVEL-MAS
VERSION=

#===========================================================


APPENDIX=""

FILENAME=eH_NOVEL_xband_loop60_static_inc40ns_zcw144_fp
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""


OUTPUTDIR=dnpsoup_analytics/outputs
OUTPUTPATH=$OUTPUTDIR/$CATEGORYNAME/$FILENAME.result$VERSION$APPENDIX
echo "./plot_results.py $OUTPUTPATH"
./plot_results.py $OUTPUTPATH

