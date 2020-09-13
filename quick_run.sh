#!/usr/bin/env bash
export OMP_NUM_THREADS=1

# =====================================
#   FILENAME | CATECORYNAME | VERSION
# =====================================

#CATEGORYNAME=Top-DNP
#CATEGORYNAME=ISE_SSE
CATEGORYNAME=CE/Hu_JCP_2008/H_orientations/seed_files
#CATEGORYNAME=CE/buildup_scan1ds
#CATEGORYNAME=CE/buildup_powder_options
#CATEGORYNAME=CE/CE_SE_fp_comparison
#CATEGORYNAME=NOVEL
#CATEGORYNAME=NOVEL-MAS
VERSION=

#===========================================================

APPENDIX=""

#INPUTDIR=dnpsoup_analytics/inputs
#echo $CATEGORYNAME
#for filename in $INPUTDIR/$CATEGORYNAME/*.json; do
#  FULLNAME=$(basename $filename)
#  FILENAME=${FULLNAME%.json}
#  echo $FILENAME
#  ./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
#  echo ""
#done

FILENAME=eeH_ee7p4A_4groups_a0b0g0_zcw21_inc40ns_fp
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""

#OUTPUTDIR=dnpsoup_analytics/outputs
#OUTPUTPATH=$OUTPUTDIR/$CATEGORYNAME/$FILENAME.result$VERSION$APPENDIX
#echo "./plot_results.py $OUTPUTPATH"
#./plot_results.py $OUTPUTPATH

