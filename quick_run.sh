#!/usr/bin/env bash
export OMP_NUM_THREADS=1

# =====================================
#   FILENAME | CATECORYNAME | VERSION
# =====================================

#CATEGORYNAME=Top-DNP
#CATEGORYNAME=ISE_SSE
CATEGORYNAME=CE/xtal_buildups
#CATEGORYNAME=NOVEL
#CATEGORYNAME=NOVEL-MAS
VERSION=

#===========================================================

APPENDIX=""

FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b0g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b15g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b30g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b45g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b60g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b75g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b90g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b105g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b120g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b135g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b150g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b165g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""
FILENAME=eeH_test_coord_e2_a0b60g150_H_xn2y2z1_p1ms_a0b180g0_9p393T_inc40ns_buildup
./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION
echo ""

#OUTPUTDIR=dnpsoup_analytics/outputs
#OUTPUTPATH=$OUTPUTDIR/$CATEGORYNAME/$FILENAME.result$VERSION$APPENDIX
#echo "./plot_results.py $OUTPUTPATH"
#./plot_results.py $OUTPUTPATH

