#!/usr/bin/env bash

# =====================================
#   FILENAME | CATECORYNAME | VERSION
# =====================================

# SE Powder Line Shape
#FILENAME=eeH_e_T1_300us_T2_1us_H_T1_4s_T2_200us_1hf_p3ms_1MHz_static_inc40ns_a0b0g0_9p385T_buildup
#FILENAME=eeH_e_T1_300us_T2_1us_H_T1_4s_T2_200us_1hf_p3ms_1MHz_static_inc40ns_a0b0g0_zcw233_9p36T_to_9p41T_fp

# cross effect
#FILENAME=eeH_1hf_p5ms_1MHz_inc40ns_a0b90g0_MAS8kHz_buildup
#FILENAME=eeH_1hf_p5ms_0MHz_inc40ns_an80bn141gn320_MAS8kHz_buildup
#FILENAME=eeH_1hf_p150us_1MHz_inc40ns_an80bn141gn320_MAS8kHz_9p394T_eigenvals
#FILENAME=eeH_1hf_p5ms_1MHz_inc20ns_a0b0g0_MAS8kHz_buildup
#FILENAME=eeH_1hf_p2ms_1MHz_inc40ns_a0b0g0_hemi21_9p36_to_9p41_step0p001_fp
#FILENAME=eeH_1hf_p150us_1MHz_inc40ns_a0b0g0_a0b90g0_MAS8kHz_scan1d_p
#FILENAME=eeH_1hf_p2ms_1MHz_inc40ns_a0b0g0_static_hemi144_9p36_to_9p41_step0p001_fp
#FILENAME=eeH_1hf_p2ms_1MHz_inc40ns_a0b0g0_mas8kHz_a0b90g0_9p36_to_9p41_step0p001_fp
#FILENAME=eeH_1hf_p1p6ms_1MHz_inc100ns_a0b0g0_mas8kHz_hemi34_9p36_to_9p41_step0p001_fp

#FILENAME=eeH_1hf_p2ms_1MHz_inc100ns_a0b120g0_mas8kHz_a0b0g0_9p36_to_9p41_step0p001_fp
#FILENAME=eeH_1hf_p2ms_1MHz_inc100ns_a0b135g0_mas8kHz_a0b0g0_9p36_to_9p41_step0p001_fp
#FILENAME=eeH_1hf_p500us_1MHz_inc40ns_a0b0g0_mas8kHz_a0b0g0_9p363T_intensity

# Trouble shooting around 9.385T no enhancement reasons
#FILENAME=eeH_1hf_p500us_1MHz_inc20ns_a0b45g0_9p385T_buildup
#FILENAME=eeH_e_T1_300us_T2_1us_H_T1_4s_T2_200us_1hf_p3ms_1MHz_inc40ns_a0b0g0_9p38T_zcw34_buildup
#FILENAME=eeH_1hf_p300us_1MHz_inc40ns_an80bn141gn320_MAS8kHz_9p394T_buildup
FILENAME=eeH_test_coord_e2_a0b50g40_H_xn2y2z1_a0b0g0_zcw34_9p36T_to_9p41T_fp
CATEGORYNAME=CW_CrossEffect/v3
VERSION=2

#FILENAME=eH_NOVEL_xband_loop100k_fp
#CATEGORYNAME=NOVEL
#===========================================================


APPENDIX=""

./take_one_input.sh $FILENAME $CATEGORYNAME $VERSION

OUTPUTPATH=results/v2/$CATEGORYNAME/$FILENAME.result$VERSION$APPENDIX
echo "./plot_results.py $OUTPUTPATH"
./plot_results.py $OUTPUTPATH

