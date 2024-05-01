
#Set path to diffusion folder (this needs to be a folder containing two things: 1. the transform T1_Regto_DTI_fsl.txt 2. the folder data.bedpostx 
DiffFolder=

#Set directory for outputs
OutDir=

STN_Motor_Region_L=seeds_to_motor_T1_L_GM_dil_top25_percentile_thr.nii.gz
STN_Motor_Region_R=seeds_to_motor_T1_R_GM_dil_top25_percentile_thr.nii.gz

avoid_mask_L=avoid_L.nii.gz
avoid_mask_R=avoid_R.nii.gz

motor_cortex_L=motor_L_barb.nii.gz
motor_cortex_R=motor_R_barb.nii.gz


#Run these commands from the folder containing the above files (and/or change their definition to the full path)


#version including the avoid mask
probtrackx2_gpu -x $STN_Motor_Region_L -l --pd --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 100000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --xfm=$DiffFolder/T1_Regto_DTI_fsl.txt --forcedir --avoid=$avoid_mask_L --stop=$motor_cortex_L --opd -s $DiffFolder/data.bedpostX/merged -m $DiffFolder/data.bedpostX/nodif_brain_mask.nii.gz --dir=$OutDir/Left --targetmasks=$motor_cortex_L
probtrackx2_gpu -x $STN_Motor_Region_R -l --pd --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 100000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --xfm=$DiffFolder/T1_Regto_DTI_fsl.txt --forcedir --avoid=$avoid_mask_R --stop=$motor_cortex_R --opd -s $DiffFolder/data.bedpostX/merged -m $DiffFolder/data.bedpostX/nodif_brain_mask.nii.gz --dir=$OutDir/Right --targetmasks=$motor_cortex_R



#version with no avoid mask
# probtrackx2_gpu -x $STN_Motor_Region_L -l --pd --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 100000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --xfm=$DiffFolder/T1_Regto_DTI_fsl.txt --forcedir --stop=$motor_cortex_L --opd -s $DiffFolder/data.bedpostX/merged -m $DiffFolder/data.bedpostX/nodif_brain_mask.nii.gz --dir=$OutDir/Left --targetmasks=$motor_cortex_L
# probtrackx2_gpu -x $STN_Motor_Region_R -l --pd --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 100000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --xfm=$DiffFolder/T1_Regto_DTI_fsl.txt --forcedir --stop=$motor_cortex_R --opd -s $DiffFolder/data.bedpostX/merged -m $DiffFolder/data.bedpostX/nodif_brain_mask.nii.gz --dir=$OutDir/Right --targetmasks=$motor_cortex_R

