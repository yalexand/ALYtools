tips for CentOS7, SLURM, Matlab 2018b

"run_AI_Powered_2D_SMLM_reconstruct_CentOS7.sh" is the file generated by Matlab's "deploytool"

the corresponding compiled linux app is named "AI_Powered_2D_SMLM_reconstruct_CentOS7" (not included)

the corresponding main m-file for "deploytool" is "AI_Powered_2d_SMLM_reconstruct.m" (in ICUtilities)

m-file gets 4 parameters - all paths are full
1) source image - file, tiff or ome.tiff 
2) destination directory
3) setting file name
4) src_filename_must_have_token - hopefully can help resolving multiple-filed ome.tiff issue

output of the procedure is a directory within desination directory, 
named as a source file name (without extension), 
containing SMLM X,Y,sigma outputs, and super-resolved image.