#!/bin/sh

src=/home/camp/alexany/working/alexany/SMLM_TEST_DATA/minimal_SMLM_image_set/source
dst=/home/camp/alexany/working/alexany/SMLM_TEST_DATA/minimal_SMLM_image_set/output
settings_file=/home/camp/alexany/working/alexany/SMLM_TEST_DATA/minimal_SMLM_image_set/ALYtools_data_settings.xml
src_filename_must_have_token=tif

ARGS="$src $dst $settings_file $src_filename_must_have_token" 

./smlm_ALYtools_LAUNCH.sh $ARGS
