#!/bin/bash

# Clean build directory
rm -rf ../build/*
cd ../build
cmake  ..
make

TAB_DEL_FILE=/Users/kot4or/cchmc/geep/testing_data/set_2/annotation_tab_del
BAM_FILE=/Users/kot4or/cchmc/geep/testing_data/set_2/reads.bam
TEST_RESULT=/Users/kot4or/cchmc/geep/testing_data/set_2/correct_results.txt
TEST_OUTPUT=/Users/kot4or/cchmc/geep/testing_data/set_2/output_results.txt

#TAB_DEL_FILE=`zenity --file-selection --title="Select tab-delimited file"`
#
#case $? in
#         0)
#                echo "\"TAB_DEL_FILE\" selected." to ${TAB_DEL_FILE};;
#         1)
#                echo "No file selected."
#                exit 1;;
#        -1)
#                echo "An unexpected error has occurred."
#                exit 1;;
#esac
#
#BAM_FILE=`zenity --file-selection --title="Select bam file"`
#
#case $? in
#         0)
#                echo "\"BAM_FILE\" selected." to ${BAM_FILE};;
#         1)
#                echo "No file selected."
#                exit 1;;
#        -1)
#                echo "An unexpected error has occurred."
#                exit 1;;
#esac
#
#
#TEST_RESULT=`zenity --file-selection --title="Select achieved result file"`
#
#case $? in
#         0)
#                echo "\"TEST_RESULT\" selected." to ${TEST_RESULT};;
#         1)
#                echo "No file selected."
#                exit 1;;
#        -1)
#                echo "An unexpected error has occurred."
#                exit 1;;
#esac

# Running with parameters
export DYLD_LIBRARY_PATH=/Users/kot4or/workspaces/geep_ws/geep/lib/
../bin/geep $BAM_FILE $TAB_DEL_FILE --test $TEST_OUTPUT > /dev/null 2>&1
diff --ignore-blank-lines --ignore-space-change ${TEST_RESULT} ${TEST_OUTPUT}
if [[ $? == "0" ]]
then
  echo "CORRECT"
else
  echo "ERROR"
fi


