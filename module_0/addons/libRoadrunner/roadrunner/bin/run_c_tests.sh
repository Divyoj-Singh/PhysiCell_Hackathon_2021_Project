#!/bin/sh

BITNESS=64
OS=win
ARCHITECTURE=amd64
COMPILER_VERSION=
COMPILER_ID=MSVC
BUILD_TYPE=Release

usage () { 
	echo "run_c_tests \'sbml_test_suite_directory\' output_dir"
}

if [ -z "$1" ]
then
    echo "No argument supplied"
	usage
	exit 0
fi

if [ ! -e $1 ]
then
	echo "error, the sbml test suite location $1 does not exist"
	usage
	exit 0
fi

if [ ! -e $2 ]
then
	echo "error, the output dir $2 does not exist"
	usage
	exit 0
fi


if [ ! -e "$1/cases" ]
then
	echo "error, the sbml test suite location $1 does not appear valid, is missing a \'cases\' subdirectory"
	usage
	exit 0
fi

if [ ! -e "$1/cases" ]
then
	echo "error, the sbml test suite location $1 does not appear valid, is missing a \'cases\' subdirectory"
	usage
	exit 0
fi

# grab the first line of a git status
GIT_STATUS=`git --git-dir=$2/.git --work-tree=$2 status | head -1`

# pick out the branch we're on
GIT_BRANCH=`echo ${GIT_STATUS} | sed -E 's/.*branch[ \t]+(.*)/\1/'`

# format the datetime as part of the filename
FILE_DATETIME=`date +%Y-%m-%d"_"%H-%M`

echo "git output: $GIT_STATUS"
echo "git branch: $GIT_BRANCH"


# could be on a new architecture, so make the dir if it does not exist
mkdir -p "$2/reports/${OS}_${ARCHITECTURE}"

REPORT_FILE="$2/reports/${OS}_${ARCHITECTURE}/c_tests_${COMPILER_ID}_${COMPILER_VERSION}_${BUILD_TYPE}_${GIT_BRANCH}_${FILE_DATETIME}.xml"

# set up the compiler, 
# if LLVM is enabled, LLVM is the default compiler
#if [ "$OS" != "win" ]
#then
#	COMPILER="gcc"
#else
#	COMPILER="tcc"
#fi

echo "run_test_suite -m $1 -r ${REPORT_FILE} -c ${COMPILER}"

# run the test
run_test_suite -m $1 -r ${REPORT_FILE} -c ${COMPILER}

echo "done, report file: ${REPORT_FILE}"

# display the summary, this is the 3'rd line of the report xml
sed -n '3p' ${REPORT_FILE}
