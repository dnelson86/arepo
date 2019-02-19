#!/bin/bash            # this line only there to enable syntax highlighting in this file
## shell script for code testing

#################################
## perform predefined examples ##
#################################

## Number of cores to compile and run the problem on
## use NUMBER_OF_TASKS=1 for 1d test problems!
NUMBER_OF_TASKS=1

## choose your tests
TESTS=""
## available 1d test cases
TESTS+="wave_1d "
TESTS+="shocktube_1d "

## loop over all tests
for TEST in $TESTS
do
  DIR="./examples/"${TEST}"/"
  RUNDIR="./run/examples/"${TEST}"/"
  
  ## clean up 
  rm -rf ./run

  ## create run directory
  mkdir ./run
  mkdir ./run/examples
  mkdir ${RUNDIR}
  
  ## create ICs in run directory
  echo ${DIR}
  python ${DIR}/create.py ${RUNDIR}
  ((return_value=$?))    ## get return value
  if [ $return_value != 0 ]    ## check return value
  then echo "ERROR: test.sh:\t" $DIR "\t python create.py failed!"
  exit $return_value
  fi
  
  ## copy Config and parameter file to run directory
  cp ${DIR}/Config.sh ${RUNDIR}
  cp ${DIR}/param.txt ${RUNDIR}

  ## compile Arepo
  make -j ${NUMBER_OF_TASKS} CONFIG=${RUNDIR}/Config.sh BUILD_DIR=${RUNDIR}/build EXEC=${RUNDIR}/Arepo
  ((return_value=$?))    ## get return value
  if [ $return_value != 0 ]    ## check return value
  then echo "ERROR: test.sh:\t" $DIR "\t make failed!"
  exit $return_value
  fi
  
  ## execute test simulation
  mpiexec -np ${NUMBER_OF_TASKS} ${RUNDIR}/Arepo ${RUNDIR}/param.txt
  ((return_value=$?))    ## get return value
  if [ $return_value != 0 ]    ## check return value
  then echo "ERROR: test.sh:\t" $DIR "\t execution failed!"
  exit $return_value
  fi
  
  ## check result in example directory, this also creates some check plots
  python ${DIR}/check.py ${RUNDIR}
  ((return_value=$?))    ## get return value
  if [ $return_value != 0 ]    ## check return value
  then echo "ERROR: test.sh: test failed!"
  exit $return_value
  else echo "test.sh:\t" $DIR "\t test passed!"
  echo "cleaning up..."
  fi
  
  ## clean up
  rm -rf ./run
done

echo "tests"
echo ${TESTS}
echo "passed!"

exit ${return_value}
