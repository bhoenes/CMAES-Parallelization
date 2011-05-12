#!/bin/bash

usage() {
  echo "usage: $0 [-p NUM_PROCESSORS] [-l lamhosts] [-e executable] Optional [-i iterations] [-v verbose] [-d num dipoles] [-m num samples]"
  exit 1
}

while getopts "p:l:d:m:e:i:v" OPTIONS
  do
    case "$OPTIONS" in
      p) NUM_PROCESSORS=$OPTARG;;
      i) ITER=$OPTARG;;
      d) NUM_DIPOLES=$OPTARG;;
      e) EXEC=$OPTARG;;
      m) NUM_SAMPLES=$OPTARG;;
      v) VERBOSE="-v";;
      l) LAMHOSTS=$OPTARG;;
	  \?) usage
         exit 1;;
    esac
  done

lamboot $VERBOSE $LAMHOSTS
#mpiexec -n 2 -machinefile ../tinyhosts cycleprunner
SUM=0
for i in { 0 .. $ITER }
do
	mpiexec -n $NUM_PROCESSORS -machinefile $LAMHOSTS $EXEC $NUM_DIPOLES $NUM_SAMPLES
done	
lamclean $VERBOSE 
lamhalt 
lamwipe $VERBOSE $LAMHOSTS
