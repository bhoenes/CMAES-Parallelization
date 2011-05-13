#!/bin/bash

usage() {
  echo "usage: $0 [-p NUM_PROCESSORS] [-l lamhosts] [-e executable] Optional [-i iterations] [-v verbose] [-d num dipoles] [-m num samples] [-t num threads] [-c num_chunks]"
  exit 1
}

ITER=1

while getopts "p:l:d:m:e:i:t:c:v" OPTIONS
  do
    case "$OPTIONS" in
      p) NUM_PROCESSORS=$OPTARG;;
      i) ITER=$OPTARG;;
      d) NUM_DIPOLES=$OPTARG;;
      e) EXEC=$OPTARG;;
      m) NUM_SAMPLES=$OPTARG;;
      v) VERBOSE="-v";;
      t) NUM_THREADS=$OPTARG;;
      c) NUM_CHUNKS=$OPTARG;;
      l) LAMHOSTS=$OPTARG;;
	  \?) usage
         exit 1;;
    esac
  done

lamboot $VERBOSE $LAMHOSTS

COUNTER=0
while [  $COUNTER -lt $ITER ]; do
	mpiexec -n $NUM_PROCESSORS -machinefile $LAMHOSTS $EXEC $NUM_DIPOLES $NUM_SAMPLES $NUM_THREADS $NUM_CHUNKS
	let COUNTER=COUNTER+1 
done	
lamclean $VERBOSE 
lamhalt 
lamwipe $VERBOSE $LAMHOSTS
