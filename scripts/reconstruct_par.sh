#!/bin/sh

FRAMES_OUTPUT_DIR=""
BATCH=0
START=0
END=0
NP=20
STEP=0
CLEAN=0

usage() {
  echo "usage: reconstruct_par [[-i sim-path] | [-h]]"
}

while [ "$1" != "" ]; do
  case $1 in
  -o | --output)
    shift
    FRAMES_OUTPUT_DIR=$(realpath "$1")
    ;;
  -s)
    shift
    START=$1
    ;;
  -e)
    shift
    END=$1
    ;;
  --step)
    shift
    STEP=$1
    ;;
  --clean)
    CLEAN=1
    ;;
  -np)
    shift
    NP=$1
    ;;
  --batch)
    BATCH=1
    ;;
  -h | --help)
    usage
    exit
    ;;
  *)
    usage
    exit 1
    ;;
  esac
  shift
done

if [[ ! -d $FRAMES_OUTPUT_DIR ]]; then
  mkdir -p $FRAMES_OUTPUT_DIR
fi

if [[ $BATCH -eq 1 ]]; then
  param="$START:$END"
  mpirun -np $NP redistributePar -reconstruct -parallel -time ':100'
fi

for t in `seq $START $STEP $END`; do 
  time=$(echo "$t" | sed -r -e 's/\.([0-9]*)0$/.\1/')
  time=$(echo "$time" | sed -r -e 's/\.0*$//')
  if [[ -z $time ]]; then
    continue
  fi
  if [[ $BATCH -eq 0 ]]; then
    if [[ -d processor0/$time ]]; then
      echo "mpirun -np $NP redistributePar -reconstruct -parallel -time $time"
      mpirun -np $NP redistributePar -reconstruct -parallel -time $time
    fi
  fi
  if [[ -d $time ]]; then
    if [[ $CLEAN -eq 1 ]]; then
      for i in $( seq 0 $NP ); do 
        rm -r processor$i/$time
      done
    fi
    mv $time $FRAMES_OUTPUT_DIR/$time
  fi
done

