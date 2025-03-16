#!/bin/bash
# This script runs a series of simulations for method analysis.
# The same scene is used under many different set of configurations.

if [ -z "$PSA_ANIM_TOOLS" ]; then
  echo "error: Nadare enviroment variables (PSA_ANIM_*) not set. Please run"
  echo "       'source PATH/TO/PSA_ANIM/BUILD/nadare_variables.sh'"
  exit
fi

OUTPUT_DIR=$(pwd)/vdb_tests
mkdir -p $OUTPUT_DIR

resolutions=(10)
skew_angles=(90)

id=0
for angle in "${skew_angles[@]}"; do
  for res in "${resolutions[@]}"; do
    filename="${angle}_${res}"
    $PSA_ANIM_TOOLS/foam2vdb -i "" -o $OUTPUT_DIR/${filename}.vdb --test \
      -s "0.2" \
      --method $method \
      --test-angle $angle \
      --test-res $res
    ((id=id+1))
  done
done
