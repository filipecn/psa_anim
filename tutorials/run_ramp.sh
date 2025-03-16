#!/bin/bash
# This script runs the simple slope tutorial: an avalanche descending an incline
# with constant slope with the following scheme
#                                 \
#     z                    \        \  1km
#    |    y                  \        \
#    |  /                      \      / 100m
#    ------x                30 ( \  /
# where
#   - slope angle is 30 degrees
#   - the slope extension is 1km
#   - the slope wide (y direction) is 100m

ramp_w=1000
ramp_d=150
ramp_h=60
ramp_a=25

while [ "$1" != "" ]; do
  case $1 in
  --length)
    shift
    ramp_w=$1
    ;;
  --width)
    shift
    ramp_d=$1
    ;;
  *)
    shift
    ;;
  esac
  shift
done

{
  echo "Running Ramp Simulation Case"
  echo "================================================================="
  echo "Simulation Parameters:"
  echo "  DOMAIN_WIDTH:        $ramp_d"
  echo "  DOMAIN_LENGTH:       $ramp_w"
  echo "  DOMAIN_HEIGHT:       $ramp_h"
  echo "  DOMAIN_ANGLE:        $ramp_a"
  echo "================================================================="
} | tee "{$LOG_FILE}_ramp"

################################################################################
#	RUN DSL SIMULATION FUNCTION                                                  #
#	  $1 - simulation path                                                       #
################################################################################
run_dsl() {
  {
    if [ $USE_DSL_SIM -eq 0 ]; then
 
      setup_dsl blockMesh

      # setup simulation geometry (terrain)
      $PYTHON "$TOOLS_DIR"/ramp_block_mesh_gen.py \
        -o "$DSL_DIR"/system                      \
        -w "$ramp_w"                              \
        --height 1                                \
        -d "$ramp_d"                              \
        -a $ramp_a                                \
        -fb walls                                 \
        --fb-type $fb_type                        \
        --b-type $b_type                          \
        --t-type $t_type                          \
        --l-type $l_type                          \
        --r-type $r_type

      # run simulation by calling ./Allrun-parallel in the sim folder
      run_sim "$DSL_DIR"
    else
      echo "using DSL simulation data from $DSL_DIR"
    fi
  } 2>&1 | tee "{$DSL_LOG_FILE}_ramp"
}

################################################################################
#	SETUP 2D PSL FUNCTION                                                        #
################################################################################
setup_2d() {

  setup_psl_2d blockMesh
 
  # setup simulation geometry (terrain)
  $PYTHON "$TOOLS_DIR"/ramp_block_mesh_gen.py \
    -o "$PSL_DIR"/system                      \
    -w "$ramp_w"                              \
    --height $ramp_h                          \
    -d "$ramp_d"                              \
    -a $ramp_a                                \
    -fb walls                                 \
    -b terrain                                \
    -cx $cell_size                            \
    --b-type $b_type                          \
    --t-type $t_type                          \
    --l-type $l_type                          \
    --r-type $r_type  --D2
  }
################################################################################
#	SETUP 3D PSL FUNCTION                                                        #
################################################################################
setup_3d() {
  
  setup_psl_3d blockMesh

  # setup simulation geometry (terrain)
  $PYTHON "$TOOLS_DIR"/ramp_block_mesh_gen.py \
    -o "$PSL_DIR"/system                      \
    -w "$ramp_w"                              \
    --height $ramp_h                          \
    -d "$ramp_d"                              \
    -a $ramp_a                                \
    -fb walls                                 \
    -b terrain                                \
    -cx $cell_size                            \
    --fb-type $fb_type                        \
    --b-type $b_type                          \
    --t-type $t_type                          \
    --l-type $l_type                          \
    --r-type $r_type
}

##################################################################
#	RUN PSL RELEASE FUNCTION                                       #
# uses all globals                                               #
##################################################################
run_release_2d() {
  # log
  echo "Running PSL release with parameters:"
  echo "  mesh type: $PSL_MESH_TYPE"
  echo "  cell size: $cell_size"
  echo "  height: $height"
  echo "  sim dir: $PSL_DIR"
  # init openfoam project directory
  $PYTHON "$TOOLS_DIR"/init_psl.py \
    -o "$PSL_DIR" \
    --template-dir "$TOOLS_DIR"/assets/psl \
    --mesh-generator "$PSL_MESH_TYPE" \
    --write-format $write_format                \
    --snow-density "$snow_density"  \
    --domain-type $boundary_condition \
    --dt $delta_t \
    --wt "$write_time" \
    --et "$duration" \
    --set-fields \
    --run-2d
  # setup simulation geometry (terrain)
  $PYTHON "$TOOLS_DIR"/ramp_block_mesh_gen.py \
    -o "$PSL_DIR"/system                      \
    -w "$ramp_w"                              \
    --height $ramp_h                          \
    -d "$ramp_d"                              \
    -a $ramp_a                                \
    -fb walls                                 \
    -b terrain                                \
    -cx $cell_size                            \
    --b-type $b_type                          \
    --t-type $t_type                          \
    --l-type $l_type                          \
    --r-type $r_type  --D2
  cd $PSL_DIR || exit
  blockMesh
  cd $WOKING_DIR || exit
  # gen empty dsl
  $PYTHON "$TOOLS_DIR"/gen_dsl.py               \
    -i "$PSL_DIR"                               \
    -o "$PSL_DIR"/constant/boundaryData         \
    -p "terrain"                                \
    --dt "$write_time"                              \
    --duration "$duration"                      
  # setup release volume
  $PYTHON "$TOOLS_DIR"/set_field.py \
    -o "$PSL_DIR"/system \
    -b 1 0 -60 15 80 60 185
  # run
  run_sim "$PSL_DIR" 
}
###########################################################
#	RUN DSL + PSL FUNCTION                                  #
# uses all globals                                        #
###########################################################
run_full() {
  # log
  {
    echo "Running PSL full sim with parameters:"
    echo "  mesh type: $PSL_MESH_TYPE"
    echo "  cell size: $cell_size"
    echo "  height: $height"
    echo "  sim dir: $PSL_DIR"
  } | tee "$PSL_LOG_FILE"

  run_dsl

  if [ $ONLY_DSL -eq 1 ]; then
    exit
  fi

  {

  if [ $RUN_2D -eq 1 ]; then
    setup_2d
  else
    setup_3d
  fi

  cd $PSL_DIR || exit
  blockMesh
  cd $WOKING_DIR || exit

  dsl_to_psl

  run_sim "$PSL_DIR" --no-clean

  } 2>&1 | tee "{$PSL_LOG_FILE}_ramp"
}

################################################################################
#  RUN                                                                         #
################################################################################
if [ $RELEASE -eq 1 ]; then
  echo "initiating release psl sim..."
  if [ $RUN_2D -eq 1 ]; then
    run_release_2d 2>&1 | tee "$PSL_LOG_FILE"
  else
    setup_3d
  fi
else
  echo "initiating full psl sim..."
  run_full
fi
