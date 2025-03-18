#!/bin/bash
# This script runs the wolfsgrube case found in RAUTTER 2018 for the DSL
# Available at https://develop.openfoam.com/Community/avalanche
# The case contains a large terrain of 3km described by a GIS file

ASSETS_DIR=$PSA_ANIM_TUTORIALS/assets/wolfsgrube
# options
DOMAIN_WIDTH=1000
DOMAIN_LENGTH=2500

# setup parameters
height=100
cell_size=3
axis_a_x=1755
axis_a_y=40
axis_b_x=600
axis_b_y=2040

while [ "$1" != "" ]; do
  case $1 in
  -W | --width)
    shift
    DOMAIN_WIDTH=$1
    ;;
  -L | --length)
    shift
    DOMAIN_LENGTH=$1
    ;;
  *)
    shift
    ;;
  esac
  shift
done

# compute domain width for the PSL simulation
# for that, lets compute the unit distance vector between point a and b of the avalanche axis
# and multiply by the domain width
# first we compute the original distance
axis_x=$(bc -l <<<"$axis_b_x - $axis_a_x;")
axis_y=$(bc -l <<<"$axis_b_y - $axis_a_y;")
# normalize axis
norm=$(bc -l <<<"sqrt($axis_x * $axis_x + $axis_y * $axis_y);")
axis_x=$(bc -l <<<"$axis_x / $norm;")
axis_y=$(bc -l <<<"$axis_y / $norm;")
# recompute axis_b
axis_b_x=$(bc -l <<<"$axis_a_x + $axis_x * $DOMAIN_LENGTH;")
axis_b_x=${axis_b_x%.*}
axis_b_y=$(bc -l <<<"$axis_a_y + $axis_y * $DOMAIN_LENGTH;")
axis_b_y=${axis_b_y%.*}

{
  echo "================================================================="
  echo "Running Wolfsgrube Simulation Case"
  echo "================================================================="
  echo "================================================================="
  echo "Generate simulation for wolfsgrube terrain with reference points:"
  echo "  ($axis_a_x, $axis_a_y) -> ($axis_b_x, $axis_b_y)"
  echo "================================================================="
  echo "Simulation Parameters:"
  echo "  DOMAIN_WIDTH:         $DOMAIN_WIDTH"
  echo "  DOMAIN_LENGTH:        $DOMAIN_LENGTH"
  echo "================================================================="
} | tee "{$LOG_FILE}_wolfsgrube"

##########################################################################
#	SETUP PSL MESH FUNCTION                                                #
# uses all globals                                                       #
##########################################################################
setup_psl_mesh() {
  # log
  echo "setup mesh with parameters:"
  echo "  mesh type: $PSL_MESH_TYPE"
  echo "  cell size: $cell_size"
  echo "  height: $height"
  echo "  sim dir: $PSL_DIR"

  # if the mesh type is block mesh then we need to generate
  if [ "$PSL_MESH_TYPE" = "blockMesh" ]; then
    echo "Perform mesh2quad conversion"
    $PYTHON "$TOOLS_DIR"/mesh2quad.py \
      -i "$DSL_DIR"/constant/surface.stl \
      -o "$PSL_DIR"/constant/surface.obj \
      --axis-a $axis_a_x $axis_a_y \
      --axis-b "$axis_b_x" "$axis_b_y" \
      --width "$DOMAIN_WIDTH" \
      --cell-size "$cell_size"
    # prepare block mesh dict
    r=$(bc -l <<<"($height / ($cell_size * 1.25))/1;")
    r=${r%.*}
    $PYTHON "$TOOLS_DIR"/quad2block_mesh_desc.py \
      -i "$PSL_DIR"/constant/surface.obj \
      -o "$PSL_DIR"/system/blockMeshDict \
      --grading 1 1 4 \
      --res 1 1 "$r" \
      --height $height \
      --inclined-top \
      --axis-a $axis_a_x $axis_a_y \
      --axis-b "$axis_b_x" "$axis_b_y" \
      --domain-type tunnel
  else
    # setup simulation geometry (terrain)
    "$C_TOOLS_DIR"/stl2foam \
      -p terrain \
      --inclined-top \
      --slope-axis-a $axis_a_x,$axis_a_y \
      --slope-axis-b "$axis_b_x","$axis_b_y" \
      --height $height \
      -i "$DSL_DIR"/constant/surface.stl \
      -o "$PSL_DIR"/constant/surface.stl
  fi
}

###########################################################
#	RUN DSL + PSL FUNCTION                                  #
# uses all globals                                        #
###########################################################
# log
{
  echo "Running PSL full sim with parameters:"
  echo "  mesh type: $PSL_MESH_TYPE"
  echo "  cell size: $cell_size"
  echo "  height: $height"
  echo "  sim dir: $PSL_DIR"
} | tee "$PSL_LOG_FILE"

{
  if [ $USE_DSL_SIM -eq 0 ]; then
    # the dsl simulation for the wolfsgrube case comes from the paper...
    # so we just copy their files and run the simulation
    [ ! -d "$DSL_DIR" ] && mkdir "$DSL_DIR"
    cp -r "$ASSETS_DIR"/* "$DSL_DIR"

    setup_dsl pMesh

    # copy assets
    #      echo "copying references from $ASSETS_DIR to $DSL_DIR"
    #cp -r "$ASSETS_DIR"/dsl/* "$DSL_DIR"

    run_sim "$DSL_DIR"
  else
    echo "using DSL simulation data from $DSL_DIR"
  fi
} | tee "{$DSL_LOG_FILE}_wolfsgrube"

# for the psl, we need to setup the mesh and input data
{
  setup_psl_3d $PSL_MESH_TYPE
  # prepare mesh
  setup_psl_mesh
  # prepare psl input data
  dsl_to_psl

  run_sim "$PSL_DIR"

} 2>&1 | tee -a "{$PSL_LOG_FILE}_wolfsgrube"


