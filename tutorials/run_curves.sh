#!/bin/bash
# This script runs a series of tests for various complex ramp inclines and avalanche parameters
# The DSL data can be simulated or procedurally generated.
# Zero dsl examples serve to test the psl propagation in a lock-exchange configuration
# The simulations are generated for diverse types of meshes
# The output directories are:
# <TERRAIN_TYPE>
#    |
#    |- slope<ANGLE>
#         	|
#          	|- <TEST_TYPE>

# list of inclines
slopes=(25)

# complex terrains
terrain_functions=(0)

DELETE_ROOT_BEFORE=0
ROOT_DIR="curves"
RENDER=0
REPORT=0
RUN_2D_FLAGS=""
if [ $RUN_2D -eq 1 ]; then
  RUN_2D_FLAGS="--d2"
fi

while [ "$1" != "" ]; do
  case $1 in
  -o | --output)
    shift
    ROOT_DIR=$1
    ;;
  -r | --render)
    RENDER=1
    ;;
  -R | --report)
    REPORT=1
    ;;
  --only-render)
    ONLY_RENDER=1
    ;;
  *)
    shift
    ;;
  esac
  shift
done

########################################################################################################################
# PARAMETERS                                                                                                           #
########################################################################################################################
# mesh size
height=100
width=1000
depth=250

echo "Running Complex Terrain Tests"
echo "================================================================="
echo "Parameters:"
echo "  ROOT_DIR:                        $ROOT_DIR"
echo "  WORKING_DIR:                     $WORKING_DIR"
echo "  DOMAIN_WIDTH:                    $width"
echo "  DOMAIN_LENGTH:                   $depth"
echo "  DELETE_ROOT_BEFORE:              $DELETE_ROOT_BEFORE"
echo "  RENDER:                          $RENDER"
echo "  REPORT:                          $REPORT"
echo "================================================================="

#######################################################################################################################
#	CREATE PATH FUNCTION                                                                                                #
# $1 path                                                                                                             #
#######################################################################################################################
create_path() {
  # create dirs if necessary
  if [ ! -d "$1" ]; then
    mkdir -p "$1"
  fi
}

#######################################################################################################################
#	SETUP MESH FUNCTION
#######################################################################################################################
setup_mesh() {
  RUN_2D_FLAG=$1
  # log
  echo "setup mesh with parameters:"
  echo "  tri terrain: $tri_input_mesh"
  echo "  quad terrain: $quad_input_mesh"
  echo "  mesh type: $PSL_MESH_TYPE"
  echo "  cell size: $cell_size"
  echo "  height: $height"
  echo "  sim dir: $sim_dir"

  $PYTHON "$TOOLS_DIR"/mesh2quad.py \
    -i $tri_input_mesh \
    -o $sim_dir/constant/quad_surface.obj \
    --terrain-aligned \
    --cell-size "$cell_size" \
    $RUN_2D_FLAG

  # resolution depends on cell size
  r=$(bc -l <<<"($height / $cell_size)/1;")
  r=${r%.*}
  $PYTHON $PSA_ANIM_SCRIPTS/quad2block_mesh_desc.py \
    -i $sim_dir/constant/quad_surface.obj \
    -o "$sim_dir"/system/blockMeshDict \
    --grading 1 1 3 \
    --res 1 1 $r \
    --height $height \
    --inclined-top \
    $RUN_2D_FLAG
}

#######################################################################################################################
#	RUN FULL SIMULATION FUNCTION
# uses globals:
# 	duration
#   width
#   slope
#######################################################################################################################
run_full() {
  full_dir=$1
  dsl_dir="$full_dir"/dsl
  psl_dir="$full_dir"/psl
  stats_dir=$full_dir/stats

  create_path $stats_dir

  # log
  echo "starting full sim with parameters:"
  echo "  duration: $duration"
  echo "  width: $width"
  echo "  slope: $slope"
  echo "  dsl path: $dsl_dir"
  echo "  psl path: $psl_dir"

  # INIT DSL
  #####################################################################################################################
  set_dsl_dir $dsl_dir
  if [ $ONLY_RENDER -eq 0 ]; then
    setup_dsl $PSL_MESH_TYPE
  fi

  # INIT MESH DICT FILE
  #####################################################################################################################
  sim_dir=$dsl_dir
  if [ $ONLY_RENDER -eq 0 ]; then
    setup_mesh
  fi

  # DSL RELEASE AREA
  #####################################################################################################################
  # here we choose a rectangle of 30% of the width and 50% of depth
  d_4=$(bc -l <<<"$depth / 6;")
  d_4=20
  w_10=$(bc -l <<<"$width * 0.05 * c($slope * 4 * a(1) / 180);")
  min_x=$(bc -l <<<"$w_10 * 2")
  max_x=$(bc -l <<<"$w_10 * 6")
  min_y=$(bc -l <<<"-$d_4 * 0.7")
  max_y=$(bc -l <<<"$d_4 * 0.7")

  max_x_e=$(bc -l <<<"$w_10 * 20")
  min_y_e=$(bc -l <<<"-$depth / 4;")
  max_y_e=$(bc -l <<<"$depth / 4;")


  # here we choose a rectangle with front in a fixed x 
  # and variable length to keep same initial mass
  front_x=200
  reference_length=100
  body_length=$(bc -l <<<"$reference_length * c($slope * 4 * a(1) / 180);")
  min_x=$(bc -l <<<"$front_x - $body_length")
  max_x=$(bc -l <<<"$front_x")
  min_y=$(bc -l <<<"-$d_4 * 0.7")
  max_y=$(bc -l <<<"$d_4 * 0.7")

  $PYTHON "$TOOLS_DIR"/set_release.py -o "$dsl_dir"/constant \
    --area 1.5 "$min_x" "$min_y" "$max_x" "$max_y" --entrain 0.4 "$min_x" "$min_y_e" "$max_x_e" "$max_y_e"

  # SIMULATE DSL
  #####################################################################################################################
  sim_dir=$dsl_dir
  run_sim dsl $sim_dir --no-clean

  # INIT PSL
  #####################################################################################################################
  set_psl_dir $psl_dir
  if [ $RUN_2D -eq 1 ]; then
    setup_psl_2d $PSL_MESH_TYPE
  else
    setup_psl_3d $PSL_MESH_TYPE
  fi

  # INIT MESH DICT FILE
  #####################################################################################################################
  sim_dir=$psl_dir
  setup_mesh $RUN_2D_FLAGS

  # DSL -> PSL
  #####################################################################################################################
  # run post process dsl
  $PYTHON $PSA_ANIM_SCRIPTS/dsl2psl.py \
    -i "$dsl_dir" -f "$dsl_dir" -o "$psl_dir"/constant/boundaryData -p terrain -op terrain -j --verbose
  # SIMULATE PSL
  #####################################################################################################################
  sim_dir=$psl_dir
  run_sim psl $sim_dir --no-clean

  # stats
  $PSA_ANIM_TOOLS/stats -i $psl_dir -o $stats_dir --dt $write_time --dx $cell_size --renumber-frames

  # vdb
  render_sim $dsl_dir $psl_dir $full_dir/vdb $full_dir $full_dir/surf

  notify "finished full sim $(pwd)"
}


# for each angle generate simulations for each different type of mesh
# for each different type of mesh run simulations
for terrain_function in "${terrain_functions[@]}"; do
  terrain_dir="$ROOT_DIR"/$terrain_function
  for slope in "${slopes[@]}"; do
    slope_dir="$terrain_dir/slope$slope"
    cd $WORKING_DIR
    create_path "$slope_dir"

    ####################
    # generate terrain #
    ####################
    quad_input_mesh="$slope_dir"/terrain.obj
    tri_input_mesh="$slope_dir"/terrain.stl
    $PYTHON $PSA_ANIM_SCRIPTS/gen_complex_terrain.py \
      --stl "$tri_input_mesh" \
      --obj "$quad_input_mesh" \
      -a "$slope" -w $width -d $depth -f "$terrain_function" -c $cell_size

    run_full $slope_dir
  done
done

#######################################################################################################################
#	REPORT                                                                                                              #
#######################################################################################################################

if [ $REPORT -eq 1 ]; then
  $PYTHON $PSA_ANIM_SCRIPTS/gen_report.py -i "$ROOT_DIR" -o "$ROOT_DIR" -n stats
fi
