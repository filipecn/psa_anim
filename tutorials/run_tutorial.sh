#!/bin/bash
# This script runs a given tutorial

if [ -z "$PSA_ANIM_SCRIPTS" ] || [ -z "$PSA_ANIM_TOOLS" ] || [ -z "$PSA_ANIM_RENDER_TOOLS" ]; then
  echo "error: PSA_ANIM enviroment variables (PSA_ANIM_*) not set. Please run"
  echo "       'source PATH/TO/PSA_ANIM/BUILD/psa_anim_variables.sh'"
  exit
fi

# check openfoam enviroment
if ! command -v blockMesh &> /dev/null; then
  echo "error: Openfoam enviroment not active."
  exit
fi

# COMMAND
BASE_COMMAND=$0
BASE_ARGS=$args

# tools
export PYTHON=python3
export TOOLS_DIR=$PSA_ANIM_SCRIPTS
export C_TOOLS_DIR=$PSA_ANIM_TOOLS
export RENDER_TOOLS=$PSA_ANIM_RENDER_TOOLS
export FOAM2VDB="$C_TOOLS_DIR"/foam2vdb
# paths
export SIM_DIR=$(pwd)
export OUTPUT_DIR=""
export WORKING_DIR=$(pwd)
export DSL_DIR="$SIM_DIR"/dsl
export PSL_DIR=$SIM_DIR/psl
export TERRAIN_SURFACE_DIR=$OUTPUT_DIR
export VDB_DIR=$OUTPUT_DIR/vdb
export DSL_SURFACE_DIR="$OUTPUT_DIR"/dsl_surf
# logging
export verbose=""
export VERBOSE=0
export LOG_FILE="$SIM_DIR"/log.run_main
export PSL_LOG_FILE="$SIM_DIR"/log.psl
export DSL_LOG_FILE="$SIM_DIR"/log.dsl
export RENDER_LOG_FILE="$OUTPUT_DIR"/log.render
# options
export RUN_DSL=1
export RUN_PSL=1
export USE_DSL_SIM=0
export RUN_PARALLEL=1
export RENDER_SIM=1
export ONLY_RENDER=0
export RUN_2D=0
export GEN_DEBUG_DATA=0
export PSL_MESH_TYPE="blockMesh"
# parallel decomposition
export PARALLEL_X=5
export PARALLEL_Y=4
export PARALLEL_Z=1
# setup parameters
export duration=50
export delta_t=0.0001
export write_time=0.05
export reconstruct_step=0
export cell_size=5
export snow_density="1.4"
export snow_viscosity="1e-04"
export dab="2e-04"
export write_format="ascii"
export boundary_condition="rev-full-wind-open"
export entrainment_method=80
export n_correctors=2
export n_northo_correctors=2
export max_co="0.2"
export entrainment_u_factor="1.6"
export entrainment_a_factor="0.1"
export entrainment_f_factor="80"
export erosion_energy="50"
# openfoam boundary patch types
export fb_type="patch"
export b_type="patch"
export l_type="patch"
export r_type="patch"
export t_type="patch"

usage() {
  cat <<-USAGE_TXT

  Usage: tutorial [OPTION]...
  Run a full PSA_ANIM case (DSL + PSL) workflow (simulation + vdb):
    - run a DSL simulation
    - convert DSL output into PSL input
    - run the PSL simulation
    - produce final animation assets (psl vdb and dsl surface mesh files)


  Example: tutorial --name wolfsgrube

  Default values are listed along with the argument types or at the end of the
  of each description between brackets (ex: [float|DEFAULT VALUE]).

  Required:
          --name [NAME]                  Case name (ex: wolfsgrube).

  Options:
          --no-render                    Skip vdb generation.
          --only-render                  No simulation, only vdb generation.
          --only-dsl                     Skip PSL simulation.
          --only-psl                     Skip DSL simulation.
          --single-thread                Run sequential (default is parallel).
          --2d                           Run 2D version (drop y axis).
          --parallel-x [N|5]             N thread divisions in x direction.
          --parallel-y [N|4]             N thread divisions in y direction.
          --parallel-z [N|1]             N thread divisions in z direction.
          --binary                       Set OpenFOAM's output files as binary.
          --verbose                    

  Simulation:
          --write-time [float|0.05]      Simulation output file time.           
          --dt [float|0.0001]            Simulation timestep (in seconds).
          --duration [float|50]          Simulation duration (in seconds).
          --cell-size [float|5]          Mesh's base cell size.
  PSL parameters:
          --u-factor [float|1.6]
          --a-factor [float|0.1]
          --f-factor [float|80]
          --snow-density [float|1.7]
          --snow-viscosity [float|1e-04]
          --dab [float|2e-04]


  Paths:
      --use-dsl-sim [PATH]               Use dsl simulation from given PATH.
      -o, --output [PATH]                Simulation output dir.
                                         [folder with case's name in current dir]
      -O [PATH]                          Result dir (vdb, dsl_surf).
                                         [simulation output dir]
      -P, --python [PATH|python3]        Python binary path 
USAGE_TXT
}

tutorial_name=""
args=$@

while [ "$1" != "" ]; do
  case $1 in
  -P | --python)
    shift
    PYTHON=$1
    ;;
  --name)
    shift
    tutorial_name=$1
    ;;
  -o | --output)
    shift
    SIM_DIR=$(realpath "$1")
    ;;
  -O)
    shift
    OUTPUT_DIR=$(realpath "$1")
    ;;
  --binary)
    write_format="binary"
    ;;
  --output-step)
    shift
    reconstruct_step=$1
    ;;
  --write-time)
    shift
    write_time=$1
    ;;
  --use-dsl-sim)
    shift
    USE_DSL_SIM=1
    DSL_DIR=$(realpath "$1")
    ;;
  --dt)
    shift
    delta_t=$1 
    ;;
  --duration)
    shift
    duration=$1
    ;;
  --cell-size)
    shift
    cell_size=$1
    ;;
  --snow-density)
    shift
    snow_density=$1
    ;;
  --snow-viscosity)
    shift
    snow_viscosity=$1
    ;;
  --dab)
    shift
    dab=$1
    ;;
  --boundary-condition)
    shift
    boundary_condition=$1
    ;;
  --entrainment-method)
    shift
    entrainment_method=$1
    ;;
  --u-factor)
    shift
    entrainment_u_factor=$1
    ;;
  --a-factor)
    shift
    entrainment_a_factor=$1
    ;;
  --f-factor)
    shift
    entrainment_f_factor=$1
    ;;
  --max-co)
    shift
    max_co=$1
    ;;
  --single-thread)
    RUN_PARALLEL=0
    ;;
  --parallel-x)
    shift
    PARALLEL_X=$1
    ;;
  --parallel-y)
    shift
    PARALLEL_Y=$1
    ;;
  --parallel-z)
    shift
    PARALLEL_Z=$1
    ;;
  --no-render)
    RENDER_SIM=0
    ;;
  --2d)
    RUN_2D=1
    ;;
  --only-dsl)
    RUN_PSL=0
    ;;
  --only-psl)
    RUN_DSL=0
    ;;
  --only-render)
    ONLY_RENDER=1
    ;;
  --open-bc)
    boundary_condition="open"
    entrainment_u_factor="1.0"
    entrainment_a_factor="0.4"
    snow_density="7.0"
    ;;
  --optimized)
    max_co="0.7"
    n_correctors=1
    n_northo_correctors=1
    ;;
  --verbose)
    VERBOSE=1
    verbose="--verbose"
    ;;
  -h | --help)
    usage
    exit
    ;;
  *)
    shift
    ;;
  esac
  shift
done

################################################################### input check
if [ "$tutorial_name" == "" ]; then
  echo "tutorial name is required"
  usage
  exit
fi

################################################################# prepare paths
SIM_DIR=$SIM_DIR/$tutorial_name
if [ $USE_DSL_SIM -eq 0 ]; then
  DSL_DIR="$SIM_DIR"/dsl
fi
PSL_DIR=$SIM_DIR/psl
if [ $RUN_2D -eq 1 ]; then
  # RUN_PARALLEL=0
  PSL_DIR=$SIM_DIR/psl
fi

########################################################### render/output paths
if [ "$OUTPUT_DIR" == "" ]; then
  OUTPUT_DIR=$SIM_DIR
fi

TERRAIN_SURFACE_DIR=$OUTPUT_DIR
VDB_DIR=$OUTPUT_DIR/vdb
DSL_SURFACE_DIR="$OUTPUT_DIR"/dsl_surf

######################################################################### tools
FOAM2VDB="$C_TOOLS_DIR"/foam2vdb

####################################################################### logging
LOG_FILE="$SIM_DIR"/log.run_${tutorial_name}
PSL_LOG_FILE="$SIM_DIR"/log.psl
DSL_LOG_FILE="$SIM_DIR"/log.dsl
RENDER_LOG_FILE="$OUTPUT_DIR"/log.render

################################################################### create dirs
[ ! -d "$SIM_DIR" ] && mkdir -p "$SIM_DIR"
[ ! -d "$OUTPUT_DIR" ] && mkdir -p "$SIM_DIR"


################################################################################
#	PARAMETERS LOG FUNCTION                                                      #
################################################################################
log_parameters() {
  {
    echo "Running Simulation Case $name"
    echo "$0"
    echo "$args"
    echo "================================================================="
    echo "Using paths:"
    echo "  PYTHON CMD:           $PYTHON"
    echo "  TOOLS_DIR:            $TOOLS_DIR"
    echo "  C_TOOLS_DIR:          $C_TOOLS_DIR"
    echo "  SIM_DIR:              $SIM_DIR"
    echo "  OUTPUT_DIR:           $OUTPUT_DIR"
    echo "  WORKING_DIR:          $WORKING_DIR"
    echo "  RENDER_TOOLS:         $RENDER_TOOLS"
    echo "  FOAM2VDB:             $FOAM2VDB"
    echo "================================================================="
    echo "Setting paths:"
    echo "  DSL_DIR:              $DSL_DIR"
    echo "  PSL_DIR:              $PSL_DIR"
    echo "  VDB_DIR:              $VDB_DIR"
    echo "  TERRAIN_SURFACE_DIR:  $TERRAIN_SURFACE_DIR"
    echo "  DSL_SURFACE_DIR:      $DSL_SURFACE_DIR"
    echo "Log files:"
    echo "  LOG_FILE:             $LOG_FILE"
    echo "  PSL_LOG_FILE:         $PSL_LOG_FILE"
    echo "  DSL_LOG_FILE:         $DSL_LOG_FILE"
    echo "  RENDER_LOG_FILE:      $RENDER_LOG_FILE"
    echo "================================================================="
    echo "Options:"
    echo "  RUN_2D:               $RUN_2D"
    echo "  RUN_DSL:              $RUN_DSL"
    echo "  RUN_PSL:              $RUN_PSL"
    echo "  RUN_PARALLEL:         $RUN_PARALLEL ($PARALLEL_X $PARALLEL_Y $PARALLEL_Z)"
    echo "  RENDER_SIM:           $RENDER_SIM"
    echo "  USE_DSL_SIM:          $USE_DSL_SIM"
    echo "  MESH TYPE:            $PSL_MESH_TYPE"
    echo "  WRITE_TIME:           $write_time"
    echo "  WRITE_FORMAT:         $write_format"
    echo "  RECONSTRUCT_STEP:     $reconstruct_step"
    echo "================================================================="
    echo "Simulation Parameters:"
    echo "  DURATION:             $duration"
    echo "  CELL_SIZE:            $cell_size"
    echo "  SNOW_DENSITY:         $snow_density"
    echo "  SNOW_VISCOSITY:       $snow_viscosity"
    echo "  DAB:                  $dab"
    echo "  MAX_CO:               $max_co"
    echo "  EROSION_ENERGY:       $erosion_energy"
    echo "  BOUNDARY_CONDITION:   $boundary_condition"
    echo "  ENTRAINMENT_METHOD:   $entrainment_method"
    echo "  ENTRAINMENT_U_FACTOR: $entrainment_u_factor"
    echo "  ENTRAINMENT_A_FACTOR: $entrainment_a_factor"
    echo "  N_CORRECTORS:         $n_correctors"
    echo "  N_NORTHO_CORRECTORS:  $n_northo_correctors"
    echo "================================================================="
  } | tee "$LOG_FILE"
}

set_dsl_dir() {
  DSL_DIR=$1
}

set_psl_dir() {
  PSL_DIR=$1
}

################################################################################
#	NOTIFY FUNCTION                                                              #
#   $1 - message                                                               #
################################################################################
notify() {
  echo "$1"
  # paplay /usr/share/sounds/gnome/default/alerts/drip.ogg
}
################################################################################
#	RUN SIMULATION FUNCTION                                                      #
#	  $1 - simulation path                                                       #
################################################################################
run_sim() {
  $PSA_ANIM_SCRIPTS/monitor_sim.sh -i $2 $3 &
  ID=$!
  cd "$2" || exit
  chmod 777 Allclean
  chmod 777 Allrun-parallel
  chmod 777 Allrun
  echo "launch monitor with pid $ID"
  if [[ ! -d "$OUTPUT_DIR/$1/frames" ]]; then
    echo "creating directory: $OUTPUT_DIR/$1/frames"
    mkdir -p $OUTPUT_DIR/$1/frames
  fi
  if [ $RUN_PARALLEL -eq 1 ]; then
    ./Allrun-parallel
    kill $ID
    if [[ "$1" == "psl" ]]; then
      $PSA_ANIM_SCRIPTS/reconstruct_par.sh -o $OUTPUT_DIR/$1/frames -e $duration --step $reconstruct_step --clean
    else
      $PSA_ANIM_SCRIPTS/reconstruct_par.sh -o $OUTPUT_DIR/$1/frames -e $duration --step $write_time --clean -np 8 --batch
    fi
  else
    ./Allrun
    kill $ID
    $PSA_ANIM_SCRIPTS/monitor_sim.sh -i $2 --move-frames
    mv -r frames/* $OUTPUT_DIR/$1/frames
  fi
  notify "... simulation complete!"
  cd "$WORKING_DIR" || exit
}
################################################################################
#	SETUP DSL FUNCTION                                                           #
################################################################################
setup_dsl() {
  if [ $USE_DSL_SIM -eq 0 ]; then
    generator=$1
    if [ -z "$2" ]; then
      $PYTHON "$TOOLS_DIR"/init_dsl.py            \
        -o "$DSL_DIR"                             \
        --template-dir "$TOOLS_DIR"/assets/dsl    \
        --et "$duration"                          \
        --wt "$write_time"                        \
        --eb $erosion_energy                      \
        --mesh-generator $generator               \
        $verbose
    else
      echo "setup_dsl with dem"
      $PYTHON "$TOOLS_DIR"/init_dsl.py            \
        -o "$DSL_DIR"                             \
        --template-dir "$TOOLS_DIR"/assets/dsl    \
        --et "$duration"                          \
        --wt "$write_time"                        \
        --eb $erosion_energy                      \
        --mesh-generator $generator               \
        --dem                                     \
        $verbose
    fi
  fi
}
################################################################################
#	SETUP 2D PSL FUNCTION                                                        #
################################################################################
setup_psl_2d() {
  generator=$1
  $PYTHON "$TOOLS_DIR"/init_psl.py              \
    -o "$PSL_DIR"                               \
    --template-dir "$TOOLS_DIR"/assets/psl      \
    --mesh-generator $generator                 \
    --domain-type $boundary_condition           \
    --write-format $write_format                \
    --dt $delta_t                               \
    --wt "$write_time"                          \
    --et "$duration"                            \
    --snow-density "$snow_density"              \
    --snow-viscosity "$snow_viscosity"          \
    --entrainment-method "$entrainment_method"  \
    --dab $dab                                  \
    --max-co $max_co                            \
    --u-factor "$entrainment_u_factor"          \
    --a-factor "$entrainment_a_factor"          \
    --f-factor "$entrainment_f_factor"          \
    --n-correctors "$n_correctors"              \
    --n-northo-correctors "$n_northo_correctors"\
    --profile                                   \
    --run-2d                                    \
    $verbose
}
################################################################################
#	SETUP 3D PSL FUNCTION                                                        #
################################################################################
setup_psl_3d() {
  generator=$1
  $PYTHON "$TOOLS_DIR"/init_psl.py                \
    -o "$PSL_DIR"                                 \
    --template-dir "$TOOLS_DIR"/assets/psl        \
    --mesh-generator $generator                   \
    --domain-type $boundary_condition             \
    --write-format $write_format                  \
    --snow-density "$snow_density"                \
    --snow-viscosity "$snow_viscosity"            \
    --entrainment-method "$entrainment_method"    \
    --dab   $dab                                  \
    --max-co $max_co                              \
    --u-factor "$entrainment_u_factor"            \
    --a-factor "$entrainment_a_factor"            \
    --f-factor "$entrainment_f_factor"            \
    --px $PARALLEL_X                              \
    --py $PARALLEL_Y                              \
    --pz $PARALLEL_Z                              \
    --dt $delta_t                                 \
    --wt "$write_time"                            \
    --et "$duration"                              \
    --n-correctors "$n_correctors"                \
    --n-northo-correctors "$n_northo_correctors"  \
    $verbose
}
################################################################################
#	PREPARE PSL DATA FUNCTION                                                    #
################################################################################
dsl_to_psl()
{
  if [ $GEN_DEBUG_DATA -eq 1 ]; then
    [ ! -d "$PSL_DIR"/debug_data ] && mkdir "$PSL_DIR"/debug_data
  fi

  # prepare psl input data
  $PYTHON "$TOOLS_DIR"/dsl2psl.py       \
    -i "$DSL_DIR"                       \
    -f "$OUTPUT_DIR/dsl/"               \
    -o "$PSL_DIR"/constant/boundaryData \
    -p terrain                          \
    -op terrain                         \
    $verbose
}
################################################################################
#	WRITE DSL INPUT                                                              #
################################################################################
dsl_release()
{
  release_files=$(ls $1/*release.obj)

  ls $1/*cover.obj > /dev/null 2> /dev/null
  if [[ $? == 0 ]]; then
    cover_files=$(ls $1/*cover.obj)
    $PYTHON "$TOOLS_DIR"/set_release.py   \
      -o "$DSL_DIR"/constant \
      --area-files $release_files \
      --entrain-files $cover_files 
    else
    $PYTHON "$TOOLS_DIR"/set_release.py   \
      -o "$DSL_DIR"/constant \
      --area-files $release_files 
  fi
}
################################################################################
#	RENDER SIM                                                                   #
################################################################################
render_sim()
{
  if [ -z $1 ]; then
    dsl_sim_dir=$DSL_DIR
  else
    dsl_sim_dir=$1
  fi

  if [ -z $2 ]; then
    psl_sim_dir=$PSL_DIR
  else
    psl_sim_dir=$2
  fi

  if [ -z $3 ]; then
    vdb_output_dir=$VDB_DIR
  else
    vdb_output_dir=$3
  fi

  if [ -z $4 ]; then
    terrain_output_dir=$TERRAIN_SURFACE_DIR
  else
    terrain_output_dir=$4
  fi

  if [ -z $5 ]; then
    dsl_output_dir=$DSL_SURFACE_DIR
  else
    dsl_output_dir=$5
  fi

  {
    echo "Render Paths"
    echo "============"
    echo "DSL SIM DIR:       $dsl_sim_dir"
    echo "PSL SIM DIR:       $psl_sim_dir"
    echo "VDB OUT DIR:       $vdb_output_dir"
    echo "TERRAIN OUT DIR:   $terrain_output_dir"
    echo "DSL OUT DIR:       $dsl_output_dir"
    # export terrain surface
    $PYTHON "$TOOLS_DIR"/foam2obj.py                \
      -i "$psl_sim_dir"                             \
      -p terrain                                    \
      -o "$terrain_output_dir"                      \
      $verbose
  
    # export terrain base
    $PYTHON "$TOOLS_DIR"/terrain_base.py            \
      -i "$terrain_output_dir"/terrain.obj          \
      -o "$terrain_output_dir"/terrain_base.obj     \
      $verbose
 
    if [ $RUN_DSL -eq 1 ]; then
      $PYTHON "$TOOLS_DIR"/dsl_surface.py \
        -i $dsl_sim_dir \
        -s $OUTPUT_DIR/dsl \
        -p terrain \
        -o "$dsl_output_dir" \
        --only-missing-frames \
        --with-height-map \
        --renumber-frames \
        $verbose
        #--with-distance-map \
    fi

    # convert psl output into openvdb volumes
    [ ! -d "$vdb_output_dir" ] && mkdir "$vdb_output_dir"

    if [ $RUN_PSL -eq 1 ]; then
      if [ $RUN_2D -eq 1 ]; then
        $FOAM2VDB -i "$psl_sim_dir"                            \
          -o "$vdb_output_dir"                                 \
          --frames-dir $OUTPUT_DIR/psl                         \
          -s 0.5                                               \
          --only-missing-frames                                \
          --renumber-frames                                    \
          --interpolation cell                                 \
          --color-by density                                   \
          --only-cells                                         \
          --ptu-size 10                                        \
          --every-n 1                                          \
          $verbose
          #--min 0.0001 \
          #--iso 0.0001 \
          # -2d  \
      else
        $FOAM2VDB -i "$psl_sim_dir"                            \
          -o "$vdb_output_dir"                                 \
          --frames-dir $OUTPUT_DIR/psl                         \
          -s 1.0                                               \
          --only-missing-frames                                \
          --renumber-frames                                    \
          --interpolation hrbf                                 \
          --ptu-size 10                                        \
          --every-n 1                                          \
          $verbose
          #--use-ptu                                            \
          #-f 600
      fi
    fi
  
  } #2>&1 | tee "$RENDER_LOG_FILE"
}

########################################################################### run

export -f log_parameters
export -f set_dsl_dir
export -f set_psl_dir
export -f setup_dsl
export -f setup_psl_2d
export -f setup_psl_3d
export -f dsl_to_psl
export -f notify
export -f run_sim
export -f dsl_release
export -f render_sim


if [ $ONLY_RENDER -eq 0 ]; then
  bash $PSA_ANIM_TUTORIALS/run_$tutorial_name.sh $args
fi

if [ $RENDER_SIM -eq 1 ]; then
  render_sim 
fi
