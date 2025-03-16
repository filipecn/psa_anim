#!/bin/bash
# This script runs a series of simulations for method analysis.
# The same scene is used under many different set of configurations.

if [ -z "$PSA_ANIM_TUTORIALS" ]; then
  echo "error: Nadare enviroment variables (PSA_ANIM_*) not set. Please run"
  echo "       'source PATH/TO/PSA_ANIM/BUILD/nadare_variables.sh'"
  exit
fi

RENDER_WITH_PARAVIEW_PROGRAM="$PSA_ANIM_SCRIPTS"/render_with_paraview.py
TUTORIAL=ramp
TUTORIAL_RUN=$PSA_ANIM_TUTORIALS/run_tutorial.sh
MATRIX_OUTPUT_DIR=$(pwd)/param_matrix
MATRIX_FILE=$MATRIX_OUTPUT_DIR/"$TUTORIAL"/parameters.txt

if [ -f $MATRIX_FILE ]; then
  rm $MATRIX_FILE
fi
mkdir -p $MATRIX_OUTPUT_DIR

# MESH OPTIONS
duration=30
# METHOD OPTIONS 
boundary_conditions=("rev-full-wind-open")
# PSL PARAMETERS
entrainment_methods=(8)
snow_densities=("1.4" "2.5" "5.0" "7.0")
snow_viscosities=("1e-06" "1e-04" "1e-01" "1")
u_factors=("1.0" "2.0" "4.0" "8.0")
a_factors=("0.01" "0.1" "0.5" "0.9")
f_factors=(10 20 40 80)
dabs=("2e-01" "2e-02" "2e-04" "2e-06")
ebs=(25 50 75 100)
cell_sizes=("0.5" "1" "5" "10" "0.1")

snow_densities=("1.4")
u_factors=("1.6")
a_factors=("0.1")
f_factors=(40)
ebs=(50)
dabs=("2e-04")
snow_viscosities=("1e-04")
entrainment_methods=(80)
cell_sizes=(2)

cell_sizes=("0.5" "1" "5" "10" "0.1")

# for each DSL variation we spawn a set of PSL sims numbered by an <ID> 
# where the directory is named as PSL_<ID> 
# A file PSL_parameters.txt contains all options for each PSL_<ID>

# First run the DSL 
#bash $TUTORIAL_RUN --name $TUTORIAL --only-dsl --duration $duration -o $MATRIX_OUTPUT_DIR --no-render --run-parallel --duration 50
mkdir param_matrix/ramp

psl_id=0
for cell_size in "${cell_sizes[@]}"; do
  for bc_option in "${boundary_conditions[@]}"; do
    for snow_density in "${snow_densities[@]}"; do
      for entrainment_method in "${entrainment_methods[@]}"; do
        for u_factor in "${u_factors[@]}"; do
          for a_factor in "${a_factors[@]}"; do
            for f_factor in "${f_factors[@]}"; do
              for dab in "${dabs[@]}"; do
                for eb in "${ebs[@]}"; do
                  for snow_viscosity in "${snow_viscosities[@]}"; do
                    output_dir=$MATRIX_OUTPUT_DIR/"$TUTORIAL"/psl_"$psl_id"
                    # use_dsl_dir=$MATRIX_OUTPUT_DIR/"$TUTORIAL"/dsl
                    use_dsl_dir=/home/filipecn/dev/nadare-sim/dsl
                    $TUTORIAL_RUN --name $TUTORIAL                      \
                      --use-dsl-sim  $use_dsl_dir                       \
                      --duration $duration                              \
                      --boundary-condition $bc_option                   \
                      --cell-size $cell_size                            \
                      --snow-density $snow_density                      \
                      --snow-viscosity $snow_viscosity                  \
                      --entrainment-method $entrainment_method          \
                      --u-factor $u_factor                              \
                      --a-factor $a_factor                              \
                      --f-factor $f_factor                              \
                      --dab $dab                                        \
                      --erosion-energy $eb                              \
                      -o $output_dir                                    \
                      --2d                                              
                      #--run-parallel                                    \
                      #--parallel-x 20 \
                      #--parallel-y 1 \
                    

                    stats_dir=$output_dir/${TUTORIAL}/stats
                    mkdir $stats_dir

                    output_dir=${output_dir}/${TUTORIAL}/psl2
                    $PSA_ANIM_TOOLS/stats -i $output_dir -o $stats_dir --dx $cell_size --renumber-frames
                    
                    current_dir=$(pwd)
                    cd $output_dir
                    # if [ ! -d "VTK" ]; then
                    #   foamToVTK
                    # fi
                    cd $current_dir
                    ((psl_id=psl_id+1))
                    echo "$output_dir $psl_id cell size $cell_size bc $bc_option entrainment_method $entrainment_method" u_factor $u_factor a_factor $a_factor f_factor $f_factor dab $dab eb $eb snow_viscosity $snow_viscosity >> $MATRIX_FILE
                  done
                done
              done
            done
          done
        done
      done
    done
  done
done

python3 $PSA_ANIM_SCRIPTS/gen_report.py -i "$MATRIX_OUTPUT_DIR" -o "$MATRIX_OUTPUT_DIR" -n stats

exit 0

# generate videos
pvbatch $RENDER_WITH_PARAVIEW_PROGRAM -i $MATRIX_OUTPUT_DIR/${TUTORIAL} -o $MATRIX_OUTPUT_DIR/${TUTORIAL}

video_files=($(ls -d $MATRIX_OUTPUT_DIR/${TUTORIAL}/*.avi))
video_count=${#video_files[@]}
group_count=$(($video_count/4))
for ((i=0; i<=$group_count; i++))
do
  input0=${video_files[$(($i * 4 + 0))]}
  input1=${video_files[$(($i * 4 + 1))]}
  input2=${video_files[$(($i * 4 + 2))]}
  input3=${video_files[$(($i * 4 + 3))]}

  output="${MATRIX_OUTPUT_DIR}/${TUTORIAL}/final"$i".avi"

  if [ ! -z $input3 ]; then
    ffmpeg  -i $input0 -i $input1 -i $input2 -i $input3                                         \
            -filter_complex "[0:v][1:v][2:v][3:v]xstack=inputs=4:layout=0_0|w0_0|0_h0|w0_h0[v]" \
            -map "[v]" $output
  elif [ ! -z $input2 ]; then
    ffmpeg  -i $input0 -i $input1 -i $input2                    \
            -filter_complex "[0:v][1:v][2:v]hstack=inputs=3[v]" \
            -map "[v]" $output
  elif [ ! -z $input1 ]; then
    ffmpeg  -i $input0 -i $input1            \
            -filter_complex hstack=inputs=2  \
            $output 
  else
    cp $input0 $output
  fi
done
