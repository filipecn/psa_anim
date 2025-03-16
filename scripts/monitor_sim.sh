#!/bin/bash
# This script monitors the progress of a simulation

WORKING_DIR=$(pwd)
SIM_DIR=$WORKING_DIR
CLEAN=0
MOVE_FRAMES=0
CLEAN_FRAMES=0

usage() {
  echo "usage: monitor_sim [[-i sim-path] | [-h]]"
}

while [ "$1" != "" ]; do
  case $1 in
  -i | --sim-path)
    shift
    SIM_DIR=$1
    ;;
  --clean)
    CLEAN=1
    ;;
  --no-clean)
    CLEAN=0
    ;;
  --move-frames)
    MOVE_FRAMES=1
    ;;
  --clean-frames)
    CLEAN_FRAMES=1
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

notify_progress() {
  paplay /usr/share/sounds/gnome/default/alerts/drip.ogg
}

#############
# variables #
#############
# try for log.pslFoam or log.faSavageHutterFoam
log_files=("log.faSavageHutterFoam" "log.pslFoam" "log.pslFoam3")

log_file=""

echo "sim path: $SIM_DIR"

get_log_file() {
  if [ -z "$log_file" ]; then
    return
  fi

  for l in "${log_files[@]}"; do
    log_file="$SIM_DIR"/$l
    echo "looking for $log_file in $SIM_DIR"
    if [ -f "$log_file" ]; then
      break
    fi
  done
  
  if [ ! -f "$log_file" ]; then
    echo "log file not found: $log_file"
  else
    echo "... $log_file found"
  fi
}

############################
# grep stats from log file #
# checks in $log_file      #
############################
get_log() {
  if [ ! -f "$log_file" ]; then
    return
  fi
  info=("deltaT" "max(mag(U))" "Min(alpha.snow)" "Courant" "Time")
  for i in "${info[@]}"; do
    tail "$log_file" -n 100 | grep "$i" | tail -1
  done

  echo "++++++++++++++++++++++++++++++++++++++++++++"
}

##########################
# gets last stored frame #
##########################
get_last_frame() {
  cd "$SIM_DIR" || return
  if [ -d "processor0" ]; then
    cd processor0 || return
  fi
  # shellcheck disable=SC2012
  ls -d -- [0-9]*.*[0-9]* > /dev/null 2> /dev/null
  if [[ $? == 0 ]]; then
    last_frame=$(ls -d -- [0-9]*.*[0-9]* | sort --version-sort | tail -1)
    echo "$(date) last frame: $last_frame"
  fi
}

##############################
# remove not necessary files #
##############################

remove_list=("alpha_inj" "dsl_divU" "dsl_slope" "dsl_U" "Uinj" "alpha_inj_U" "alpha.snow_0" "dsl_h" "dsl_W" "Us_0" "Us_mag.asc" "Us_x.asc" "Us_y.asc" "Us_z.asc" "W.asc" "divUs.asc" "Us" "fa_hentrain" "hentrain" "hentrain.asc" "h.asc" "hentrain_0" "W" "phi" "rhoPhi" "p_rgh" "p" "dsl_front_sdf")
#remove_list=("dsl_front_sdf" "Uinj" "alpha_inj_U" "alpha.snow_0" "dsl_h")

clean_dir() {
  cd "$1" || return
  ls -d -- [0-9]*.*[0-9]* > /dev/null 2> /dev/null
  if [[ $? == 0 ]]; then
    times=($(ls -d -- [0-9]*.*[0-9]* | sort -g))
    array_size=${#times[@]}
    if (( $array_size > 1 )); then
      # delete files from at most n - 1 times 
      # i.e. keep the 5 most recent ones
      last_i=$(( $array_size - 1 ))
      for i in "${!times[@]}"; do
        if (( $i < $last_i )); then
          # remove unecessary files
          for rem in ${remove_list[@]}; do
            file=${times[$i]}/${rem}
            if [[ -f $file ]]; then
              echo "rm $1 / $(du -h $file)"
              rm -r $file
            fi
          done
        fi
      done
    fi
  fi
}

clean_par() {
  cd $SIM_DIR
  if [ ! -d "processor0" ]; then
    return
  fi
  folders=($(ls -d processor*))
  ls -d -- [0-9]*.*[0-9]* > /dev/null 2> /dev/null
  if [[ $? == 0 ]]; then
    times=($(ls -d -- [0-9]*.*[0-9]* | sort -g))
    array_size=${#times[@]}
    if (( $array_size > 5 )); then
      # delete files from at most n - 5 times 
      # i.e. keep the 5 most recent ones
      last_i=$(( $array_size - 5 ))
      for i in "${!times[@]}"; do
        if (( $i < $last_i )); then
          for f in "${folders[@]}"; do
            if [ -d $f/${times[$i]} ]; then
              echo "remove $f/${times[$i]}"
              # rm -r  $f/${times[$i]}
            fi
          done
        fi
      done
    fi
  fi
}

clean() {
  folders=($SIM_DIR)
  if [ -d "${SIM_DIR}/processor0" ]; then
    folders+=($(ls -d ${SIM_DIR}/processor*))
  fi
  for folder in ${folders[@]}; do
    cd $WORKING_DIR
    clean_dir $folder
  done
}

move_frames() {
  if [[ ! -d frames ]]; then
    mkdir frames
  fi
  echo "$(pwd)"
  ls -d -- [0-9]*.*[0-9]* > /dev/null 2> /dev/null
  if [[ $? == 0 ]]; then
    times=($(ls | egrep "^[0-9]+\.?[0-9]*$"| sort -g))
    array_size=${#times[@]}
    if (( $array_size > 0 )); then
      for i in "${!times[@]}"; do
        echo "storing frame " ${times[$i]}
        mv ${times[$i]} frames/${times[$i]}
      done
    fi
  fi
}

remove_par_frames() {
  cd $SIM_DIR
  if [ ! -d "processor0" ]; then
    return
  fi
  folders=($(ls -d processor*))
  for f in "${folders[@]}"; do
    cd $f 
    ls -d -- [0-9]*.*[0-9]* > /dev/null 2> /dev/null
    if [[ $? == 0 ]]; then
      times=($(ls | grep -E "^[0-9]+\.?[0-9]*$"| sort -g))
      array_size=${#times[@]}
      if (( $array_size > 1 )); then
        echo $f $array_size
        last_i=$(( $array_size - 1 ))
        for i in "${!times[@]}"; do
          if (( $i < $last_i )); then
            if [[ "${times[$i]}" != "0" ]]; then
              rm -r ${times[$i]}
              echo "remove $f/${times[$i]}"
            fi
          fi
        done
      fi
    fi
    cd .. 
  done
}

# main loop
###########
if [ $MOVE_FRAMES -eq 1 ]; then 
  cd "$SIM_DIR"
  move_frames
  exit
fi

while true; do
  if [ ! -d "$SIM_DIR" ]; then
    pwd
    echo "invalid simulation path: $SIM_DIR"
  else
    cd "$WORKING_DIR"
    get_last_frame
    cd "$WORKING_DIR"
    get_log_file
    cd "$WORKING_DIR"
    get_log
    if [ $CLEAN -eq 1 ]; then
      cd "$WORKING_DIR"
      clean
    fi
    if [ $CLEAN_FRAMES -eq 1 ]; then
      cd "$WORKING_DIR"
      remove_par_frames
    fi
    cd "$WORKING_DIR"
  fi
  sleep 10
done
