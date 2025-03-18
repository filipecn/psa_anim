#!/bin/bash

export ASSETS_DIR=$PSA_ANIM_TUTORIALS/assets/wolfsgrube

export DOMAIN_WIDTH=1000
export DOMAIN_LENGTH=2500
export TERRAIN_ALIGNED=0
export height=100
export cell_size=3
export axis_a_x=1755
export axis_a_y=40
export axis_b_x=600
export axis_b_y=2040

bash $PSA_ANIM_TUTORIALS/run_gis.sh 
