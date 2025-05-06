#!/bin/bash

export ASSETS_DIR=$PSA_ANIM_TUTORIALS/assets/tantalo

export DOMAIN_WIDTH=600
export DOMAIN_LENGTH=1600
export TERRAIN_ALIGNED=0
export height=100
export cell_size=2
export axis_a_x=-375
export axis_a_y=-460
export axis_b_x=440
export axis_b_y=991

bash $PSA_ANIM_TUTORIALS/run_gis.sh 
