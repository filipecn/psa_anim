#!/bin/bash
# Another example of running a full case of DSL+PSL for a given input mesh.

export ASSETS_DIR=$PSA_ANIM_TUTORIALS/assets/logo

# setup parameters
export height=200
export cell_size=3
export TERRAIN_ALIGNED=0
export axis_a_x=0
export axis_a_y=0
export axis_b_x=1100
export axis_b_y=0
export DOMAIN_WIDTH=300
export DOMAIN_LENGTH=1170

bash $PSA_ANIM_TUTORIALS/run_gis.sh 
