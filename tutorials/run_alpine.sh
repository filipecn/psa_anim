#!/bin/bash
# Another example of running a full case of DSL+PSL for a given input mesh.

export ASSETS_DIR=$PSA_ANIM_TUTORIALS/assets/alpine

# setup parameters
export height=200
export cell_size=2
export TERRAIN_ALIGNED=1
export DOMAIN_WIDTH=1500
export DOMAIN_LENGTH=1000

bash $PSA_ANIM_TUTORIALS/run_gis.sh 
