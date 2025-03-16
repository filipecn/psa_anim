#!/bin/bash
# Another example of running a full case of DSL+PSL for a given input mesh.

{
  echo "================================================================="
  echo "Running GIS Simulation Case"
  echo "================================================================="
  echo "================================================================="
  if [ $TERRAIN_ALIGNED -eq 0 ]; then
    echo "Generate simulation for terrain with reference points:"
    echo "  ($axis_a_x, $axis_a_y) -> ($axis_b_x, $axis_b_y)"
  fi
  echo "  height                $height"
  echo "================================================================="
  echo "Simulation Parameters:"
  echo "  DOMAIN_WIDTH:         $DOMAIN_WIDTH"
  echo "  DOMAIN_LENGTH:        $DOMAIN_LENGTH"
  echo "  TERRAIN_ALIGNED:      $TERRAIN_ALIGNED"
  echo "================================================================="
} | tee "$LOG_FILE"

# compute domain width for the PSL simulation
# for that, lets compute the unit distance vector between point a and b of the avalanche axis
# and multiply by the domain width
# first we compute the original distance
if [ $TERRAIN_ALIGNED -eq 0 ]; then
  axis_x=$(bc -l <<<"$axis_b_x - $axis_a_x;")
  axis_y=$(bc -l <<<"$axis_b_y - $axis_a_y;")
  echo "($axis_a_x, $axis_a_y) and ($axis_b_x, $axis_b_y)"
  echo "($axis_x, $axis_y)"

  # normalize axis
  norm=$(bc -l <<<"sqrt($axis_x * $axis_x + $axis_y * $axis_y);")
  axis_x=$(bc -l <<<"$axis_x / $norm;")
  axis_y=$(bc -l <<<"$axis_y / $norm;")
  
  # recompute axis_b
  axis_b_x=$(bc -l <<<"$axis_a_x + $axis_x * $DOMAIN_LENGTH;")
  axis_b_x=${axis_b_x%.*}
  axis_b_y=$(bc -l <<<"$axis_a_y + $axis_y * $DOMAIN_LENGTH;")
  axis_b_y=${axis_b_y%.*}
fi

################################################################
#	SETUP PSL MESH FUNCTION                                      #
# uses all globals                                             #
################################################################
setup_psl_mesh() {
  # log
  echo "setup mesh with parameters:"
  echo "  mesh type: $PSL_MESH_TYPE"
  echo "  cell size: $cell_size"
  echo "  height: $height"
  echo "  sim dir: $PSL_DIR"

  # if the mesh type is block mesh then we need to generate
  if [ "$PSL_MESH_TYPE" = "blockMesh" ]; then
    # prepare block mesh dict
    r=$(bc -l <<<"($height / $cell_size)/2;")
    r=${r%.*}
    if [ $TERRAIN_ALIGNED -eq 0 ]; then
      $PYTHON "$TOOLS_DIR"/mesh2quad.py \
        -i "$DSL_DIR"/constant/surface.stl \
        -o "$PSL_DIR"/constant/surface.obj \
        --axis-a $axis_a_x $axis_a_y \
        --axis-b "$axis_b_x" "$axis_b_y" \
        --width "$DOMAIN_WIDTH" \
        --cell-size "$cell_size"
      $PYTHON "$TOOLS_DIR"/quad2block_mesh_desc.py \
        -i "$PSL_DIR"/constant/surface.obj \
        -o "$PSL_DIR"/system/blockMeshDict \
        --grading 1 1 4 \
        --res 1 1 "$r" \
        --height $height \
        --inclined-top \
        --axis-a $axis_a_x $axis_a_y \
        --axis-b "$axis_b_x" "$axis_b_y" \
        --domain-type $boundary_condition
    else
      # use terrain geometry to align boundaries
      $PYTHON "$TOOLS_DIR"/mesh2quad.py \
        -i "$DSL_DIR"/constant/surface.stl \
        -o "$PSL_DIR"/constant/surface.obj \
        --terrain-aligned \
        --cell-size "$cell_size" 
      $PYTHON "$TOOLS_DIR"/quad2block_mesh_desc.py \
        -i "$PSL_DIR"/constant/surface.obj \
        -o "$PSL_DIR"/system/blockMeshDict \
        --domain-type $boundary_condition \
        --grading 1 1 3 \
        --res 1 1 "$r" \
        --height $height \
        --inclined-top 
    fi
  else
    if [ $TERRAIN_ALIGNED -eq 0 ]; then
      # setup simulation geometry (terrain)
      "$C_TOOLS_DIR"/stl2foam \
        -p terrain \
        --slope-axis-a $axis_a_x,$axis_a_y \
        --slope-axis-b "$axis_b_x","$axis_b_y" \
        --height $height \
        -i "$DSL_DIR"/constant/surface.stl \
        -o "$PSL_DIR"/constant/surface.stl
    else
      "$C_TOOLS_DIR"/stl2foam \
        -p terrain \
        -i "$DSL_DIR"/constant/surface.stl \
        -o "$PSL_DIR"/constant/surface.stl \
        --height $height \
        --align-faces
        # --inclined-top \
    fi
  fi
}

##############################################################
#	RUN DSL + PSL FUNCTION                                     #
# uses all globals                                           #
##############################################################
run_full() {
  # log
  {
    echo "Running PSL full sim with parameters:"
    echo "  mesh type: $PSL_MESH_TYPE"
    echo "  cell size: $cell_size"
    echo "  height: $height"
    echo "  sim dir: $PSL_DIR"
    echo "  asset dir: $ASSETS_DIR"
  } | tee "$PSL_LOG_FILE"

  {
    if [ $USE_DSL_SIM -eq 0 ]; then

      if [[ -f "$ASSETS_DIR"/dem.asc ]]; then
        setup_dsl pMesh dem
      else
        setup_dsl pMesh
      fi
  
      # prepare release data
      [[ ! -d "$DSL_DIR"/constant/gisdata ]] && mkdir "$DSL_DIR"/constant/gisdata

      if [[ -f ${ASSETS_DIR}/releaseArea ]]; then
        #cp "$ASSETS_DIR"/release.dbf "$DSL_DIR"/constant/gisdata/release.dbf
        #cp "$ASSETS_DIR"/release.shp "$DSL_DIR"/constant/gisdata/release.shp
        #cp "$ASSETS_DIR"/release.shx "$DSL_DIR"/constant/gisdata/release.shx
        #
        #cp "$ASSETS_DIR"/aoi.dbf "$DSL_DIR"/constant/gisdata/aoi.dbf
        #cp "$ASSETS_DIR"/aoi.shp "$DSL_DIR"/constant/gisdata/aoi.shp
        #cp "$ASSETS_DIR"/aoi.shx "$DSL_DIR"/constant/gisdata/aoi.shx

        #cp "$ASSETS_DIR"/dem.asc "$DSL_DIR"/constant/gisdata/dem.asc
        
        cp "$ASSETS_DIR"/releaseArea "$DSL_DIR"/constant/releaseArea
      else
        dsl_release $ASSETS_DIR
      fi

        # setup simulation geometry (terrain)
      if [ -f $ASSETS_DIR/surface.stl ]; then
        if [ $TERRAIN_ALIGNED -eq 0 ]; then
          "$C_TOOLS_DIR"/stl2foam \
            -p terrain \
            -i "$ASSETS_DIR"/surface.stl \
            --inclined-top \
            --height 200 \
            -o "$DSL_DIR"/constant/surface.stl
        else
          "$C_TOOLS_DIR"/stl2foam \
            -p terrain \
            -i "$ASSETS_DIR"/surface.stl \
            --inclined-top \
            --height 200 \
            -o "$DSL_DIR"/constant/surface.stl \
            --align-faces
        fi
      fi

      run_sim "$DSL_DIR"
    else
      echo "using DSL simulation data from $DSL_DIR"
    fi
  } | tee "$DSL_LOG_FILE"

  # for the psl, we need to setup the mesh and input data
  {
    if [ $ONLY_DSL -eq 0 ]; then
      setup_psl_3d $PSL_MESH_TYPE
      # prepare mesh
      setup_psl_mesh
      # prepare psl input data
      dsl_to_psl

      run_sim "$PSL_DIR"
    fi

  } 2>&1 | tee -a "$PSL_LOG_FILE"
}

echo "initiating full psl sim..."
run_full
