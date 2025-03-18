// profile
const int do_profile = readLabel(runTime.controlDict().lookup("pslProfile"));
const word profile_file = word(runTime.controlDict().lookup("pslProfileFile"));
// terrain patch
const word terrain_patch_name = word(runTime.controlDict().lookup("pslTerrainPatchName"));
// entrainment
const int psl_entrainment_method = readLabel(runTime.controlDict().lookup("pslEntrainmentMethod"));
const scalar psl_snow_cover_density = readScalar(runTime.controlDict().lookup("pslSnowCoverDensity"));

const scalar psl_u_entrainment_factor = readScalar(runTime.controlDict().lookup("pslUEntrainmentFactor"));
const scalar psl_alpha_entrainment_factor = readScalar(runTime.controlDict().lookup("pslAEntrainmentFactor"));
const scalar psl_front_length_factor = readScalar(runTime.controlDict().lookup("pslFrontLengthFactor"));

// Find out which surface patch id corresponds to the terrain surface
int surface_patch_id = 0;
forAll(mesh.boundaryMesh(), patchi) {
  Info<< patchi << " " << mesh.boundaryMesh()[patchi].name() << endl;
  if(mesh.boundaryMesh()[patchi].name() == terrain_patch_name) {
    surface_patch_id = patchi;
  }
}

// get mountain surface patch
Info<<"picking " << surface_patch_id << " as terrain patch id" << endl;
