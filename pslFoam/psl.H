// compute snow injection into the powder cloud
// the mass increment is stored as snow concentration (alpha1) by the scalar
// field alpha_inj

/******************************************************************************/
/*                                                                      INPUT */
/******************************************************************************/
// DSL surface patch fields
auto dslH_surface_field = dslH.boundaryField()[surface_patch_id];
auto dslU_surface_field = dslU.boundaryField()[surface_patch_id];
auto sdf_field_ref = dslFrontSDF.boundaryField()[surface_patch_id];
auto slope_angle_field = dslSlope.boundaryField()[surface_patch_id];

// Most methods work on the idea that entrained snow comes from a snow
// cover sitting on the mountain surface
// The density of the snow cover is
dimensionedScalar snow_cover_density("snow_cover_density",
                                     dimensionSet(1, -3, 0, 0, 0, 0, 0),
                                     scalar(psl_snow_cover_density));

/******************************************************************************/
/*                                                                      DEBUG */
/******************************************************************************/
// DEBUG INFO (per terrain face)
struct PSLDebugData {
  int global_face_id{0};
  int local_face_id{0};
  Foam::Vector<double> fluid_velocity;
  Foam::Vector<double> dsl_velocity;
  double dsl_area{0};
  double snow_cover_mass{0};
  double dm{0};
  double entrainment_mass{0};
  double front_weight{0};
  double dsl_velocity_mag{0};
  double dsl_height{0};
  double ce{0};
  double we{0};
  double slope_angle{0};
  double entrained_depth{0};
  double acc_entrained_depth{0};
  double local_snow_cover_height{0};
};
// set of debug values
std::vector<PSLDebugData> psl_debug_data;

// auxiliary macros for debug data collect
#define INIT_CURRENT_DEBUG_DATA PSLDebugData current_debug_data
#define END_CURRENT_DEBUG_DATA psl_debug_data.emplace_back(current_debug_data)
#define DEBUG_DATA_FIELD(A) current_debug_data.A = A
#define DEBUG_DIM_DATA_FIELD(A) current_debug_data.A = A.value()

/******************************************************************************/
/*                                                                       MAIN */
/******************************************************************************/
// Entrainment happens close to mountain's surface
// Here we will consider only close cells from mountain surface, so let's
// iterate only over terrain boundary faces
forAll(U.boundaryField()[surface_patch_id], local_face_id) {

  // for each boundary face id in local index space, we need its owner cell id
  // in global index space The owner cell id can be retrieved from the global
  // face id
  auto global_face_id =
      local_face_id + mesh.boundaryMesh()[surface_patch_id].start();
  // The boundary cell id
  auto cell_id = mesh.owner()[global_face_id];
  // All computation is based on the the closest boundary face
  // The local index space face id:
  auto closest_boundary_mesh_face = local_face_id;

  // reset output fields
  alpha_inj_U.boundaryFieldRef()[surface_patch_id][closest_boundary_mesh_face] =
      0;
  alpha_inj[cell_id] = 0;
  alpha_inj.boundaryFieldRef()[surface_patch_id][closest_boundary_mesh_face] =
      0;
  Uinj.boundaryFieldRef()[surface_patch_id][closest_boundary_mesh_face] =
      vector(0, 0, 0);

  // only faces near front are considered here
  // remember that the input field marks faces outside the DSL with distance -1
  if (sdf_field_ref[local_face_id] < 0 || sdf_field_ref[local_face_id] > 1000) {
    auto acc_entrained_depth = acc_entrained_height[local_face_id];
    auto local_snow_cover_height = snow_cover_height[local_face_id];
    INIT_CURRENT_DEBUG_DATA;
    DEBUG_DATA_FIELD(local_face_id);
    DEBUG_DATA_FIELD(global_face_id);
    DEBUG_DATA_FIELD(acc_entrained_depth);
    DEBUG_DATA_FIELD(local_snow_cover_height);
    END_CURRENT_DEBUG_DATA;
    continue;
  }
  /****************************************************************************/
  /*                                                        ENTRAINMENT INPUT */
  /* 	For each boundary cell, entrainment may be computed based on: */
  /*    1  slope angle                                                        */
  /*    2  psl velocity                                                       */
  /*    3  cell size                                                          */
  /*    4  snow cover                                                         */
  /*    5  dsl                                                                */
  /****************************************************************************/
  /*                                 1 - Slope angle at closest boundary face */
  /* 	slope_angle      -> slope angle */
  /****************************************************************************/
  auto slope_angle = slope_angle_field[closest_boundary_mesh_face];
  /****************************************************************************/
  /*                                                         2 - PSL velocity */
  /* 	fluid_velocity    -> PSL velocity */
  /****************************************************************************/
  auto fluid_velocity = U[cell_id];
  /****************************************************************************/
  /*                                                            3 - cell size */
  /*  dsl_area 	       -> cell face area in m^2 */
  /*  cell_volume      -> cell volume in m^3 */
  /****************************************************************************/
  // Area of the closest DSL face
  dimensionedScalar dsl_area(
      "dsl_area", dimensionSet(0, 2, 0, 0, 0, 0, 0),
      scalar(mesh.magSf().boundaryField()[surface_patch_id]
                                         [closest_boundary_mesh_face]));
  // Cell volume
  dimensionedScalar cell_volume("cell_volume",
                                dimensionSet(0, 3, 0, 0, 0, 0, 0),
                                scalar(mesh.V()[cell_id]));
  /****************************************************************************/
  /*                                                           4 - snow cover */
  /*  snow_cover_height  -> snow cover height in meters */
  /*  snow_cover_mass    -> snow cover mass in kg */
  /****************************************************************************/
  dimensionedScalar local_snow_cover_height(
      "snow_cover_height", dimensionSet(0, 1, 0, 0, 0, 0, 0),
      // scalar(snow_cover_height[local_face_id])
      scalar(1.0));
  // Snow cover quantity
  auto snow_cover_mass =
      snow_cover_density * local_snow_cover_height * dsl_area;
  /****************************************************************************/
  /*                                                                  5 - DSL */
  /* 	dsl_velocity   	  -> DSL velocity */
  /*	dsl_velocity_mag  -> DSL velocity magnitude */
  /*  dsl_height        -> DSL height */
  /****************************************************************************/
  // DSL velocity at surface
  auto dsl_velocity =
      dslU.boundaryField()[surface_patch_id][closest_boundary_mesh_face];
  // DSL velocity magnitude
  dimensionedScalar dsl_velocity_mag("dsl_v_mag",
                                     dimensionSet(0, 1, -1, 0, 0, 0, 0),
                                     scalar(mag(dsl_velocity)));
  // DSL height in meters
  dimensionedScalar dsl_height(
      "dsl_height", dimensionSet(0, 1, 0, 0, 0, 0, 0),
      scalar(
          dslH.boundaryField()[surface_patch_id][closest_boundary_mesh_face]));
  /****************************************************************************/
  /*                                               ENTRAINMENT METHODS (MASS) */
  /* Here we compute the total mass entrained over the last time step through
   * the current surface face. The output will serve as boundary condition in a
   * later stage, so mass needs to be converted into volume concentration. */
  /*  The output is */
  /*		- ce: the volume concentration representing the entrained mass
   * quantity                                       */
  /****************************************************************************/
  // injection mass
  auto ce = 0.0;
  /****************************************************************************/
  /*                                 ENTRAINMENT METHODS (INJECTION VELOCITY) */
  /*  Now the actual values for mass entrainment are calculated based on */
  /*  the entrained mass quantity defined above. The output will serve as */
  /*  boundary conditions in a later stage, so mass needs to be converted */
  /*  into volume concentration. The output values are stored in */
  /*    - we:  the injection velocity */
  /****************************************************************************/
  // injection velocity
  dimensionedScalar we("we", dimensionSet(0, 1, -1, 0, 0, 0, 0), scalar(0));
  /****************************************************************************/
  /*                                                             FRONT WEIGHT */
  /*  The intensity of the injection is considered to happen only close to the
   * front. Threfore its influence decays as the distance d increases by a
   * factor of                                         1 */
  /*                                              ---------- */
  /*                                              exp(d / L) */
  /****************************************************************************/
  // compute front distance decay factor
  auto front_weight = 0.0;

  if (sdf_field_ref[local_face_id] >= 0 &&
      sdf_field_ref[local_face_id] < 10000) {
    auto weight =
        std::exp(sdf_field_ref[local_face_id] / psl_front_length_factor);
    front_weight = 1.0 / weight;
  }

  // auxiliary information
  auto face_center =
      mesh.Cf().boundaryField()[surface_patch_id][closest_boundary_mesh_face];
  auto face_normal =
      mesh.Sf().boundaryField()[surface_patch_id][closest_boundary_mesh_face];
  face_normal = -face_normal.normalise();

  switch (psl_entrainment_method) {
  case 80: {
    // WE NOISE + CE NOISE
    // we -> dsl U modulated by the noise

    float velocity_trigger = std::min(
        1.0,
        std::max(0.0, 1.0 / (1.0 + std::exp(-dsl_velocity_mag.value() + 10))));

    // w_h is the entrainment velocity
    double entrainment_velocity =
        2.0 * abs(dslW.boundaryFieldRef()[surface_patch_id]
                                         [closest_boundary_mesh_face]);
    auto entrainment_mass =
        entrainment_velocity * dsl_area * 1000.0 * runTime.deltaT();

    auto u_noise = fractalAllicatorNoise(face_center, 0.5);
    auto a_noise = fractalAllicatorNoise(face_center, 1.0);

    we.value() = u_noise * dsl_velocity_mag.value();
    we.value() += entrainment_velocity;
    we.value() *= psl_u_entrainment_factor; // gamma u
    we.value() *= front_weight;
    we.value() *= velocity_trigger;

    // ce comes from entrainment/deposition from DSL
    // total cell mass capacity
    auto full_cell_mass = rho1 * cell_volume;
    ce = (entrainment_mass / full_cell_mass).value();

    ce *= psl_alpha_entrainment_factor;
    ce *= a_noise;
    ce *= front_weight;
    ce *= velocity_trigger;

    // clamp
    ce = std::min(1.0, ce);
    ce = std::max(ce, 0.005);
  } break;
  default: {
    // compute we
    Info << "INVALID ENTRAINMENT METHOD" << endl;
  }
  }

  /****************************************************************************/
  /*                                                       ENTRAINMENT OUTPUT */
  /*  - alpha_inj:    the volume concentration representing the entrained mass
   * quantity                               */
  /*  - alpha_inj_U:  the injection velocity */
  /****************************************************************************/
  alpha_inj_U.boundaryFieldRef()[surface_patch_id][closest_boundary_mesh_face] =
      we.value();
  alpha_inj[cell_id] = ce;
  // here we also fix boundary conditions to sync with cell increment
  alpha_inj.boundaryFieldRef()[surface_patch_id][closest_boundary_mesh_face] =
      alpha_inj[cell_id];
  Uinj.boundaryFieldRef()[surface_patch_id][closest_boundary_mesh_face] =
      we.value() * face_normal;

  // decrement entrained height from snow cover
  double entrained_depth =
      we.value() * runTime.deltaT().value() * dsl_area.value();
  acc_entrained_height[local_face_id] += entrained_depth;
  auto acc_entrained_depth = acc_entrained_height[local_face_id];

  /****************************************************************************/
  /* 	                                                           DEBUG DATA */
  /****************************************************************************/
  auto lambda = 1.0;
  auto dm = snow_cover_density *
            (local_snow_cover_height.value() * lambda * std::cos(slope_angle) *
             dsl_area) *
            dsl_velocity_mag;
  auto entrainment_mass = dm * runTime.deltaT();
  INIT_CURRENT_DEBUG_DATA;
  DEBUG_DATA_FIELD(local_face_id);
  DEBUG_DATA_FIELD(global_face_id);
  DEBUG_DATA_FIELD(fluid_velocity);
  DEBUG_DATA_FIELD(dsl_velocity);
  DEBUG_DIM_DATA_FIELD(dsl_velocity_mag);
  DEBUG_DIM_DATA_FIELD(dsl_height);
  DEBUG_DIM_DATA_FIELD(dsl_area);
  DEBUG_DIM_DATA_FIELD(snow_cover_mass);
  DEBUG_DIM_DATA_FIELD(dm);
  DEBUG_DIM_DATA_FIELD(entrainment_mass);
  DEBUG_DATA_FIELD(ce);
  DEBUG_DIM_DATA_FIELD(we);
  DEBUG_DATA_FIELD(front_weight);
  DEBUG_DATA_FIELD(slope_angle);
  DEBUG_DATA_FIELD(entrained_depth);
  DEBUG_DATA_FIELD(acc_entrained_depth);
  DEBUG_DIM_DATA_FIELD(local_snow_cover_height);
  END_CURRENT_DEBUG_DATA;
}

/******************************************************************************/
/*                                                                      DEBUG */
/******************************************************************************/
// append data to the respective processor profile file
if (Pstream::myProcNo() == 0) {
  // write data labels
  FILE *fp = fopen("debug_data/labels", "w+");
  if (fp != nullptr) {
    fprintf(fp, "global_face_id "
                "local_face_id "
                "fluid_velocity.x "
                "fluid_velocity.y "
                "fluid_velocity.z "
                "dsl_velocity.x "
                "dsl_velocity.y "
                "dsl_velocity.z "
                "dsl_area "
                "snow_cover_mass "
                "dm "
                "entrainment_mass "
                "front_weight "
                "dsl_velocity_mag "
                "dsl_height "
                "ce "
                "we "
                "slope_angle "
                "entrained_depth "
                "acc_entrained_depth "
                "local_snow_cover_height ");
    fclose(fp);
  }
}

if (runTime.outputTime()) {
  alpha_inj.write();
  alpha_inj_U.write();
  // create debug variables file
  FILE *fp =
      fopen((std::string("debug_data/variables_") + runTime.timeName()).c_str(),
            "w+");
  if (fp) {
    for (const auto &psl_d : psl_debug_data) {
      fprintf(fp,
              "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
              "%lf %lf %lf %lf %lf\n",
              psl_d.global_face_id, psl_d.local_face_id,
              psl_d.fluid_velocity.x(), psl_d.fluid_velocity.y(),
              psl_d.fluid_velocity.z(), psl_d.dsl_velocity.x(),
              psl_d.dsl_velocity.y(), psl_d.dsl_velocity.z(), psl_d.dsl_area,
              psl_d.snow_cover_mass, psl_d.dm, psl_d.entrainment_mass,
              psl_d.front_weight, psl_d.dsl_velocity_mag, psl_d.dsl_height,
              psl_d.ce, psl_d.we, psl_d.slope_angle, psl_d.entrained_depth,
              psl_d.acc_entrained_depth, psl_d.local_snow_cover_height);
    }
    fclose(fp);
  }
}

#undef DEBUG_DATA_FIELD
#undef DEBUG_DIM_DATA_FIELD
#undef INIT_CURRENT_DEBUG_DATA
#undef END_CURRENT_DEBUG_DATA
