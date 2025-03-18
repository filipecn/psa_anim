// compute the various statistics for each time step
// The folling data is computed:
// 	1 -	mass
// 	2 - kinetic energy
// 	3 - entrainment
// 	4 - plumes

/******************************************************************************/
/*                                                                       MASS */
/******************************************************************************/
scalar total_mass = 0;
scalar dsl_total_mass = 0;
scalar total_entrained_mass = 0;
scalar max_alpha_inj = max(alpha_inj.boundaryField()[surface_patch_id]);
/******************************************************************************/
/*                                                                   VELOCITY */
/******************************************************************************/
scalar min_mag_u = min(mag(U)).value();
scalar max_mag_u = max(mag(U)).value();
scalar max_mag_alpha_inj_u =
    max(mag(alpha_inj_U.boundaryField()[surface_patch_id]));

std::map<std::string, double> max_mag_u_at_boundaries;
forAll(mesh.boundaryMesh(), patchi) {
  max_mag_u_at_boundaries[mesh.boundaryMesh()[patchi].name()] =
      max(mag(U.boundaryField()[patchi]));
}

/******************************************************************************/
/*                                                                     ENERGY */
/******************************************************************************/
scalar total_potential_energy = 0;
scalar total_kinetic_energy = 0;
scalar total_entrained_kinetic_energy = 0;
/******************************************************************************/
/*                                                                       SIZE */
/******************************************************************************/
// TODO

/******************************************************************************/
/*                                                        INTERNAL FIELD DATA */
/******************************************************************************/

forAll(mesh.V(), cell_id) {
  /****************************************************************************/
  // 	                                              					           mass
  /****************************************************************************/
  scalar cell_mass = alpha1[cell_id] * rho1.value() * mesh.V()[cell_id];
  total_mass += cell_mass;
  /****************************************************************************/
  // 	                                              						        energy
  /****************************************************************************/
  auto cell_velocity = mag(U[cell_id]);
  total_kinetic_energy += 0.5 * cell_mass * cell_velocity * cell_velocity;
  total_potential_energy += cell_mass * mesh.C()[cell_id].z() * 9.81;
}

/******************************************************************************/
//                                                         BOUNDARY FIELD DATA
/******************************************************************************/

forAll(U.boundaryField()[surface_patch_id], local_face_id) {
  // The owner cell id can be retrieved from the global face id
  auto global_face_id =
      local_face_id + mesh.boundaryMesh()[surface_patch_id].start();
  // The boundary cell id
  auto cell_id = mesh.owner()[global_face_id];
  // The local index space face id:
  auto closest_boundary_mesh_face = local_face_id;

  auto surface_area =
      mesh.magSf()
          .boundaryField()[surface_patch_id][closest_boundary_mesh_face];
  auto entrainment_velocity =
      mag(U.boundaryFieldRef()[surface_patch_id][closest_boundary_mesh_face]);
  /****************************************************************************/
  //         							entrainment mass
  /****************************************************************************/

  // dt Flux = dt * (rho * u) * (A * n)
  // all mutiplied by ce
  scalar entrained_mass =
      runTime.deltaT().value() * // dt
      alpha1.boundaryFieldRef()[surface_patch_id]
                               [closest_boundary_mesh_face] * // ce
      rho1.value() *                                          // rho
      entrainment_velocity *                                  // u
      surface_area;                                           // A
  // * mesh.V()[cell_id] *
  total_entrained_mass += entrained_mass;

  // dsl mass
  scalar dsl_mass =
      scalar(dslH.boundaryFieldRef()[surface_patch_id]
                                    [closest_boundary_mesh_face]) *
      surface_area * 1000;
  dsl_total_mass += dsl_mass;
  /****************************************************************************/
  //	      entrainment energy
  /****************************************************************************/
  total_entrained_kinetic_energy +=
      0.5 * entrained_mass * entrainment_velocity * entrainment_velocity;
}

/******************************************************************************/
//                                                                  OUTPUT
/******************************************************************************/
// append data to the respective processor profile file
if (Pstream::myProcNo() == 0) {
  // write data labels
  FILE *fp = fopen((profile_file + "_labels").c_str(), "w+");
  if (fp != nullptr) {
    fprintf(fp,
            "                          time "   // 0
            "                    total_mass "   // 1
            "          total_entrained_mass "   // 2
            "                     min_mag_u "   // 3
            "                     max_mag_u "   // 4
            "          total_kinetic_energy "   // 5
            "total_entrained_kinetic_energy "   // 6
            "        total_potential_energy "   // 7
            "                 max_alpha_inj "   // 8
            "           max_mag_alpha_inj_U "   // 9
            "                dsl_total_mass "); // 1
    // for(const auto& patch_name : max_mag_u_at_boundaries)
    // 	fprintf(fp, "%30s ", (patch_name.first + "_max_u").c_str());
    fclose(fp);
  }
}

FILE *fp =
    fopen((profile_file + std::to_string(Pstream::myProcNo())).c_str(), "a+");
if (fp != nullptr) {
  //                0      1      2      3      4      5      6      7      8 9
  //                10
  // fprintf(fp, "%26.4f %26.4f %26.4f %26.4f %26.4f %26.4f %26.4f %26.4f %26.4f
  // %26.4f % ",
  fprintf(fp, "%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f",
          runTime.value() - runTime.deltaT().value(), // 0
          total_mass,                                 // 1
          total_entrained_mass,                       // 2
          min_mag_u,                                  // 3
          max_mag_u,                                  // 4
          total_kinetic_energy,                       // 5
          total_entrained_kinetic_energy,             // 6
          total_potential_energy,                     // 7
          max_alpha_inj,                              // 8
          max_mag_alpha_inj_u,                        // 9
          dsl_total_mass);                            // 10
  // for(const auto& patch_max : max_mag_u_at_boundaries)
  // 	fprintf(fp, "%26.4f ", patch_max.second);
  fprintf(fp, "\n");
  fclose(fp);
}
