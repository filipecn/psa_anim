/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    twoLiquidMixingFoam

Group
    grpMultiphaseSolvers

Description
    Solver for mixing two incompressible fluids.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "MULES.H"

#include "subCycle.H"

#include "incompressibleTwoPhaseMixture.H"

#include "turbulentTransportModel.H"

#include "pimpleControl.H"

#include "cellDistFuncs.H"

#include "fvMesh.H"

#include "noise.H"

#include "stopWatch.H"

#include <cmath>
// #include <limits>
// #include <queue>
// #include <string>
// #include <unordered_map>

// profiling objects
#define MAX_PROCESS_COUNT 30
StopWatch totalTime[MAX_PROCESS_COUNT];
StopWatch initTime[MAX_PROCESS_COUNT];
StopWatch mainLoopTime[MAX_PROCESS_COUNT];
StopWatch mainLoopIterationTime[MAX_PROCESS_COUNT];
StopWatch PIMPLELoopTime[MAX_PROCESS_COUNT];
StopWatch PISOLoopTime[MAX_PROCESS_COUNT];
StopWatch boundaryUpdateTime[MAX_PROCESS_COUNT];
StopWatch alphaUpdateTime[MAX_PROCESS_COUNT];
StopWatch momentumPredictionTime[MAX_PROCESS_COUNT];
StopWatch pressureEquationTime[MAX_PROCESS_COUNT];
StopWatch pressureSetupTime[MAX_PROCESS_COUNT];
StopWatch correctNonOrthogonalLoopTime[MAX_PROCESS_COUNT];
StopWatch correctNonOrthogonalLoopIterationTime[MAX_PROCESS_COUNT];
StopWatch constrainPressureTime[MAX_PROCESS_COUNT];
StopWatch adjustPhiTime[MAX_PROCESS_COUNT];
StopWatch pressureSolveTime[MAX_PROCESS_COUNT];

#define PSA_ANIM_TICK(WATCH) WATCH[Pstream::myProcNo()].start()
#define PSA_ANIM_TACK(WATCH) WATCH[Pstream::myProcNo()].stop()
#define PSA_ANIM_TIME(WATCH) WATCH[Pstream::myProcNo()].getTotalTime()
#define PSA_ANIM_COUNT(WATCH) WATCH[Pstream::myProcNo()].getCount()

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[]) {

  argList::addNote("Solver for mixing two incompressible fluids");

  Info << "Initing PSL simulation in OpenFOAM" << endl;

  PSA_ANIM_TICK(totalTime);
  PSA_ANIM_TICK(initTime);

#include "postProcess.H"

#include "addCheckCaseOptions.H"

#include "setRootCaseLists.H"

#include "createTime.H"

#include "createMesh.H"

#include "createControl.H"

#include "initContinuityErrs.H"

// read terrain_patch_name
#include "pslControls.H"

// create extra fields
// Input fields contain DSL data:
// dsl_h
// dsl_U
// dsl_front_sdf
// dsl_slope
// Output fields contain PSL data:
// alpha_inj
// alpha_inj_U
#include "createFields.H"

  // create snow cover height fields
  std::vector<double> snow_cover_height, acc_entrained_height;
  if (mesh.boundaryMesh()[surface_patch_id].size()) {
    snow_cover_height.resize(mesh.boundaryMesh()[surface_patch_id].size(), 20);
    acc_entrained_height.resize(mesh.boundaryMesh()[surface_patch_id].size(),
                                0);
  }

#include "createTimeControls.H"

#include "CourantNo.H"

#include "setInitialDeltaT.H"

  turbulence->validate();

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info << "\nStarting time loop\n" << endl;

  // Initialize DSL input fields
  dslH.correctBoundaryConditions();
  dslU.correctBoundaryConditions();
  dslW.correctBoundaryConditions();
  dslFrontSDF.correctBoundaryConditions();
  dslSlope.correctBoundaryConditions();

  PSA_ANIM_TACK(initTime);
  PSA_ANIM_TICK(mainLoopTime);
  while (runTime.run()) {

    PSA_ANIM_TICK(mainLoopIterationTime);
    PSA_ANIM_TICK(boundaryUpdateTime);

    // force boundary condition values update
    // to update DSL data
    dslH.correctBoundaryConditions();
    dslU.correctBoundaryConditions();
    dslW.correctBoundaryConditions();
    dslFrontSDF.correctBoundaryConditions();
    dslSlope.correctBoundaryConditions();
// compute mass injection
#include "psl.H"

    PSA_ANIM_TACK(boundaryUpdateTime);

#include "readTimeControls.H"

#include "CourantNo.H"

#include "alphaCourantNo.H"

#include "setDeltaT.H"

    ++runTime;

    PSA_ANIM_TICK(alphaUpdateTime);
    mixture.correct();

#include "alphaEqnSubCycle.H"

#include "alphaDiffusionEqn.H"
    PSA_ANIM_TACK(alphaUpdateTime);

    // compute grad alpha field
    gradAlpha1 = fvc::grad(alpha1);

    // --- Pressure-velocity PIMPLE corrector loop
    PSA_ANIM_TICK(PIMPLELoopTime);
    while (pimple.loop()) {
      PSA_ANIM_TICK(momentumPredictionTime);
#include "UEqn.H"
      PSA_ANIM_TACK(momentumPredictionTime);

      // --- Pressure corrector loop
      PSA_ANIM_TICK(PISOLoopTime);
      while (pimple.correct()) {
        PSA_ANIM_TICK(pressureEquationTime);
#include "pEqn.H"
        PSA_ANIM_TACK(pressureEquationTime);
      }
      PSA_ANIM_TACK(PISOLoopTime);

      if (pimple.turbCorr()) {
        turbulence->correct();
      }
    }
    PSA_ANIM_TACK(PIMPLELoopTime);
    PSA_ANIM_TACK(mainLoopIterationTime);

    if (do_profile) {
#include "profile.H"
    }

    runTime.write();

    runTime.printExecutionTime(Info);
  }
  PSA_ANIM_TACK(mainLoopTime);
  PSA_ANIM_TACK(totalTime);

  FILE *fp = fopen(
      ("profile_times_" + std::to_string(Pstream::myProcNo())).c_str(), "a+");
  if (fp != nullptr) {
#define PRINT_TIME(WATCH)                                                      \
  fprintf(fp, "%s %lf %d\n", #WATCH, PSA_ANIM_TIME(WATCH),                     \
          PSA_ANIM_COUNT(WATCH))

    PRINT_TIME(totalTime);
    PRINT_TIME(initTime);
    PRINT_TIME(mainLoopTime);
    PRINT_TIME(mainLoopIterationTime);
    PRINT_TIME(PIMPLELoopTime);
    PRINT_TIME(PISOLoopTime);
    PRINT_TIME(boundaryUpdateTime);
    PRINT_TIME(alphaUpdateTime);
    PRINT_TIME(momentumPredictionTime);
    PRINT_TIME(pressureEquationTime);
    PRINT_TIME(pressureSetupTime);
    PRINT_TIME(correctNonOrthogonalLoopTime);
    PRINT_TIME(correctNonOrthogonalLoopIterationTime);
    PRINT_TIME(constrainPressureTime);
    PRINT_TIME(adjustPhiTime);
    PRINT_TIME(pressureSolveTime);

    fclose(fp);
  }

  Info << "End\n" << endl;

  return 0;
}

// ************************************************************************* //
