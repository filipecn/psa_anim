{
  PSA_ANIM_TICK(pressureSetupTime);
  volScalarField rAU("rAU", 1.0 / UEqn.A());
  surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
  volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p_rgh));
  surfaceScalarField phiHbyA("phiHbyA",
                             fvc::flux(HbyA) + fvc::interpolate(rho * rAU) *
                                                   fvc::ddtCorr(U, phi));
  PSA_ANIM_TACK(pressureSetupTime);
  PSA_ANIM_TICK(adjustPhiTime);
  adjustPhi(phiHbyA, U, p_rgh);
  PSA_ANIM_TACK(adjustPhiTime);
  PSA_ANIM_TICK(constrainPressureTime);
  surfaceScalarField phig(-ghf * fvc::snGrad(rho) * rAUf * mesh.magSf());

  phiHbyA += phig;

  // Update the pressure BCs to ensure flux consistency
  constrainPressure(p_rgh, U, phiHbyA, rAUf);
  PSA_ANIM_TACK(constrainPressureTime);

  PSA_ANIM_TICK(correctNonOrthogonalLoopTime);
  while (pimple.correctNonOrthogonal()) {
    PSA_ANIM_TICK(correctNonOrthogonalLoopIterationTime);
    fvScalarMatrix p_rghEqn(fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA));

    p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

    PSA_ANIM_TICK(pressureSolveTime);
    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
    PSA_ANIM_TACK(pressureSolveTime);
    if (pimple.finalNonOrthogonalIter()) {
      phi = phiHbyA - p_rghEqn.flux();

      U = HbyA + rAU * fvc::reconstruct((phig - p_rghEqn.flux()) / rAUf);
      U.correctBoundaryConditions();
    }
    PSA_ANIM_TACK(correctNonOrthogonalLoopIterationTime);
  }
  PSA_ANIM_TACK(correctNonOrthogonalLoopTime);

#include "continuityErrs.H"

  p == p_rgh + rho *gh;

  if (p_rgh.needReference()) {
    p += dimensionedScalar("p", p.dimensions(),
                           pRefValue - getRefCellValue(p, pRefCell));
    p_rgh = p - rho * gh;
  }
}
