fvVectorMatrix UEqn(fvm::ddt(rho, U) + fvm::div(rhoPhi, U) -
                    fvm::laplacian(rho * turbulence->nuEff(), U));

UEqn.relax();

if (pimple.momentumPredictor()) {
  solve(UEqn ==
        fvc::reconstruct((-ghf * fvc::snGrad(rho) - fvc::snGrad(p_rgh)) *
                         mesh.magSf()));
}
