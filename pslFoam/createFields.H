Info << "Reading field p_rgh\n" << endl;
volScalarField p_rgh(IOobject("p_rgh", runTime.timeName(), mesh,
                              IOobject::MUST_READ, IOobject::AUTO_WRITE),
                     mesh);

Info << "Reading field U\n" << endl;
volVectorField U(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ,
                          IOobject::AUTO_WRITE),
                 mesh);

// #include "createPhi.H"
Info << "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi(IOobject("phi", runTime.timeName(), mesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::AUTO_WRITE),
                       fvc::flux(U));

Info << "Reading transportProperties\n" << endl;
incompressibleTwoPhaseMixture mixture(U, phi);

volScalarField &alpha1(mixture.alpha1());
volScalarField &alpha2(mixture.alpha2());

const dimensionedScalar &rho1 = mixture.rho1();
const dimensionedScalar &rho2 = mixture.rho2();

dimensionedScalar Dab("Dab", dimViscosity, mixture);

// Read the reciprocal of the turbulent Schmidt number
dimensionedScalar alphatab("alphatab", dimless, mixture);

// Need to store rho for ddt(rho, U)
volScalarField rho("rho", alpha1 *rho1 + alpha2 * rho2);
rho.oldTime();

// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi(IOobject("rhoPhi", runTime.timeName(), mesh,
                                   IOobject::NO_READ, IOobject::AUTO_WRITE),
                          rho1 *phi);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel>
    turbulence(incompressible::turbulenceModel::New(U, phi, mixture));

#include "readGravitationalAcceleration.H"
#include "readhRef.H"

#include "gh.H"

volScalarField p(IOobject("p", runTime.timeName(), mesh, IOobject::NO_READ,
                          IOobject::AUTO_WRITE),
                 p_rgh + rho * gh);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, p_rgh, pimple.dict(), pRefCell, pRefValue);

if (p_rgh.needReference()) {
  p += dimensionedScalar("p", p.dimensions(),
                         pRefValue - getRefCellValue(p, pRefCell));
  p_rgh = p - rho * gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// create surface field to store dsl height
volScalarField dslH(IOobject("dsl_h", runTime.timeName(), mesh,
                             IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                    mesh);

// create surface vector field to store dsl velocity
volVectorField dslU(IOobject("dsl_U", runTime.timeName(), mesh,
                             IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                    mesh);

// create surface field to store dsl W
volScalarField dslW(IOobject("dsl_W", runTime.timeName(), mesh,
                             IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                    mesh);

volScalarField dslDivU(IOobject("dsl_divU", runTime.timeName(), mesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::AUTO_WRITE),
                       mesh);

volVectorField Uinj(IOobject("Uinj", runTime.timeName(), mesh,
                             IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                    mesh);

// alpha gradient
volVectorField gradAlpha1(IOobject("gradAlpha1", runTime.timeName(), mesh,
                                   IOobject::NO_READ, IOobject::AUTO_WRITE),
                          mesh);

// create front distance field
volScalarField dslFrontSDF(IOobject("dsl_front_sdf", runTime.timeName(), mesh,
                                    IOobject::READ_IF_PRESENT,
                                    IOobject::AUTO_WRITE),
                           mesh);

// create terrain slope field
volScalarField dslSlope(IOobject("dsl_slope", runTime.timeName(), mesh,
                                 IOobject::READ_IF_PRESENT,
                                 IOobject::AUTO_WRITE),
                        mesh);

// create mass injection field
volScalarField alpha_inj(IOobject("alpha_inj", runTime.timeName(), mesh,
                                  IOobject::NO_READ, IOobject::AUTO_WRITE),
                         mesh);

// create mass injection velocity field
volScalarField alpha_inj_U(IOobject("alpha_inj_U", runTime.timeName(), mesh,
                                    IOobject::NO_READ, IOobject::AUTO_WRITE),
                           mesh);
