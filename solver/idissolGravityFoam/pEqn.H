
volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(HbyA)
      +rAUf*fvc::ddtCorr(rEps, U, phiByEpsf)  
      + phig
    )
);

constrainPressure(p_rgh, U, phiHbyA, rAUf);

		volVectorField gradEps = fvc::grad(eps);
		surfaceVectorField gradEpsf = fvc::interpolate(gradEps);
		surfaceVectorField nEpsv = -gradEpsf/(mag(gradEpsf) + deltaN);
		nEpsf = nEpsv & mesh.Sf();
		surfaceScalarField phiR = reactionRate*fvc::interpolate(Yeq-Y)*nEpsf;
		surfaceScalarField phiREps =
            fvc::flux
            (
                phiR,
                eps,
                "div(phiR,eps)"
            );
		mdotY = eps*fvc::div(phiR) - fvc::div(phiREps);

	mdotSolid = mdotY*beta;
	rhof = rho0*rhok;
	while (pimple.correctNonOrthogonal())
	{
		fvScalarMatrix p_rghEqn
		(
		  fvm::laplacian(rAUf, p_rgh)
		==
		  fvc::div(phiHbyA)
		- fvm::Su(mdotSolid*(1.0/rhof-1.0/rhos),p_rgh)		
		);

		p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

		p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

		if (pimple.finalNonOrthogonalIter())
		{
			// Calculate the conservative fluxes
			phi = phiHbyA -  p_rghEqn.flux(); 

			// Explicitly relax pressure for momentum corrector
			p_rgh.relax();

			// Correct the momentum source with the pressure gradient flux
			// calculated from the relaxed pressure
			U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
			U.correctBoundaryConditions();
			fvOptions.correct(U);
		}
	}
    #include "continuityErrs.H"

    p = p_rgh + rhok*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }
