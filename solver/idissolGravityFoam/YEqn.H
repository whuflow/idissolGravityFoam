{
    Eps2Df = eps*eps*Df;

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

    fvScalarMatrix YEqn
    (
		  fvm::ddt(eps, Y)
        + fvm::div(phi,Y)
        - fvm::laplacian(Eps2Df,Y)
        ==
		  fvm::Su(mdotY,Y)
        + fvOptions(Y)
    );

    YEqn.relax();

    fvOptions.constrain(YEqn);

    YEqn.solve();

    fvOptions.correct(Y);
    Y.max(0.0);
    Y.min(Yeq);	
	
    rhok = 1 + beta0*(Y-Y0);
	
}