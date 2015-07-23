{

    //Cp = thermo.Cp();
    //Info<< "Cp max/min : " << max(Cp).value() << " " << min(Cp).value() << endl;
    
    //kappa = thermo.kappa();
    //Info<< "kappa max/min : " << max(kappa).value() << " " << min(kappa).value() << endl;

    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
       //-fvm::laplacian(mu_h, h)      
      - fvm::laplacian(turbulence->alphat(), h)
     ==
      - fvc::div(phi, 0.5*magSqr(U), "div(phi,K)")
      + fvc::laplacian(kappa, T)
    );

    hEqn.relax();

    hEqn.solve();

    h.relax();

    //h.correctBoundaryConditions();

    Info<< "h max/min : " << max(h).value() << " " << min(h).value() << endl;
    Info<< "T max/min : " << max(T).value() << " " << min(T).value() << endl;   
    Info<< "mu max/min : " << max(mu).value() << " " << min(mu).value() << endl; 

    thermo.correct();
    rho = thermo.rho();
    rho = max(rho, rhoMin);   
    rho = min(rho, rhoMax);
    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

    kappa = thermo.kappa();
    kappa = max(kappa, kappaMin);
    kappa = min(kappa, kappaMax);
    Info<< "kappa max/min : " << max(kappa).value() << " " << min(kappa).value() << endl;
    //rho = thermo.rho();
    //rho = max(rho, rhoMin);
    //rho = min(rho, rhoMax);
}
