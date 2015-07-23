    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicRhoThermo> rThermo
    (
        basicRhoThermo::New(mesh)
    );
    basicRhoThermo& thermo = rThermo();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();

    Info<< "h max/min : " << max(h).value() << " " << min(h).value() << endl;
    const volScalarField& psi = thermo.psi();
    const volScalarField& T = thermo.T();
    const volScalarField& mu = thermo.mu();
    volScalarField Cp = thermo.Cp();
    Info<< "Cp max/min : " << max(Cp).value() << " " << min(Cp).value() << endl;
    volScalarField kappa = thermo.kappa();
    Info<< "kappa max/min : " << max(kappa).value() << " " << min(kappa).value() << endl;
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "compressibleCreatePhi.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);

    dimensionedScalar rhoMax(simple.dict().lookup("rhoMax"));
    dimensionedScalar rhoMin(simple.dict().lookup("rhoMin"));
    dimensionedScalar kappaMax(simple.dict().lookup("kappaMax"));
    dimensionedScalar kappaMin(simple.dict().lookup("kappaMin"));

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);
