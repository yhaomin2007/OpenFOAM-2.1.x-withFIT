/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "hRhoThermo.H"
#include "FITprop.cpp"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MixtureType>
void Foam::hRhoThermo<MixtureType>::calculate()
{
    const scalarField& hCells = this->h_.internalField();
    const scalarField& pCells = this->p_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& rhoCells = this->rho_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();
    scalarField& xCells = this->x_.internalField();

    forAll(TCells, celli)
    {
        /*const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.TH(hCells[celli], TCells[celli]);
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(TCells[celli]);
        alphaCells[celli] = mixture_.alpha(TCells[celli]);*/

	FITph_all(&pCells[celli],&hCells[celli],&TCells[celli],&rhoCells[celli],&psiCells[celli],&alphaCells[celli],&muCells[celli],&xCells[celli]);

    }

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];
        fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];

        fvPatchScalarField& ph = this->h_.boundaryField()[patchi];

        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];
	fvPatchScalarField& px = this->x_.boundaryField()[patchi];

        /*if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.H(pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pT[facei]);
                palpha[facei] = mixture_.alpha(pT[facei]);

		//ph[facei] = FITtp_h(pT[facei],pp[facei]);
		//FITph_all(&pp[facei],&ph[facei],&pT[facei],&prho[facei],&ppsi[facei],&palpha[facei],&pmu[facei]);
		FITtp_all(&pT[facei],&pp[facei],&ph[facei],&prho[facei],&ppsi[facei],&palpha[facei],&pmu[facei]);  
		px[facei] = 1.0;      
	    }
        }
        else
        { */
            forAll(pT, facei)
            {   
                /*const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.TH(ph[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pT[facei]);
                palpha[facei] = mixture_.alpha(pT[facei]);*/

		FITph_all(&pp[facei],&ph[facei],&pT[facei],&prho[facei],&ppsi[facei],&palpha[facei],&pmu[facei],&px[facei]);
            }
        //}
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hRhoThermo<MixtureType>::hRhoThermo(const fvMesh& mesh)
:
    basicRhoThermo(mesh),
    MixtureType(*this, mesh),
    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    x_
    (
        IOobject
        (
            "x",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hRhoThermo<MixtureType>::~hRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::hRhoThermo<MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering hRhoThermo<MixtureType>::correct()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting hRhoThermo<MixtureType>::correct()" << endl;
    }
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hRhoThermo<MixtureType>::h
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    const scalarField& pCells = this->p_.internalField();

    forAll(T, celli)
    {
	h[celli] = FITtp_h(T[celli],pCells[celli]);  
        //h[celli] = this->cellMixture(cells[celli]).H(T[celli]);
    }

    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hRhoThermo<MixtureType>::h
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];

    forAll(T, facei)
    {
        h[facei] = FITtp_h(T[facei],pp[facei]);
        //h[facei] = this->patchFaceMixture(patchi, facei).H(T[facei]);
    }

    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hRhoThermo<MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
    const fvPatchScalarField& ph = this->h_.boundaryField()[patchi];
    
    forAll(T, facei)
    {
	cp[facei]=FITph_cp(pp[facei],ph[facei]);         
	//cp[facei] = this->patchFaceMixture(patchi, facei).Cp(T[facei]);
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::hRhoThermo<MixtureType>::Cp() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cp = tCp();
    
    const scalarField& hCells = this->h_.internalField();
    const scalarField& pCells = this->p_.internalField();
    
    forAll(this->T_, celli)
    {
        cp[celli]=FITph_cp(pCells[celli],hCells[celli]);
	//cp[celli] = this->cellMixture(celli).Cp(this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cp.boundaryField()[patchi];
	const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& ph = h_.boundaryField()[patchi];

        forAll(pT, facei)
        {
	    pCp[facei] = FITph_cp(pp[facei],ph[facei]);   
            //pCp[facei] = this->patchFaceMixture(patchi, facei).Cp(pT[facei]);
        }
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hRhoThermo<MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& cv = tCv();

    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
    const fvPatchScalarField& ph = this->h_.boundaryField()[patchi];

    forAll(T, facei)
    {
	cv[facei] = FITph_cv(pp[facei],ph[facei]);
        //cv[facei] = this->patchFaceMixture(patchi, facei).Cv(T[facei]);
    }

    return tCv;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::hRhoThermo<MixtureType>::Cv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cv = tCv();

    const scalarField& hCells = this->h_.internalField();
    const scalarField& pCells = this->p_.internalField();

    forAll(this->T_, celli)
    {
	cv[celli]= FITph_cv(pCells[celli],hCells[celli]);     
        //cv[celli] = this->cellMixture(celli).Cv(this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
	fvPatchScalarField& pCv = cv.boundaryField()[patchi];

	const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& ph = h_.boundaryField()[patchi];
	forAll(pp, facei)
        {
	pCv[facei]=FITph_cv(pp[facei],ph[facei]);
	}
        //cv.boundaryField()[patchi] =
        //    Cv(this->T_.boundaryField()[patchi], patchi);
    }

    return tCv;
}

// calculate thermal conductivity
template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::hRhoThermo<MixtureType>::kappa() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tKappa
    (
        new volScalarField
        (
            IOobject
            (
                "kappa",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/(dimTime*dimLength*dimTemperature) // may have problem here
        )
    );

    volScalarField& kappa = tKappa();
    
    const scalarField& hCells = this->h_.internalField();
    const scalarField& pCells = this->p_.internalField();
    
    forAll(this->T_, celli)
    {
        kappa[celli]=FITph_kappa(pCells[celli],hCells[celli]);
	//cp[celli] = this->cellMixture(celli).Cp(this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pKappa = kappa.boundaryField()[patchi];
	const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& ph = h_.boundaryField()[patchi];

        forAll(pT, facei)
        {
	    pKappa[facei] = FITph_kappa(pp[facei],ph[facei]);   
            //pCp[facei] = this->patchFaceMixture(patchi, facei).Cp(pT[facei]);
        }
    }

    return tKappa;
}

template<class MixtureType>
bool Foam::hRhoThermo<MixtureType>::read()
{
    if (basicRhoThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
