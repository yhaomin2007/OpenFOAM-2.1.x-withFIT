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

Class
    Foam::FITprop

Description
	call FIT property function    
	
SourceFiles
    FITprop.cpp

\*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <iostream>

extern"C" {void fittp_(double *temp,double *pres,double *enth,double *dens,double *cv,double *cp,double *cond,double *visc);
void fitph_(double *pres,double *enth,double *temp,double *dens,double *cv,double *cp,double *cond,double *visc,double *qual);}

inline void FITph_all(const double *pres_given,const double *enth_given,double *temp_need,double *dens_need,double *psi_need,double *alpha_need,double *visc_need, double *qual_need)
{	
	// alpha is thermal diffusivity
	double temp, pres, enth, dens, cv, cp,cond,visc,qual;
	
	enth = *enth_given/1000; // convert from j/kg to kj/kg
	pres = *pres_given/1000; // convert from Pa to kPa

	fitph_(&pres, &enth,&temp,&dens,&cv, &cp, &cond, &visc,&qual);

	*temp_need = temp;  // K
	*dens_need = dens;  // kg/m^3
	*psi_need  = dens/ *pres_given; // kg/m^3 / Pa
	*alpha_need = cond/(cp*1000);  // alpha is thermal diffusivity
	//*alpha_need = 1e-5;	
	*visc_need = visc;
	*qual_need = qual;

}

inline void FITtp_all(double *temp_given,double *pres_given,double *enth_need,double *dens_need,double *psi_need,double *alpha_need,double *visc_need)
{
	double temp, pres, enth, dens, cv, cp,cond,visc;
	
	temp = *temp_given;
	//enth = *enth_given/1000; // convert from j/kg to kj/kg
	pres = *pres_given/1000; // convert from Pa to kPa

	fittp_(&temp, &pres, &enth, &dens,&cv, &cp, &cond, &visc);

	*enth_need = enth*1000;
	*dens_need = dens;  // kg/m^3
	*psi_need  = dens/ *pres_given; // kg/m^3 / Pa
	*alpha_need = cond/(cp*1000); // alpha is thermal diffusivity
	//*alpha_need = 1e-5;	
	*visc_need = visc;

}

inline double FITtp_h(double temp_given,double pres_given)
{
	double temp, pres, enth, dens, cv, cp,cond,visc;
	
	temp = temp_given;
	pres = pres_given/1000; // convert from Pa to kPa

	fittp_(&temp, &pres, &enth, &dens,&cv, &cp, &cond, &visc);

	return enth*1000; // convert from kj/kg to j/kg
}

inline double FITph_cp(double pres_given,double enth_given)
{
	double temp, pres, enth, dens, cv, cp,cond,visc,qual;
	
	enth = enth_given/1000; // convert from j/kg to kj/kg
	pres = pres_given/1000; // convert from Pa to kPa

	fitph_(&pres, &enth,&temp,&dens,&cv, &cp, &cond, &visc,&qual);

	return cp*1000; // convert from kj/kg-K to j/kg-K
}

inline double FITph_cv(double pres_given,double enth_given)
{
	double temp, pres, enth, dens, cv, cp,cond,visc,qual;
	
	enth = enth_given/1000; // convert from j/kg to kj/kg
	pres = pres_given/1000; // convert from Pa to kPa

	fitph_(&pres, &enth,&temp,&dens,&cv, &cp, &cond, &visc,&qual);

	return cv*1000; // convert from kj/kg-K to j/kg-K
}

inline double FITph_kappa(double pres_given,double enth_given)
{
	double temp, pres, enth, dens, cv, cp,cond,visc,qual;
	
	enth = enth_given/1000; // convert from j/kg to kj/kg
	pres = pres_given/1000; // convert from Pa to kPa

	fitph_(&pres, &enth,&temp,&dens,&cv, &cp, &cond, &visc,&qual);

	return cond; 
}
