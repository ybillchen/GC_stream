/** BSD 3-Clause License
    Copyright (c) 2023 Yingtian Chen
    All rights reserved.
*/

#include "potential_varying.h"

namespace potential {

PotentialVarying::PotentialVarying(
	const std::vector<double> &times,
    const std::vector<PtrPotential> &pots) :
    time_seq(times), pot_seq(pots) {};


double PotentialVarying::enclosedMass(const double radius) const
    { return pot_seq.back()->enclosedMass(radius); }

void PotentialVarying::evalCyl(const coord::PosCyl &pos,
    double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2, 
    double time) const 
{

    // use the boundary values if out of time range
    if (time < time_seq.front())
        pot_seq.front()->eval(pos, potential, deriv, deriv2);
    if (time >= time_seq.back())
        pot_seq.back()->eval(pos, potential, deriv, deriv2);

    for (int i = 0; i < time_seq.size() - 1; ++i)
        if ((time >= time_seq[i]) && (time < time_seq[i+1]))
        {
    		double r; // ratio within time intervals, in [0, 1)
            r = (time - time_seq[i]) / (time_seq[i+1] - time_seq[i]);

            double potential_0, potential_1;
            coord::GradCyl deriv_0, deriv_1;
            coord::HessCyl deriv2_0, deriv2_1;

            pot_seq[i]->eval(pos, 
            	&potential_0, &deriv_0, &deriv2_0);
            pot_seq[i+1]->eval(pos, 
                &potential_1, &deriv_1, &deriv2_1);

            // linear interpolation
    		if(potential)
            	*potential = potential_0*(1-r) + potential_1*r;
		    if(deriv) 
		    {
		        deriv->dR = deriv_0.dR*(1-r) + deriv_1.dR*r;
		        deriv->dz = deriv_0.dz*(1-r) + deriv_1.dz*r;
		        deriv->dphi = deriv_0.dphi*(1-r) + deriv_1.dphi*r;
		    }
		    if(deriv2) 
		    {
		        deriv2->dR2 = deriv2_0.dR2*(1-r) + deriv2_1.dR2*r;
		        deriv2->dz2 = deriv2_0.dz2*(1-r) + deriv2_1.dz2*r;
		        deriv2->dphi2 = deriv2_0.dphi2*(1-r) + deriv2_1.dphi2*r;
		        deriv2->dRdz = deriv2_0.dRdz*(1-r) + deriv2_1.dRdz*r;
		        deriv2->dRdphi = deriv2_0.dRdphi*(1-r) + deriv2_1.dRdphi*r;
		        deriv2->dzdphi = deriv2_0.dzdphi*(1-r) + deriv2_1.dzdphi*r;
		    }

            break; // no need to continue
        }
}

double PotentialVarying::densityCyl(const coord::PosCyl &pos, 
	double time) const
{
    // use the boundary values if out of time range
    if (time < time_seq.front())
        return pot_seq.front()->density(pos);
    if (time >= time_seq.back())
        return pot_seq.back()->density(pos);

    for (int i = 0; i < time_seq.size() - 1; ++i)
        if ((time >= time_seq[i]) && (time < time_seq[i+1]))
        {
    		double r; // ratio within time intervals, in [0, 1)
            r = (time - time_seq[i]) / (time_seq[i+1] - time_seq[i]);

            double dens_0 = pot_seq[i]->density(pos);
            double dens_1 = pot_seq[i+1]->density(pos);

            return dens_0*r + dens_1*(1-r);
        }

    return 0; // this should never happen
}

}  // namespace