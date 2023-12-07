/** BSD 3-Clause License
    Copyright (c) 2023 Yingtian Chen
    All rights reserved.

    The time varying potential interpolates the values at a sequence 
    of snapshots to get the values in between.

    Currently only implemented linear interpolation.
*/

#pragma once
#include <vector>

#include "smart.h"
#include "coord.h"
#include "potential_base.h"

namespace potential {

/// Time-varying potentials
class PotentialVarying: public BasePotentialCyl{
public:
    /** Construct the potential from a time sequence and a set of 
        correpsonding Potential pointers.
        \param[in]  times  is the time sequence (ascending).
        \param[in]  pots   is the correpsonding Potential pointers.
    */
    PotentialVarying(const std::vector<double> &times,
        const std::vector<PtrPotential> &pots);

    virtual coord::SymmetryType symmetry() const
        { return pot_seq.back()->symmetry(); }
    virtual const char* name() const
        { return myName(); }
    static const char* myName() 
        { static const char* my_name = "Varying"; return my_name; }

    virtual double enclosedMass(const double radius) const;

private:
    const std::vector<double> time_seq;
    const std::vector<PtrPotential> pot_seq;

    virtual void evalCyl(const coord::PosCyl &pos,
        double* potential, coord::GradCyl* deriv, 
        coord::HessCyl* deriv2, double time) const;

    virtual double densityCyl(const coord::PosCyl &pos, 
        double time) const;
};

}  // namespace