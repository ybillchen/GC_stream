/** BSD 3-Clause License
    Copyright (c) 2023 Yingtian Chen
    All rights reserved.
*/

#include <iostream>
#include <cmath>
#include <vector>
// #include <memory>

#include "units.h"
#include "smart.h"
#include "potential_analytic.h"

#include "potential_varying.h"

int main(int argc, char const *argv[])
{
    units::InternalUnits unit = units::InternalUnits(1e3*units::pc, 1e6*units::yr); // kpc and Myr

	std::vector<double> times = {0., 1.};

	potential::PtrPotential pot0(new potential::NFW(1e12*unit.from_Msun, 10.));
	potential::PtrPotential pot1(new potential::NFW(2e12*unit.from_Msun, 10.));

	std::vector<potential::PtrPotential> pots;
	pots.push_back(pot0);
	pots.push_back(pot1);

	potential::PtrPotential potV(new potential::PotentialVarying(times, pots));

    double energyUnit = unit.to_kms * unit.to_kms;
    double p0 = potV->value(coord::PosCar(0.,0.,0.), 0.) * energyUnit;
    double p1 = potV->value(coord::PosCar(0.,0.,0.), 1.) * energyUnit;
    double p_between = potV->value(coord::PosCar(0.,0.,0.), 0.5) * energyUnit;

    std::cout.precision(4);
    std::cout << "center potential at t = 0: " << p0 << " km2/s2" << std::endl;
    std::cout << "center potential at t = 0.5: " << p_between << " km2/s2" << std::endl;
    std::cout << "center potential at t = 1: " << p1 << " km2/s2" << std::endl;

    bool passed = (fabs(0.5*p0 + 0.5*p1 - p_between) / fabs(p_between) < 1e-6);

    if(passed)
        std::cout << "\033[1;32mALL TESTS PASSED\033[0m\n";
    else
        std::cout << "\033[1;31mSOME TESTS FAILED\033[0m\n";

    return 0;
}