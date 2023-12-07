/** BSD 3-Clause License
    Copyright (c) 2023 Yingtian Chen
    All rights reserved.
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "orbit.h"
#include "potential_analytic.h"
#include "potential_utils.h"
#include "units.h"
#include "utils.h"
#include "debug_utils.h"
#include "math_random.h"

#include "potential_varying.h"

// const double epsE   = 1e-6;  // accuracy of energy conservation
const double epsL   = 1e-6;  // accuracy of AM conservation
const double epsCS  = 1e-3;  // accuracy of comparison of orbits in different coordinate systems
const double epsrot = 1e-3;  // accuracy of comparison between inertial and rotating frames
const double Omega  = 2.718; // rotation frequency (arbitrary)

const int NUMPOINTS=7;
const double posvel_car[NUMPOINTS][6] = {
    {0,-1e3, 2e3, 0, 200, -50},   // point in y-z plane
    {2e3, 0,-1e3, 0, -30, 40},   // point in x-z plane
    {0, 0, 1e3, 0, 0, 50},   // point along z axis and only vz!=0 (1d orbit)
    {0, 0, 1e3, 20, 30, 40},   // point along z axis, but with all three velocity components !=0
    {1e-3, 0,-1, 200, 1e-1,50},   // point nearly on z axis but with nonzero Lz
    {1e-3, 0, 0, 200, 30, 60},  // point at origin with nonzero velocity
    {10e3, 2e3, 3e3, 0, -100, -50}};   // ordinary point

inline double difposvel(const coord::PosVelCar& a, const coord::PosVelCar& b) {
    return sqrt(
        (pow_2(a.x-b.x) + pow_2(a.y-b.y) + pow_2(a.z-b.z)) / (pow_2(b.x) + pow_2(b.y) + pow_2(b.z)) +
        (pow_2(a.vx-b.vx) + pow_2(a.vy-b.vy) + pow_2(a.vz-b.vz)) / (pow_2(b.vx) + pow_2(b.vy) + pow_2(b.vz)));
}

template<typename CoordT>
bool test_coordsys(const potential::BasePotential& potential,
    const coord::PosVelCar& initial_conditions, double total_time,
    std::vector< std::pair<coord::PosVelCar, double> > &traj)
{
    std::string name = CoordT::name();
    name.resize(12, ' ');
    std::cout << name;
    double timestep = 2e-4 * total_time;
    double init_time = -0.01 * total_time;
    std::vector< std::pair<coord::PosVelCar, double> > trajFull, trajRot;
    orbit::OrbitIntParams params(/*accuracy*/ 1e-11, /*maxNumSteps*/10000000);
    orbit::OrbitIntegrator<CoordT> orbint(potential, /*Omega*/0, params);
    // record the orbit at regular intervals of time
    orbint.addRuntimeFnc(orbit::PtrRuntimeFnc(new orbit::RuntimeTrajectory(
        orbint, timestep, /*output*/ traj)));
    // record the orbit at each timestep of the ODE integrator
    orbint.addRuntimeFnc(orbit::PtrRuntimeFnc(new orbit::RuntimeTrajectory(
        orbint, 0, trajFull)));
    // run the orbit
    orbint.init(initial_conditions, init_time);
    orbint.run(total_time);
    // whether compare the result of orbit integration in rotating and inertial frames
    bool checkRot = isAxisymmetric(potential);
    if(checkRot) {
        orbit::OrbitIntegrator<CoordT> orbrot(potential, Omega, params);
        orbrot.addRuntimeFnc(orbit::PtrRuntimeFnc(new orbit::RuntimeTrajectory(
            orbrot, timestep, trajRot)));
        orbrot.init(initial_conditions, init_time);
        // orbrot.run(total_time/2);  // check that a continuation of the orbit integration
        // orbrot.run(total_time/2);  // gives the same result as a single run for a longer time
        orbrot.run(total_time);  // The check above is not necessary
        if(trajRot.size() != traj.size()) {
            std::cout << "\033[1;34mOrbit in the rotating frame has different length\033[0m: " << 
                trajRot.size() << " vs "<< traj.size() << "\n";
            return false;
        }
    }
    math::Averager avgE, avgL;
    double difrot = 0;
    for(size_t i=0; i<traj.size(); i++) {
        avgE.add(totalEnergy(potential, traj[i].first));
        avgL.add(Lz(traj[i].first));
        if(checkRot) {
            double angle = (traj[i].second-init_time)*Omega, cosa = cos(angle), sina = sin(angle);
            coord::PosVelCar& rxv = trajRot[i].first;  // cartesian coords in the rotating frame
            coord::PosVelCar ixv(  // converted to the inertial frame
                rxv.x  * cosa - rxv.y  * sina, rxv.y  * cosa + rxv.x  * sina, rxv.z,
                rxv.vx * cosa - rxv.vy * sina, rxv.vy * cosa + rxv.vx * sina, rxv.vz);
            difrot = fmax(difrot, difposvel(ixv, traj[i].first));
        }
    }
    bool ok = traj.back().second > 0.999999 * (total_time + init_time);
    if(ok)
        std::cout << trajFull.size()-1 << " steps" <<
            " steps at time " << trajFull.back().second;
    else {
        std::cout << "\033[1;31mCRASHED\033[0m after " << trajFull.size()-1 <<
            " steps at time " << trajFull.back().second;
    }
    std::cout << ", E=" << avgE.mean() << " +- " << sqrt(avgE.disp());
    // Energy conservation is not required for time-varying potential
    // if(avgE.disp() > pow_2(epsE*avgE.mean())) {
    //     std::cout << " \033[1;31m**\033[0m";
    //     ok = false;
    // }
    if(isAxisymmetric(potential)) {
        std::cout << ", Lz=" << avgL.mean() << " +- " << sqrt(avgL.disp());
        if(avgL.disp() > pow_2(epsL*avgL.mean()) && (fabs(Lz(initial_conditions)) < 1e-6 && fabs(avgL.mean()) > 1e-6)) {
            std::cout << " \033[1;35m**\033[0m";
            ok = false;
        }
    }
    if(checkRot) {
        std::cout << ", |rot-nonrot|=" << difrot;
        if(difrot > epsrot) {
            std::cout << " \033[1;34m**\033[0m";
            ok = false;
        }
    }
    std::cout << "\n";
    return ok;
}

bool test_potential(const potential::BasePotential& potential,
    const coord::PosVelCar& initial_conditions)
{
    // double total_time = 10.0 * T_circ(potential, totalEnergy(potential, initial_conditions));
    // if(!isFinite(total_time))
    //     total_time = 100.0;
    double total_time = 1e4; // 10 Gyr
    std::cout << "\033[1;37m" << potential.name() << "\033[0m, ic=(" << initial_conditions <<
        "), Torb=" << T_circ(potential, totalEnergy(potential, initial_conditions)) << "\n";

    std::vector< std::pair<coord::PosVelCar, double> > trajCar, trajCyl, trajSph;
    bool ok = true;
    ok &= test_coordsys<coord::Car>(potential, initial_conditions, total_time, trajCar);
    ok &= test_coordsys<coord::Cyl>(potential, initial_conditions, total_time, trajCyl);
    ok &= test_coordsys<coord::Sph>(potential, initial_conditions, total_time, trajSph);
    if(trajCar.size() == trajCyl.size() && trajCar.size() == trajSph.size()) {
        double maxdifcyl=0, maxdifsph=0;
        for(size_t i=0; i<trajCar.size(); i++) {
            maxdifcyl = fmax(maxdifcyl, difposvel(trajCar[i].first, trajCyl[i].first));
            maxdifsph = fmax(maxdifsph, difposvel(trajCar[i].first, trajSph[i].first));
        }
        std::cout << "|Car-Cyl|=" << maxdifcyl;
        if(maxdifcyl > epsCS) {
            std::cout << " \033[1;31m**\033[0m";
            ok = false;
        }
        std::cout << ", |Car-Sph|=" << maxdifsph;
        if(maxdifsph > epsCS) {
            std::cout << " \033[1;31m**\033[0m";
            ok = false;
        }
        std::cout << "\n";
    } else {
        std::cout << "\033[1;33mOrbit lengths differ between coordinate systems\033[0m\n";
        ok = false;
    }

    std::ofstream outfile;
    outfile.open("save.txt");//std::ios_base::app
    for(size_t i=0; i<trajCar.size(); i++) {
        outfile << trajCar[i].second << "\t" << 
            trajCar[i].first.x << "\t" <<
            trajCar[i].first.y << "\t" <<
            trajCar[i].first.z << "\n"; 
    }
    outfile.close();
    return ok;
}

int main() {
    units::InternalUnits unit = units::InternalUnits(units::pc, 1e6*units::yr); // pc and Myr, V ~ 1 km/s

    // std::vector<double> times = {0, 1e3, 1e4};

    // potential::PtrPotential pot0(new potential::NFW(1e11*unit.from_Msun, 1e4));
    // potential::PtrPotential pot1(new potential::NFW(1e12*unit.from_Msun, 1e4));
    // potential::PtrPotential pot2(new potential::NFW(2e12*unit.from_Msun, 1e4));

    std::vector<double> times = {0, 9e3, 1e4};

    potential::PtrPotential pot0(new potential::NFW(2e12*unit.from_Msun, 2e3));
    potential::PtrPotential pot1(new potential::NFW(1e12*unit.from_Msun, 8e3));
    potential::PtrPotential pot2(new potential::NFW(1e11*unit.from_Msun, 1e4));

    std::vector<potential::PtrPotential> pots;
    pots.push_back(pot0);
    pots.push_back(pot1);
    pots.push_back(pot2);

    potential::PtrPotential potV(new potential::PotentialVarying(times, pots));

    bool allok = true;

    for(int ic=0; ic<NUMPOINTS; ic++)
        allok &= test_potential(*potV, coord::PosVelCar(posvel_car[ic]));

    if(allok)
        std::cout << "\033[1;32mALL TESTS PASSED\033[0m\n";
    else
        std::cout << "\033[1;31mSOME TESTS FAILED\033[0m\n";

    return 0;
}