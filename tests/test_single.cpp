/** BSD 3-Clause License
    Copyright (c) 2023 Yingtian Chen
    All rights reserved.
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <dirent.h>
// #include <memory>

#include "orbit.h"
#include "units.h"
#include "smart.h"
#include "potential_multipole.h"
#include "potential_composite.h"
#include "potential_factory.h"
#include "debug_utils.h"
#include "math_spline.h"

#include "potential_varying.h"

// // load particles: 519311
// const int NUMPOINTS=5;
// const double posvel_car[NUMPOINTS][6] = {
//     {-22.398896, 5.244870, 10.457926, -44.358045, -125.459707, -27.125288},
//     {5.634460, -21.808264, 12.040468, -123.880807, -34.653069, -29.051565},
//     {-16.128471, -17.151681, 10.107484, -95.652157, 104.160769, 18.551908},
//     {12.114484, -19.822092, 10.342809, -117.011988, -52.317583, -52.191099},
//     {16.196563, 17.084948, 10.926397, -84.062415, 86.687878, 45.398332}};

// // load particles: Neil's halo
// const int NUMPOINTS=15;
// const double posvel_car[NUMPOINTS][6] = {
//     {-3.618246107853337 ,6.030165013813528 ,-5.833312459099942 ,-43.740066528320284 ,-49.51930236816405 ,111.4039611816406 ,},
//     {0.23067631920639536 ,-0.6538858508647535 ,-3.9184619067450503 ,19.576240539550795 ,-74.96347045898435 ,85.33363342285153 ,},
//     {-2.621206551975774 ,-0.9470247307381212 ,1.4946985525766647 ,-86.78588867187497 ,57.04235504567623 ,62.45512390136716 ,},
//     {-9.239878107597177 ,-5.734296824799457 ,-1.1699624932703043 ,-31.870979309082017 ,81.82342910766602 ,47.05665588378905 ,},
//     {-14.93694870875333 ,38.960647251624316 ,10.13558831964792 ,-56.02555847167966 ,34.49814796447754 ,-69.32206213474274 ,},
//     {-10.0414303958496 ,9.827418450346157 ,-28.120918311145935 ,3.8800811767578267 ,-0.3858489990234304 ,3.0782394409179545 ,},
//     {-7.699450018411879 ,-24.786599620036206 ,17.88296227863884 ,-1.034156799316392 ,1.7171287536621165 ,-6.3139915466308665 ,},
//     {-7.897986343083174 ,11.890601846114803 ,-29.251725468700435 ,-0.3485565185546733 ,1.8461761474609446 ,2.363563537597642 ,},
//     {8.44561705653541 ,-21.02171773149166 ,23.331297821285713 ,1.3712768554687642 ,6.179698944091804 ,-13.159801483154304 ,},
//     {-19.06114585747127 ,-19.41029601854643 ,15.87291669116075 ,6.121383666992202 ,5.2511253356933665 ,-7.52788162231446 ,},
//     {-9.608388308319261 ,24.390794105262096 ,-19.357639999901036 ,3.3097610473632955 ,0.19513320922852273 ,1.7011337280273295 ,},
//     {8.199841546833339 ,-7.5595045729514805 ,30.657973415727604 ,7.616592407226577 ,-6.750991821289055 ,8.565246582031236 ,},
//     {-6.08675735741417 ,15.03419188302723 ,-28.641345562157305 ,0.7259063720703267 ,-1.9004554748535085 ,3.1934967041015483 ,},
//     {18.3419055022714 ,21.23044509779356 ,15.09175071385425 ,-2.499053955078111 ,-7.3015670776367045 ,-10.033638000488288 ,},
//     {-2.064243723347317 ,-18.574536294572678 ,27.067532760775514 ,1.3982009887695455 ,0.34297943115235086 ,-13.582866668701179 ,},};


int main(int argc, char const *argv[])
{
    units::InternalUnits unit = units::InternalUnits(units::Kpc, units::Gyr);

    // load pots: 519311
    const double h100 = 1.0;
    int first_snap = 10;
    units::ExternalUnits extunit = units::ExternalUnits(unit, units::Kpc/h100, units::kms, units::Msun);
	std::vector<potential::PtrPotential> pots;
    std::string base_path = "/Users/ybchen/Documents/GC-model/GC-model/model_v1/results/pot_analogues/pot_total_id519311_snap_";
    for (int snap=first_snap; snap<=99; ++snap) {
        std::string file_path = base_path + std::to_string(snap) + "_align_at_last.pot";
        pots.push_back(potential::readPotential(file_path, extunit));
    }

    // // load pots: Neil's halo
    // const double h100 = 0.6774;
    // int first_snap = 1;
    // units::ExternalUnits extunit = units::ExternalUnits(unit, units::Kpc/h100, units::kms, units::Msun);
    // std::vector<potential::PtrPotential> pots;
    // std::string base_path = "/Users/ybchen/Downloads/GrNr_1425_coeffs_phys/snapshot_";
    // for (int snap=first_snap; snap<=99; ++snap) {
    //     std::string file_path = base_path + std::to_string(snap) + ".txt";
    //     pots.push_back(potential::readPotential(file_path, extunit));
    // }

    // load times
    std::vector<double> times;
    std::ifstream infile("/Users/ybchen/Documents/GC-model/GC-model/data/test/TNG_t_list.txt");
    double time;
    while (infile >> time) {
        times.push_back(time*unit.from_Gyr);
    }
    infile.close();
    times.erase(times.begin(), times.begin() + first_snap);


    // create time-varying pot
	// potential::PtrPotential potVSatic(new potential::PotentialVarying(times, pots));
    potential::PtrPotential potVSatic(new potential::Evolving(times, pots, true));

    // load acceleration
    std::vector<double> axs, ays, azs;
    double a;
    std::ifstream infilex("/Users/ybchen/Documents/GC-model/GC-model/model_v1/results/id519311_ax_mpb.txt");
    while (infilex >> a) {
        axs.push_back(-1.0*a*unit.from_kms/unit.from_Gyr);
    }
    infilex.close();
    std::ifstream infiley("/Users/ybchen/Documents/GC-model/GC-model/model_v1/results/id519311_ay_mpb.txt");
    while (infiley >> a) {
        ays.push_back(-1.0*a*unit.from_kms/unit.from_Gyr);
    }
    infiley.close();
    std::ifstream infilez("/Users/ybchen/Documents/GC-model/GC-model/model_v1/results/id519311_az_mpb.txt");
    while (infilez >> a) {
        azs.push_back(-1.0*a*unit.from_kms/unit.from_Gyr);
    }
    infilez.close();
    axs.erase(axs.begin(), axs.begin() + first_snap);
    ays.erase(ays.begin(), ays.begin() + first_snap);
    azs.erase(azs.begin(), azs.begin() + first_snap);

    // create time-varying acceleration
    potential::PtrPotential accField(new potential::UniformAcceleration(
        math::CubicSpline(times, axs), math::CubicSpline(times, ays), math::CubicSpline(times, azs)));

    std::vector<potential::PtrPotential> pots_with_acc = {potVSatic, accField};
    potential::PtrPotential potV_raw(new potential::Composite(pots_with_acc));

    // load offset
    std::ifstream infileoffset("/Users/ybchen/Documents/GC-model/GC-model/model_v1/results/hpos_relative_to_cm_id519311.txt");
    std::vector<double> dxs, dys, dzs;
    double temp_shap, dx, dy, dz;
    while (infileoffset >> temp_shap) {
        infileoffset >> dx; // dx = 0;
        dxs.push_back(1.0*dx*unit.from_Kpc);
        infileoffset >> dy; // dy = 0;
        dys.push_back(1.0*dy*unit.from_Kpc);
        infileoffset >> dz; // dz = 0;
        dzs.push_back(1.0*dz*unit.from_Kpc);
    }
    infileoffset.close();
    math::CubicSpline dx_sp = math::CubicSpline(times, dxs);
    math::CubicSpline dy_sp = math::CubicSpline(times, dys);
    math::CubicSpline dz_sp = math::CubicSpline(times, dzs);

    // create time-varying offset
    potential::PtrPotential potV(new potential::Shifted(potV_raw, dx_sp, dy_sp, dz_sp));

    // double potential;
    // coord::GradCyl deriv;
    // coord::HessCyl deriv2;
    // potV->eval(coord::PosCyl(50.*unit.from_Kpc,0.,0.), &potential, &deriv, &deriv2, times.back());
    // std::cout << deriv << std::endl;
    // potV->eval(coord::PosCyl(50.*unit.from_Kpc,0.,1.), &potential, &deriv, &deriv2, times.back());
    // std::cout << deriv << std::endl;

    // // load the particles:
    // // std::string filename = "/Users/ybchen/Documents/GC-model/GC-model/model_v1/results/selected_particles_id519311_snap_99_iom.txt";  // Replace with your file name
    // std::string filename = "/Users/ybchen/Documents/GC-model/GC-model/model_v1/results/pid4463640962_id519311_snap_99.txt";  // Replace with your file name
    // std::ifstream file(filename);
    // // First pass: count the number of lines in the file
    // std::string line;
    // int NUMPOINTS = 0;
    // while (std::getline(file, line)) {NUMPOINTS++;}
    // double (*posvel_car)[6] = new double[NUMPOINTS][6];
    // // Rewind the file to the beginning
    // file.clear();
    // file.seekg(0, std::ios::beg);
    // // Second pass: read the data
    // int i = 0;
    // while (std::getline(file, line)) {
    //     std::istringstream iss(line);
    //     for (int j = 0; j < 6; ++j) {
    //         if (!(iss >> posvel_car[i][j])) {
    //             std::cerr << "Error reading value." << std::endl;
    //             return 1;
    //         }
    //     }
    //     int dummy;  // For the seventh integer value
    //     iss >> dummy;
    //     i++;
    // }
    // file.close();

    // load the stream stars:
    std::string filename = "outputs/sample_gc/save_ic.txt";  // Replace with your file name
    std::ifstream file(filename);
    // First pass: count the number of lines in the file
    std::string line;
    int NUMPOINTS = 0;
    while (std::getline(file, line)) {NUMPOINTS++;}
    double (*posvel_car)[6] = new double[NUMPOINTS][6];
    std::vector<double> times_ic;
    double time_ic;
    // Rewind the file to the beginning
    file.clear();
    file.seekg(0, std::ios::beg);
    // Second pass: read the data
    int i = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> time_ic;
        times_ic.push_back(time_ic * unit.from_Gyr);
        for (int j = 0; j < 6; ++j) {
            if (!(iss >> posvel_car[i][j])) {
                std::cerr << "Error reading value." << std::endl;
                return 1;
            }
        }
        i++;
    }
    file.close();

    for (int ic=0; ic<NUMPOINTS; ++ic) {
        std::cout << ic << std::endl;

        // double total_time = -6 * unit.from_Gyr; // backward, for test particles
        // double timestep = -10 * unit.from_Myr;
        // double init_time = times.back();

        double init_time = times_ic[ic]; // forward, for gc streams
        double total_time = 14 * unit.from_Gyr - init_time;
        double timestep = 10 * unit.from_Myr;

        std::vector< std::pair<coord::PosVelCar, double> > traj;
        orbit::OrbitIntParams params(/*accuracy*/ 1e-10, /*maxNumSteps*/10000000);
        orbit::OrbitIntegrator<coord::Car> orbint(*potV, /*Omega*/0, params);
        // record the orbit at regular intervals of time
        orbint.addRuntimeFnc(orbit::PtrRuntimeFnc(new orbit::RuntimeTrajectory(
            orbint, timestep, /*output*/ traj)));
        // run the orbit
        double posvel_single[6] = {
            // (posvel_car[ic][0]*unit.from_Kpc + dx)/h100, // backward
            // (posvel_car[ic][1]*unit.from_Kpc + dy)/h100,
            // (posvel_car[ic][2]*unit.from_Kpc + dz)/h100,
            (posvel_car[ic][0]*unit.from_Kpc)/h100, // forward
            (posvel_car[ic][1]*unit.from_Kpc)/h100,
            (posvel_car[ic][2]*unit.from_Kpc)/h100,
            (posvel_car[ic][3]) * unit.from_kms,
            (posvel_car[ic][4]) * unit.from_kms,
            (posvel_car[ic][5]) * unit.from_kms};
        orbint.init(coord::PosVelCar(posvel_single), init_time);
        orbint.run(total_time);

        std::ofstream outfile;
        std::string file_path = "outputs/sample_gc/save" + std::to_string(ic) + ".txt";
        outfile.open(file_path, std::ios_base::trunc);
        for(size_t i=0; i<traj.size(); i++) {
            outfile << traj[i].second << "\t" << 
                traj[i].first.x*unit.to_Kpc*h100 << "\t" <<
                traj[i].first.y*unit.to_Kpc*h100 << "\t" <<
                traj[i].first.z*unit.to_Kpc*h100 << "\t" <<
                traj[i].first.vx*unit.to_kms << "\t" <<
                traj[i].first.vy*unit.to_kms << "\t" <<
                traj[i].first.vz*unit.to_kms << "\n"; 
        }
        outfile.close();
    }
    delete[] posvel_car;
    return 0;
}