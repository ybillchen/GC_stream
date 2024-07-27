"""
BSD 3-Clause License
Copyright (c) 2024 Yingtian Chen
All rights reserved.
"""

"""
Excecude the code with:

$ python main.py params.json

"""

import sys
import json

import numpy as np

import agama
agama.setUnits(mass=1, length=1, velocity=1) # Msun, kpc, km/s

from pot_evolving import Evolving, M17Pot
from streams_generator import Orbit, MassHistory, \
    F15StreamsGenerator, SphStreamsGenerator, R24StreamsGenerator

def main(params):
    # random number generator
    rng = np.random.default_rng(params["seed"])

    # load potential
    if params["pot_type"] == "M17": # McMillan+2017
        pot = M17Pot()
    elif params["pot_type"] == "Evolving":
        pot_base_path = params["pot_base_path"]
        pots = []
        for snap in range(10,100):
            pots.append(agama.Potential(
                pot_base_path + "%d_align_at_last.pot"%snap))
        times = np.loadtxt(params["times_path"])[10:]
        pot = Evolving(pots, times)
    else:
        print("Potential type",
            "%s not recogonized."%params["pot_type"], 
            "Using default: M17")
        pot = M17Pot()

    if params["streams_generator"] == "F15":
        Rapo = params["Rapo"]
        Rperi = params["Rperi"]
        ft = params["ft"]
        gala = bool(params["gala"])
        generator = F15StreamsGenerator(pot, Rapo, Rperi, ft, gala)
    elif params["streams_generator"] == "R24":
        fe = np.array(params["fe"])
        eps = np.array(params["eps"])
        generator = R24StreamsGenerator(pot, fe, eps)
    elif params["streams_generator"] == "Sph":
        mean = np.array(params["mean"])
        cov = np.array(params["cov"])
        orbit_dependent = bool(params["orbit_dependent"])
        generator = SphStreamsGenerator(pot, mean, cov, orbit_dependent)
    else:
        print("Streams generator",
            "%s not recogonized."%params["streams_generator"], 
            "Using default: F15")
        generator = F15StreamsGenerator(pot)

    # load gc orbit
    gc_orbit = np.loadtxt(params["gc_orbit_path"])
    gc_orbit = Orbit(gc_orbit[:,1:7], gc_orbit[:,0])

    # load mass history
    mass_history = np.loadtxt(params["mass_history_path"])
    mass_history = MassHistory(mass_history[:,1], mass_history[:,0])

    time_sample, posvel_ej = generator.sample_stream(
        params["t_begin_in_Gyr"], params["t_end_in_Gyr"], 
        gc_orbit, mass_history, params["dm_in_Msun"], rng)

    out = np.column_stack((time_sample, posvel_ej))
    np.savetxt(params["save_path_base"], out)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Params filename not provided.",
            "Using default: params.json")
        filename = "params.json"
    else:
        filename = sys.argv[1]

    with open(filename) as f:
        params = json.load(f)

    main(params)
