# GC_stream

## Intro

## Prerequisite

`python >= 3.6`
```
numpy
h5py
agama
```

## Related papers

### Potential models

- doi:10.3847/1538-4357/ac93fb (Arora et al. 2022, On the Stability of Tidal Streams in Action Space)
- doi:10.1051/0004-6361/202245117 (Ishchenko et al. 2023, TNG-TVP)

### Stream models

- doi:10.1093/mnras/stu1986 (Gibbons et al. 2014, mLCS model)
- doi:10.1088/0004-637X/795/1/94 (Bonaca et al. 2014, stream model)
- doi:10.1093/mnras/stv1198 (Fardal et al. 2015, stream model)

## Pipeline

1. Backward orbit integration of GCs (survived and dead) until
  1. Some fixed end time such as t_lb = 10 Gyr
  2. The formation time of clusters
2. Forward evolving
  1. Distribute tracer particles mass according to GC mass history
  2. Distribute tracer particles positions and velocities according to some stream sampling algorithm 
  3. Intergrate the orbits of tracers until t_lb = 0

Input data structure:
- A table for time-dependent galaxy potential parameters
- A catalogue of GC present-day positions and velocities
- A table for each GC's mass history

Output data structure
- GC
  - GC orbits from backward intergation
- Tracer
  - Tracer orbits from forward intergation
  - Distribution time of tracers
  - Host GC ID of tracers

## Usage

## Author

- [@Yingtian (Bill) Chen](https://github.com/ybillchen)