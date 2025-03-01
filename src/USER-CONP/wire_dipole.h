/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "boundary_correction.h"

namespace LAMMPS_NS {

class WireDipole : public BoundaryCorrection {
 public:
  WireDipole(LAMMPS *);
  void vector_corr(bigint *, double *);
  void matrix_corr(bigint *, double **);
  void compute_corr(double, int, int, double &, double *);
  void setup(double);
};
}  // namespace LAMMPS_NS
