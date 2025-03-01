/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_BOUNDARY_CORRECTION_H
#define LMP_BOUNDARY_CORRECTION_H

#include "pointers.h"

namespace LAMMPS_NS {

class BoundaryCorrection : protected Pointers {
 public:
  BoundaryCorrection(LAMMPS *);
  virtual void vector_corr(bigint *, double *){};
  virtual void matrix_corr(bigint *, double **){};
  virtual void compute_corr(double, int, int, double &, double *){};
  void setup(double, double, double);
  void setup(double, double, double, double);

 protected:
  double area;
  double volume;
  double xprd_wire;
  double yprd_wire;
  double zprd_slab;
  double qqrd2e;
  double scale;
  double g_ewald;

  std::vector<bigint> gather_jmat(bigint *);
  std::vector<int> gather_recvcounts(int);
  std::vector<int> gather_displs(std::vector<int>);
};
}  // namespace LAMMPS_NS
#endif
