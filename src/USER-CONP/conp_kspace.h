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

#ifndef LMP_CONP_KSPACE_H
#define LMP_CONP_KSPACE_H

#include "lmptype.h"

namespace LAMMPS_NS {
class ConpKspace {
 public:
  virtual void compute_vector(bigint *, double *) = 0;
  virtual void compute_vector_corr(bigint *, double *) = 0;
  virtual void compute_matrix(bigint *, double **) = 0;
  virtual void compute_matrix_corr(bigint *, double **) = 0;
};
}  // namespace LAMMPS_NS

#endif
