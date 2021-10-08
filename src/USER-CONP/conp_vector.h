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

#include "conp_kspace.h"
#include "pointers.h"

namespace LAMMPS_NS {

class ConpVector : protected Pointers {
 public:
  ConpVector(class LAMMPS *, int, double);
  ~ConpVector();
  void setup(std::vector<int>);
  void compute_vector();
  double *vector;
  int igroup;

 private:
  int groupbit;
  bigint ngroup;
  double **cutsq;
  double g_ewald, eta;
  std::vector<int> tag_to_iele;
  std::vector<bigint> mpos;
  class Pair *pair;
  // class NeighList *list;
  class ConpKspace *conp_kspace;

  void update_mpos();

  void pair_contribution();
  double calc_erfc(double);

  double setup_time_total;
  double reduce_time_total;
  double kspace_time_total;
  double pair_time_total;
  double boundary_time_total;
  double b_time_total;
  double alloc_time_total;
  double mpos_time_total;
};

}  // namespace LAMMPS_NS

