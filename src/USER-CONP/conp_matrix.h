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

class ConpMatrix : protected Pointers {
 public:
  ConpMatrix(class LAMMPS *, int, double);
  ~ConpMatrix();
  void setup(std::vector<int>);
  void compute_array();
  double **array;
  int igroup;

 private:
  int  groupbit;
  bigint ngroup;
  double **cutsq;
  double g_ewald, eta;
  std::vector<int> tag_to_iele;
  bool assigned;
  std::vector<bigint> mpos;
  class Pair *pair;
  // class NeighList *list;
  class ConpKspace *conp_kspace;

  void update_mpos();
  void pair_contribution();
  void self_contribution();
  double calc_erfc(double);
};

}  // namespace LAMMPS_NS

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute group/group group ID does not exist

Self-explanatory.

E: Compute group/group molecule requires molecule IDs

UNDOCUMENTED

E: No pair style defined for compute group/group

Cannot calculate group interactions without a pair style defined.

E: Pair style does not support compute group/group

The pair_style does not have a single() function, so it cannot be
invoked by the compute group/group command.

E: No Kspace style defined for compute group/group

Self-explanatory.

E: Kspace style does not support compute group/group

Self-explanatory.

W: Both groups in compute group/group have a net charge; the Kspace boundary
correction to energy will be non-zero

Self-explanatory.

*/
