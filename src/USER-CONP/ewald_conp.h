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

#ifdef KSPACE_CLASS

KSpaceStyle(ewald/conp, EwaldConp)

#else

#ifndef LMP_EWALDCONP_H
#define LMP_EWALDCONP_H

#include "conp_kspace.h"
#include "ewald.h"

namespace LAMMPS_NS {

class EwaldConp : public Ewald, public ConpKspace {
 public:
  EwaldConp(class LAMMPS *);
  virtual ~EwaldConp();
  void init();
  void setup();
  virtual void compute(int, int);
  void compute_group_group(int, int, int);

  // k-space part of coulomb matrix computation
  void compute_vector(bigint *, double *);
  void compute_vector_corr(bigint *, double *);
  void compute_matrix(bigint *, double **);
  void compute_matrix_corr(bigint *, double **);

 protected:
  class BoundaryCorrection *boundcorr;
  double area;
  void coeffs();
  virtual void eik_dot_r();

 private:
  int eikr_step;
  void update_eikr(bool);
};

}  // namespace LAMMPS_NS

#endif
#endif

    /* ERROR/WARNING messages:

    E: Illegal ... command

    Self-explanatory.  Check the input script syntax and compare to the
    documentation for the command.  You can use -echo screen as a
    command-line option when running LAMMPS to see the offending line.

    E: Cannot use Ewald with 2d simulation

    The kspace style ewald cannot be used in 2d simulations.  You can use
    2d Ewald in a 3d simulation; see the kspace_modify command.

    E: Kspace style requires atom attribute q

    The atom style defined does not have these attributes.

    E: Cannot use non-periodic boundaries with Ewald

    For kspace style ewald, all 3 dimensions must have periodic boundaries
    unless you use the kspace_modify command to define a 2d slab with a
    non-periodic z dimension.

    E: Incorrect boundaries with slab Ewald

    Must have periodic x,y dimensions and non-periodic z dimension to use
    2d slab option with Ewald.

    E: Cannot (yet) use Ewald with triclinic box and slab correction

    This feature is not yet supported.

    E: KSpace style is incompatible with Pair style

    Setting a kspace style requires that a pair style with matching
    long-range Coulombic or dispersion components be used.

    E: KSpace accuracy must be > 0

    The kspace accuracy designated in the input must be greater than zero.

    E: Must use 'kspace_modify gewald' for uncharged system

    UNDOCUMENTED

    E: Cannot (yet) use K-space slab correction with compute group/group for
    triclinic systems

    This option is not yet supported.

    */
