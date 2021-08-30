/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
     K-space terms added by Stan Moore (BYU)
------------------------------------------------------------------------- */

#include "conp_matrix.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "neigh_list.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

/* ---------------------------------------------------------------------- */

ConpMatrix::ConpMatrix(LAMMPS *lmp, int electrode_group, double eta)
    : Pointers(lmp) {
  igroup = electrode_group;  // group of all electrode atoms
  groupbit = group->bitmask[igroup];
  ngroup = group->count(igroup);
  array = new double *[ngroup];
  for (bigint i = 0; i < ngroup; i++) array[i] = new double[ngroup];
  this->eta = eta;
}

/* ---------------------------------------------------------------------- */

ConpMatrix::~ConpMatrix() {
  for (bigint i = 0; i < ngroup; i++) delete[] array[i];
  delete[] array;
}

/* ---------------------------------------------------------------------- */

void ConpMatrix::setup(std::vector<int> tag_ids) {
  if (force->kspace == nullptr)
    error->all(FLERR, "No Kspace style defined for compute conp matrix");

  // check if coul pair style is active, no need for single() since done
  // explicitly
  int itmp;
  double *p_cutoff = (double *)force->pair->extract("cut_coul", itmp);
  if (p_cutoff == nullptr)
    error->all(FLERR, "conp matrix is incompatible with Pair style");
  pair = force->pair;
  cutsq = force->pair->cutsq;

  conp_kspace = dynamic_cast<ConpKspace *>(force->kspace);
  if (conp_kspace == nullptr)
    error->all(FLERR, "Kspace does not implement ConpKspace");
  g_ewald = force->kspace->g_ewald;

  tag_to_iele = tag_ids;
}

/* ---------------------------------------------------------------------- */

void ConpMatrix::compute_array() {
  // setting all entries of coulomb matrix to zero
  size_t nbytes = sizeof(double) * ngroup;
  if (nbytes)
    for (int i = 0; i < ngroup; i++) memset(&array[i][0], 0, nbytes);

  MPI_Barrier(world);
  double kspace_time = MPI_Wtime();
  update_mpos();
  conp_kspace->compute_matrix(&mpos[0], array);
  MPI_Barrier(world);
  if (comm->me == 0)
    utils::logmesg(lmp,
                   fmt::format("Kspace time: {}\n", MPI_Wtime() - kspace_time));
  pair_contribution();
  self_contribution();
  conp_kspace->compute_matrix_corr(&mpos[0], array);

  // reduce coulomb matrix with contributions from all procs
  // all procs need to know full matrix for matrix inversion
  for (int i = 0; i < ngroup; i++) {
    MPI_Allreduce(MPI_IN_PLACE, &array[i][0], ngroup, MPI_DOUBLE, MPI_SUM,
                  world);
  }
}

/* ---------------------------------------------------------------------- */

void ConpMatrix::pair_contribution() {
  int inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double r, rinv, rsq, grij, etarij, expm2, t, erfc, aij;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  double etaij =
      eta * eta / sqrt(2.0 * eta * eta);  // see mw ewald theory eq. (29)-(30)

  // neighbor list will be ready because called from post_neighbor
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I,J are not in 2 groups

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    // skip if atom I is not in either group
    if (!(mask[i] & groupbit)) continue;

    bigint const ipos = mpos[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // real-space part of matrix is symmetric
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];  // neighlists take care of pbc
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1.0 / r;
        aij = rinv;

        grij = g_ewald * r;
        expm2 = exp(-grij * grij);
        t = 1.0 / (1.0 + EWALD_P * grij);
        erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
        aij *= erfc;

        etarij = etaij * r;
        expm2 = exp(-etarij * etarij);
        t = 1.0 / (1.0 + EWALD_P * etarij);
        erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
        aij -= erfc * rinv;

        // newton on or off?
        if (!(newton_pair || j < nlocal)) aij *= 0.5;
        bigint jpos = tag_to_iele[tag[j]];
        array[ipos][jpos] += aij;
        array[jpos][ipos] += aij;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ConpMatrix::self_contribution() {
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  const double selfint = 2.0 / MY_PIS * g_ewald;
  const double preta = MY_SQRT2 / MY_PIS;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      array[mpos[i]][mpos[i]] += preta * eta - selfint;
    }
}

/* ---------------------------------------------------------------------- */

void ConpMatrix::update_mpos() {
  int const nall = atom->nlocal + atom->nghost;
  int *tag = atom->tag;
  int *mask = atom->mask;
  mpos = std::vector<bigint>(nall, -1);

  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit)
      mpos[i] = tag_to_iele[tag[i]];
    else
      mpos[i] = -1;
  }
}
