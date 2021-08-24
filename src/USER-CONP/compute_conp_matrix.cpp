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

#include "compute_conp_matrix.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

enum { OFF, INTER, INTRA };

/* ---------------------------------------------------------------------- */

ComputeConpMatrix::ComputeConpMatrix(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), mpos(nullptr), fp(nullptr), fp_inv(nullptr) {
  if (narg < 4) error->all(FLERR, "Illegal compute coul/matrix command");

  array_flag = 1;
  size_array_cols = 0;
  size_array_rows = 0;
  size_array_rows_variable = 0;
  extarray = 0;

  fp = nullptr;
  fp_inv = nullptr;
  array = nullptr;
  mpos = nullptr;

  pairflag = 1;
  kspaceflag = 1;
  boundaryflag = 1;  // include infite boundary correction term
  selfflag = 1;
  gaussians = 1;
  recalc_every = 0;
  overwrite = 1;
  assigned = 0;

  g_ewald = 0.0;

  // TODO infer from pair_style!
  eta = utils::numeric(FLERR, arg[3], false, lmp);

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "pair") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/matrix command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        pairflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        pairflag = 0;
      else
        error->all(FLERR, "Illegal compute coul/matrix command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "kspace") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/matrix command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        kspaceflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        kspaceflag = 0;
      else
        error->all(FLERR, "Illegal compute coul/matrix command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "self") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/matrix command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        selfflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        selfflag = 0;
      else
        error->all(FLERR, "Illegal compute coul/matrix command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "boundary") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/matrix command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        boundaryflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        boundaryflag = 0;
      else
        error->all(FLERR, "Illegal compute coul/matrix command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "overwrite") ==
               0) {  // TODO  if matrix is recalculated overwrite or append
                     // output
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/matrix command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        overwrite = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        overwrite = 0;
      else
        error->all(FLERR, "Illegal compute coul/matrix command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "file") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/matrix command");
      if (comm->me == 0) {
        fp = fopen(arg[iarg + 1], "w");
        if (fp == nullptr)
          error->one(FLERR,
                     fmt::format("Cannot open compute coul/matrix file {}: {}",
                                 arg[iarg + 1], utils::getsyserror()));
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "file_inv") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/matrix command");
      if (comm->me == 0) {
        fp_inv = fopen(arg[iarg + 1], "w");
        if (fp_inv == nullptr)
          error->one(FLERR,
                     fmt::format("Cannot open compute coul/matrix file {}: {}",
                                 arg[iarg + 1], utils::getsyserror()));
      }
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute coul/matrix command");
  }

  // open file handler for matrix output
  for (FILE *f : {fp, fp_inv}) {
    if (f && comm->me == 0) {
      clearerr(f);
      if (ferror(f)) error->one(FLERR, "Error writing file header");
      filepos = ftell(f);
    }
  }
}

/* ---------------------------------------------------------------------- */

ComputeConpMatrix::~ComputeConpMatrix() {
  deallocate();

  if (assigned) memory->destroy(mpos);

  for (FILE *f : {fp, fp_inv}) {
    if (f && comm->me == 0) fclose(f);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeConpMatrix::init() {
  // if non-hybrid, then error if single_enable = 0
  // if hybrid, let hybrid determine if sub-style sets single_enable = 0

  // error if Kspace style does not compute coul/matrix interactions

  if ((boundaryflag || kspaceflag) && force->kspace == nullptr)
    error->all(FLERR, "No Kspace style defined for compute coul/matrix");

  // TODO need another flag since we don't use compute_group_group()
  if (kspaceflag && force->kspace->group_group_enable == 0)
    error->all(FLERR, "Kspace style does not support compute coul/matrix");

  // check if coul pair style is active, no need for single() since done
  // explicitly

  if (pairflag) {
    int itmp;
    double *p_cutoff = (double *)force->pair->extract("cut_coul", itmp);
    if (p_cutoff == nullptr)
      error->all(FLERR, "compute coul/matrix is incompatible with Pair style");
    pair = force->pair;
    cutsq = force->pair->cutsq;
  } else
    pair = nullptr;

  if (boundaryflag || kspaceflag) {
    conp_kspace = dynamic_cast<ConpKspace *>(force->kspace);
    if (conp_kspace == nullptr) 
      error->all(FLERR, "Kspace does not implement ConpKspace");
    g_ewald = force->kspace->g_ewald;
  } else
    conp_kspace = nullptr;

  // need an occasional half neighbor list

  if (pairflag) {
    int irequest = neighbor->request(this, instance_me);
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->compute = 1;
    neighbor->requests[irequest]->occasional = 1;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeConpMatrix::setup() {
  ngroup = group->count(igroup);

  // TODO could be useful to assign homogenously all atoms in both groups to
  // all procs for calculating matrix to distribute evenly the workload
  allocate();

  // assign atom tags to matrix locations and vice versa
  matrix_assignment();
}

/* ---------------------------------------------------------------------- */

void ComputeConpMatrix::init_list(int /*id*/, NeighList *ptr) { list = ptr; }

/* ---------------------------------------------------------------------- */

void ComputeConpMatrix::compute_array() {
  // setting all entries of coulomb matrix to zero
  size_t nbytes = sizeof(double) * ngroup;
  if (nbytes)
    for (int i = 0; i < ngroup; i++) memset(&array[i][0], 0, nbytes);

  MPI_Barrier(world);
  double kspace_time = MPI_Wtime();
  if (kspaceflag) conp_kspace->compute_matrix(mpos, array);
  MPI_Barrier(world);
  if (comm->me == 0)
    utils::logmesg(lmp,
                   fmt::format("Kspace time: {}\n", MPI_Wtime() - kspace_time));
  if (pairflag) pair_contribution();
  if (selfflag) self_contribution();
  if (boundaryflag) conp_kspace->compute_matrix_corr(mpos, array);

  // reduce coulomb matrix with contributions from all procs
  // all procs need to know full matrix for matrix inversion
  for (int i = 0; i < ngroup; i++) {
    MPI_Allreduce(MPI_IN_PLACE, &array[i][0], ngroup, MPI_DOUBLE, MPI_SUM,
                  world);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeConpMatrix::pair_contribution() {
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

  bigint jpos;

  double etaij =
      eta * eta / sqrt(2.0 * eta * eta);  // see mw ewald theory eq. (29)-(30)

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

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
      // skip if atom J is not in either group
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

        // kspace solver?

        if (kspaceflag || boundaryflag) {
          grij = g_ewald * r;
          expm2 = exp(-grij * grij);
          t = 1.0 / (1.0 + EWALD_P * grij);
          erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;

          aij *= erfc;
        }

        // real-space gaussians?

        if (gaussians) {
          // TODO infer eta from coeffs of pair coul/long/gauss
          etarij = etaij * r;
          expm2 = exp(-etarij * etarij);
          t = 1.0 / (1.0 + EWALD_P * etarij);
          erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;

          aij -= erfc * rinv;
        }

        // TODO we don't assign ghost atoms to matrix positions - which would
        // be a non-trivial task in matrix_assignment() - so I'm using this
        // rather slow approach here... It seems that lammps stores also tags of
        // ghost atoms
        for (jpos = 0; jpos < ngroup; jpos++)
          if (mat2tag[jpos] == tag[j]) break;

        // newton on or off?

        if (!(newton_pair || j < nlocal)) aij *= 0.5;
        array[mpos[i]][jpos] += aij;
        array[jpos][mpos[i]] += aij;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeConpMatrix::self_contribution() {
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  const double selfint = 2.0 / MY_PIS * g_ewald;
  const double preta = MY_SQRT2 / MY_PIS;

  // TODO infer eta from pair_coeffs
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      array[mpos[i]][mpos[i]] += preta * eta - selfint;
    }
}

/* ----------------------------------------------------------------------
   looks up to which proc each atom in each group belongs and creates a
   local array which locates the position of each local atom in the global
   matrix. entries are sorted: first A then B. need to be so complex here
   b/c atom tags might not be be consecutive or sorted in any way.
------------------------------------------------------------------------- */

void ComputeConpMatrix::matrix_assignment() {
  // assign local matrix indices to local atoms on each proc

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nprocs = comm->nprocs;
  tagint *tag = atom->tag;

  tagint *itaglist, *itaglist_local;
  int igroupnum_local;
  int *igroupnum_list;
  int *idispls;

  igroupnum_local = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) igroupnum_local++;
  }

  memory->create(idispls, nprocs, "coul/matrix:idispls");
  memory->create(igroupnum_list, nprocs, "coul/matrix:ilist");

  MPI_Allgather(&igroupnum_local, 1, MPI_INT, igroupnum_list, 1, MPI_INT,
                world);

  idispls[0] = 0;
  for (int i = 1; i < nprocs; i++) {
    idispls[i] = idispls[i - 1] + igroupnum_list[i - 1];
  }

  memory->create(itaglist_local, igroupnum_local, "coul/matrix:itaglist_local");

  igroupnum_local = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      itaglist_local[igroupnum_local] = tag[i];
      igroupnum_local++;
    }
  }

  memory->create(itaglist, ngroup, "coul/matrix:itaglist");
  MPI_Allgatherv(itaglist_local, igroupnum_local, MPI_LMP_TAGINT, itaglist,
                 igroupnum_list, idispls, MPI_LMP_TAGINT, world);
  // must be sorted for compatibility with fix_charge_update
  std::sort(itaglist, itaglist + ngroup);

  // if local+ghost matrix assignment already created, recreate

  if (assigned) memory->destroy(mpos);
  memory->create(mpos, nlocal, "coul/matrix:mpos");

  assigned = 1;

  // local+ghost non-matrix atoms are -1 in mpos

  size_t nbytes = sizeof(bigint) * nlocal;
  if (nbytes) memset(mpos, -1, nbytes);

  // store which tag represents value in matrix
  for (bigint i = 0; i < ngroup; i++) mat2tag[i] = itaglist[i];

  // create global matrix indices for local+ghost atoms
  for (bigint ii = 0; ii < ngroup; ii++)
    for (int i = 0; i < nlocal; i++)
      if (mat2tag[ii] == tag[i]) mpos[i] = ii;

  memory->destroy(igroupnum_list);
  memory->destroy(idispls);
  memory->destroy(itaglist_local);
  memory->destroy(itaglist);
}

/* ---------------------------------------------------------------------- */

void ComputeConpMatrix::allocate() {
  memory->create(mat2tag, ngroup, "coul/matrix:mat2tag");
  array = new double *[ngroup];
  for (bigint i = 0; i < ngroup; i++) array[i] = new double[ngroup];
}

/* ---------------------------------------------------------------------- */

void ComputeConpMatrix::deallocate() {
  memory->destroy(mat2tag);
  for (bigint i = 0; i < ngroup; i++) delete[] array[i];
  delete[] array;
}

