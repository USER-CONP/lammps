
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

#include "compute_conp_vector.h"

#include "assert.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"

using namespace LAMMPS_NS;

#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

ComputeConpVector::ComputeConpVector(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), fp(nullptr) {
  if (narg < 4) error->all(FLERR, "Illegal compute coul/vector command");

  vector_flag = 1;
  size_array_cols = 0;
  size_array_rows = 0;
  size_array_rows_variable = 0;
  extarray = 0;

  fp = nullptr;
  vector = nullptr;

  pairflag = 1;
  kspaceflag = 1;
  boundaryflag = 1;  // include infite boundary correction term
  gaussians = 1;
  recalc_every = 0;
  overwrite = 1;

  g_ewald = 0.0;

  eta =
      utils::numeric(FLERR, arg[3], false, lmp);  // TODO infer from pair_style!

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "pair") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/vector command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        pairflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        pairflag = 0;
      else
        error->all(FLERR, "Illegal compute coul/vector command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "kspace") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/vector command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        kspaceflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        kspaceflag = 0;
      else
        error->all(FLERR, "Illegal compute coul/vector command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "boundary") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/vector command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        boundaryflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        boundaryflag = 0;
      else
        error->all(FLERR, "Illegal compute coul/vector command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "overwrite") ==
               0) {  // TODO  if vector is recalculated overwrite or append
                     // output
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/vector command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        overwrite = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        overwrite = 0;
      else
        error->all(FLERR, "Illegal compute coul/vector command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "file") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute coul/vector command");
      if (comm->me == 0) {
        fp = fopen(arg[iarg + 1], "w");
        if (fp == nullptr)
          error->one(FLERR,
                     fmt::format("Cannot open compute coul/vector file {}: {}",
                                 arg[iarg + 1], utils::getsyserror()));
      }
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute coul/vector command");
  }

  // print file comment lines

  if (fp && comm->me == 0) {
    clearerr(fp);
    fprintf(fp, "# Constant potential coulomb vector\n");
    if (ferror(fp)) error->one(FLERR, "Error writing file header");
    filepos = ftell(fp);
  }
  setup_time_total = 0;
  reduce_time_total = 0;
  kspace_time_total = 0;
  pair_time_total = 0;
  boundary_time_total = 0;
  b_time_total = 0;

  alloc_time_total = 0;
  mpos_time_total = 0;
}

/* ---------------------------------------------------------------------- */

ComputeConpVector::~ComputeConpVector() {
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("B time: {}\n", b_time_total));
    utils::logmesg(lmp, fmt::format("B kspace time: {}\n", kspace_time_total));
    utils::logmesg(lmp, fmt::format("B pair time: {}\n", pair_time_total));
    utils::logmesg(lmp,
                   fmt::format("B boundary time: {}\n", boundary_time_total));
    utils::logmesg(lmp, fmt::format("B setup time: {}\n", setup_time_total));
    utils::logmesg(lmp, fmt::format("B reduce time: {}\n", reduce_time_total));
    utils::logmesg(lmp, fmt::format("B alloc time: {}\n", alloc_time_total));
    utils::logmesg(lmp, fmt::format("B mpos time: {}\n", mpos_time_total));
  }
  delete[] vector;
  if (fp && comm->me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void ComputeConpVector::init_list(int /*id*/, NeighList *ptr) { list = ptr; }

/* ---------------------------------------------------------------------- */

void ComputeConpVector::init() {
  // if non-hybrid, then error if single_enable = 0
  // if hybrid, let hybrid determine if sub-style sets single_enable = 0

  // error if Kspace style does not compute coul/vectr interactions

  if ((boundaryflag || kspaceflag) && force->kspace == nullptr)
    error->all(FLERR, "No Kspace style defined for compute coul/vectr");

  // check if coul pair style is active, no need for single() since done
  // explicitly

  if (pairflag) {
    int itmp;
    double *p_cutoff = (double *)force->pair->extract("cut_coul", itmp);
    if (p_cutoff == nullptr)
      error->all(FLERR, "compute coul/vector is incompatible with Pair style");
    pair = force->pair;
    cutsq = force->pair->cutsq;
  } else
    pair = nullptr;

  if (boundaryflag || kspaceflag) {
    //kspace = force->kspace;
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

void ComputeConpVector::setup() {
  ngroup = group->count(igroup);

  // TODO could be useful to assign homogenously all atoms in both groups to
  // all procs for calculating vector to distribute evenly the workload
  vector = new double[ngroup]();  // init to zero
  // assign atom tags to vector locations and vice versa
  create_taglist();
}
/* ---------------------------------------------------------------------- */

void ComputeConpVector::compute_vector() {
  MPI_Barrier(world);
  double start_time = MPI_Wtime();
  // setup
  double setup_start_time = MPI_Wtime();
  update_mpos();
  for (int i = 0; i < ngroup; i++) vector[i] = 0.;
  MPI_Barrier(world);
  setup_time_total += MPI_Wtime() - setup_start_time;
  // pair
  double pair_start_time = MPI_Wtime();
  if (pairflag) pair_contribution();
  MPI_Barrier(world);
  pair_time_total += MPI_Wtime() - pair_start_time;
  // kspace
  double kspace_start_time = MPI_Wtime();
  if (kspaceflag) conp_kspace->compute_vector(&mpos[0], vector);
  MPI_Barrier(world);
  kspace_time_total += MPI_Wtime() - kspace_start_time;
  // boundary
  double boundary_start_time = MPI_Wtime();
  if (boundaryflag) conp_kspace->compute_vector_corr(&mpos[0], vector);
  MPI_Barrier(world);
  boundary_time_total += MPI_Wtime() - boundary_start_time;
  // reduce
  double reduce_start_time = MPI_Wtime();
  MPI_Allreduce(MPI_IN_PLACE, vector, ngroup, MPI_DOUBLE, MPI_SUM, world);
  MPI_Barrier(world);
  reduce_time_total += MPI_Wtime() - reduce_start_time;
  b_time_total += MPI_Wtime() - start_time;
}

/* ---------------------------------------------------------------------- */

void ComputeConpVector::pair_contribution() {
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int *mask = atom->mask;
  neighbor->build_one(list);
  int const nlocal = atom->nlocal;
  int const inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  for (int ii = 0; ii < inum; ii++) {
    int const i = ilist[ii];
    bool const i_in_electrode = (mask[i] & groupbit);
    double const xtmp = x[i][0];
    double const ytmp = x[i][1];
    double const ztmp = x[i][2];
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int const j = jlist[jj] & NEIGHMASK;
      bool const j_in_electrode = (mask[j] & groupbit);
      if (i_in_electrode == j_in_electrode) continue;
      double const delx = xtmp - x[j][0];  // neighlists take care of pbc
      double const dely = ytmp - x[j][1];
      double const delz = ztmp - x[j][2];
      double const rsq = delx * delx + dely * dely + delz * delz;
      int jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;
      double const r = sqrt(rsq);
      double const rinv = 1.0 / r;
      double aij = rinv;
      if (kspaceflag || boundaryflag) {
        aij *= calc_erfc(g_ewald * r);
        // TODO real-space gaussians?
        if (gaussians) {
          // TODO infer eta from coeffs of pair coul/long/gauss
          aij -= calc_erfc(eta * r) * rinv;
        }
      }
      if (!(newton_pair || j < nlocal)) aij *= 0.5;
      if (i_in_electrode) {
        vector[mpos[i]] += aij * q[j];
      } else if (j_in_electrode) {
        vector[mpos[j]] += aij * q[i];
      }
    }
  }
}
/* ---------------------------------------------------------------------- */

void ComputeConpVector::create_taglist() {
  // assign a tag to each matrix index

  int *mask = atom->mask;
  int const nlocal = atom->nlocal;
  int const nprocs = comm->nprocs;
  tagint *tag = atom->tag;

  std::vector<int> taglist_local;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) taglist_local.push_back(tag[i]);
  }
  int igroupnum_local = taglist_local.size();

  std::vector<int> idispls(nprocs);
  std::vector<int> igroupnum_list(nprocs);
  MPI_Allgather(&igroupnum_local, 1, MPI_INT, &igroupnum_list.front(), 1,
                MPI_INT, world);
  idispls[0] = 0;
  for (int i = 1; i < nprocs; i++) {
    idispls[i] = idispls[i - 1] + igroupnum_list[i - 1];
  }

  std::vector<int> taglist = std::vector<int>(ngroup);
  MPI_Allgatherv(&taglist_local.front(), igroupnum_local, MPI_LMP_TAGINT,
                 &taglist.front(), &igroupnum_list.front(), &idispls.front(),
                 MPI_LMP_TAGINT, world);
  // must be sorted for compatibility with fix_charge_update
  std::sort(taglist.begin(), taglist.end());

  int const tag_max = taglist[ngroup - 1];
  tag_to_iele = std::vector<int>(tag_max + 1, -1);
  for (size_t i = 0; i < taglist.size(); i++) {
    tag_to_iele[taglist[i]] = i;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeConpVector::update_mpos() {
  MPI_Barrier(world);
  double alloc_start = MPI_Wtime();
  int const nall = atom->nlocal + atom->nghost;
  int *tag = atom->tag;
  int *mask = atom->mask;
  mpos = std::vector<bigint>(nall, -1);

  MPI_Barrier(world);
  alloc_time_total += MPI_Wtime() - alloc_start;
  double mpos_start = MPI_Wtime();
  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit)
      mpos[i] = tag_to_iele[tag[i]];
    else
      mpos[i] = -1;
  }
  MPI_Barrier(world);
  mpos_time_total += MPI_Wtime() - mpos_start;
}

/* ---------------------------------------------------------------------- */

double ComputeConpVector::calc_erfc(double x) {
  double expm2 = exp(-x * x);
  double t = 1.0 / (1.0 + EWALD_P * x);
  return t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
}

