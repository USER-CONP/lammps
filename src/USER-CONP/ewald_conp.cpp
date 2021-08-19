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
   Contributing authors: Roy Pollock (LLNL), Paul Crozier (SNL)
     per-atom energy/virial added by German Samolyuk (ORNL), Stan Moore (BYU)
     group/group energy/force added by Stan Moore (BYU)
     triclinic added by Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "ewald_conp.h"

#include <cmath>
#include <iostream>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "pair.h"
#include "slab_2d.h"
#include "slab_dipole.h"
#include "update.h"
#include "wire_dipole.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

EwaldConp::EwaldConp(LAMMPS *lmp)
    : KSpace(lmp),
      kxvecs(nullptr),
      kyvecs(nullptr),
      kzvecs(nullptr),
      ug(nullptr),
      eg(nullptr),
      vg(nullptr),
      ek(nullptr),
      sfacrl(nullptr),
      sfacim(nullptr),
      sfacrl_all(nullptr),
      sfacim_all(nullptr),
      cs(nullptr),
      sn(nullptr),
      sfacrl_A(nullptr),
      sfacim_A(nullptr),
      sfacrl_A_all(nullptr),
      sfacim_A_all(nullptr),
      sfacrl_B(nullptr),
      sfacim_B(nullptr),
      sfacrl_B_all(nullptr),
      sfacim_B_all(nullptr) {
  group_allocate_flag = 0;
  kmax_created = 0;
  ewaldflag = 1;
  group_group_enable = 1;

  accuracy_relative = 0.0;

  kmax = 0;
  kxvecs = kyvecs = kzvecs = nullptr;

  ug = nullptr;
  eg = vg = nullptr;
  sfacrl = sfacim = sfacrl_all = sfacim_all = nullptr;

  nmax = 0;
  ek = nullptr;
  cs = sn = nullptr;

  kcount = 0;
  eikr_step = -1;
}

void EwaldConp::settings(int narg, char **arg) {
  if (narg != 1) error->all(FLERR, "Illegal kspace_style ewald command");

  accuracy_relative = fabs(utils::numeric(FLERR, arg[0], false, lmp));
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

EwaldConp::~EwaldConp() {
  delete boundcorr;
  deallocate();
  if (group_allocate_flag) deallocate_groups();
  memory->destroy(ek);
  memory->destroy3d_offset(cs, -kmax_created);
  memory->destroy3d_offset(sn, -kmax_created);
}

/* ---------------------------------------------------------------------- */

void EwaldConp::init() {
  if (comm->me == 0) utils::logmesg(lmp, "Ewald initialization ...\n");

  // error check
  triclinic_check();
  if (domain->dimension == 2)
    error->all(FLERR, "Cannot use Ewald with 2d simulation");

  if (!atom->q_flag)
    error->all(FLERR, "Kspace style requires atom attribute q");

  if (slabflag == 0 && wireflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR, "Cannot use non-periodic boundaries with Ewald");
  if (slabflag) {
    if (wireflag)
      error->all(FLERR, "Cannot use slab and wire corrections together");
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR, "Incorrect boundaries with slab Ewald");
    if (domain->triclinic)
      error->all(FLERR,
                 "Cannot (yet) use Ewald with triclinic box "
                 "and slab correction");
  } else if (wireflag) {
    if (domain->zperiodic != 1 || domain->boundary[0][0] != 1 ||
        domain->boundary[0][1] != 1 || domain->boundary[1][0] != 1 ||
        domain->boundary[1][1] != 1)
      error->all(FLERR, "Incorrect boundaries with wire Ewald");
    if (domain->triclinic)
      error->all(FLERR,
                 "Cannot (yet) use Ewald with triclinic box "
                 "and wire correction");
  }

  if (slabflag == 1) {
    // EW3Dc dipole correction
    boundcorr = new SlabDipole(lmp);
  } else if (slabflag == 3) {
    boundcorr = new Slab2d(lmp);
    cout << "WARNING WIP" << endl;
  } else if (wireflag == 1) {
    // EW3Dc wire correction
    boundcorr = new WireDipole(lmp);
  } else {
    error->all(FLERR, "pppm conp with dipole corrections, only");
  }

  // compute two charge force

  two_charge();

  // extract short-range Coulombic cutoff from pair style

  triclinic = domain->triclinic;
  pair_check();

  int itmp;
  double *p_cutoff = (double *)force->pair->extract("cut_coul", itmp);
  if (p_cutoff == nullptr)
    error->all(FLERR, "KSpace style is incompatible with Pair style");
  double cutoff = *p_cutoff;

  // compute qsum & qsqsum and warn if not charge-neutral

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  qsum_qsq();

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0)
    accuracy = accuracy_absolute;
  else
    accuracy = accuracy_relative * two_charge_force;

  // setup K-space resolution

  bigint natoms = atom->natoms;

  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab Ewald
  // 3d Ewald just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;

  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  if (!gewaldflag) {
    if (accuracy <= 0.0) error->all(FLERR, "KSpace accuracy must be > 0");
    if (q2 == 0.0)
      error->all(FLERR, "Must use 'kspace_modify gewald' for uncharged system");
    g_ewald =
        accuracy * sqrt(natoms * cutoff * xprd * yprd * zprd) / (2.0 * q2);
    if (g_ewald >= 1.0)
      g_ewald = (1.35 - 0.15 * log(accuracy)) / cutoff;
    else
      g_ewald = sqrt(-log(g_ewald)) / cutoff;
  }

  // setup Ewald coefficients so can print stats

  setup();

  // final RMS accuracy

  double lprx = rms(kxmax_orig, xprd_wire, natoms, q2);
  double lpry = rms(kymax_orig, yprd_wire, natoms, q2);
  double lprz = rms(kzmax_orig, zprd_slab, natoms, q2);
  double lpr = sqrt(lprx * lprx + lpry * lpry + lprz * lprz) / sqrt(3.0);
  double q2_over_sqrt =
      q2 / sqrt(natoms * cutoff * xprd_wire * yprd_wire * zprd_slab);
  double spr = 2.0 * q2_over_sqrt * exp(-g_ewald * g_ewald * cutoff * cutoff);
  double tpr = estimate_table_accuracy(q2_over_sqrt, spr);
  double estimated_accuracy = sqrt(lpr * lpr + spr * spr + tpr * tpr);

  // stats

  if (comm->me == 0) {
    std::string mesg =
        fmt::format("  G vector (1/distance) = {:.8g}\n", g_ewald);
    mesg += fmt::format("  estimated absolute RMS force accuracy = {:.8g}\n",
                        estimated_accuracy);
    mesg += fmt::format("  estimated relative force accuracy = {:.8g}\n",
                        estimated_accuracy / two_charge_force);
    mesg += fmt::format("  KSpace vectors: actual max1d max3d = {} {} {}\n",
                        kcount, kmax, kmax3d);
    mesg += fmt::format("                  kxmax kymax kzmax  = {} {} {}\n",
                        kxmax, kymax, kzmax);
    utils::logmesg(lmp, mesg);
  }
}

/* ----------------------------------------------------------------------
   adjust Ewald coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void EwaldConp::setup() {
  // volume-dependent factors

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // adjustment of z dimension for 2d slab Ewald
  // 3d Ewald just uses zprd since slab_volfactor = 1.0

  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;
  volume = xprd_wire * yprd_wire * zprd_slab;

  area = xprd_wire * yprd_wire;

  unitk[0] = 2.0 * MY_PI / xprd_wire;
  unitk[1] = 2.0 * MY_PI / yprd_wire;
  unitk[2] = 2.0 * MY_PI / zprd_slab;

  int kmax_old = kmax;

  if (kewaldflag == 0) {
    // determine kmax
    // function of current box size, accuracy, G_ewald (short-range cutoff)

    bigint natoms = atom->natoms;
    double err;
    kxmax = 1;
    kymax = 1;
    kzmax = 1;

    err = rms(kxmax, xprd_wire, natoms, q2);
    while (err > accuracy) {
      kxmax++;
      err = rms(kxmax, xprd_wire, natoms, q2);
    }

    err = rms(kymax, yprd_wire, natoms, q2);
    while (err > accuracy) {
      kymax++;
      err = rms(kymax, yprd_wire, natoms, q2);
    }

    err = rms(kzmax, zprd_slab, natoms, q2);
    while (err > accuracy) {
      kzmax++;
      err = rms(kzmax, zprd_slab, natoms, q2);
    }

    kmax = MAX(kxmax, kymax);
    kmax = MAX(kmax, kzmax);
    kmax3d = 4 * kmax * kmax * kmax + 6 * kmax * kmax + 3 * kmax;

    double gsqxmx = unitk[0] * unitk[0] * kxmax * kxmax;
    double gsqymx = unitk[1] * unitk[1] * kymax * kymax;
    double gsqzmx = unitk[2] * unitk[2] * kzmax * kzmax;
    gsqmx = MAX(gsqxmx, gsqymx);
    gsqmx = MAX(gsqmx, gsqzmx);

    kxmax_orig = kxmax;
    kymax_orig = kymax;
    kzmax_orig = kzmax;

    // scale lattice vectors for triclinic skew

    if (triclinic) {
      double tmp[3];
      tmp[0] = kxmax / xprd;
      tmp[1] = kymax / yprd;
      tmp[2] = kzmax / zprd;
      lamda2xT(&tmp[0], &tmp[0]);
      kxmax = MAX(1, static_cast<int>(tmp[0]));
      kymax = MAX(1, static_cast<int>(tmp[1]));
      kzmax = MAX(1, static_cast<int>(tmp[2]));

      kmax = MAX(kxmax, kymax);
      kmax = MAX(kmax, kzmax);
      kmax3d = 4 * kmax * kmax * kmax + 6 * kmax * kmax + 3 * kmax;
    }

  } else {
    kxmax = kx_ewald;
    kymax = ky_ewald;
    kzmax = kz_ewald;

    kxmax_orig = kxmax;
    kymax_orig = kymax;
    kzmax_orig = kzmax;

    kmax = MAX(kxmax, kymax);
    kmax = MAX(kmax, kzmax);
    kmax3d = 4 * kmax * kmax * kmax + 6 * kmax * kmax + 3 * kmax;

    double gsqxmx = unitk[0] * unitk[0] * kxmax * kxmax;
    double gsqymx = unitk[1] * unitk[1] * kymax * kymax;
    double gsqzmx = unitk[2] * unitk[2] * kzmax * kzmax;
    gsqmx = MAX(gsqxmx, gsqymx);
    gsqmx = MAX(gsqmx, gsqzmx);
  }

  gsqmx *= 1.00001;

  // if size has grown, reallocate k-dependent and nlocal-dependent arrays

  if (kmax > kmax_old) {
    deallocate();
    allocate();
    group_allocate_flag = 0;

    memory->destroy(ek);
    memory->destroy3d_offset(cs, -kmax_created);
    memory->destroy3d_offset(sn, -kmax_created);
    nmax = atom->nmax;
    memory->create(ek, nmax, 3, "ewald/conp:ek");
    memory->create3d_offset(cs, -kmax, kmax, 3, nmax, "ewald/conp:cs");
    memory->create3d_offset(sn, -kmax, kmax, 3, nmax, "ewald/conp:sn");
    kmax_created = kmax;
  }
  boundcorr->setup(xprd_wire, yprd_wire, zprd_slab, g_ewald);

  // pre-compute Ewald coefficients

  if (triclinic == 0)
    coeffs();
  else
    coeffs_triclinic();
}

/* ----------------------------------------------------------------------
   compute RMS accuracy for a dimension
------------------------------------------------------------------------- */

double EwaldConp::rms(int km, double prd, bigint natoms, double q2) {
  if (natoms == 0) natoms = 1;  // avoid division by zero
  double value =
      2.0 * q2 * g_ewald / prd * sqrt(1.0 / (MY_PI * km * natoms)) *
      exp(-MY_PI * MY_PI * km * km / (g_ewald * g_ewald * prd * prd));

  return value;
}

/* ----------------------------------------------------------------------
   compute the Ewald long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldConp::compute(int eflag, int vflag) {
  // set energy/virial flags
  ev_init(eflag, vflag);

  qsum_qsq(0);

  // return if there are no charges
  if (qsqsum == 0.0) return;

  update_eikr(true);

  MPI_Allreduce(sfacrl, sfacrl_all, kcount, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(sfacim, sfacim_all, kcount, MPI_DOUBLE, MPI_SUM, world);

  // K-space portion of electric field
  // double loop over K-vectors and local atoms
  // perform per-atom calculations if needed

  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  int kx, ky, kz;
  double cypz, sypz, exprl, expim, partial, partial_peratom;

  for (int i = 0; i < nlocal; i++) {
    ek[i][0] = 0.0;
    ek[i][1] = 0.0;
    ek[i][2] = 0.0;
  }

  for (int k = 0; k < kcount; k++) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];

    for (int i = 0; i < nlocal; i++) {
      cypz = cs[ky][1][i] * cs[kz][2][i] - sn[ky][1][i] * sn[kz][2][i];
      sypz = sn[ky][1][i] * cs[kz][2][i] + cs[ky][1][i] * sn[kz][2][i];
      exprl = cs[kx][0][i] * cypz - sn[kx][0][i] * sypz;
      expim = sn[kx][0][i] * cypz + cs[kx][0][i] * sypz;
      partial = expim * sfacrl_all[k] - exprl * sfacim_all[k];
      ek[i][0] += partial * eg[k][0];
      ek[i][1] += partial * eg[k][1];
      ek[i][2] += partial * eg[k][2];

      if (evflag_atom) {
        partial_peratom = exprl * sfacrl_all[k] + expim * sfacim_all[k];
        if (eflag_atom) eatom[i] += q[i] * ug[k] * partial_peratom;
        if (vflag_atom)
          for (int j = 0; j < 6; j++)
            vatom[i][j] += ug[k] * vg[k][j] * partial_peratom;
      }
    }
  }

  // convert E-field to force

  const double qscale = qqrd2e * scale;

  for (int i = 0; i < nlocal; i++) {
    if (wireflag != 2) {
      f[i][0] += qscale * q[i] * ek[i][0];
      f[i][1] += qscale * q[i] * ek[i][1];
    }
    if (slabflag != 2) {
      f[i][2] += qscale * q[i] * ek[i][2];
    }
  }

  // sum global energy across Kspace vevs and add in volume-dependent term

  if (eflag_global) {
    for (int k = 0; k < kcount; k++)
      energy += ug[k] *
                (sfacrl_all[k] * sfacrl_all[k] + sfacim_all[k] * sfacim_all[k]);

    energy -= g_ewald * qsqsum / MY_PIS +
              MY_PI2 * qsum * qsum / (g_ewald * g_ewald * volume);
    energy *= qscale;
  }

  // global virial

  if (vflag_global) {
    double uk;
    for (int k = 0; k < kcount; k++) {
      uk = ug[k] *
           (sfacrl_all[k] * sfacrl_all[k] + sfacim_all[k] * sfacim_all[k]);
      for (int j = 0; j < 6; j++) virial[j] += uk * vg[k][j];
    }
    for (int j = 0; j < 6; j++) virial[j] *= qscale;
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    if (eflag_atom) {
      for (int i = 0; i < nlocal; i++) {
        eatom[i] -= g_ewald * q[i] * q[i] / MY_PIS +
                    MY_PI2 * q[i] * qsum / (g_ewald * g_ewald * volume);
        eatom[i] *= qscale;
      }
    }

    if (vflag_atom)
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < 6; j++) vatom[i][j] *= q[i] * qscale;
  }

  boundcorr->compute_corr(qsum, eflag_atom, eflag_global, energy, eatom);
}

/* ---------------------------------------------------------------------- */

void EwaldConp::eik_dot_r() {
  int i, k, l, m, n, ic;
  double cstr1, sstr1, cstr2, sstr2, cstr3, sstr3, cstr4, sstr4;
  double sqk, clpm, slpm;

  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  n = 0;

  // (k,0,0), (0,l,0), (0,0,m)

  for (ic = 0; ic < 3; ic++) {
    sqk = unitk[ic] * unitk[ic];
    if (sqk <= gsqmx) {
      cstr1 = 0.0;
      sstr1 = 0.0;
      for (i = 0; i < nlocal; i++) {
        cs[0][ic][i] = 1.0;
        sn[0][ic][i] = 0.0;
        cs[1][ic][i] = cos(unitk[ic] * x[i][ic]);
        sn[1][ic][i] = sin(unitk[ic] * x[i][ic]);
        cs[-1][ic][i] = cs[1][ic][i];
        sn[-1][ic][i] = -sn[1][ic][i];
        cstr1 += q[i] * cs[1][ic][i];
        sstr1 += q[i] * sn[1][ic][i];
      }
      if (slabflag != 3 || ic < 2) {  // skip (0, 0, m) for ew2d
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
      }
    }
  }

  for (m = 2; m <= kmax; m++) {
    for (ic = 0; ic < 3; ic++) {
      sqk = m * unitk[ic] * m * unitk[ic];
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cs[m][ic][i] =
              cs[m - 1][ic][i] * cs[1][ic][i] - sn[m - 1][ic][i] * sn[1][ic][i];
          sn[m][ic][i] =
              sn[m - 1][ic][i] * cs[1][ic][i] + cs[m - 1][ic][i] * sn[1][ic][i];
          cs[-m][ic][i] = cs[m][ic][i];
          sn[-m][ic][i] = -sn[m][ic][i];
          cstr1 += q[i] * cs[m][ic][i];
          sstr1 += q[i] * sn[m][ic][i];
        }
        if (slabflag != 3 || ic < 2) {  // skip (0, 0, m) for ew2d
          sfacrl[n] = cstr1;
          sfacim[n++] = sstr1;
        }
      }
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (k * unitk[0] * k * unitk[0]) + (l * unitk[1] * l * unitk[1]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cstr1 +=
              q[i] * (cs[k][0][i] * cs[l][1][i] - sn[k][0][i] * sn[l][1][i]);
          sstr1 +=
              q[i] * (sn[k][0][i] * cs[l][1][i] + cs[k][0][i] * sn[l][1][i]);
          cstr2 +=
              q[i] * (cs[k][0][i] * cs[l][1][i] + sn[k][0][i] * sn[l][1][i]);
          sstr2 +=
              q[i] * (sn[k][0][i] * cs[l][1][i] - cs[k][0][i] * sn[l][1][i]);
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (l * unitk[1] * l * unitk[1]) + (m * unitk[2] * m * unitk[2]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cstr1 +=
              q[i] * (cs[l][1][i] * cs[m][2][i] - sn[l][1][i] * sn[m][2][i]);
          sstr1 +=
              q[i] * (sn[l][1][i] * cs[m][2][i] + cs[l][1][i] * sn[m][2][i]);
          cstr2 +=
              q[i] * (cs[l][1][i] * cs[m][2][i] + sn[l][1][i] * sn[m][2][i]);
          sstr2 +=
              q[i] * (sn[l][1][i] * cs[m][2][i] - cs[l][1][i] * sn[m][2][i]);
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (k * unitk[0] * k * unitk[0]) + (m * unitk[2] * m * unitk[2]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cstr1 +=
              q[i] * (cs[k][0][i] * cs[m][2][i] - sn[k][0][i] * sn[m][2][i]);
          sstr1 +=
              q[i] * (sn[k][0][i] * cs[m][2][i] + cs[k][0][i] * sn[m][2][i]);
          cstr2 +=
              q[i] * (cs[k][0][i] * cs[m][2][i] + sn[k][0][i] * sn[m][2][i]);
          sstr2 +=
              q[i] * (sn[k][0][i] * cs[m][2][i] - cs[k][0][i] * sn[m][2][i]);
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (k * unitk[0] * k * unitk[0]) + (l * unitk[1] * l * unitk[1]) +
              (m * unitk[2] * m * unitk[2]);
        if (sqk <= gsqmx) {
          cstr1 = 0.0;
          sstr1 = 0.0;
          cstr2 = 0.0;
          sstr2 = 0.0;
          cstr3 = 0.0;
          sstr3 = 0.0;
          cstr4 = 0.0;
          sstr4 = 0.0;
          for (i = 0; i < nlocal; i++) {
            clpm = cs[l][1][i] * cs[m][2][i] - sn[l][1][i] * sn[m][2][i];
            slpm = sn[l][1][i] * cs[m][2][i] + cs[l][1][i] * sn[m][2][i];
            cstr1 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
            sstr1 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);

            clpm = cs[l][1][i] * cs[m][2][i] + sn[l][1][i] * sn[m][2][i];
            slpm = -sn[l][1][i] * cs[m][2][i] + cs[l][1][i] * sn[m][2][i];
            cstr2 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
            sstr2 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);

            clpm = cs[l][1][i] * cs[m][2][i] + sn[l][1][i] * sn[m][2][i];
            slpm = sn[l][1][i] * cs[m][2][i] - cs[l][1][i] * sn[m][2][i];
            cstr3 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
            sstr3 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);

            clpm = cs[l][1][i] * cs[m][2][i] - sn[l][1][i] * sn[m][2][i];
            slpm = -sn[l][1][i] * cs[m][2][i] - cs[l][1][i] * sn[m][2][i];
            cstr4 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
            sstr4 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);
          }
          sfacrl[n] = cstr1;
          sfacim[n++] = sstr1;
          sfacrl[n] = cstr2;
          sfacim[n++] = sstr2;
          sfacrl[n] = cstr3;
          sfacim[n++] = sstr3;
          sfacrl[n] = cstr4;
          sfacim[n++] = sstr4;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void EwaldConp::eik_dot_r_triclinic() {
  int i, k, l, m, n, ic;
  double cstr1, sstr1;
  double sqk, clpm, slpm;

  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  double unitk_lamda[3];

  double max_kvecs[3];
  max_kvecs[0] = kxmax;
  max_kvecs[1] = kymax;
  max_kvecs[2] = kzmax;

  // (k,0,0), (0,l,0), (0,0,m)

  for (ic = 0; ic < 3; ic++) {
    unitk_lamda[0] = 0.0;
    unitk_lamda[1] = 0.0;
    unitk_lamda[2] = 0.0;
    unitk_lamda[ic] = 2.0 * MY_PI;
    x2lamdaT(&unitk_lamda[0], &unitk_lamda[0]);
    sqk = unitk_lamda[ic] * unitk_lamda[ic];
    if (sqk <= gsqmx) {
      for (i = 0; i < nlocal; i++) {
        cs[0][ic][i] = 1.0;
        sn[0][ic][i] = 0.0;
        cs[1][ic][i] = cos(unitk_lamda[0] * x[i][0] + unitk_lamda[1] * x[i][1] +
                           unitk_lamda[2] * x[i][2]);
        sn[1][ic][i] = sin(unitk_lamda[0] * x[i][0] + unitk_lamda[1] * x[i][1] +
                           unitk_lamda[2] * x[i][2]);
        cs[-1][ic][i] = cs[1][ic][i];
        sn[-1][ic][i] = -sn[1][ic][i];
      }
    }
  }

  for (ic = 0; ic < 3; ic++) {
    for (m = 2; m <= max_kvecs[ic]; m++) {
      unitk_lamda[0] = 0.0;
      unitk_lamda[1] = 0.0;
      unitk_lamda[2] = 0.0;
      unitk_lamda[ic] = 2.0 * MY_PI * m;
      x2lamdaT(&unitk_lamda[0], &unitk_lamda[0]);
      sqk = unitk_lamda[ic] * unitk_lamda[ic];
      for (i = 0; i < nlocal; i++) {
        cs[m][ic][i] =
            cs[m - 1][ic][i] * cs[1][ic][i] - sn[m - 1][ic][i] * sn[1][ic][i];
        sn[m][ic][i] =
            sn[m - 1][ic][i] * cs[1][ic][i] + cs[m - 1][ic][i] * sn[1][ic][i];
        cs[-m][ic][i] = cs[m][ic][i];
        sn[-m][ic][i] = -sn[m][ic][i];
      }
    }
  }

  for (n = 0; n < kcount; n++) {
    k = kxvecs[n];
    l = kyvecs[n];
    m = kzvecs[n];
    cstr1 = 0.0;
    sstr1 = 0.0;
    for (i = 0; i < nlocal; i++) {
      clpm = cs[l][1][i] * cs[m][2][i] - sn[l][1][i] * sn[m][2][i];
      slpm = sn[l][1][i] * cs[m][2][i] + cs[l][1][i] * sn[m][2][i];
      cstr1 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
      sstr1 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);
    }
    sfacrl[n] = cstr1;
    sfacim[n] = sstr1;
  }
}

/* ----------------------------------------------------------------------
   pre-compute coefficients for each Ewald K-vector
------------------------------------------------------------------------- */

void EwaldConp::coeffs() {
  int k, l, m;
  double sqk, vterm;

  double g_ewald_sq_inv = 1.0 / (g_ewald * g_ewald);
  double preu = 4.0 * MY_PI / volume;

  kcount = 0;

  // (k,0,0), (0,l,0), (0,0,m), skip (0,0) in case of EW2D (slabflag == 3)
  for (m = 1; m <= kmax; m++) {
    sqk = (m * unitk[0]) * (m * unitk[0]);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = m;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = 0;
      ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
      eg[kcount][0] = 2.0 * unitk[0] * m * ug[kcount];
      eg[kcount][1] = 0.0;
      eg[kcount][2] = 0.0;
      vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
      vg[kcount][0] = 1.0 + vterm * (unitk[0] * m) * (unitk[0] * m);
      vg[kcount][1] = 1.0;
      vg[kcount][2] = 1.0;
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
    sqk = (m * unitk[1]) * (m * unitk[1]);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = m;
      kzvecs[kcount] = 0;
      ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
      eg[kcount][0] = 0.0;
      eg[kcount][1] = 2.0 * unitk[1] * m * ug[kcount];
      eg[kcount][2] = 0.0;
      vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
      vg[kcount][0] = 1.0;
      vg[kcount][1] = 1.0 + vterm * (unitk[1] * m) * (unitk[1] * m);
      vg[kcount][2] = 1.0;
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
    sqk = (m * unitk[2]) * (m * unitk[2]);
    if (sqk <= gsqmx && slabflag != 3) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = m;
      ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
      eg[kcount][0] = 0.0;
      eg[kcount][1] = 0.0;
      eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
      vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
      vg[kcount][0] = 1.0;
      vg[kcount][1] = 1.0;
      vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (unitk[0] * k) * (unitk[0] * k) + (unitk[1] * l) * (unitk[1] * l);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = k;
        kyvecs[kcount] = l;
        kzvecs[kcount] = 0;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
        eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
        eg[kcount][2] = 0.0;
        vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
        vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
        vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
        vg[kcount][2] = 1.0;
        vg[kcount][3] = vterm * unitk[0] * k * unitk[1] * l;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = 0.0;
        kcount++;

        kxvecs[kcount] = k;
        kyvecs[kcount] = -l;
        kzvecs[kcount] = 0;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
        eg[kcount][1] = -2.0 * unitk[1] * l * ug[kcount];
        eg[kcount][2] = 0.0;
        vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
        vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
        vg[kcount][2] = 1.0;
        vg[kcount][3] = -vterm * unitk[0] * k * unitk[1] * l;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = 0.0;
        kcount++;
        ;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (unitk[1] * l) * (unitk[1] * l) + (unitk[2] * m) * (unitk[2] * m);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 0.0;
        eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
        eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
        vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
        vg[kcount][0] = 1.0;
        vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
        vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
        vg[kcount][3] = 0.0;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = vterm * unitk[1] * l * unitk[2] * m;
        kcount++;

        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = -m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 0.0;
        eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
        eg[kcount][2] = -2.0 * unitk[2] * m * ug[kcount];
        vg[kcount][0] = 1.0;
        vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
        vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
        vg[kcount][3] = 0.0;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = -vterm * unitk[1] * l * unitk[2] * m;
        kcount++;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (unitk[0] * k) * (unitk[0] * k) + (unitk[2] * m) * (unitk[2] * m);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = k;
        kyvecs[kcount] = 0;
        kzvecs[kcount] = m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
        eg[kcount][1] = 0.0;
        eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
        vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
        vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
        vg[kcount][1] = 1.0;
        vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
        vg[kcount][3] = 0.0;
        vg[kcount][4] = vterm * unitk[0] * k * unitk[2] * m;
        vg[kcount][5] = 0.0;
        kcount++;

        kxvecs[kcount] = k;
        kyvecs[kcount] = 0;
        kzvecs[kcount] = -m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
        eg[kcount][1] = 0.0;
        eg[kcount][2] = -2.0 * unitk[2] * m * ug[kcount];
        vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
        vg[kcount][1] = 1.0;
        vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
        vg[kcount][3] = 0.0;
        vg[kcount][4] = -vterm * unitk[0] * k * unitk[2] * m;
        vg[kcount][5] = 0.0;
        kcount++;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (unitk[0] * k) * (unitk[0] * k) +
              (unitk[1] * l) * (unitk[1] * l) + (unitk[2] * m) * (unitk[2] * m);
        if (sqk <= gsqmx) {
          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
          eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
          eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
          vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
          vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
          vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
          vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
          vg[kcount][3] = vterm * unitk[0] * k * unitk[1] * l;
          vg[kcount][4] = vterm * unitk[0] * k * unitk[2] * m;
          vg[kcount][5] = vterm * unitk[1] * l * unitk[2] * m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = -l;
          kzvecs[kcount] = m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
          eg[kcount][1] = -2.0 * unitk[1] * l * ug[kcount];
          eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
          vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
          vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
          vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
          vg[kcount][3] = -vterm * unitk[0] * k * unitk[1] * l;
          vg[kcount][4] = vterm * unitk[0] * k * unitk[2] * m;
          vg[kcount][5] = -vterm * unitk[1] * l * unitk[2] * m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = -m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
          eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
          eg[kcount][2] = -2.0 * unitk[2] * m * ug[kcount];
          vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
          vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
          vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
          vg[kcount][3] = vterm * unitk[0] * k * unitk[1] * l;
          vg[kcount][4] = -vterm * unitk[0] * k * unitk[2] * m;
          vg[kcount][5] = -vterm * unitk[1] * l * unitk[2] * m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = -l;
          kzvecs[kcount] = -m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
          eg[kcount][1] = -2.0 * unitk[1] * l * ug[kcount];
          eg[kcount][2] = -2.0 * unitk[2] * m * ug[kcount];
          vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
          vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
          vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
          vg[kcount][3] = -vterm * unitk[0] * k * unitk[1] * l;
          vg[kcount][4] = -vterm * unitk[0] * k * unitk[2] * m;
          vg[kcount][5] = vterm * unitk[1] * l * unitk[2] * m;
          kcount++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pre-compute coefficients for each Ewald K-vector for a triclinic
   system
------------------------------------------------------------------------- */

void EwaldConp::coeffs_triclinic() {
  int k, l, m;
  double sqk, vterm;

  double g_ewald_sq_inv = 1.0 / (g_ewald * g_ewald);
  double preu = 4.0 * MY_PI / volume;

  double unitk_lamda[3];

  kcount = 0;

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = -kymax; l <= kymax; l++) {
      for (m = -kzmax; m <= kzmax; m++) {
        unitk_lamda[0] = 2.0 * MY_PI * k;
        unitk_lamda[1] = 2.0 * MY_PI * l;
        unitk_lamda[2] = 2.0 * MY_PI * m;
        x2lamdaT(&unitk_lamda[0], &unitk_lamda[0]);
        sqk = unitk_lamda[0] * unitk_lamda[0] +
              unitk_lamda[1] * unitk_lamda[1] + unitk_lamda[2] * unitk_lamda[2];
        if (sqk <= gsqmx) {
          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk_lamda[0] * ug[kcount];
          eg[kcount][1] = 2.0 * unitk_lamda[1] * ug[kcount];
          eg[kcount][2] = 2.0 * unitk_lamda[2] * ug[kcount];
          vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
          vg[kcount][0] = 1.0 + vterm * unitk_lamda[0] * unitk_lamda[0];
          vg[kcount][1] = 1.0 + vterm * unitk_lamda[1] * unitk_lamda[1];
          vg[kcount][2] = 1.0 + vterm * unitk_lamda[2] * unitk_lamda[2];
          vg[kcount][3] = vterm * unitk_lamda[0] * unitk_lamda[1];
          vg[kcount][4] = vterm * unitk_lamda[0] * unitk_lamda[2];
          vg[kcount][5] = vterm * unitk_lamda[1] * unitk_lamda[2];
          kcount++;
        }
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = -kzmax; m <= kzmax; m++) {
      unitk_lamda[0] = 0.0;
      unitk_lamda[1] = 2.0 * MY_PI * l;
      unitk_lamda[2] = 2.0 * MY_PI * m;
      x2lamdaT(&unitk_lamda[0], &unitk_lamda[0]);
      sqk = unitk_lamda[1] * unitk_lamda[1] + unitk_lamda[2] * unitk_lamda[2];
      if (sqk <= gsqmx) {
        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 0.0;
        eg[kcount][1] = 2.0 * unitk_lamda[1] * ug[kcount];
        eg[kcount][2] = 2.0 * unitk_lamda[2] * ug[kcount];
        vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
        vg[kcount][0] = 1.0;
        vg[kcount][1] = 1.0 + vterm * unitk_lamda[1] * unitk_lamda[1];
        vg[kcount][2] = 1.0 + vterm * unitk_lamda[2] * unitk_lamda[2];
        vg[kcount][3] = 0.0;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = vterm * unitk_lamda[1] * unitk_lamda[2];
        kcount++;
      }
    }
  }

  // (0,0,m)

  for (m = 1; m <= kmax; m++) {
    unitk_lamda[0] = 0.0;
    unitk_lamda[1] = 0.0;
    unitk_lamda[2] = 2.0 * MY_PI * m;
    x2lamdaT(&unitk_lamda[0], &unitk_lamda[0]);
    sqk = unitk_lamda[2] * unitk_lamda[2];
    if (sqk <= gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = m;
      ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
      eg[kcount][0] = 0.0;
      eg[kcount][1] = 0.0;
      eg[kcount][2] = 2.0 * unitk_lamda[2] * ug[kcount];
      vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
      vg[kcount][0] = 1.0;
      vg[kcount][1] = 1.0;
      vg[kcount][2] = 1.0 + vterm * unitk_lamda[2] * unitk_lamda[2];
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
  }
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors
------------------------------------------------------------------------- */

void EwaldConp::allocate() {
  kxvecs = new int[kmax3d];
  kyvecs = new int[kmax3d];
  kzvecs = new int[kmax3d];

  ug = new double[kmax3d];
  memory->create(eg, kmax3d, 3, "ewald/conp:eg");
  memory->create(vg, kmax3d, 6, "ewald/conp:vg");

  sfacrl = new double[kmax3d];
  sfacim = new double[kmax3d];
  sfacrl_all = new double[kmax3d];
  sfacim_all = new double[kmax3d];
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors
------------------------------------------------------------------------- */

void EwaldConp::deallocate() {
  delete[] kxvecs;
  delete[] kyvecs;
  delete[] kzvecs;

  delete[] ug;
  memory->destroy(eg);
  memory->destroy(vg);

  delete[] sfacrl;
  delete[] sfacim;
  delete[] sfacrl_all;
  delete[] sfacim_all;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double EwaldConp::memory_usage() {
  double bytes = 3 * kmax3d * sizeof(int);
  bytes += (1 + 3 + 6) * kmax3d * sizeof(double);
  bytes += 4 * kmax3d * sizeof(double);
  bytes += nmax * 3 * sizeof(double);
  bytes += 2 * (2 * kmax + 1) * 3 * nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   group-group interactions
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   compute the Ewald total long-range force and energy for groups A and B
 ------------------------------------------------------------------------- */

void EwaldConp::compute_group_group(int groupbit_A, int groupbit_B,
                                    int AA_flag) {
  if (slabflag && triclinic)
    error->all(FLERR,
               "Cannot (yet) use K-space slab "
               "correction with compute group/group for triclinic systems");

  int i, k;

  if (!group_allocate_flag) {
    allocate_groups();
    group_allocate_flag = 1;
  }

  e2group = 0.0;     // energy
  f2group[0] = 0.0;  // force in x-direction
  f2group[1] = 0.0;  // force in y-direction
  f2group[2] = 0.0;  // force in z-direction

  // partial and total structure factors for groups A and B

  for (k = 0; k < kcount; k++) {
    // group A

    sfacrl_A[k] = 0.0;
    sfacim_A[k] = 0.0;
    sfacrl_A_all[k] = 0.0;
    sfacim_A_all[k] = 0.0;

    // group B

    sfacrl_B[k] = 0.0;
    sfacim_B[k] = 0.0;
    sfacrl_B_all[k] = 0.0;
    sfacim_B_all[k] = 0.0;
  }

  double *q = atom->q;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  int kx, ky, kz;
  double cypz, sypz, exprl, expim;

  // partial structure factors for groups A and B on each processor

  for (k = 0; k < kcount; k++) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];

    for (i = 0; i < nlocal; i++) {
      if (!((mask[i] & groupbit_A) && (mask[i] & groupbit_B)))
        if (AA_flag) continue;

      if ((mask[i] & groupbit_A) || (mask[i] & groupbit_B)) {
        cypz = cs[ky][1][i] * cs[kz][2][i] - sn[ky][1][i] * sn[kz][2][i];
        sypz = sn[ky][1][i] * cs[kz][2][i] + cs[ky][1][i] * sn[kz][2][i];
        exprl = cs[kx][0][i] * cypz - sn[kx][0][i] * sypz;
        expim = sn[kx][0][i] * cypz + cs[kx][0][i] * sypz;

        // group A

        if (mask[i] & groupbit_A) {
          sfacrl_A[k] += q[i] * exprl;
          sfacim_A[k] += q[i] * expim;
        }

        // group B

        if (mask[i] & groupbit_B) {
          sfacrl_B[k] += q[i] * exprl;
          sfacim_B[k] += q[i] * expim;
        }
      }
    }
  }

  // total structure factor by summing over procs

  MPI_Allreduce(sfacrl_A, sfacrl_A_all, kcount, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(sfacim_A, sfacim_A_all, kcount, MPI_DOUBLE, MPI_SUM, world);

  MPI_Allreduce(sfacrl_B, sfacrl_B_all, kcount, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(sfacim_B, sfacim_B_all, kcount, MPI_DOUBLE, MPI_SUM, world);

  const double qscale = qqrd2e * scale;
  double partial_group;

  // total group A <--> group B energy
  // self and boundary correction terms are in compute_group_group.cpp

  for (k = 0; k < kcount; k++) {
    partial_group =
        sfacrl_A_all[k] * sfacrl_B_all[k] + sfacim_A_all[k] * sfacim_B_all[k];
    e2group += ug[k] * partial_group;
  }

  e2group *= qscale;

  // total group A <--> group B force

  for (k = 0; k < kcount; k++) {
    partial_group =
        sfacim_A_all[k] * sfacrl_B_all[k] - sfacrl_A_all[k] * sfacim_B_all[k];
    if (wireflag != 2) {
      f2group[0] += eg[k][0] * partial_group;
      f2group[1] += eg[k][1] * partial_group;
    }
    if (slabflag != 2) f2group[2] += eg[k][2] * partial_group;
  }

  f2group[0] *= qscale;
  f2group[1] *= qscale;
  f2group[2] *= qscale;

  // 2d slab correction

  if (slabflag == 1 || slabflag == 3)
    slabcorr_groups(groupbit_A, groupbit_B, AA_flag);
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void EwaldConp::slabcorr_groups(int groupbit_A, int groupbit_B, int AA_flag) {
  if (slabflag == 3)
    error->all(FLERR,
               "Cannot (yet) use EW2D correction with compute group/group");

  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double zprd = domain->zprd;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double qsum_A = 0.0;
  double qsum_B = 0.0;
  double dipole_A = 0.0;
  double dipole_B = 0.0;
  double dipole_r2_A = 0.0;
  double dipole_r2_B = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (!((mask[i] & groupbit_A) && (mask[i] & groupbit_B)))
      if (AA_flag) continue;

    if (mask[i] & groupbit_A) {
      qsum_A += q[i];
      dipole_A += q[i] * x[i][2];
      dipole_r2_A += q[i] * x[i][2] * x[i][2];
    }

    if (mask[i] & groupbit_B) {
      qsum_B += q[i];
      dipole_B += q[i] * x[i][2];
      dipole_r2_B += q[i] * x[i][2] * x[i][2];
    }
  }

  // sum local contributions to get total charge and global dipole moment
  //  for each group

  double tmp;
  MPI_Allreduce(&qsum_A, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  qsum_A = tmp;

  MPI_Allreduce(&qsum_B, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  qsum_B = tmp;

  MPI_Allreduce(&dipole_A, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  dipole_A = tmp;

  MPI_Allreduce(&dipole_B, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  dipole_B = tmp;

  MPI_Allreduce(&dipole_r2_A, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  dipole_r2_A = tmp;

  MPI_Allreduce(&dipole_r2_B, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  dipole_r2_B = tmp;

  // compute corrections

  const double qscale = qqrd2e * scale;
  const double efact = qscale * MY_2PI / volume;

  e2group += efact * (dipole_A * dipole_B -
                      0.5 * (qsum_A * dipole_r2_B + qsum_B * dipole_r2_A) -
                      qsum_A * qsum_B * zprd * zprd / 12.0);

  // add on force corrections

  const double ffact = qscale * (-4.0 * MY_PI / volume);
  f2group[2] += ffact * (qsum_A * dipole_B - qsum_B * dipole_A);
}

/* ----------------------------------------------------------------------
   Wire-geometry correction term to dampen inter-wire interactions between
   periodically repeating wires.  Yields good approximation to 1D Ewald if
   adequate empty space is left between repeating wires (J. Mol. Struct.
   704, 101). x and y are non-periodic.
------------------------------------------------------------------------- */

void EwaldConp::wirecorr_groups(int groupbit_A, int groupbit_B, int AA_flag) {
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double qsum_A = 0.0;
  double qsum_B = 0.0;
  double xdipole_A = 0.0;
  double xdipole_B = 0.0;
  double xdipole_r2_A = 0.0;
  double xdipole_r2_B = 0.0;
  double ydipole_A = 0.0;
  double ydipole_B = 0.0;
  double ydipole_r2_A = 0.0;
  double ydipole_r2_B = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (!((mask[i] & groupbit_A) && (mask[i] & groupbit_B)))
      if (AA_flag) continue;

    if (mask[i] & groupbit_A) {
      qsum_A += q[i];
      xdipole_A += q[i] * x[i][0];
      xdipole_r2_A += q[i] * x[i][0] * x[i][0];
      ydipole_A += q[i] * x[i][1];
      ydipole_r2_A += q[i] * x[i][1] * x[i][1];
    }

    if (mask[i] & groupbit_B) {
      qsum_B += q[i];
      xdipole_B += q[i] * x[i][0];
      xdipole_r2_B += q[i] * x[i][0] * x[i][0];
      ydipole_B += q[i] * x[i][1];
      ydipole_r2_B += q[i] * x[i][1] * x[i][1];
    }
  }

  // sum local contributions to get total charge and global dipole moment
  //  for each group

  double tmp;
  MPI_Allreduce(&qsum_A, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  qsum_A = tmp;

  MPI_Allreduce(&qsum_B, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  qsum_B = tmp;

  if (fabs(qsum_A + qsum_B) > SMALL)
    error->all(FLERR,
               "Cannot (yet) use K-space wire "
               "correction with compute group/group for charged groups");

  MPI_Allreduce(&xdipole_A, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  xdipole_A = tmp;
  MPI_Allreduce(&ydipole_A, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  ydipole_A = tmp;

  MPI_Allreduce(&xdipole_B, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  xdipole_B = tmp;
  MPI_Allreduce(&ydipole_B, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  ydipole_B = tmp;

  MPI_Allreduce(&xdipole_r2_A, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  xdipole_r2_A = tmp;
  MPI_Allreduce(&ydipole_r2_A, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  ydipole_r2_A = tmp;

  MPI_Allreduce(&xdipole_r2_B, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  xdipole_r2_B = tmp;
  MPI_Allreduce(&ydipole_r2_B, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
  ydipole_r2_B = tmp;

  // compute corrections

  const double qscale = qqrd2e * scale;
  const double efact = qscale * MY_PI / volume;

  // TODO do math for non-neutral wire geometries of Ballenegger et al.

  e2group += efact * (xdipole_A * xdipole_B + ydipole_A * ydipole_B -
                      0.5 * (qsum_A * xdipole_r2_B + qsum_A * ydipole_r2_B +
                             qsum_B * xdipole_r2_A + qsum_B * ydipole_r2_A) -
                      qsum_A * qsum_B * xprd * xprd / 12.0 -
                      qsum_A * qsum_B * yprd * yprd / 12.0);

  // add on force corrections

  const double ffact = qscale * (-2.0 * MY_PI / volume);
  f2group[1] += ffact * (qsum_A * xdipole_B + qsum_A * ydipole_B -
                         qsum_B * xdipole_A - qsum_B * ydipole_A);
}

/* ----------------------------------------------------------------------
   allocate group-group memory that depends on # of K-vectors
------------------------------------------------------------------------- */

void EwaldConp::allocate_groups() {
  // group A

  sfacrl_A = new double[kmax3d];
  sfacim_A = new double[kmax3d];
  sfacrl_A_all = new double[kmax3d];
  sfacim_A_all = new double[kmax3d];

  // group B

  sfacrl_B = new double[kmax3d];
  sfacim_B = new double[kmax3d];
  sfacrl_B_all = new double[kmax3d];
  sfacim_B_all = new double[kmax3d];
}

/* ----------------------------------------------------------------------
   deallocate group-group memory that depends on # of K-vectors
------------------------------------------------------------------------- */

void EwaldConp::deallocate_groups() {
  // group A

  delete[] sfacrl_A;
  delete[] sfacim_A;
  delete[] sfacrl_A_all;
  delete[] sfacim_A_all;

  // group B

  delete[] sfacrl_B;
  delete[] sfacim_B;
  delete[] sfacrl_B_all;
  delete[] sfacim_B_all;
}

/* ----------------------------------------------------------------------
   compute b-vector of constant potential approach
 ------------------------------------------------------------------------- */

void EwaldConp::compute_vector(bigint *imat, double *vector) {
  update_eikr(false);

  int const nlocal = atom->nlocal;
  double *q = atom->q;
  std::vector<double> q_cos(kcount);
  std::vector<double> q_sin(kcount);

  for (int k = 0; k < kcount; k++) {
    int const kx = kxvecs[k];
    int const ky = kyvecs[k];
    int const kz = kzvecs[k];
    double q_cos_k = 0;
    double q_sin_k = 0;
    for (int i = 0; i < nlocal; i++) {
      if (imat[i] >= 0) continue;  // only electrode atoms
      double const cos_kxky =
          cs[kx][0][i] * cs[ky][1][i] - sn[kx][0][i] * sn[ky][1][i];
      double const sin_kxky =
          sn[kx][0][i] * cs[ky][1][i] + cs[kx][0][i] * sn[ky][1][i];
      double const cos_kr = cos_kxky * cs[kz][2][i] - sin_kxky * sn[kz][2][i];
      double const sin_kr = sin_kxky * cs[kz][2][i] + cos_kxky * sn[kz][2][i];

      q_cos_k += q[i] * cos_kr;
      q_sin_k += q[i] * sin_kr;
    }
    q_cos[k] = q_cos_k;
    q_sin[k] = q_sin_k;
  }

  MPI_Allreduce(MPI_IN_PLACE, &q_cos.front(), kcount, MPI_DOUBLE, MPI_SUM,
                world);
  MPI_Allreduce(MPI_IN_PLACE, &q_sin.front(), kcount, MPI_DOUBLE, MPI_SUM,
                world);

  for (int i = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;
    double bi = 0;
    for (int k = 0; k < kcount; k++) {
      int const kx = kxvecs[k];
      int const ky = kyvecs[k];
      int const kz = kzvecs[k];
      double const cos_kxky =
          cs[kx][0][i] * cs[ky][1][i] - sn[kx][0][i] * sn[ky][1][i];
      double const sin_kxky =
          sn[kx][0][i] * cs[ky][1][i] + cs[kx][0][i] * sn[ky][1][i];
      double const cos_kr = cos_kxky * cs[kz][2][i] - sin_kxky * sn[kz][2][i];
      double const sin_kr = sin_kxky * cs[kz][2][i] + cos_kxky * sn[kz][2][i];
      bi += 2 * ug[k] * (cos_kr * q_cos[k] + sin_kr * q_sin[k]);
      // different sign than fix_conp for now
    }
    vector[imat[i]] += bi;
  }
}

/* ----------------------------------------------------------------------
   compute b-vector EW3DC correction of constant potential approach
 ------------------------------------------------------------------------- */

void EwaldConp::compute_vector_corr(bigint *imat, double *vec) {
  update_eikr(false);
  boundcorr->vector_corr(imat, vec);
}

/* ----------------------------------------------------------------------
   compute individual interactions between all pairs of atoms in group A
   and B. see lammps_gather_atoms_concat() on how all sn and cs have been
   obtained.
 ------------------------------------------------------------------------- */

void EwaldConp::compute_matrix(bigint *imat, double **matrix) {
  update_eikr(false);
  int nlocal = atom->nlocal;
  int nprocs = comm->nprocs;

  double *csx, *csy, *csz, *snx, *sny, *snz;
  double *csx_all, *csy_all, *csz_all;
  double *snx_all, *sny_all, *snz_all;
  bigint *jmat, *jmat_local;
  // how many local group atoms owns each proc and how many in total
  bigint ngroup = 0;
  int ngrouplocal =
      std::count_if(&imat[0], &imat[nlocal], [](int i) { return i >= 0; });
  MPI_Allreduce(&ngrouplocal, &ngroup, 1, MPI_INT, MPI_SUM, world);

  // gather only subset of local sn and cs on each proc

  memory->create(csx, ngrouplocal * (kxmax + 1), "ewald/conp:csx");
  memory->create(snx, ngrouplocal * (kxmax + 1), "ewald/conp:snx");
  memory->create(csy, ngrouplocal * (kymax + 1), "ewald/conp:csy");
  memory->create(sny, ngrouplocal * (kymax + 1), "ewald/conp:sny");
  memory->create(snz, ngrouplocal * (kzmax + 1), "ewald/conp:snz");
  memory->create(csz, ngrouplocal * (kzmax + 1), "ewald/conp:csz");

  memory->create(jmat_local, ngrouplocal, "ewald/conp:jmat_local");

  // copy subsets of local sn and cn to new local group arrays
  // beeing as memory efficient as one can possibly be ...

  for (int i = 0, n = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;

    for (int k = 0; k <= kxmax; k++) {
      csx[k + n * (kxmax + 1)] = cs[k][0][i];
      snx[k + n * (kxmax + 1)] = sn[k][0][i];
    }
    for (int k = 0; k <= kymax; k++) {
      csy[k + n * (kymax + 1)] = cs[k][1][i];
      sny[k + n * (kymax + 1)] = sn[k][1][i];
    }
    for (int k = 0; k <= kzmax; k++) {
      csz[k + n * (kzmax + 1)] = cs[k][2][i];
      snz[k + n * (kzmax + 1)] = sn[k][2][i];
    }
    jmat_local[n] = imat[i];
    n++;
  }

  // TODO check if ((bigint) kxmax+1)*ngroup overflows ...

  memory->create(csx_all, ((bigint)kxmax + 1) * ngroup, "ewald/conp:csx_all");
  memory->create(snx_all, ((bigint)kxmax + 1) * ngroup, "ewald/conp:snx_all");
  memory->create(csy_all, ((bigint)kymax + 1) * ngroup, "ewald/conp:csy_all");
  memory->create(sny_all, ((bigint)kymax + 1) * ngroup, "ewald/conp:sny_all");
  memory->create(csz_all, ((bigint)kzmax + 1) * ngroup, "ewald/conp:csz_all");
  memory->create(snz_all, ((bigint)kzmax + 1) * ngroup, "ewald/conp:snz_all");

  memory->create(jmat, ngroup, "ewald/conp:jmat");

  int *recvcounts, *displs;  // TODO allgather requires int for displs but
                             // displs might overflow!
  memory->create(recvcounts, nprocs, "ewald/conp:recvcounts");
  memory->create(displs, nprocs, "ewald/conp:displs");

  // gather subsets global cs and sn
  int n = (kxmax + 1) * ngrouplocal;
  // TODO check if (kxmax+1)*ngrouplocal, etc.
  // overflows int n! typically kxmax small

  MPI_Allgather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++)
    displs[i] = displs[i - 1] + recvcounts[i - 1];
  MPI_Allgatherv(csx, n, MPI_DOUBLE, csx_all, recvcounts, displs, MPI_DOUBLE,
                 world);
  MPI_Allgatherv(&snx[0], n, MPI_DOUBLE, snx_all, recvcounts, displs,
                 MPI_DOUBLE, world);
  n = (kymax + 1) * ngrouplocal;
  MPI_Allgather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++)
    displs[i] = displs[i - 1] + recvcounts[i - 1];
  MPI_Allgatherv(&csy[0], n, MPI_DOUBLE, csy_all, recvcounts, displs,
                 MPI_DOUBLE, world);
  MPI_Allgatherv(&sny[0], n, MPI_DOUBLE, sny_all, recvcounts, displs,
                 MPI_DOUBLE, world);

  n = (kzmax + 1) * ngrouplocal;
  MPI_Allgather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++)
    displs[i] = displs[i - 1] + recvcounts[i - 1];
  MPI_Allgatherv(&csz[0], n, MPI_DOUBLE, csz_all, recvcounts, displs,
                 MPI_DOUBLE, world);
  MPI_Allgatherv(&snz[0], n, MPI_DOUBLE, snz_all, recvcounts, displs,
                 MPI_DOUBLE, world);

  // gather subsets global matrix indexing

  MPI_Allgather(&ngrouplocal, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++)
    displs[i] = displs[i - 1] + recvcounts[i - 1];
  MPI_Allgatherv(&jmat_local[0], ngrouplocal, MPI_LMP_BIGINT, jmat, recvcounts,
                 displs, MPI_LMP_BIGINT, world);

  memory->destroy(displs);
  memory->destroy(recvcounts);

  memory->destroy(jmat_local);

  // aij for each atom pair in groups; first loop over i,j then over k to
  // reduce memory access
  for (int i = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;

    for (bigint j = 0; j < ngroup; j++) {
      // matrix is symmetric, skip upper triangular matrix
      if (jmat[j] > imat[i]) continue;

      double aij = 0.0;

      for (int k = 0; k < kcount; k++) {
        // local  indexing  cs[k_idim][idim][i]       <>
        // csx_all[i+k*ngrouplocal+displs[comm->me]]]

        // anyway, use local sn and cs for simplicity

        int const kx = kxvecs[k];
        int const ky = kyvecs[k];
        int const kz = kzvecs[k];
        int const sign_ky = (ky > 0) - (ky < 0);
        int const sign_kz = (kz > 0) - (kz < 0);

        double cos_kxky =
            cs[kx][0][i] * cs[ky][1][i] - sn[kx][0][i] * sn[ky][1][i];
        double sin_kxky =
            sn[kx][0][i] * cs[ky][1][i] + cs[kx][0][i] * sn[ky][1][i];

        double const cos_kxkykz_i =
            cos_kxky * cs[kz][2][i] - sin_kxky * sn[kz][2][i];
        double const sin_kxkykz_i =
            sin_kxky * cs[kz][2][i] + cos_kxky * sn[kz][2][i];

        // global indexing  csx_all[kx+j*(kxmax+1)]  <>  csx_all[kx][j]

        int const kxj = kx + j * (kxmax + 1);
        int const kyj = abs(ky) + j * (kymax + 1);
        int const kzj = abs(kz) + j * (kzmax + 1);

        cos_kxky =
            csx_all[kxj] * csy_all[kyj] - snx_all[kxj] * sny_all[kyj] * sign_ky;
        sin_kxky =
            snx_all[kxj] * csy_all[kyj] + csx_all[kxj] * sny_all[kyj] * sign_ky;

        double const cos_kxkykz_j =
            cos_kxky * csz_all[kzj] - sin_kxky * snz_all[kzj] * sign_kz;
        double const sin_kxkykz_j =
            sin_kxky * csz_all[kzj] + cos_kxky * snz_all[kzj] * sign_kz;

        aij += 2.0 * ug[k] *
               (cos_kxkykz_i * cos_kxkykz_j + sin_kxkykz_i * sin_kxkykz_j);
      }
      matrix[imat[i]][jmat[j]] += aij;
      if (imat[i] != jmat[j]) matrix[jmat[j]][imat[i]] += aij;
    }

    if ((i + 1) % 500 == 0) printf("(%d/%d) on %d\n", i + 1, nlocal, comm->me);
  }

  memory->destroy(jmat);
  memory->destroy(csx_all);
  memory->destroy(snx_all);
  memory->destroy(csy_all);
  memory->destroy(sny_all);
  memory->destroy(csz_all);
  memory->destroy(snz_all);
  memory->destroy(csx);
  memory->destroy(snx);
  memory->destroy(csy);
  memory->destroy(sny);
  memory->destroy(csz);
  memory->destroy(snz);
}

/* ----------------------------------------------------------------------
   compute individual corrections between all pairs of atoms in group A
   and B. see lammps_gather_atoms_concat() on how all sn and cs have been
   obtained.
 ------------------------------------------------------------------------- */

void EwaldConp::compute_matrix_corr(bigint *imat, double **matrix) {
  update_eikr(false);
  boundcorr->matrix_corr(imat, matrix);
}

/* ---------------------------------------------------------------------- */

void EwaldConp::update_eikr() { update_eikr(false); }

/* ---------------------------------------------------------------------- */

void EwaldConp::update_eikr(bool force_update) {
  if (eikr_step < update->ntimestep || force_update) {
    // extend size of per-atom arrays if necessary
    if (atom->nmax > nmax) {
      memory->destroy(ek);
      memory->destroy3d_offset(cs, -kmax_created);
      memory->destroy3d_offset(sn, -kmax_created);
      nmax = atom->nmax;
      memory->create(ek, nmax, 3, "ewald/conp:ek");
      memory->create3d_offset(cs, -kmax, kmax, 3, nmax, "ewald/conp:cs");
      memory->create3d_offset(sn, -kmax, kmax, 3, nmax, "ewald/conp:sn");
      kmax_created = kmax;
    }
    eikr_step = update->ntimestep;
    if (triclinic == 0)
      eik_dot_r();
    else
      eik_dot_r_triclinic();
  }
}

