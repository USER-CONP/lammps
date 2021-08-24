
#ifdef FIX_CLASS

FixStyle(charge_update, FixChargeUpdate)

#else

#ifndef LMP_FIX_CHARGE_UPDATE_H
#define LMP_FIX_CHARGE_UPDATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixChargeUpdate : public Fix {
 public:
  FixChargeUpdate(class LAMMPS *, int, char **);
  ~FixChargeUpdate();
  int setmask();
  // void post_constructor();
  void init();
  // void setup(int);
  void setup_post_neighbor();
  // void setup_pre_force(int);
  void setup_pre_reverse(int, int);
  void pre_force(int);
  void pre_reverse(int, int);
  double compute_scalar();
  // void post_run();
  // void setup_pre_force_respa(int,int);
  // void pre_force_respa(int,int,int);
  // void set_arrays(int);
  // void write_restart(FILE *);
  // void restart(char *);

 private:
  FILE *f_inv, *f_vec, *f_mat;  // files for capacitance, eleastance and vector
  std::string input_file_inv, input_file_mat;
  class ComputeConpMatrix *array_compute;
  class ComputeConpVector *vector_compute;
  std::vector<int> groups, group_bits;
  std::vector<double> group_psi;
  std::vector<double> psi;
  std::vector<std::vector<double> > elastance;
  bigint ngroup;
  std::vector<tagint> taglist, taglist_bygroup, group_idx;
  std::vector<int> tag_to_iele;
  bool read_inv, read_mat;
  bool symm;  // symmetrize elastance for charge neutrality
  double eta;
  double update_time, mult_time;
  void create_taglist();
  void invert(std::vector<std::vector<double> >);
  void symmetrize();
  double gausscorr(int, bool);
  void update_charges();
  double potential_energy(int, std::vector<int>);
  double self_energy(int);
  std::vector<int> local_to_matrix();
  void write_to_file(FILE *, std::vector<tagint>,
                     std::vector<std::vector<double> >);
  std::vector<std::vector<double> > read_from_file(std::string input_file);
};

}  // namespace LAMMPS_NS

#endif
#endif
