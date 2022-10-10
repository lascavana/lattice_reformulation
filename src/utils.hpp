#include <vector>

// SCIP
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

// NTL library
#include <NTL/ZZ.h>
#include <NTL/LLL.h>
#include <NTL/mat_ZZ.h>

NTL::ZZ double2zz(double value);

int double2int(double value);

/* get new variable bounds */
SCIP_RETCODE get_new_varbounds(
  std::vector<std::vector<int>> basis,
  std::vector<int> lhs,
  std::vector<int> rhs,
  std::vector<int> &upper_bounds,
  std::vector<int> &lower_bounds
);

void print_reformulation(
  SCIP* scip,
  const char *filename,
  std::vector<std::vector<int>> basis,
  std::vector<int> x0,
  std::vector<double> upper,
  std::vector<double> lower,
  std::vector<double> objfun,
  bool maximization
);

void print_for_ls(
  std::string filename,
  std::vector<std::vector<int>> basis,
  std::vector<int> lhs,
  std::vector<int> rhs,
  std::vector<int> upper,
  std::vector<int> lower
);






