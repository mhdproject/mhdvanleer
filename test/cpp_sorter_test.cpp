

#include "riemann.h"
#include "gtest/gtest.h"

TEST(sgn_test, int_arr_sort
) {
  double num = -7;
  double result = sgn(num);
  EXPECT_EQ(result,
            -1);
}

TEST(riemann_test, int_arr_sort
) {
  double lhs[] = {1, 0, 0, 0, 5, 1, 1, 1};
  double rhs[] = {1, 0, 0, 0, 5, 1, 1, 1};
  double flux[7];
  double res_state[7];
  double exp_flux[] = {0, 0, 0, 0, 0, 0, 0, 0};
  int idir = 1;
  Riemann riem;
  double result = riem.solver(lhs, rhs, flux, res_state, idir);
  assert(result == 0);
  EXPECT_FLOAT_EQ(exp_flux[0], flux[0]);
}

TEST(hlld_test, int_arr_sort
) {
  double lhs[] = {1, 0, 0, 0, 1, 1, 1, 1};
  double rhs[] = {1, 0, 0, 0, 1, 1, 1, 1};
  double flux[7];
  double exp_flux[] = {0, 0, 0, 0, 0, 0, 0, 0};
  Hlld hlld;
  double result = hlld.solver(lhs, rhs, flux);
  assert(result == 0);
  EXPECT_FLOAT_EQ(exp_flux[0], flux[0]);

}

