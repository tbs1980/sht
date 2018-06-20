#include <gtest/gtest.h>
#include "../include/sht/legendre.h"
#include <vector>
#include <cmath>

namespace sht
{
  TEST(legendre_p, n_plms)
  {
    unsigned int const l_max = 10;
    unsigned int const m_max = 8;
    double const cos_theta = 0.5;
    std::vector<double> plms = sht::legendre_p<double>(l_max, m_max, cos_theta);

    std::size_t n_plms = 0;
    for (unsigned m = 0; m <= m_max; ++m)
    {
      for (unsigned l = m; l <= l_max; ++l)
      {
        ++n_plms;
      }
    }
    EXPECT_EQ(n_plms, plms.size());
  }

  TEST(legendre_p, p00)
  {
    double const pi(3.141592653589793238462643383279502884197);
    unsigned int const l_max = 0;
    unsigned int const m_max = 0;
    double const cos_theta = 0.5;
    std::vector<double> plms = sht::legendre_p<double>(l_max, m_max, cos_theta);
    EXPECT_EQ(plms.size(), 1);
    ASSERT_DOUBLE_EQ(plms[0], 1. / std::sqrt(4.*pi));
  }

  TEST(legendre_p, p_10_p_11)
  {
    double const pi(3.141592653589793238462643383279502884197);
    unsigned int const l_max = 2;
    unsigned int const m_max = 2;
    double const cos_theta = 0.5;
    double const sin_theta = std::sqrt(1. - cos_theta * cos_theta);
    std::vector<double> plms = sht::legendre_p<double>(l_max, m_max, cos_theta);
    EXPECT_EQ(plms.size(), 6);
    ASSERT_DOUBLE_EQ(plms[0], 1. / std::sqrt(4.*pi)); // p00
    ASSERT_DOUBLE_EQ(plms[1], std::sqrt(3.) * cos_theta / std::sqrt(4. * pi)); // p10
    ASSERT_DOUBLE_EQ(plms[2], std::sqrt(5.) * (3. * cos_theta * cos_theta - 1.) / std::sqrt(16. * pi)); // p20
    ASSERT_DOUBLE_EQ(plms[3], std::sqrt(3.) * sin_theta / std::sqrt(4. * 2. * pi)); //p11
    ASSERT_DOUBLE_EQ(plms[4], std::sqrt(15.) * sin_theta * cos_theta / std::sqrt(4. * 2. * pi)); // p21
    ASSERT_DOUBLE_EQ(plms[5], std::sqrt(15.) * sin_theta * sin_theta / std::sqrt(16. * 2. * pi)); // p22
  }
}