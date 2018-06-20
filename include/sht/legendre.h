#pragma once

#include <cassert>
#include <cmath>
#include <vector>
#include <cstddef>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


namespace sht
{
  template<class RealScalar>
  std::vector<RealScalar> legendre_p(unsigned int const l_max, unsigned int const m_max, RealScalar const cos_theta)
  {
    assert(l_max >= 0);
    assert(m_max <= l_max);
    assert((cos_theta >= -1.) && (cos_theta <= 1.));

    using boost::math::constants::pi;
    RealScalar const sin_theta = sqrt(RealScalar(1. - cos_theta * cos_theta));

    std::size_t const n_plms = ((m_max + 1) * (m_max + 2)) / 2 + (m_max + 1) * (l_max - m_max);
    std::vector<RealScalar> p_lm(n_plms);

    std::size_t ind(0);
    RealScalar p_mm1mm1 = sqrt(RealScalar(1. / (4. * pi<RealScalar>())));
    for (unsigned m = 0; m <= m_max; ++m)
    {
      RealScalar p_mm = p_mm1mm1;
      if (m > 0)
      {
        p_mm = sin_theta * sqrt(RealScalar((2. * m + 1.) / (2. * m))) * p_mm1mm1;
        p_mm1mm1 = p_mm;
      }
      RealScalar alpha_lm1m(1);
      for (unsigned l = m; l <= l_max; ++l)
      {
        if (l == m)
        {
          p_lm[ind] = p_mm;
          ++ind;
        }
        else
        {
          RealScalar const alpha_lm = sqrt(RealScalar((4. * l * l - 1.) / (l * l - m * m)));
          if (l == m + 1)
          {
            RealScalar const p_mmp1 = alpha_lm * cos_theta * p_mm;
            p_lm[ind] = p_mmp1;
            ++ind;
          }
          else
          {
            RealScalar const beta_lm = -alpha_lm / alpha_lm1m;
            p_lm[ind] = alpha_lm * cos_theta * p_lm[ind - 1] + beta_lm * p_lm[ind - 2];
            ++ind;
          }
          alpha_lm1m = alpha_lm;
        }
      }
    }
    return p_lm;
  }
}