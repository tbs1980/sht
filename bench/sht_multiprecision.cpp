#include <fstream>
#include <iomanip>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "../include/sht/legendre.h"

void print_ylm_values_to_file(unsigned const l_max, unsigned const m_max, double const cos_theta)
{
  using boost::multiprecision::cpp_dec_float_50;
  std::vector<double> plms = sht::legendre_p<double>(l_max, m_max, cos_theta);
  std::vector<cpp_dec_float_50> plms_mp = sht::legendre_p<cpp_dec_float_50>(l_max, m_max, cos_theta);

  std::ofstream out_file("sht_mp_comparison.txt", std::ios::trunc);
  if (out_file.is_open())
  {
    out_file << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    std::size_t ind_sht(0);
    for (unsigned int m = 0; m <= m_max; ++m)
    {
      for (unsigned int l = m; l <= l_max; ++l)
      {
        out_file << plms[ind_sht] << "," << plms_mp[ind_sht] << std::endl;
        ++ind_sht;
      }
    }
  }
}

int main()
{
  unsigned const l_max = 15;
  unsigned const m_max = 15;
  double const cos_theta = 0.5;
  print_ylm_values_to_file(l_max, m_max, cos_theta);
  return 0;
}
