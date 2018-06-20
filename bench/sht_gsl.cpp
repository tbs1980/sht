#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <gsl/gsl_sf_legendre.h>
#include "../include/sht/legendre.h"

void print_ylm_values_to_file(unsigned const l_max, unsigned const m_max, double const cos_theta)
{
  std::vector<double> plms = sht::legendre_p<double>(l_max, m_max, cos_theta);
  std::size_t const legendre_array_size = gsl_sf_legendre_array_n(l_max);
  std::vector<double> legendre_array(legendre_array_size);
  int res = gsl_sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, l_max, cos_theta, &legendre_array[0]);

  std::ofstream out_file("sht_gsl_comparison.txt", std::ios::trunc);
  if (out_file.is_open())
  {
    out_file << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    std::size_t ind_sht(0);
    std::size_t ind_gsl(0);
    for (unsigned int m = 0; m <= m_max; ++m)
    {
      for (unsigned int l = m; l <= l_max; ++l)
      {
        ind_gsl = gsl_sf_legendre_array_index(l, m);
        out_file << plms[ind_sht] << "," << legendre_array[ind_gsl] << "," << plms[ind_sht] - legendre_array[ind_gsl] << std::endl;
        ++ind_sht;
      }
    }
  }
}

int main()
{
  unsigned const l_max = 128;
  unsigned const m_max = 128;
  double const cos_theta = 0.5;
  print_ylm_values_to_file(l_max, m_max, cos_theta);
  return 0;
}
