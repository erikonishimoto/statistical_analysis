# -*- coding: utf-8 -*-

require 'narray'
require 'gsl'
module Statistic

module_function

  def two_sample_t_val( xmean, ymean, xvariance, yvariance, nx, ny )
    bunbo1=(((nx-1)*xvariance + (ny-1)*yvariance)/(nx+ny-2))**0.5
    bunbo2=nx**(-1)+ny**(-1)
    t_val=(xmean-ymean)*(bunbo1*bunbo2)**(-1)

    return t_val
  end

  def two_sample_t_test( xmean, ymean, xvariance, yvariance, nx, ny, p=0.05, i=2 )
    t_val=two_sample_t_val( xmean, ymean, xvariance, yvariance, nx, ny )

    nu=nx+ny-2
    t_a=GSL::Cdf.tdist_Qinv(p/i.to_f,nu)

    if t_val.abs > t_a
      state=true
    elsif t_val.abs <= t_a
      state=false
    end

    return state, t_val
  end

end
