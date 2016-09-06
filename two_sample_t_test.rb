# -*- coding: utf-8 -*-

=begin
= Statistical_Analysis
==Index
* ((<module Statistic>))
  * ((<rcov_sample>))

=module Statistic
    statistic functions made by eriko, modified by koshiro.
    in $HOME/lib/ruby/1.8 or $HOME/lib/ruby/1.9

==Module Functions
---kurtosis_sample( fx ) => kurtosis_s
    標本尖度を計算する。

    ARGUMENTS
    *  fx : NArray or NArrayMiss
    RETURN VALUES
    * kurtosis_s : float, 標本尖度
=end

require 'narray'
require 'gsl'
module Statistic

module_function

  def two_sample_t_val( xmean, ymean, xvariance, yvariance, nx, ny )
    sp2=((nx-1).to_f*xvariance.to_f + (ny-1).to_f*yvariance.to_f)/(nx+ny-2).to_f
    bunbo2=1.0/nx.to_f+1.0/ny.to_f
    t_val=(xmean.to_f-ymean.to_f)/(Math::sqrt(sp2)*Math::sqrt(bunbo2))

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
