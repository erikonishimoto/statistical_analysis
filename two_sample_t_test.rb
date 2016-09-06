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

  #-- 等分散の差の検定
  def calc_f0( xvariance, yvariance )
    f0=xvariance/yvariance

    return f0
  end

  def two_sample_f_test( xvariance, yvariance, nx, ny, p=0.05, i=2 )
    f0=calc_f0( xvariance, yvariance )

    f_val1=GSL::Cdf.fdist_Qinv( p/i.to_f, nx-1, ny-1 )
    f_val2=GSL::Cdf.fdist_Qinv( 1.0-p/i.to_f, nx-1, ny-1 )

    #------ 2組の標本の母分散に差がない(帰無仮説)
    state=true
    #-- 片側検定
    if i==1
      state=false if f0<f_val1   #--2組の標本の母分散に差がある
    #-- 両側検定
    elsif i==2
      state=false if f0<f_val1 or f_val2<f0   #--2組の標本の母分散に差がある
    end

    return state, f0, f_val1, f_val2
  end

  #-- 母平均の差の検定 (分散が大きく違わない場合)
  def calc_two_sample_t_val( xmean, ymean, xvariance, yvariance, nx, ny )
    sp2=((nx-1).to_f*xvariance.to_f + (ny-1).to_f*yvariance.to_f)/(nx+ny-2).to_f
    bunbo2=1.0/nx.to_f+1.0/ny.to_f
    t0=(xmean.to_f-ymean.to_f)/(Math::sqrt(sp2)*Math::sqrt(bunbo2))

    return t0
  end

  def two_sample_t_test( xmean, ymean, xvariance, yvariance, nx, ny, p=0.05, i=2 )
    t0=calc_two_sample_t_val( xmean, ymean, xvariance, yvariance, nx, ny )

    #---- i=2; 両側検定
    nu=nx+ny-2
    t_val=GSL::Cdf.tdist_Qinv(p/i.to_f,nu)

    #---- 帰無仮説: 2組の標本の平均に差がない
    state=true
    state=false if t0.abs > t_val   #--- 棄却、差がある

    return state, t0, t_val
  end

  #-- 母平均の差の検定 (分散が大きく違う場合)
  def welch_t_test( xmean, ymean, xvariance, yvariance, nx, ny, p=0.05, i=2 )
    t0=( xmean.to_f - ymean.to_f )/Math::sqrt( xvariance/nx.to_f + yvariance/ny.to_f )

    #-- 自由度k
    c=( xvariance/nx.to_f )/( xvariance/nx.to_f + yvariance/ny.to_f )
    k=1.0/( c**2/(nx-1).to_f + (1.0-c**2)/(ny-1).to_f )

    #-- t0 が自由度k のt分布に近似的に従う
    t_val=GSL::Cdf.tdist_Qinv( p/i.to_f, k )

    #---- 帰無仮説: 2組の標本の平均に差がない
    state=true
    state=false if t0.abs > t_val   #--- 棄却、差がある

    return state, t0, t_val
  end
end
