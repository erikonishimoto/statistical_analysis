# -*- coding: utf-8 -*-

=begin
= Statistical_Analysis
==Index
* ((<module Statistic>))
  * ((<rcov_population>))
  * ((<rcov>))
  * ((<rcor>))
  * ((<student_test>))
  * ((<significance_level>))
  * ((<weighted_regression>))
  * ((<normal_distribution>))
  * ((<expectation>))
  * ((<skewness>))
  * ((<kurtosis>))

=module Statistic
    statistic functions made by eriko, modified by koshiro.
    in $HOME/lib/ruby/1.8 or $HOME/lib/ruby/1.9.1

==Module Functions
---rcov_sample(fx, fy) => cov
    標本共分散を返す

    ARGUMENTS
    * fx : Array or NArray
    * fy : Array or NArray
    RETURN VALUES
    * cov: Real, 共分散

---rcov(fx, fy) => cov
    不偏共分散を返す

    ARGUMENTS
    * fx : Array or NArray
    * fy : Array or NArray
    RETURN VALUES
    * cov: Real, 共分散

---rcor(fx, fy) => a,b,r
    相関解析を行なう. 

    ARGUMENTS
    * fx : Array or NArray
    * fy : Array or NArray
    RETURN VALUES
    * a  : Real
    * b  : Real
      * y = a*x + b
    * r  : Real, 相関係数

---student_test(rcor,n,p=0.05,i=2) => t0,st,state
    t検定を行なう. 

    ARGUMENTS
    * rcor : 相関係数
    * n : データ数
    * p : 有意水準(default 0.05)
    * i : 片側検定の場合 1, 両側検定の場合 2 (default 2)
    RETURN VALUES
    * t0 : 相関係数rcor, データ数nでのt値
    * st : 自由度n-2, 有意水準pでのt値
    * state : |t0|>=st で有意だった場合true, |t0|<st で有意ではなかった場合false

---significance_level(n,p=0.05,i=2) => rcor
    ある標本数において有意水準を満す相関係数を求める

    ARGUMENTS
    * n : データ数
    * p : 有意水準 (default 0.05)
    * i : 片側検定の場合 1, 両側検定の場合 2 (default 2)
    RETURN VALUES
    * rcor : 相関係数

---weighted_regression( x, y, w ) => a,b
    重み付き回帰直線

    ARGUMENTS
    * fx : NArray or NArrayMiss
    * fy : NArray or NArrayMiss
    *  w : weighted function
    RETURN VALUES
    * a  : Real
    * b  : Real
      * y = a*x + b

---normal_distribution( x, mu, sigma ) => w
    正規分布N(mu,sigma^2)を返す。平均値mu=0、標準偏差sigma=1の場合は標準正規分布N(0,1)を返す。

    ARGUMENTS
    *   x : NArray or NArrayMiss
    *  mu : Float, mean value ( default 0 )
    * sigma : Float, standard deviation ( default 1 )
    RETURN VALUES
    *   w : NArray, Normal distribution

---expectation( fx, pdf ) => gx
    期待値を計算する。

    ARGUMENTS
    *  fx : NArray or NarrayMiss, 確率変数
    * pdf : NArray or NarrayMiss, 確率密度関数
    RETURN VALUES
    *  gx : float, 期待値

---skewness( fx, pdf ) => skew
    歪度を計算する。

    ARGUMENTS
    *  fx : NArray or NArrayMiss, 確率変数
    * pdf : NArray or NArrayMiss, 確率密度関数
    RETURN VALUES
    * skew : float, 歪度

---kurtosis( fx, pdf ) => kurtosis
    尖度を計算する。

    ARGUMENTS
    *  fx : NArray or NArrayMiss, 確率変数
    * pdf : NArray or NArrayMiss, 確率密度関数
    RETURN VALUES
    * kurtosis : float, 尖度

=end

require 'narray'
require 'gsl'
module Statistic

  module_function

  #<<< sample covariance >>>
  def rcov_sample(x,y)
    if !x.is_a?(NArray)
      if !x.is_a?(Array)
        raise 'array must be NArray or Array'
      end
    end
    if !y.is_a?(NArray)
      if !y.is_a?(Array)
        raise 'array must be NArray or Array'
      end
    end
    if x.length != y.length
      raise 'Must be same length'
    end

    cov = (x*y).mean - (x.mean)*(y.mean)
    return cov
  end

  #<<< unbiased covariance >>>
  def rcov(x,y)
    if !x.is_a?(NArray)
      if !x.is_a?(Array)
        raise 'array must be NArray or Array'
      end
    end
    if !y.is_a?(NArray)
      if !y.is_a?(Array)
        raise 'array must be NArray or Array'
      end
    end
    if x.length != y.length
      raise 'Must be same length'
    end

    cov = NMath::covariance(x,y)
    return cov
  end

  #<<< linear regression analysis >>>
  def rcor(x,y)
    if !x.is_a?(NArray)
      if !x.is_a?(Array)
        raise 'array must be NArray or Array'
      end
    end
    if !y.is_a?(NArray)
      if !y.is_a?(Array)
        raise 'array must be NArray or Array'
      end
    end

    sxy = NMath::covariance(x,y)
    sx2 = NMath::covariance(x,x)
    sy2 = NMath::covariance(y,y)

    #-- y = a*x + b
    a = sxy/sx2
    b = y.mean - a * x.mean
    r = sxy/(sqrt(sx2) * sqrt(sy2))
    
    return a,b,r
  end

  #<<< student t test >>>
  # t検定を行なう
  def student_test(rcor,n,p=0.05,i=2)
    nu=n-2
    t0 = (rcor*Math::sqrt(nu/(1.0-rcor*rcor))).abs
    st = GSL::Cdf.tdist_Qinv(p/i.to_f,nu)

    if t0>=st
      state = true
    else
      state = false
    end

    return t0,st,state
  end

  #<<< significance level >>>
  # 標本数nに対して有意水準pを満す相関係数を求める
  def significance_level(n,p=0.05,i=2)
    nu=n-2
    st = GSL::Cdf.tdist_Qinv(p/i.to_f,nu)
    r = st/Math::sqrt(nu+st**2)

    return r
  end

  #<<< weighted regression >>>
  #-- 重み付き回帰直線
  def weighted_regression( fx, fy, w )
    wx = (w*fx).sum
    wy = (w*fy).sum
    wxx = (w*fx*fx).sum
    wxy = (w*fx*fy).sum
    wsum = w.sum

    a = ( wxy/wx - wy/wsum )/( wxx/wx - wx/wsum )
    b = ( wy - a*wx ) / wsum

    return a,b
  end

  #<<< Normal distribution / Gaussian distribution >>>
  #-- 正規分布
  def normal_distribution( x, mu=0.0, sigma=1.0 )
    pi = NMath::PI

    xx = -( (x-mu)**2.0 )/( 2.0*sigma**2 )
    exp_xx = NMath::exp( xx )

    w = 1.0/NMath::sqrt( 2.0*pi*sigma**2 ) * exp_xx
    return w
  end

  #<<< expectation >>>
  #-- 期待値
  def expectation( fx, pdf )
    gx = fx.mul_add( pdf, 0 )
  end

  #<<< Skewness >>>
  #-- 歪度
  def skewness( fx, pdf )
    mu = fx.mean
    sigma = fx.stddev

    e_fx3 = expectation( fx**3, pdf )
    e_fx2 = 3.0 * mu * expectation( fx**2, pdf )
    e_fx1 = 2.0 * mu**3
    skew = ( e_fx3 + e_fx2 + e_fx1 ) / sigma**3

    return skew
  end

  #<<< Kurtosis >>>
  #-- 尖度
  def kurtosis( fx, pdf )
    mu = fx.mean
    sigma = fx.stddev

    e_fx4 = expectation( fx**4, pdf )
    e_fx3 =  -4.0 * mu * expectation( fx**3, pdf )
    e_fx2 =   6.0 * mu**2 * expectation( fx**2, pdf )
    e_fx1 =   3.0 * mu**4
    kurtosis = ( e_fx4 + e_fx3 + e_fx2 + e_fx1 )/sigma**4

    return kurtosis
  end

end

