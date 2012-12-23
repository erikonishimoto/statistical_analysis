#Statistical_Analysis
#====================
#
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
    * p : 有意水準(default 0.05)
    * i : 片側検定の場合 1, 両側検定の場合 2 (default 2)
    RETURN VALUES
    * rcor : 相関係数

---weighted_regression( x, y, w )
    重み付き回帰直線

