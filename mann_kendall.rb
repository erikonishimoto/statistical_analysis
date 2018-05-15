require "narray"
require "numru/gphys"
include NumRu

module Statistic

module_function

  #<<< mann kendall test >>>
  # see 松山・谷本(2008) p.46
  def mann_kendall( val )
    nlength=val.length

    #---- calc rank static n_i
    icount=NArray.sfloat(nlength)
    nlength.times{|i|
      if i==nlength-1
        icount[i]=0
      else
        vali=val[i]
        icount[i]=(val[i+1..-1].ge(vali)).where.length
      end
    }#- nlength.times{|i|
    #p "icount",icount

    #---- calc tau
    tau=4.0*(icount/(nlength*(nlength-1))).sum - 1.0
    #p "tau",tau

    #---- calc tau_g
    t_g=1.960   #-- for 5% !!!!!
    tau_g=t_g*((4*nlength+10.0)/(9*nlength*(nlength-1)))**0.5

    return tau,tau_g
  end
  
  def detrend_linear_trend( gphys )
    load  "/home/eriko/db/git/ruby/statistical_analysis/statistic.rb"
    fx=gphys.coord(0).val-gphys.coord(0).mean
    fy=gphys.val.to_na
    a,b,r=Statistic::rcor(fx, fy)

    linear_trend=a*fx+b

    detrend=gphys.copy-linear_trend
    detrend=detrend.rename("detrend")
    
    return detrend,linear_trend,r
  end

end
