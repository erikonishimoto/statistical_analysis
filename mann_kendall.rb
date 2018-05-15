require "narray"
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

end
