#Statistical_Analysis
#====================
#
=module Statistic
    statistic functions made by eriko, modified by koshiro.
    in $HOME/lib/ruby/1.8 or $HOME/lib/ruby/1.9.1

==Module Functions
---rcov_sample(fx, fy) => cov
    $BI8K\6&J,;6$rJV$9(B

    ARGUMENTS
    * fx : Array or NArray
    * fy : Array or NArray
    RETURN VALUES
    * cov: Real, $B6&J,;6(B

---rcov(fx, fy) => cov
    $BITJP6&J,;6$rJV$9(B

    ARGUMENTS
    * fx : Array or NArray
    * fy : Array or NArray
    RETURN VALUES
    * cov: Real, $B6&J,;6(B

---rcor(fx, fy) => a,b,r
    $BAj4X2r@O$r9T$J$&(B. 

    ARGUMENTS
    * fx : Array or NArray
    * fy : Array or NArray
    RETURN VALUES
    * a  : Real
    * b  : Real
      * y = a*x + b
    * r  : Real, $BAj4X78?t(B

---student_test(rcor,n,p=0.05,i=2) => t0,st,state
    t$B8!Dj$r9T$J$&(B. 

    ARGUMENTS
    * rcor : $BAj4X78?t(B
    * n : $B%G!<%??t(B
    * p : $BM-0U?e=`(B(default 0.05)
    * i : $BJRB&8!Dj$N>l9g(B 1, $BN>B&8!Dj$N>l9g(B 2 (default 2)
    RETURN VALUES
    * t0 : $BAj4X78?t(Brcor, $B%G!<%??t(Bn$B$G$N(Bt$BCM(B
    * st : $B<+M3EY(Bn-2, $BM-0U?e=`(Bp$B$G$N(Bt$BCM(B
    * state : |t0|>=st $B$GM-0U$@$C$?>l9g(Btrue, |t0|<st $B$GM-0U$G$O$J$+$C$?>l9g(Bfalse

---significance_level(n,p=0.05,i=2) => rcor
    $B$"$kI8K\?t$K$*$$$FM-0U?e=`$rK~$9Aj4X78?t$r5a$a$k(B

    ARGUMENTS
    * n : $B%G!<%??t(B
    * p : $BM-0U?e=`(B(default 0.05)
    * i : $BJRB&8!Dj$N>l9g(B 1, $BN>B&8!Dj$N>l9g(B 2 (default 2)
    RETURN VALUES
    * rcor : $BAj4X78?t(B

---weighted_regression( x, y, w )
    $B=E$_IU$-2s5"D>@~(B

