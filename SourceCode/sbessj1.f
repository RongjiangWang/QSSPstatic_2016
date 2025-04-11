      double precision function sbessj1(x)
      implicit none
      double precision x
c
      sbessj1=(dsin(x)/x-dcos(x))/x
      return
      end