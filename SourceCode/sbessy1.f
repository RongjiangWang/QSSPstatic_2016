      double precision function sbessy1(x)
      implicit none
      double precision x
c
      sbessy1=(-dcos(x)/x-dsin(x))/x
      return
      end