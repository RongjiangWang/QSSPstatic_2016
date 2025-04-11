      program qsspst
      implicit none
c
      include 'qpsglobal.h'
c
c     work space
c
      integer ig,ierr,runtime
      integer time
      character*80 inputfile
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#      QQQ    SSSS   SSSS  PPPP    SSSS  TTTTT       #'
      print *,'#     Q   Q  S      S      P   P  S        T         #'
      print *,'#     Q   Q   SSS    SSS   PPPP    SSS     T         #'
      print *,'#     Q  QQ      S      S  P          S    T         #'
      print *,'#      QQQQ  SSSS   SSSS   P      SSSS     T         #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#              Co- and postseismic deformation       #'
      print *,'#                      based on                      #'
      print *,'#          a spherically symmetric earth model       #'
      print *,'#                                                    #'
      print *,'#                  (Version 2016)                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#             Last modified: Feb 2017                #'
      print *,'######################################################'
      print *,'                                                      '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      open(10,file=inputfile,status='old')
      call qpsgetinp(10)
      close(10)
c
      call qpssublayer(ierr)
c
      do ig=1,ngrn
        if(grnsel(ig).eq.1)then
          lys=lygrn(ig)
          call qpsgrnspec(ig)
        endif
      enddo
      call qpswvint(ierr)
      call qpsfftinv(ierr)
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with qsspst      #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
      stop
      end
