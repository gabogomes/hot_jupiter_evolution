      subroutine sim_setup(A0,E0,O0,zms,rs,aps,zmp,rp,app,f_p,gms,gmp)
      implicit real*8 (a-h,o-z)
      CHARACTER*50 INFILE001
      CHARACTER*100 OUTFILE001,OUTFILE002
      CHARACTER*11 datadir

      INFILE001="../cases/sim.init"
      datadir="../results/"

      OPEN(2,FILE=INFILE001,STATUS='unknown')
      
      read(2,*)OUTFILE001
      read(2,*)A0,E0,O0
      read(2,*)zms,rs,aps,gms
      read(2,*)zmp,rp,app,gmp
      read(2,*)f_p

      OUTFILE002=trim(OUTFILE001)//"_stellar_rotation.txt"
      OUTFILE001=trim(OUTFILE001)//"_spin_orbit.txt"
      OUTFILE001=trim(datadir//OUTFILE001)
      OUTFILE002=trim(datadir//OUTFILE002)

      OPEN(1,FILE=OUTFILE001,STATUS='unknown')
      OPEN(5,FILE=OUTFILE002,STATUS='unknown')

      return

      end