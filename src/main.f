	include 'radau15.f'
      include 'force_sb.f'
      include 'sim_setup.f'
      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DE0(-7:7),DE2(-7:7)
      DIMENSION X(3),V(3),F(3)			
	CHARACTER*50 INFILE001
      CHARACTER*100 OUTFILE001,OUTFILE002
      CHARACTER*30 datadir

!	COMMON BLOCKS

      common/star/zms,rs,alphas,gammas
      common/planet/zmp,rp,alphap,gammap
      common/gravity/g
      common/constants/pi
      common/wind/brak,wsat,rat
      common/roche/a_roche
      common/planet_rotation/znup

!     FUNDAMENTAL CONSTANTS

      PI=DACOS(-1.0D0)
      day=86400.0
      yr=365.25
      brak=2.7d40                        
      earth_mass=5.972d24
      zjup_mass=1.898d27
      solar_mass=1.988d30
      earth_radius=6371.d3
      zjup_radius=69911.d3
      solar_radius=696340.d3
      solar_rot=2.0*pi/(28.0)                   
      au=1.495978707d11
      G=6.67408D-11/(au**3)*(solar_mass*day**2)				
      brak=brak/(solar_mass*au*au*day)
      
!	SPIN, ORBIT AND SYSTEM PARAMETERS

      call sim_setup(A0,E0,O0,zms,rs,aps,zmp,rp,app,f_p,gms,gmp)
      
!     Interpolating saturation frequency based on stellar mass. See scripts file, routine "rotation_rate_interpolation.py"
      
      wsat=zms**2*26.6666667-zms*18.0+5.33333334 
      wsat=wsat*solar_rot
      
      zmp=zmp*zjup_mass/solar_mass 
      rs=rs*solar_radius/au
      rp=1.0*zjup_radius/au                                  
      gammas=gms*day                                       
      gammap=gmp*day
      alphas=aps
      alphap=app
      OMEGA0=2.0*pi/O0
      
      rat=dsqrt((rs/(solar_radius/au))*(1.0/zms))
      a_roche=2.8*(zms/zmp)**(1.0/3.0)*rp
      brak=brak*f_p
      
      NV=3
      LL=12
      XL=0.001d0
      NCLASS=1
      NOR=15        

!	INTEGRATOR PARAMETERS

      dn=5.0d10*yr 	       					! TOTAL INTEGRATION TIME
      tf=1.0d5*yr    						      ! OUTPUT TIMESTEP
      j=dn/tf								! HOW MUCH TIMES I WILL CALL RADAU15

      X(1)=E0*E0
      X(2)=A0
      X(3)=OMEGA0

      destruction_flag=0

      do i=1,j      
      di=dfloat(i)	
      call RA15(x,v,tf,xl,ll,nv,nclass,nor)
      t_tot=tf*(di)
      period_rot=2.0*pi/X(3)
      e2=X(1)
      e=dsqrt(e2)
      a=X(2)
      omegas=X(3)
      dnsat=DSQRT((G*(zms+zmp))/(A**3.0D0))
      period_orb=2.0*pi/dnsat
      p_ratio=period_orb/period_rot
      if(a.ge.(1.01*a_roche)) then
      destruction_flag=0
      pg=a*(1.0-e)
      ap=a*(1.0+e)
      write(1,*)t_tot/365.25d9,a,e,period_rot,period_orb
      else
      destruction_flag=1
      endif
      write(5,*)t_tot/365.25d9,period_rot,destruction_flag
      enddo

      close(1)
      close(2)
      close(5)

      END	