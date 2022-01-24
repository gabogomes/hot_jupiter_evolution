      include 'cayley_coefficients.f'

	subroutine force(X,V,TM,F)    
	IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(3),V(3),F(3),DE0(-7:7),DE2(-7:7)

      common/star/zms,rs,alphas,gammas
      common/planet/zmp,rp,alphap,gammap
      common/gravity/g
      common/constants/pi
      common/wind/brak,wsat,rat
      common/roche/a_roche
      common/planet_rotation/znup

      e2=X(1)
      e=dsqrt(e2)
      a=X(2)
      omegas=X(3)

      call CAYLEY0(E,DE0)
      call CAYLEY2(E,DE2)

      dnsat=dsqrt(g*(zms+zmp)/a**3)
      znus=2.0*omegas-2.0*dnsat
      znup=1.0+15.0/2.0*e*e+45.0/8.0*e**4+5.0/16.0*e**6            ! Eq. (42) from Hut 1981
      znup=znup/((1.0+3.0*e*e+3.0/8.0*e**4)*(1.0-e*e)**(1.5))      ! Eq. (42) from Hut 1981
      znup=znup*dnsat
      znup=2.0*znup-2.0*dnsat
      
      eprop=15.0/4.0*zms/zmp*(rp/a)**3
      epros=15.0/4.0*zmp/zms*(rs/a)**3
      zin_s=alphas*zms*rs*rs

      f_exc_p=-3.0*dnsat*rp*rp*eprop/(5.0*a**2)*alphap/0.4          ! Correction for non-homogeneous planet (FM2020)
      f_exc_s=-3.0*dnsat*rs*rs*epros/(5.0*a**2)*alphas/0.4          ! Correction for non-homogeneous star (FM2020)

      f_a_p=dnsat*rp**2*eprop/(2.0*a)*alphap
      f_a_s=dnsat*rs**2*epros/(2.0*a)*alphas                    

      f_om_s=-3.0*g*zmp*epros/(2.0*a*a*a)

      zl_orb=zms*zmp*dsqrt(g*a*(1.0-e*e)/(zms+zmp))
      if(omegas.lt.wsat) then
      zl_wind_dot=-brak*omegas*rat*omegas**2
      else
      zl_wind_dot=-brak*omegas*rat*wsat**2
      endif

      dedt=0.0
      dadt=0.0
      domsdt=0.0

      ! PSEUDO SYNCHRONOUS SOLUTION FOR PLANETARY TIDAL DA/DT AND DE/DT

      f_p_tot_exc=-10.5*dnsat*alphap*rp**2*eprop*e/(a**2)
      f_p_tot_exc=f_p_tot_exc*gammap*dnsat/(gammap**2+dnsat**2)
      f_p_tot_exc=2.0*e*f_p_tot_exc

      f_p_tot_a=-21.0*dnsat*alphap*rp**2*eprop*e**2/a
      f_p_tot_a=f_p_tot_a*gammap*dnsat/(gammap**2+dnsat**2)

      do jj=-7,7
      
      djj=dfloat(jj)      
      pk1=2.0*dsqrt(1.0-e*e)-(2.0-djj)*(1.0-e*e)
      pk2=(1.0-e*e)*djj*djj/3.0

      f_s_exc=DE2(jj)**2*gammas*(znus+djj*dnsat)
      f_s_exc=f_s_exc/(gammas**2+(znus+djj*dnsat)**2)*pk1
      f_s_exc2=pk2*gammas*dnsat*DE0(jj)**2/(gammas**2+djj**2*dnsat**2)
      f_s_tot_exc=f_exc_s*(f_s_exc+f_s_exc2)      

      f_s_a=DE2(jj)**2*3.0*(2.0-djj)*gammas*(znus+djj*dnsat)
      f_s_a=f_s_a/(gammas**2+(znus+djj*dnsat)**2)
      f_s_a2=-djj**2*gammas*dnsat*DE0(jj)**2
      f_s_a2=f_s_a2/(gammas**2+djj**2*dnsat**2)
      f_s_tot_a=f_a_s*(f_s_a+f_s_a2)

      f_s_rot_t=DE2(jj)**2*gammas*(znus+djj*dnsat)
      f_s_rot_t=f_s_rot_t/(gammas**2+(znus+djj*dnsat)**2)

      dedt=dedt+f_s_tot_exc!+f_p_tot_exc
      dadt=dadt+f_s_tot_a!+f_p_tot_a
      domsdt=domsdt+f_om_s*f_s_rot_t

      enddo

      dadt=dadt+f_p_tot_a
      dedt=dedt+f_p_tot_exc

      domsdt=domsdt+zl_wind_dot/zin_s
      zl_rot_dot=zin_s*domsdt

      if(a.le.a_roche) then
      dadt=0.0
      a=1.0d8
      domsdt=zl_wind_dot/zin_s
      dedt=0.0
      endif

      f(1)=dedt
      f(2)=dadt
      f(3)=domsdt

      end