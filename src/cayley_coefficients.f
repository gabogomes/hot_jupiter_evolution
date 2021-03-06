      SUBROUTINE CAYLEY2(E,DE2)	
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DE2(-7:7)
      
      e2=e*e
      e3=e2*e
      e4=e3*e
      e5=e4*e
      e6=e5*e
      e7=e6*e
      
      c1=12144273.0/71680.0
      c2=73369.0/720.0
      c3=228347.0/3840.0
      c4=3071075.0/18432.0
      c5=533.0/16.0
      c6=13827.0/160.0
      c7=845.0/48.0
      c8=32525.0/768.0
      c9=208225.0/6144.0
      c10=17.0/2.0
      c11=115.0/6.0
      c12=601.0/48.0
      c13=7.0/2.0
      c14=123.0/16.0
      c15=489.0/128.0
      c16=1763.0/2048.0
      c17=1.0
      c18=5.0/2.0
      c19=13.0/16.0
      c20=35.0/288.0
      c21=1.0/2.0
      c22=1.0/16.0
      c23=5.0/384.0
      c24=143.0/18432.0     
      c25=0.0
      c26=1.0/48.0
      c27=11.0/768.0
      c28=313.0/30720.0
      c29=1.0/24.0
      c30=7.0/240.0
      c31=81.0/1280.0
      c32=81.0/2048.0
      c33=4.0/45.0
      c34=15625.0/129024.0
      
      DE2(-7)=c1*e7
      DE2(-6)=c2*e6
      DE2(-5)=c3*e5-c4*e7
      DE2(-4)=c5*e4-c6*e6
      DE2(-3)=c7*e3-c8*e5+c9*e7
      DE2(-2)=c10*e2-c11*e4+c12*e6
      DE2(-1)=c13*e-c14*e3+c15*e5-c16*e7
      DE2(0)=c17-c18*e2+c19*e4-c20*e6
      DE2(1)=-c21*e+c22*e3-c23*e5-c24*e7
      DE2(2)=c25
      DE2(3)=c26*e3+c27*e5+c28*e7
      DE2(4)=c29*e4+c30*e6
      DE2(5)=c31*e5+c32*e7
      DE2(6)=c33*e6
      DE2(7)=c34*e7
      
      return
      end

      SUBROUTINE CAYLEY0(E,DE0)	
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DE0(-7:7)
      
      e2=e*e
      e3=e2*e
      e4=e3*e
      e5=e4*e
      e6=e5*e
      e7=e6*e
      
      c1=1.0
      c2=1.0/2.0
      c3=15.0/8.0
      c4=35.0/16.0
      c5=3.0/2.0
      c6=27.0/16.0
      c7=261.0/128.0
      c8=14309.0/6144.0
      c9=9.0/4.0
      c10=7.0/4.0
      c11=141.0/64.0
      c12=53.0/16.0
      c13=393.0/256.0
      c14=24753.0/10240.0
      c15=77.0/16.0
      c16=129.0/160.0
      c17=1773.0/256.0
      c18=4987.0/6144.0
      c19=3167.0/320.0
      c20=432091.0/30720.0
      
      DE0(0)=c1+c2*e2+c3*e4-c4*e6
      DE0(1)=c5*e+c6*e3-c7*e5+c8*e7
      DE0(2)=c9*e2+c10*e4+c11*e6
      DE0(3)=c12*e3+c13*e5+c14*e7
      DE0(4)=c15*e4+c16*e6
      DE0(5)=c17*e5-c18*e7
      DE0(6)=c19*e6
      DE0(7)=c20*e7
            
      DE0(-1)=DE0(1)
      DE0(-2)=DE0(2)
      DE0(-3)=DE0(3)
      DE0(-4)=DE0(4)
      DE0(-5)=DE0(5)
      DE0(-6)=DE0(6)
      DE0(-7)=DE0(7)
      
      return
      end