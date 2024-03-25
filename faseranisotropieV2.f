!c
!c
!c#######################################################################
! c.Wrapper for UMAT for standalone testing: take out INCLUDE 'ABA_PARAM.INC
!c#######################################################################
      subroutine myumat(stress,props, dfgrd1,nprops)
      implicit none
      integer ntens
      integer nprops, nstatv
      real*8 ::  STRESS(6)
      real*8 :: STRAN(6), statev(10), DSTRAN(6)
      real*8 ::  DDSDDE(6,6),Props(nprops),STRANMAT(3,3)
      real*8  dfgrd1(3,3)
! c
      integer:: i
      character*80 cmname
      integer ndi,nshr,noel,npt
      integer layer,kspt,kstep,kinc
      real*8 theta1, theta2
      real*8  sse,spd,scd,rpl,ddsddt(1),drplde(1)
      real*8  drpldt,time(2),dtime,dfgrd0(3,3)
      real*8  temp,dtemp,predef(3),dpred(3)
      real*8  coords(3),drot(3,3),pnewdt,celent
! cCF2PY intent(inout) STRESS    
! c      
      do i=1,6
      STRAN(i)  =  0.0d0
      DSTRAN(i) =  0.0d0
      STRESS(i) =  0.0d0
      end do
      do i=1,10
      STATEV(i) = 0.0d0
      end do
! c
! c
      call umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      return
      end
! c
! c
! c
!c#######################################################################
! c.Actual UMAT Wrapper
!c#######################################################################
! c
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
! c
! c
! c        INCLUDE 'ABA_PARAM.INC'
      implicit none
! c.... ABAQUS variables
! c
! c    include 'aba_param.inc'
! c
      integer nprecd
! c      parameter (nprecd=2)
! c      
      character*80 cmname
      integer ndi,nshr,ntens,nstatv,nprops,noel,npt
      integer layer,kspt,kstep,kinc
      real*8  stress(ntens),statev(nstatv),ddsdde(ntens,ntens)
      real*8  sse,spd,scd,rpl,ddsddt(ntens),drplde(ntens)
      real*8  drpldt,stran(ntens),dstran(ntens),time(2),dtime
! c
      real*8  temp,dtemp,predef(1),dpred(1),props(nprops)
      real*8  coords(3),pnewdt,celent,dfgrd0(3,3), drot(3,3)
      real*8  dfgrd1(3,3)
!C
!c------------------------------------------------------------------------
!C.....Variables needed for elsewhere
!C
      real*8 n1ref(3), n2ref(3),n1(3),n2(3)
      real*8 RCG(3,3), RCGtemp(3,3), SPK(3,3), CaSi(3,3)
      real*8 SPKtemp(3,3), dSdE(3,3,3,3)
      real*8 detF, eps, theta1, theta2
      integer i, j, k, l
!c....temp to check C
      real*8 detRCG, RCGinv(3,3),Ccheck(3,3,3,3)
      real*8 fac, SPKcheck(3,3), Eins(3,3)
      real*8 dSdEAb(6,6), CcheckAb(6,6)
!c    Initiate variables??? Don't know if this is necessary????
      ddsddt=0.0d0
      n1ref=0.0d0
      n2ref=0.0d0
      n1=0.0d0
      n2=0.0d0
      RCG=0.0d0
      RCGtemp=0.0d0
      SPK=0.0d0
      CaSi=0.0d0
      SPKtemp=0.0d0
      dSdE=0.0d0
      detF=0.0d0
      eps=0.0d0
      theta1=0.0d0
      theta2=0.0d0
      detRCG=0.0d0
      RCGinv=0.0d0
      Ccheck=0.0d0
      fac=0.0d0
      SPKcheck=0.0d0
      Eins=0.0d0
      dSdEAb=0.0d0
      CcheckAb=0.0d0
!c------------------------------------------------------------------
!c    Find the fiber vectors
!c    (Note we do this before finding stiffness: assume angle changed
!C    is so small that it has no effect in the range of slope)
!c------------------------------------------------------------------
!C
!        write(*,*) '____________step number = ', kstep
!        write(*,*) '_______increment number = ', kinc
!        write(*,*) '_________element number = ', noel
!        write(*,*) '__________int pt number = ', npt
!C
!C    First find the fiber directions in the initial config
!C    use d to make them double precision
      theta1=props(25) * 3.1415927d0 / 180.0d0
      theta2=props(26) * 3.1415927d0 / 180.0d0
      n1ref = [dcos(theta1), dsin(theta1), 0.d0]
      n2ref = [dcos(theta2), dsin(theta2), 0.d0]

!C    Now rotate it into the current config 
 !     call multmavec3_BS(dfgrd1,n1ref, n1)
!      call multmavec3_BS(dfgrd1,n2ref, n2)
!C    Normalize the vector
!      call normVec(n1, 3)
!      call normVec(n2, 3) 
!       call printmat(0,'n1ref,','f10.4',n1ref,3,1)
!       call printmat(0,'n1,','f10.4',n1,3,1)
!c
!c------------------------------------------------------------------
!c     Get the right-cauchy-green strain and determinant of F
!c------------------------------------------------------------------
!C
      call multMaTMa(dfgrd1,dfgrd1,RCG,3,3,3)
!C
      call determinant33(dfgrd1,detF)
!C
!c------------------------------------------------------------------
!c Get the second-piola kirchhoff stress 
!c------------------------------------------------------------------
!C
      call Fasermaterial(props, RCG, n1ref, n2ref, SPK, dSdE)
!C
!       call printmat(0,'CaSi,','e10.2',CaSi,3,3)
!c------------------------------------------------------------------
!c Now do a loop to get the numerical tangent
!c------------------------------------------------------------------
!c
!      eps=1.0d-10
!C........First copy RCG
!      call copy3b3(RCG,RCGtemp)
!c
!      do i=1,3
!        do j=i,3
!C........First increment RCG (use 2eps for normal components)
!          RCGtemp(i,j) = RCGtemp(i,j) + 0.5d0 * eps
!          RCGtemp(j,i) = RCGtemp(j,i) + 0.5d0 * eps
!C
!C........Then get the spk stress for the increment
!          call Fasermaterial(props, RCGtemp, n1ref, n2ref, SPKtemp)
!C
!           write(*,*)'i', i, 'j', j
!           call printmat(0,'SPKtemp,','e10.2',SPKtemp,3,3)
!           call printmat(0,'RCGtemp,','e10.2',RCGtemp,3,3)
!C
!C........Find the difference (slope)
!C
!          do k=1,3
!            do l=k,3
!C  The 2 is for dE/dC, then the difference between normal and shear is  
!C  going from tensorial strain to engineering strain ??  
!              dSdE(k,l,i,j)=2.0d0*(SPKtemp(k,l)-SPK(k,l))/(eps)
!              dSdE(k,l,j,i)=2.0d0*(SPKtemp(k,l)-SPK(k,l))/(eps)
!              dSdE(l,k,i,j)=2.0d0*(SPKtemp(l,k)-SPK(l,k))/(eps)
!              dSdE(l,k,j,i)=2.0d0*(SPKtemp(l,k)-SPK(l,k))/(eps)
!            end do
!          end do
!C......Reset RCG
!        call copy3b3(RCG,RCGtemp)
!        end do
!      end do
!C
!C
!c------------------------------------------------------------------
!c Transfer to cauchy
!c------------------------------------------------------------------
!C
!C....Transfer stress 
      call spk_to_cauchy (detF, dfgrd1, SPK, CaSi)
!C....Transfer cauchy stress to abaqus voigt notation 
      call mat33n6_ABAQ(CaSi,stress)
!C....Transfer Slope (Also transfers it to abaqus form)
      call dSdE_to_dSigdEps(dfgrd1,detF, SPK,dSdE,ddsdde)
!C
!c------------------------------------------------------------------
!c Check Tangent against analytical solution for neo-hookean material
!c------------------------------------------------------------------
!c
!c     fourth order material tangent tensor d SPK / d E
! !c
 !      call determinant33(RCG,detRCG)
! !c
! !c     inverse of right Cauchy-Green strain
! !c
!       call move_gf(RCG,RCGinv,3,3)
!       call Inverse33(RCGinv,detRCG)
!       fac = 2.0d0*props(1)-props(2)*(detRCG-1.0d0)
! !c
!       do i = 1,3
!         do j = 1,3
!           do k = 1,3
!             do l = 1,3
! !c
!               Ccheck(i,j,k,l) = 0.0d0
! !c
!               Ccheck(i,j,k,l) = props(2)*detRCG*RCGinv(i,j)*RCGinv(k,l)
!     &                   + 0.5d0*fac*RCGinv(i,k)*RCGinv(j,l)
!     &                   + 0.5d0*fac*RCGinv(i,l)*RCGinv(j,k)
!             end do
!           end do
!         end do
!       end do
! !c
!       call vier_zweiABAQ(dSdE,dSdEAb)
!       call vier_zweiABAQ(Ccheck,CcheckAb)
!       call printmat(0,'dSdE','f10.4',dSdEAb,6,6)
!       call printmat(0,'Ccheck','f10.4',CcheckAb,6,6)
!C
!c------------------------------------------------------------------
!c Check Stress against analytical solution for neo-hookean material
!c------------------------------------------------------------------
!c
! !c
!       call idendityMatrix(Eins, 3)
! !       write(*,*)'__________Stress_____________________'
!       do i = 1,3
!         do j = 1,3
!           SPKcheck(i,j) = 0.5d0*props(2)*(detRCG-1.0d0)*RCGinv(i,j)
!     &               + props(1)*(Eins(i,j)-RCGinv(i,j))     
!         end do
!       end do
! c
!        call printmat(0,'spkCheck,','f10.8',SPKcheck,3,3)
!        call printmat(0,'SPK,','f10.8',SPK,3,3)
!c------------------------------------------------------------------
!c Check Stress and slope for 0 deg uniaxial fiber material
!c------------------------------------------------------------------
!c
!       write(*,*)'__________Stress_____________________'
!       SPKcheck(1,1)=2.0d0*props(7)*props(8)*(RCG(1,1)-1.0d0)**(props(8)-1.0d0)
!       CcheckAb(1,1)=4.0d0*props(7)*props(8)*(props(8)-1.0d0)
!     &            *(RCG(1,1)-1.0d0)**(props(8)-2.0d0)   
!c
!       write(*,*)'spk1',SPK(1,1), 'check', SPKcheck(1,1)
!       write(*,*)'dSdE',dSdE(1,1,1,1), 'check', CcheckAb(1,1)
!C
!c------------------------------------------------------------------
!c Check outputs by printing them
!c------------------------------------------------------------------
!C   
!        call printmat(0,'dfgrd1,','f10.4',dfgrd1,3,3)
!        call printmat(0,'stress,','f10.4',stress,6,1)
!c       call printmat(0,'RCG,','f10.4',RCG,3,3)
!c       call printmat(0,'RCGinv,','f10.4',RCGinv,3,3)
!        call printmat(0,'ddsdde,','f10.4',ddsdde,6,6)
!        write(*,*)'----------------'
      return
      end
!C
!c#######################################################################
! c....ANISOTROPIC ELASTICITY: Adapted from code from Bertram Stier 
! c....by Scott Stapleton
! C
! c---------------------------------------------------------------------
! Inputs
! 
! d: Matrix: Neo-Hookean 
!      1-mu (lame constant of NeaHookean Matrix)
!      2-lambda (lc of NH Matrix)
!    Matrix: Isotropic
!      3-K_iso1
!      4-K_iso2
!      5-alpha_1
!      6-alpha_2
!    Fiber1/2: Anisotropic
!      7-K_1ani1
!      8-beta_1
!      9-K_2ani1
!      10-beta_2
!      11-K_1ani2
!      12-gamma_1
!      13-K_2ani2
!      14-gamma_2
!      15-K_kop1
!      16-delta_1
!      17-K_kop2
!      18-delta_2
!      19-K_kop12
!      20-zeta
!      21-C1 delay until fiber 1 is loaded (should be 1 by default)
!      22-C2 delay until fiber 2 is loaded (should be 1 by default)
!      23-Vf_1    Volume fraction of fiber 1 
!      24-Vf_2    Volume fraction of fiber 2
!      25-theta_1  Angle of fiber 1 with on the xy plane from x (deg)
!      26-theta_2  Angle of fiber 2 with on the xy plane from x (deg)
!
! N1, N2 unit vectors of the current fiber directions
! F      deformation gradient
! 
! c---------------------------------------------------------------------
! Outputs
!
! CaSi Cauchy Stress in a 6-vector
! Cdach  is the dCaSi / d Legrange-Green Strain

!c#######################################################################
! c
      subroutine Fasermaterial(d, RCG, N1, N2, SPK, dSdE)
      implicit none
!c
      integer i, j, k, l
!C....Inputs
      real*8 N1(3), N2(3), d(*), RCG(3,3)
!C....Outputs
      real*8 SPK(3,3), dSdE(3,3,3,3)
!C....Tensors
      real*8 detRCG, RCGinv(3,3), Eins(3,3)
!c....Parameters
      real*8 mu, lambda
      real*8 K_iso1, K_iso2, alpha_1, alpha_2
      real*8 K_1ani1, beta_1,K_2ani1, beta_2, detF
      real*8 K_1ani2, gamma_1, K_2ani2, gamma_2, NR(3)
      real*8 K_kop1, delta_1, K_kop2, delta_2
      real*8 K_kop12, zeta
      real*8 Mtens1(3,3), Mtens2(3,3) 
      real*8 C1, C2, Vf_1, Vf_2, Vm, Vf_12
!c....Invarients and derivatives
      real*8 I1, I2, I3, I4, I5, I6, I7
      real*8 dI1dC(3,3), dI2dC(3,3), dI3dC(3,3), dI4dC(3,3)
      real*8 dI5dC(3,3), dI6dC(3,3), dI7dC(3,3), dI8dC(3,3)
      real*8 dI5dC1(3,3), dI5dC2(3,3), dI7dC1(3,3), dI7dC2(3,3)
      real*8 dWdI4, dWdI5, dWdI6, dWdI7, dWdI1
!c....Stress Components
      real*8 SPKnh(3,3), SPKiso(3,3), SPKani(3,3)
!c....Temp Variables
      real*8 trRCGRCG, RCGRCG(3,3), dI5dChelp(3,3), dI7dChelp(3,3)
      real*8 I7help(3), I6help(3), I5help(3), I4help(3), skalar
      real*8 hilf(3,3), gen2(3,3), gen5(3,3), gen7(3,3)
      real*8 gen22(3,3,3,3), gen33(3,3,3,3), gen55(3,3,3,3)  
      real*8 gen77(3,3,3,3)
!c
!C --------------------------------------------------------------
!C Set all of the constants to their values
!c
      mu=d(1)
      lambda=d(2)
!c
      K_iso1=d(3)
      K_iso2=d(4)
      alpha_1=d(5)
      alpha_2=d(6)
!c
      K_1ani1=d(7)
      beta_1=d(8)
      K_2ani1=d(9)
      beta_2=d(10)
      K_1ani2=d(11)
      gamma_1=d(12)
      K_2ani2=d(13)
      gamma_2=d(14)
      K_kop1=d(15)
      delta_1=d(16)
      K_kop2=d(17)
      delta_2=d(18)
      K_kop12=d(19)
      zeta=d(20)
      C1=d(21)
      C2=d(22)
      Vf_1=d(23)
      Vf_2=d(24)
!C ...Volume fraction of matrix 
      Vm=1.d0 - Vf_1 - Vf_2
      if((Vf_1+Vf_1) == 0.0d0)then
        Vf_12 = 0.0d0
      else 
        Vf_12 = (Vf_1*Vf_2)/(Vf_1 + Vf_2) 
      endif
!c
!C --------------------------------------------------------------
!C 
!c     unit tensor
      call idendityMatrix(Eins, 3)
!C
!C     Find the structural Tensors
!C
      call dyprod3_BS(n1,n1,Mtens1)
      call dyprod3_BS(n2,n2,Mtens2)
!C
!C    Some RCG stuff
!c
      call determinant33(RCG,detRCG)
!c
!     inverse of right Cauchy-Green strain
!c
      call move_gf(RCG,RCGinv,3,3)
      call determinant33(RCG,detRCG)
      call Inverse33(RCGinv,detRCG)
!c
!      call printmat(0,'Mtens1,','f10.4',Mtens1,3,3)
!      call printmat(0,'Mtens2,','f10.4',Mtens2,3,3)
!c      call printmat(0,'RCG,','f10.4',RCG,3,3)
!c------------------------------------------------------------------
!c Neo-Hookean Part (below 20% strain, matrix)
!c------------------------------------------------------------------
!c
!c     Berechnung der ableitungen Invarianten
!c
      call trace3b3(RCG,I1)
      dI1dC=Eins
!c
      call multmama(RCG,RCG,RCGRCG,3,3,3)
      call trace3b3(RCGRCG,trRCGRCG)
      I2=0.5*(I1*I1-trRCGRCG)
      dI2dC=I1*Eins-RCG
!c
      I3=detRCG
      dI3dC=I3*RCGinv
!c
!c     Second Piola-Kirchhoff  stress tensor for neo-hookean
!c      write(*,*) I1, I2, I3, lambda, mu
!c
      do i = 1,3
        do j = 1,3
          SPKnh(i,j) =   0.5d0*lambda*(I3-1.0d0)*RCGinv(i,j)
     &               + mu*(Eins(i,j)-RCGinv(i,j))     
        end do
      end do  
!c------------------------------------------------------------------
!c Isotropic Part (above 20% strain, matrix)
!c------------------------------------------------------------------
! c
!C Turn off Isotropic and coupled parts when in tension
!      if(I1.le.3.d0)then
!      K_iso1=0.d0
!      K_kop1=0.d0
!      K_kop2=0.d0
!      endif
!      if(I2.le.3.d0)then
!      K_iso2=0.d0
!      endif
!c
!c
!     Second Piola-Kirchhoff  stress tensor for Isotropic Part
!c
      do i = 1,3
        do j = 1,3
          SPKiso(i,j) =2.0d0*K_iso1*alpha_1*(I1-3.d0)**(alpha_1-1.d0)
     &       *Eins(i,j)
     &       +2.0d0*K_iso2*alpha_2*(I2-3.d0)**(alpha_2-1.d0)*dI2dC(i,j)  
        end do
      end do
!C
!c------------------------------------------------------------------
!c Anisotropic Fiber Part 
!c------------------------------------------------------------------
! c           Definieren der Konstanten
!c
!c    Invariants
!c
      call multmavec3_BS(RCG,N1,I4help)
      call multvecvec3_BS(N1,I4help,I4)
!c
      call multmavec3_BS(RCGRCG,N1,I5help)
      call multvecvec3_BS(N1,I5help,I5)
!c
      call multmavec3_BS(RCG,N2,I6help)
      call multvecvec3_BS(N2,I6help,I6)
!c
      call multmavec3_BS(RCGRCG,N2,I7help)
      call multvecvec3_BS(N2,I7help,I7)
!c
!   Hier werden die Ableitungen der Invarianten nach C berechnet
!c
!c   zwei werden für die anisotroie nicht gebraucht!
!c
      dI4dC=Mtens1
!c
      call multmavec3_BS(RCG,N1,dI5dChelp)
      call dyprod3_BS(N1,dI5dChelp,dI5dC1)
      dI5dC2=transpose(dI5dC1)
      dI5dC=dI5dC1+dI5dC2
!c
      dI6dC=Mtens2
!c
      call multmavec3_BS(RCG,N2,dI7dChelp)
      call dyprod3_BS(N2,dI7dChelp,dI7dC1)
      dI7dC2=transpose(dI7dC1)
      dI7dC=dI7dC1+dI7dC2
!c
!c
!   Now turn off the fiber components that are in tension
!c  (If there are values for K_ani1 and K_ani2, you should probably set 
!c  C1 and C2 to 1
!c
!cTook this stuff out because it makes the dSdE assymmetric and makes convergence very difficult!!!!
      if(I4.le.C1)then
!c      write(*,*),"inv4 set to 0"
      K_1ani1=0.d0
      K_kop1=0.d0
      K_kop12=0.d0
      endif
!c
      if(I5.le.1.d0)then
      K_2ani1=0.d0
      endif
!c      !if(inv6.le.1.d0)then !Took this out (SES)
!c      call printmat(0,'inv6,','f10.4',inv6,1,1)  
!c
      if(I6.le.C2)then
!c      write(*,*),"inv6 set to 0"
      K_1ani2=0.d0
      K_kop2=0.d0
      K_kop12=0.0
      endif
!c
      if(I7.le.1.d0)then
      K_2ani2=0.d0
      endif
!c
!c   Now find the derivatives of W
!c   erstmal nur gültig für beta,gamma,zeta,delta =element der natürlichen Zahlen!!
!c
      dWdI1=  delta_1*K_kop1*((I1-3.0d0)**(delta_1-1.0d0))
     &        *(I4-1.0d0)**delta_1
     &        +delta_2*K_kop2*((I1-3.0d0)**(delta_2-1.0d0))
     &        *(I6-1.0d0)**delta_2
!c
      dWdI4=  beta_1*K_1ani1*((I4-1.0d0)**(beta_1-1.0d0))
     &        +K_kop1*delta_1*((I1-3.0d0)**(delta_1))
     &        *((I4-1.0d0)**(delta_1-1.0d0))
     &        +K_kop12*zeta*((I6-1.0d0)**(zeta))
     &        *((I4-1.0d0)**(zeta-1.0d0))
!c
      dWdI5=  beta_2*K_2ani1*((I5-1.0d0)**(beta_2-1.0d0))
!c
      dWdI6=  gamma_1*K_1ani2*((I6-1.0d0)**(gamma_1-1.0d0))
     &        +K_kop2*delta_2*((I1-3.0d0)**(delta_2))
     &        *((I6-1.0d0)**(delta_2-1.0d0))
     &        +K_kop12*zeta*((I4-1.0d0)**zeta)
     &        *((I6-1.0d0)**(zeta-1.0d0))
!c
      dWdI7=  gamma_2*K_2ani2*((I7-1.0d0)**(gamma_2-1.0d0))
!c
!c    Now add all of the components together to get the stress

      SPKani=2.0d0 * (Vf_12 * dWdI1*dI1dC + (Vf_1) 
     &       * dWdI4*dI4dC + (Vf_1) * dWdI5*dI5dC
     &       +(Vf_2) * dWdI6*dI6dC + (Vf_2) * dWdI7*dI7dC)
!c 
!c------------------------------------------------------------------
!c Put all the Stresses Together
!c------------------------------------------------------------------
!c 
      SPK = Vm * SPKnh + Vm * SPKiso + SPKani
!c
!c------------------------------------------------------------------
!c Calculate the material tangent.  Maybe clean this up later
!c------------------------------------------------------------------
! c
! c.... compute gen2=I_C 1-C
      do i=1,3
        do j=1,3
          gen2(i,j)=I1*Eins(i,j)-RCG(i,j)
        end do
      end do
! ! c
! ! c.... compute gen5=C.M_1+M_1.C
      call multmama(RCG,MTens1,hilf,3,3,3)
      call multmama(MTens1,RCG,gen5,3,3,3)
      do i=1,3
        do j=1,3
          gen5(i,j)=gen5(i,j)+hilf(i,j)
        end do
      end do
! c
! c.... compute gen7=C.M_2+M_2.C
      call multmama(RCG,MTens2,hilf,3,3,3)
      call multmama(MTens2,RCG,gen7,3,3,3)
      do i=1,3
        do j=1,3
          gen7(i,j)=gen7(i,j)+hilf(i,j)
        end do
      end do
! c
! c.... compute several dyadic products and matrices which
! c.... are needed for the material tensor
! c
! c.... (1) d gen2/d C  (2) d gen3/d C  (gen3=C inverse)
! c
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
            gen22(i,j,k,l)=Eins(i,j)*Eins(k,l)
     &                  -0.5d0*Eins(i,k)*Eins(j,l)
     &                  -0.5d0*Eins(i,l)*Eins(j,k)
            gen33(i,j,k,l)=-0.5d0*RCGinv(i,k)*RCGinv(j,l)
     &                  -0.5d0*RCGinv(i,l)*RCGinv(j,k)
            gen55(i,j,k,l)=0.5d0*Eins(k,j)*MTens1(l,i)
     &                  +0.5d0*Eins(i,k)*MTens1(j,l)
     &                  +0.5d0*Eins(i,l)*MTens1(j,k)
     &                  +0.5d0*Eins(l,j)*MTens1(k,i)
            gen77(i,j,k,l)=0.5d0*Eins(k,j)*MTens2(l,i)
     &                  +0.5d0*Eins(i,k)*MTens2(j,l)
     &                  +0.5d0*Eins(i,l)*MTens2(j,k)
     &                  +0.5d0*Eins(l,j)*MTens2(k,i)
!C
            dSdE(i,j,k,l)=0.d0
            end do
          end do
        end do
      end do
! c
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
! c....Isotropic part
      dSdE(i,j,k,l)=Vm*(K_iso1*alpha_1*(alpha_1-1.d0)*
     &         (I1-3.d0)**(alpha_1-2.d0)*Eins(i,j)*Eins(k,l)
! c
     &         +K_iso2*alpha_2*(alpha_2-1.d0)*
     &         (I2-3.d0)**(alpha_2-2.d0)*gen2(i,j)*gen2(k,l)
     &         +K_iso2*alpha_2*(I2-3.d0)**(alpha_2-1.d0)
     &         *gen22(i,j,k,l)
! c
! c.... Hyperelastic Part
     &         +lambda/4.d0*I3*RCGinv(i,j)*RCGinv(k,l)
     &         +lambda/4.d0*(I3-1.d0)*gen33(i,j,k,l)
     &         -mu/2.d0 *gen33(i,j,k,l))
!        write(*,*)i,j,k,l,dSdE(i,j,k,l),gen33(i,j,k,l)
! c
! c..... Fiber 1
      dSdE(i,j,k,l)=dSdE(i,j,k,l)
     &         +Vf_1*(K_1ani1*beta_1*(beta_1-1.d0)*
     &         (I4-1.d0)**(beta_1-2.d0)*MTens1(i,j)*MTens1(k,l)
! c
     &         +K_2ani1*beta_2*(beta_2-1.d0)*
     &         (I5-1.d0)**(beta_2-2.d0)*gen5(i,j)*gen5(k,l)
     &         +K_2ani1*beta_2*(I5-1.d0)**(beta_2-1.d0)
     &         *gen55(i,j,k,l))
! c
! c.... Fiber 1 and matrix coupling
      dSdE(i,j,k,l)=dSdE(i,j,k,l)
     &         +Vf_1*(K_kop1*delta_1*(delta_1-1.d0)*(I1-3.d0)
     &         **(delta_1-2.d0)*
     &         (I4-1.d0)**delta_1*Eins(i,j)*Eins(k,l)
     &         +K_kop1*delta_1*(delta_1-1.d0)*(I4-1.d0)**(delta_1-2.d0)*
     &         (I1-3.d0)**delta_1*MTens2(i,j)*MTens2(k,l)
     &         +K_kop1*delta_1*delta_1*(I1-3.d0)**(delta_1-1.d0)*
     &         (I4-1.d0)**(delta_1-1.d0)*
     &         (MTens1(i,j)*Eins(k,l)+Eins(i,j)*MTens1(k,l)))
! c
! c..... Fiber 2
      dSdE(i,j,k,l)=dSdE(i,j,k,l)
     &         +Vf_2*(K_1ani2*gamma_1*(gamma_1-1.d0)*
     &         (I6-1.0d0)**(gamma_1-2.d0)*MTens2(i,j)*MTens2(k,l)
     &         +K_2ani2*gamma_2*(gamma_2-1.d0)*
     &         (I7-1.d0)**(gamma_2-2.d0)*gen7(i,j)*gen7(k,l)
     &         +K_2ani2*gamma_2*(I7-1.d0)**(gamma_2-1.d0)
     &         *gen77(i,j,k,l))
! c
! c.... Fiber 2 and matrix coupling
      dSdE(i,j,k,l)=dSdE(i,j,k,l)
     &         +Vf_2*(K_kop2*delta_2*(delta_2-1.d0)*(I1-3.d0)
     &         **(delta_2-2.d0)*
     &         (I6-1.d0)**delta_2*Eins(i,j)*Eins(k,l)
     &         +K_kop2*delta_2*(delta_2-1.d0)*(I6-1.d0)**(delta_2-2.d0)*
     &         (I1-3.d0)**delta_2*MTens2(i,j)*MTens2(k,l)
     &         +K_kop2*delta_2*delta_2*(I1-3.d0)**(delta_2-1.d0)*
     &         (I4-1.d0)**(delta_2-1.d0)*
     &         (MTens2(i,j)*Eins(k,l)+Eins(i,j)*MTens2(k,l)))
! c
! c.... Fiber 1 and 2 coupling
     &         +Vf_12*(K_kop12*zeta*(zeta-1.d0)*(I4-1.d0)**(zeta-2.d0)*
     &         (I6-1.d0)**zeta*MTens1(i,j)*MTens1(k,l)
     &         +K_kop12*zeta*(zeta-1.d0)*(I6-1.d0)**(zeta-2.d0)*
     &         (I4-1.d0)**zeta*MTens2(i,j)*MTens2(k,l)
     &         +K_kop12*zeta*zeta*(I4-1.d0)**(zeta-1.d0)*
     &         (I6-1.d0)**(zeta-1.d0)*
     &         (MTens2(i,j)*MTens1(k,l)+MTens1(i,j)*MTens2(k,l)))
! c)
      dSdE(i,j,k,l)=4.d0*dSdE(i,j,k,l)
!      write(*,*)i,j,k,l,dSdE(i,j,k,l)
            end do
          end do
        end do
      end do
!
      return
      end
!C
!c#######################################################################
!c#######################################################################
!C   Math Subroutines (could replace some with built-in Fortran functions)
!c#######################################################################
!c#######################################################################
      subroutine spk_to_cauchy (detF, F, spk, CaSi)
      implicit none
      real*8 F(3,3), detF, spk(3,3), CaSi(3,3)
      integer i, j, a, b
!C
      do a=1,3
        do b=1,3
          CaSi(a,b)=0.d0
          do i=1,3
            do j=1,3
              CaSi(a,b)=CaSi(a,b)+F(a,i)*F(b,j)*spk(i,j)/detF
            end do
          end do
        end do
      end do 
!C
      return
      end
!C
!c#######################################################################
      subroutine dSdE_to_dCdE (detF, F, dSdE, dCdE)
      implicit none
      real*8 F(3,3), detF, dSdE(3,3,3,3), dCdE(3,3,3,3)
      integer i, j, a, b, c, k, e, l
!C! c
!c.... push forward of material tensor L
!c
      do a=1,3
        do b=1,3
          do c=1,3
            do e=1,3
              dCdE(a,b,c,e)=0.d0
              do i=1,3
                do j=1,3
                  do k=1,3
                    do l=1,3
                    dCdE(a,b,c,e)=dCdE(a,b,c,e)
     &              +F(a,i)*F(b,j)*F(c,k)*F(e,l)*dSdE(i,j,k,l)/detF
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      return
      end
!C
!c#######################################################################
      subroutine dyprod3_BS (V1,V2,M)
      implicit none
      real*8 V1(3), V2(3), M(3,3)
      M(1,1)=V1(1)*V2(1)
      M(1,2)=V1(1)*V2(2)
      M(1,3)=V1(1)*V2(3)
      M(2,1)=V1(2)*V2(1)
      M(2,2)=V1(2)*V2(2)
      M(2,3)=V1(2)*V2(3)
      M(3,1)=V1(3)*V2(1)
      M(3,2)=V1(3)*V2(2)
      M(3,3)=V1(3)*V2(3)
      return
      end
!c#######################################################################
      subroutine multmavec3_BS (M,V,U)
      implicit none
      real*8 M(3,3), V(3), U(3)
      U(1)=M(1,1)*V(1)+M(1,2)*V(2)+M(1,3)*V(3)
      U(2)=M(2,1)*V(1)+M(2,2)*V(2)+M(2,3)*V(3)
      U(3)=M(3,1)*V(1)+M(3,2)*V(2)+M(3,3)*V(3)
      return
      end
!c#######################################################################
      subroutine multvecvec3_BS (V1,V2,skalar)
      implicit none
      real*8 V1(3), V2(3), skalar
      skalar=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
      return
      end
!#######################################################################
      subroutine mat33n6_ABAQ(b9,a6)
      implicit none
! c
! c.... transfer 3x3 matrix to 6-dim. Voigt notation 
! c
      integer i,j
      real*8  a6(6),b9(3,3)
! c
      do i=1,3
      a6(i)=b9(i,i)
      end do
      a6(4)=b9(1,2)
      a6(5)=b9(1,3)
      a6(6)=b9(2,3)
! c
      return
      end
!#######################################################################
      subroutine vier_zweiABAQ(d,dss66)
      implicit none
! c   Takes 4th order tensorial and changes it to matrix with abaqus 
!C    notation and tensorial strain definition
      integer i
      real*8  dss66(6,6),d(3,3,3,3)
! c
      do i=1,3
      dss66(i,1)=d(i,i,1,1)
      dss66(i,2)=d(i,i,2,2)
      dss66(i,3)=d(i,i,3,3)
      dss66(i,4)=d(i,i,1,2)
      dss66(i,5)=d(i,i,1,3)
      dss66(i,6)=d(i,i,2,3)
      end do
! c
      dss66(4,1)=d(1,2,1,1)
      dss66(4,2)=d(1,2,2,2)
      dss66(4,3)=d(1,2,3,3)
      dss66(4,4)=d(1,2,1,2)
      dss66(4,5)=d(1,2,1,3)
      dss66(4,6)=d(1,2,2,3)
! c
      dss66(5,1)=d(1,3,1,1)
      dss66(5,2)=d(1,3,2,2)
      dss66(5,3)=d(1,3,3,3)
      dss66(5,4)=d(1,3,1,2)
      dss66(5,5)=d(1,3,1,3)
      dss66(5,6)=d(1,3,2,3)
! c
      dss66(6,1)=d(2,3,1,1)
      dss66(6,2)=d(2,3,2,2)
      dss66(6,3)=d(2,3,3,3)
      dss66(6,4)=d(2,3,1,2)
      dss66(6,5)=d(2,3,1,3)
      dss66(6,6)=d(2,3,2,3)
! c
      return
      end
!#######################################################################
!#######################################################################
      subroutine zweiTensorial_zweiEngineering(d_Tens,d_Engineering)
      implicit none
! c   Takes 4th order tensorial and changes it to matrix with abaqus 
!C    notation and engineering strain definition
      integer i,j
      real*8  d_Tens(6,6),d_Engineering(6,6)
! c
      do i=1,6
        do j=1,6
          if (j<4)then
            d_Engineering(i,j) = d_Tens(i,j)
          else
            d_Engineering(i,j) = 2.0d0 * d_Tens(i,j)
          end if
        end do
      end do
! c
      return
      end
!#######################################################################
      subroutine printmat(un,name,fmt,ma,m,n)
      implicit none
! c
      integer, intent(in) :: un
      integer, intent(in) :: m
      integer, intent(in) :: n
      character*(*), intent(in) :: name
      character*(*), intent(in) :: fmt
      double precision, intent(in) :: ma(m*n)
! c
      integer :: lb
      integer :: ub
      integer :: i
      integer :: j
      character(len=2) nc
      character(len=100) ffmt
! c
      write(nc,'(i2)') n
! c
      ffmt = '(' // trim(adjustl(nc)) // '(' // trim(fmt) // ',x))'
! c
      write(un,*)
! c
      write(un,*) name
! c
      do i=1,m
! c
      write(un,trim(ffmt)) ( ma( (j-1)*m + i ), j=1, n )

      end do
! c
      return
      end subroutine printmat
!#######################################################################
! c
      subroutine determinant33(mat,det)
      implicit none
      real*8, intent(in) ::  mat(3,3)
      real*8, intent(inout) :: det
! c.... compute determinant
      det= mat(1,1)*mat(2,2)*mat(3,3)
     &     + mat(1,2)*mat(2,3)*mat(3,1)
     &     + mat(1,3)*mat(2,1)*mat(3,2)
     &     - mat(1,1)*mat(2,3)*mat(3,2)
     &     - mat(2,2)*mat(3,1)*mat(1,3)
     &     - mat(3,3)*mat(1,2)*mat(2,1)
      return
      end
!#######################################################################
      subroutine move_gf(a,b,na,nb)
      implicit none
! c
      integer i,j,na,nb
      real*8  a(na,nb),b(na,nb)
! c.... move array a to array b
      do i=1,na
      do j=1,nb
      b(i,j)=a(i,j)
      end do
      end do
      return
      end
!#######################################################################
      subroutine invert(a,nmax,ndm)
! 
! c      * * F E A P * * A Finite Element Analysis Program
! 
! c....  Copyright (c) 1984-2002: Regents of the University of California
! c                               All rights reserved
! c-----[--.----+----.----+----.-----------------------------------------]
! c      Purpose: Invert small square matrix
! 
! c      Inputs:
! c         a(ndm,*) - Matrix to be inverted
! c         nmax     - Size of upper submatrix to invert
! c         ndm      - Dimension of array
! 
! c      Outputs:
! c         a(ndm,*) - Submatrix replaces original terms, others not
! c                    changed
! c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,n,ndm,nmax
      real*8    d, a(ndm,*)
CF2PY intent(in,out) a(ndm,*) 

      do n = 1,nmax
      if(a(n,n).ne.0.0d0) then
      d = 1.d0/a(n,n)
      do j = 1,nmax
            a(n,j) = -a(n,j)*d
      end do

      do i = 1,nmax
            if(n.ne.i) then
            do j = 1,nmax
            if(n.ne.j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
            end do
            endif
            a(i,n) = a(i,n)*d
      end do
      a(n,n) = d
      else
      write(*,*) ' *WARNING* Zero pivot in INVERT'
      endif
      end do
      end
!#######################################################################
!
      subroutine matzero(b,n,m)
      implicit none
!
      integer i, j, n, m
      real*8 b(n,m)
!
      do j=1,m
        do i=1,n
          b(i,j)=0.d0
        end do
      end do
!
      return
      end
!#######################################################################
! c
      subroutine trace3b3(a,c)
      implicit none
! c.... trace of a 3x3 Matrix 
      real*8 a(3,3),c
! c
      c=a(1,1)+a(2,2)+a(3,3)
! c
      return
      end
!#######################################################################
      subroutine normVec(a, n)
      implicit none
!c
!c.... normalize a vector
!c
      integer i,n
      real*8 a(n), norm
      norm = 0  
      do i=1,n
      norm = norm + a(i) * a(i)
      end do
!c
      norm = sqrt(norm)
      do i=1,n
      a(i) = a(i) / norm
      end do
      return
      end
!c
!#######################################################################
! c
      subroutine copy3b3(a,c)
      implicit none
      integer i,j 
      real*8 a(3,3),c(3,3)
! c
      do i=1,3
        do j=1,3
          c(i,j)=a(i,j)
        end do
      end do
! c
      return
      end
! c
!#######################################################################
!c        
      subroutine idendityMatrix(ident, n)
      implicit none
!c
!c.... Creation of idendity Matrix
!c
      integer n, i, j
      real*8 ident(n,n)
!c      
      do i=1,n
       do j=1,n
        ident(i,j)=0.d0
       end do
       ident(i,i) = 1.d0
      end do
!c
      return
      end
!#######################################################################
! c
      subroutine multMaMa(a,b,c,l,m,n)
      implicit none
! c
! c.... Multiplication of a matrix with a matrix -> matrix
! c      
      integer l,m,n,i,j,k
      real*8 a(l,m),b(m,n),c(l,n)
! c
      do i=1,l
      do j=1,n
      c(i,j)=0.d0
      do k=1,m
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
      end do
      end do
      end do
! ! c
      return
      end
! c
!#######################################################################
! c
      subroutine multMaTMa(a,b,c,l,m,n)
      implicit none
! c
! c.... Multiplication of a matrix^T with a matrix -> matrix
! c      
      integer l,m,n,i,j,k
      real*8 a(m,l),b(m,n),c(l,n)
! c
      do i=1,l
      do j=1,n
      c(i,j)=0.d0
      do k=1,m
            c(i,j)=c(i,j)+a(k,i)*b(k,j)
      end do
      end do
      end do
! c
      return
      end
! c
!#######################################################################
!
      subroutine Inverse33(A,DetA)
!
!     Bestimmung der Inversen einer 3x3-Matrix
!
!     A = Matrix A
!
      implicit none
!
      integer i, j
      real*8 A(3,3)
!
      real*8 DetA
      real*8 InvA(3,3)
!
      DetA =    A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) )
     &        - A(1,2) * ( A(2,1) * A(3,3) - A(2,3) * A(3,1) )
     &        + A(1,3) * ( A(2,1) * A(3,2) - A(2,2) * A(3,1) )
!
      if (DetA .LE. 1.0d-8) write(*,*)'Inverse33: A singulaer'
!
      InvA(1,1) = (  A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / DetA
      InvA(2,1) = ( -A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / DetA
      InvA(3,1) = (  A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / DetA
!
      InvA(1,2) = ( -A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / DetA
      InvA(2,2) = (  A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / DetA
      InvA(3,2) = ( -A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / DetA
!
      InvA(1,3) = (  A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / DetA
      InvA(2,3) = ( -A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / DetA
      InvA(3,3) = (  A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / DetA
!
      do i=1,3
        do j=1,3
          A(i,j) = InvA(i,j)
        end do
      end do
!
      return
      end
!
!#######################################################################
      subroutine dSdE_to_dSigdEps(dfgrd1,detF,spk,dSdE,dSigdEps_Abaqus)
      implicit none
      
      double precision, intent(in) :: dfgrd1(3,3)
      double precision, intent(in) :: spk(3,3)
      double precision, intent(in) :: detF
      double precision, intent(in) :: dSdE(3,3,3,3)
      double precision, intent(out) :: dSigdEps_Abaqus(6,6)
      
      integer :: i,j,k,l,m,n
      double precision :: fp33(3,3)
      double precision :: on33(3,3)
      double precision :: A_4(3,3,3,3)
      double precision :: T_4(3,3,3,3)
      double precision :: Tsum4(3,3,3,3)
      double precision :: dSigdEps_Tens(6,6)
      

      on33      = 0.0d0
      on33(1,1) = 1.0d0
      on33(2,2) = 1.0d0
      on33(3,3) = 1.0d0
      
c
c.... Push-forward of the tangent
c
      do l=1,3
          do k=1,3
            do j=1,3
              do i=1,3
                A_4(i,j,k,l)=0.d0
                T_4(i,j,k,l)=0.d0
                Tsum4(i,j,k,l)=0.d0
              end do
            end do
          end do
        end do
c
      do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                A_4(i,j,k,l)=on33(i,k)*spk(j,l)
                do m=1,3
                  do n=1,3
                    A_4(i,j,k,l)=A_4(i,j,k,l)
     &                +dfgrd1(i,m)*dSdE(m,j,n,l)*dfgrd1(k,n)
                  end do
                end do
              end do
            end do
          end do
        end do
c
      call multmama(dfgrd1,spk,fp33,3,3,3)
c
      do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                do m=1,3
                    T_4(i,j,k,l)=T_4(i,j,k,l)+A_4(i,m,k,l)*dfgrd1(j,m)
     &                +fp33(i,m)*on33(j,k)*on33(m,l)
                end do
              end do
            end do
          end do
        end do
c
      do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                do m=1,3
                 Tsum4(i,j,k,l)=Tsum4(i,j,k,l)+(T_4(i,j,k,m)*dfgrd1(l,m)
     &                +T_4(i,j,l,m)*dfgrd1(k,m))/2.d0
                end do
                Tsum4(i,j,k,l)=Tsum4(i,j,k,l)/detF
              end do
            end do
          end do
        end do
c
c
!      call vier_zweiABAQ(Tsum4,dSigdEps_Tens)
       call vier_zweiABAQ(Tsum4,dSigdEps_Abaqus)
!      call zweiTensorial_zweiEngineering(dSigdEps_Tens,dSigdEps_Abaqus)
!        call printmat(0,'Tsum4,','f10.4',Tsum4,9,9)
!        call printmat(0,'T_4,','f10.4',T_4,9,9)
!        write(*,*)'detF', detF
!        call printmat(0,'A_4,','f10.4',A_4,9,9)
!        call printmat(0,'fp33,','f10.4',fp33,3,3)
!        call printmat(0,'on33,','f10.4',on33,3,3)
      return
      end subroutine dSdE_to_dSigdEps
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc]
