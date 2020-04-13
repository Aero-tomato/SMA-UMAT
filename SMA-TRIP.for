!DEC$ FREEFORM
!***********************************************************************************************XXX
!              3D UNIFIED SMA UMAT INCORPORATING TANSFORMATION-INDUCED PLASTICITY(TRIP)         XXX
!                                                                                               XXX
!    	    COPYRIGHT   LEI XU1, Alexandros Solomou1, Theo Baxevanis2, Dimitris Lagoudas1 	XXX
!           1 Department of Aerospace Engineering, Texas A&M University    XXX
!           2 Department of Mechanical Engineering, University of Houston  XXX
!               Version  (3-D)  13/MAY/2017   
!
! XXX******************************************************************************************************XXX
! XXX  *MATERIAL, NAME=SMA                                                                                                              XXX
! XXX  *DEPVAR                                                                                                                                         XXX
! XXX    42                                                                                                                                                  XXX 
! XXX     statev(1)           =      x                                !Martensite Vol. Frac.                                              XXX 
! XXX     statev(2)           =      xd                              !Detwinned Martensite Vol. Frac.                           XXX 
! XXX     statev(3)           =      zt                               !Accumulated Mart. Vol. Frac.                                 XXX
! XXX     statev(4)           =      ztd                             !Accumulated Detwinned Martensite Vol. Frac.     XXX
! XXX     statev(5)           =      dxd                            !Incremental  Detwinned Martensite Vol. Frac.       XXX 
! XXX     statev(6)           =      real(flag_fwd)             !FWD finish at 1 or not 0                                       XXX
! XXX     statev(7)           =      real(flag_rev)             !REV  finish at 1 or not 0                                         XXX
! XXX     statev(8)           =      real(transformation)   !Trans. DIR  0-elastic; -1 REV; 1 FWD                   XXX
! XXX     statev(9:14)      =      et(1:6)                      !Trans. Strain                                                             XXX 
! XXX     statev(15:20)    =      ep(1:6)                      !TRIP Strain                                                               XXX
! XXX     statev(21:26)    =      et_tr(1:6)                  !Trans. Strain at reverse                                           XXX
! XXX     statev(27:32)    =      Lamdat_r(1:6)            !Reverse DIR tensor                                                XXX
! XXX     statev(33:38)    =      Backstress(1:6)          ! Backstress Tensor                                                 XXX
! XXX     statev(39)         =      RPLC                                                                                                             XXX
! XXX     statev(40)         =      BOUND_REACHED                                                                                      XXX
! XXX     statev(41)         =      NR_CONVERGENCE                                                                                    XXX
! XXX     statev(42)         =      TOO                                                                                                              XXX
! XXX************************************************************************************************************XXX
! XXX  *User Material, constants=30                                                                                                                XXX
! XXX  <Ea>,              <Em>,          <v>,           <alpha>,        <CA>,            <CM>,          <Ms>,         <Mf>    XXX 
! XXX  <As>,              <Af>,          <Hmax>,         <kt>,          <SigCal>,       <n1>,           <n2>,         <n3>     XXX
! XXX  <n4>,             <Eb>,          <C1>,          <C2>,           <L1>,       <X_initial>,     <Etd>.  <Tube_Flag>  XXX
! XXX  <Tube_Axis>,  <Tube_Radius>, <Elastic>,    <Coupling>,    <Model>,      <kp>                                  XXX
! XXX*************************************************************************************************************XXX
! XXX*************************************************************************************************************XXX


SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
                            & NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
INCLUDE 'ABA_PARAM.INC'
integer NTENS, NPROPS
CHARACTER*80 CMNAME
DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
                 & PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)				 
REAL*8  PARAM(35),VAR(70)
INTEGER NPARAM, NVAR
REAL*8 variable
NPARAM=35 							!  NUMBER OF INPUT MODEL PARAMETERS
NVAR=70  								!  NUMBER OF USER DEFINED DEPENDENT VARIABLES 

CALL INITIALIZATION(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,DFGRD0,DFGRD1,PARAM,NPARAM,VAR,NVAR)
CALL ELASTIC_PREDICTOR(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,PARAM,NPARAM,VAR,NVAR)
CALL TRANS_CORRECTOR(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,PARAM,NPARAM,VAR,NVAR)
CALL STRESS_UPDATE(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,PARAM,NPARAM,VAR,NVAR,STRESS)
CALL JACOBIAN_MATRIX(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,PARAM,NPARAM,VAR,NVAR,STRESS,DDSDDE)
CALL Accumulated_zt(PARAM,NPARAM,VAR,NVAR)
CALL Backstress_Update(PARAM,NPARAM,VAR,NVAR)
CALL SDV_UPDATE(PARAM,NPARAM,VAR,NVAR,STATEV,NSTATV)
RETURN
END SUBROUTINE

!**********************************************************************
!***********************END of MAIN UMAT********************************
!**********************************************************************

SUBROUTINE INITIALIZATION(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,DFGRD0,DFGRD1,PARAM,NPARAM,VAR,NVAR)
IMPLICIT NONE



!DEFINITION OF ABAQUS PROVIDED VARIABLES
REAL*8 STRAN(6),DSTRAN(6),COORDS(3),DROT(3,3),TIME(2),DFGRD0(3,3),DFGRD1(3,3)
REAL*8 TEMP,DTEMP,DTIME
REAL*8 STATEV(NSTATV),PROPS(NPROPS)
INTEGER JSTEP(4)
INTEGER NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,NSTATV,NPROPS


!DEFINITION OF  INPUT MODEL PARAMETERSS

REAL*8   Sa,Sm,v,alpha,CA,CM,As,Af,Ms,Mf,Hmax,kt,SigCal,n1,n2,n3,n4,Eb,C1,C2,l1,X_Initial,Tube_Radius,kp
REAL*8   rdso, rduo, Dc, Yo, a1, a2, a3
INTEGER EtD,Tube_Flag,Tube_Axis,ELASTIC,COUPLING,MODEL,NLGEOM
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 TOO,RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),BACKSTRESS(6)
INTEGER transformation,flag_fwd,flag_rev,NR_CONVERGENCE,BOUND_REACHED
INTEGER NVAR
REAL*8 VAR(NVAR)



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& MATLAB GENERATED VARIABLES SPECIFIC FOR THIS SUBROUTINE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
real*8 bs1,bs2,bs3,bs4,bs5,bs6
real*8 T0,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
real*8 t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40
real*8 t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CODE START @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


! Finite Strain Model Parameters
REAL*8 FN(3,3),FN1(3,3),RLOG(3,3)

!Some Additional Features Parameters
REAL*8 X_UP_BOUND,X_LOW_BOUND,RATIO,Rcur, A0(1,7)


NLGEOM=JSTEP(3)

e  =0.0D0
eo=0.0D0
et =0.0D0
eto=0.0D0
ep  =0.0D0
epo  =0.0D0
Lamdat_r=0.0D0

x_up_bound=0.99999999_8
x_low_bound=0.00000001_8
BOUND_REACHED=0
NR_CONVERGENCE=0

Sa                  =1._8/props(1)
Sm                 =1._8/props(2)
v                    =props(3)
alpha              =props(4)
CA                  =props(5)
CM                 =props(6)
Ms                  =props(7)
Mf                  =props(8)
As                  =props(9)
Af                   =props(10)
Hmax             =props(11)

kt                   =props(12)
SigCal             =props(13)
n1                  =props(14)
n2                  =props(15)
n3                  =props(16)
n4                  =props(17)
Eb                  =props(18)
C1                  =props(19)
C2                  =props(20)
L1                  =props(21)
x_initial          =props(22)
Etd                 =nint(props(23))
Tube_Flag       =nint(props(24))
Tube_Axis       =nint(props(25))
Tube_Radius   =props(26)
ELASTIC         =nint(props(27))
Coupling         =nint(props(28))
Model             =nint(props(29))
kp                  =props(30)


! Tube Setup Codes
IF(Tube_Flag==1)THEN

     IF(Tube_Axis==1)THEN
     Rcur=sqrt(COORDS(2)**2+COORDS(3)**2)
     ELSEIF(Tube_Axis==2)THEN
     Rcur=sqrt(COORDS(1)**2+COORDS(3)**2)
     ELSEIF(Tube_Axis==3)THEN
     Rcur=sqrt(COORDS(2)**2+COORDS(1)**2)
     END IF	 
     RATIO=Rcur/Tube_Radius
     elseif(Tube_Flag/=1)then
     ratio=1._8
ENDIF



IF (JSTEP(1)==1.AND.KINC<=1) THEN
	x =x_initial                       ! Chec SMA phase
    xo=x_initial
    xd=0._8
    xdo=0._8   
	dxd=0._8
	dxdo=0._8
    zt=0._8
    zto=0._8
    ztd=0._8
    ztdo=0._8
	too=TEMP
	backstress=0.0d0
	transformation=0                            ! Assuming first increment is elastic.	
	IF (x<=x_low_bound) then ! austensite
		x=x_low_bound
		flag_fwd=0                        ! fwd allow
		flag_rev=1                        ! rev stop
	elseif(x>=x_up_bound) then ! martensite
		x=x_up_bound
		flag_fwd=1                       ! fwd stop
		flag_rev=0                         ! rev allow
	elseif(x>=x_low_bound .and. x<=x_up_bound) then ! mixed phase
		flag_fwd=0                         ! fwd allow
		flag_rev=0                         ! rev allow
	endif
		
	IF(x<=X_LOW_BOUND) THEN  ! TO AVOID SINGULARITY IN CASE X_INITIAL<=X_LOW_BOUND
       LAMDAT_R(1:6)=0.0D0
    ELSE 
       LAMDAT_R(1:6)=ET(1:6)/X
    END IF
ELSE
CALL SDV_ASSIGN(x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,transformation,et,eto,ep,epo,et_tr,Lamdat_r,backstress,rplc,BOUND_REACHED,NR_CONVERGENCE,too,STATEV,NSTATV)
ENDIF




 bs1=backstress(1)
 bs2=backstress(2)
 bs3=backstress(3)
 bs4=backstress(4)
 bs5=backstress(5)
 bs6=backstress(6)

! bs1=0
! bs2=0
! bs3=0
! bs4=0
! bs5=0
! bs6=0

! CALCULATE COMMON MATERIAL PROPERTIES such as: a1 a2 a3, rdso, rduo, Dc, Y
!DEC$ NOFREEFORM  
      t2 = SigCal+bs1-bs2
      t3 = SigCal+bs1-bs3
      t4 = bs2-bs3
      t5 = t2**2
      t6 = t5*(1.0D0/2.0D0)
      t7 = t3**2
      t8 = t7*(1.0D0/2.0D0)
      t9 = bs4**2
      t10 = t9*3.0D0
      t11 = bs5**2
      t12 = t11*3.0D0
      t13 = bs6**2
      t14 = t13*3.0D0
      t15 = t4**2
      t16 = t15*(1.0D0/2.0D0)
      t17 = t6+t8+t10+t12+t14+t16
      t18 = sqrt(t17)
      t20 = kt*t18
      t19 = exp(-t20)
      t21 = SigCal*2.0D0
      t22 = bs1*2.0D0
      t23 = 1.0D0/sqrt(t17)
      t24 = bs1-bs2
      t25 = bs1-bs3
      t26 = ca+cm
      t27 = 1.0D0/t26
      t28 = sa-sm
      t29 = SigCal*t28
      t30 = t19-1.0D0
      t31 = Hmax*t30
      t43 = bs2+bs3-t21-t22
      t32 = Hmax*SigCal*kt*t19*t23*t43*(1.0D0/2.0D0)
      t33 = t24**2
      t34 = t33*(1.0D0/2.0D0)
      t35 = t25**2
      t36 = t35*(1.0D0/2.0D0)
      t37 = t10+t12+t14+t16+t34+t36
      t38 = sqrt(t37)
      t39 = Hmax*kt*t19*t23*t38*(bs2+bs3-t21-t22)*(1.0D0/2.0D0)
      t40 = t29+t31+t32+t39
      t41 = mf-ms
      t42 = af-as
      t44 = n1+1.0D0
      t45 = 1.0D0/t44
      t46 = t45+1.0D0
      t47 = n3+1.0D0
      t48 = 1.0D0/t47
      t49 = t48+1.0D0
      t50 = Hmax*SigCal*kt*t19*t23*(bs2+bs3-t21-t22)*(1.0D0/2.0D0)
      t51 = t29+t31+t39+t50
      A0(1,1) = ca*cm*t27*t40*t41*2.0D0
      A0(1,2) = ca*cm*t27*t40*t42*(-2.0D0)
      A0(1,3) = ca*cm*t27*t40*t41*t46*(-1.0D0/2.0D0)-ca*cm*t27*t40*t42*t
     &49*(1.0D0/2.0D0)
      A0(1,4) = ca*cm*t27*t40*2.0D0
      A0(1,5) = -Hmax*t30*t38+ca*cm*t27*t40*(af+ms)
      A0(1,6) = -(t27*(t40*(ca-cm)-t26*(t32-(c1*c2)/(c2*zt+1.0D0))))/(t3
     &1+t32)
      A0(1,7) = -ca*cm*t27*t40*(af-ms)+ca*cm*t27*t41*t46*t51*(1.0D0/2.0D
     &0)+ca*cm*t27*t42*t49*t51*(1.0D0/2.0D0)
!DEC$ FREEFORM
     
a1     =A0(1,1)
a2     =A0(1,2)
a3     =A0(1,3)
rdso  =A0(1,4)
rduo  =A0(1,5)
Dc     =A0(1,6)
Yo     =A0(1,7) 

      ! write(*,*)'***************************************'
      ! write(*,*)'a1',a1
	  ! write(*,*)'a2',a2
	  ! write(*,*)'a3',a3
	  ! write(*,*)'rdso',rdso
	  ! write(*,*)'rduo',rduo
	  ! write(*,*)'Dc',Dc
	  ! write(*,*)'Yo',Yo
	  
!Hmax=-4.15e-5*zt+0.039 ! Only for NiTI Atli Simulation.  

if(NLGEOM==1)THEN ! 1-> ENABLED , 0--> DISABLED
  FN=DFGRD0
  FN1=DFGRD1
  CALL FIND_LOG(RLOG,FN,FN1)   !Find lorarithmic incrementation rotation matrix.
  CALL FIND_E(FN1,e)                  ! Find current  total lorarithmic strain  e for n+1 step.
  CALL FIND_E(FN,eo)                 ! Find previous total lorarithmic strain eo for  n step.
  CALL ROTATION_SDV(eo,et,eto,ep,epo,et_tr,Lamdat_r,RLOG,backstress,NDI,NSHR)  ! RLOG
  !CALL ROTATION_SDV(eo,et,eto,ep,epo,et_tr,Lamdat_r,DROT,NDI,NSHR)  ! DROT   
ELSEIF(NLGEOM==0)THEN
  e  =STRAN+DSTRAN
  eo=STRAN
ENDIF


!@@@@@@@@@@@@@  DEFINE COMMON MATRICES @@@@@@@@@@@@@@@@@@@@@@@@@@@
PARAM(1)  =Sa
PARAM(2)  =Sm
PARAM(3)  =v
PARAM(4)  =alpha
PARAM(5)  =Ca                                 
PARAM(6)  =Cm     
PARAM(7)  =Ms      
PARAM(8)  =Mf              
PARAM(9)  =As              
PARAM(10)=Af                       
PARAM(11)=Hmax                  
PARAM(12)=kt               
PARAM(13)=n1              
PARAM(14)=n2      
PARAM(15)=n3      
PARAM(16)=n4 
PARAM(17)=a1   
PARAM(18)=a2 
PARAM(19)=a3 
PARAM(20)=rdso                                                
PARAM(21)=rduo                                                
PARAM(22)=Dc                                                   
PARAM(23)=Yo                                                   
PARAM(24)=Eb                           
PARAM(25)=C1                           
PARAM(26)=C2 
PARAM(27)=L1 
PARAM(28)=TOO     
PARAM(29)=X_UP_BOUND        
PARAM(30)=X_LOW_BOUND
PARAM(31)=RATIO
PARAM(32)=NLGEOM
PARAM(33)=ELASTIC
PARAM(34)=COUPLING
PARAM(35)=MODEL     
PARAM(36)=kp
   

VAR(1)           =      x 
VAR(2)           =      xo
VAR(3)           =      xd
VAR(4)           =      xdo
VAR(5)           =      zt
VAR(6)           =      zto
VAR(7)           =      ztd            
VAR(8)           =      ztdo               
VAR(9)           =      dxd                            
VAR(10)         =      real(flag_fwd)             
VAR(11)         =      real(flag_rev)             
VAR(12)         =      real(transformation)    
VAR(13:18)    =      et(1:6) 
VAR(19:24)    =      eto(1:6)                       
VAR(25:30)    =      ep(1:6)
VAR(31:36)    =      epo(1:6)                        
VAR(37:42)    =      et_tr(1:6)                  
VAR(43:48)    =      Lamdat_r(1:6)            
VAR(49:54)    =      Backstress(1:6)          
VAR(55:60)    =      e(1:6)
VAR(61:66)    =      eo(1:6)
VAR(67)         =      RPLC
VAR(68)         =      real(BOUND_REACHED)
VAR(69)         =      real(NR_CONVERGENCE)
VAR(70)         =      TOO

END SUBROUTINE

SUBROUTINE ROTATION_SDV(eo,et,eto,ep,epo,et_tr,Lamdat_r,RLOG,backstress,NDI,NSHR)

IMPLICIT NONE
REAL*8 eo(6),et(6),eto(6),ep(6),epo(6),et_tr(6),Lamdat_r(6),backstress(6),RLOG(3,3)

integer nshr,ndi

!CALL ROTATION_PROCEDURE(eo,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(et,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(eto,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(ep,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(epo,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(et_tr,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(Lamdat_r,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(backstress,RLOG,1,NDI,NSHR)

end subroutine

SUBROUTINE ROTATION_PROCEDURE(A0,RLOG,INDEX1,NDI,NSHR)
IMPLICIT NONE
REAL*8 A0(NDI+NSHR),B0(NDI+NSHR),RLOG(3,3)
INTEGER INDEX1
integer nshr,ndi
! index 1 strain like variables
! index 2 stress like variables
B0=0


CALL ROTSIG(A0,RLOG,B0,INDEX1,NDI,NSHR)
A0=B0



END SUBROUTINE

SUBROUTINE SDV_ASSIGN(x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,Lamdat_r,backstress,rplc,bound_reached,NR_Convergence,too,STATEV,NSTATV)
implicit none

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 TOO,RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),BACKSTRESS(6)
INTEGER transformation,flag_fwd,flag_rev,NR_CONVERGENCE,BOUND_REACHED

integer NSTATV
real*8 STATEV(NSTATV)

     x                             =statev(1)                             
     xo                           =statev(1)                             
     xd                           =statev(2)                             
     xdo                         =statev(2)                             
     zt                            =statev(3)                             
     zto                          =statev(3)                             
     ztd                          =statev(4)                             
     ztdo                        =statev(4)                             
	 dxd                         =statev(5)                             
     dxdo                       =statev(5)                             
     flag_fwd                  =nint(statev(6))                     
     flag_rev                   =nint(statev(7))                    
     transformation         =nint(statev(8))                     
     et(1:6)                    =statev(9:14)                        
     eto(1:6)                  =statev(9:14)                        
     ep(1:6)                   =statev(15:20)                      
     epo(1:6)                 =statev(15:20)                       
	 et_tr(1:6)                =statev(21:26)                      
     Lamdat_r(1:6)         =statev(27:32)                      
	 backstress(1:6)       =statev(33:38)                       
                              
     RPLC                         =               statev(39)                                                                                      
     BOUND_REACHED      =              nint(statev(40))   	                                                                              
     NR_CONVERGENCE    =               nint(statev(41))  	                                                                            
     TOO                            =             statev(42)																				

end
                    
SUBROUTINE PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
IMPLICIT NONE
!DEFINITION OF USER DEFINED MODEL PARAMETERS
REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
INTEGER NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CODE START @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Sa                   =PARAM(1)  
Sm                  =PARAM(2)  
v                     =PARAM(3)  
Alpha               =PARAM(4)  
Ca                   =PARAM(5)             
Cm                  =PARAM(6)  
Ms                   =PARAM(7)  
Mf                   =PARAM(8)  
As                   =PARAM(9)  
Af                    =PARAM(10)
Hmax              =PARAM(11)
kt                    =PARAM(12)
n1                   =PARAM(13)
n2                   =PARAM(14)
n3                   =PARAM(15)
n4                   =PARAM(16)
a1                   =PARAM(17)
a2                   =PARAM(18)
a3                   =PARAM(19)
rdso                 =PARAM(20)                            
rduo                =PARAM(21)                           
Dc                   =PARAM(22)                             
Yo                   =PARAM(23)                           
Eb                   =PARAM(24)     
C1                   =PARAM(25)     
C2                   =PARAM(26)
L1                    =PARAM(27)
TOO                  =PARAM(28)
X_UP_BOUND    =PARAM(29)  
X_LOW_BOUND  =PARAM(30)
RATIO               =PARAM(31)
NLGEOM            =nint(PARAM(32))
ELASTIC            =nint(PARAM(33))
COUPLING         =nint(PARAM(34))
MODEL              =nint(PARAM(35)) 
kp=                  PARAM(36)             

END SUBROUTINE

SUBROUTINE VAR_ASSIGNMENTS(INDEX,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too) 
IMPLICIT NONE

!DEFINITIONS OF ARGUMENTS SPECIFIC FOR THIS SUBROUTINE
INTEGER INDEX


!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 too,RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CODE START @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

IF (INDEX==1)THEN

x                                     =           VAR(1)         
xo                                    =           VAR(2)          
xd                                    =           VAR(3)          
xdo                                  =           VAR(4)         
zt                                    =           VAR(5)         
zto                                   =           VAR(6)          
ztd                                   =           VAR(7)          
ztdo                                 =           VAR(8)         
dxd                                  =           VAR(9)          
flag_fwd                           =           nint(VAR(10))       
flag_rev                            =           nint(VAR(11))       
transformation                  =           nint(VAR(12))       
et(1:6)                             =           VAR(13:18)   
eto(1:6)                           =           VAR(19:24)  
ep(1:6)                            =           VAR(25:30)  
epo(1:6)                          =           VAR(31:36)  
et_tr(1:6)                        =           VAR(37:42)  
Lamdat_r(1:6)                  =           VAR(43:48)   
Backstress(1:6)                =           VAR(49:54)  
e(1:6)                              =           VAR(55:60)   
eo(1:6)                            =           VAR(61:66)  
RPLC                                =           VAR(67)       
BOUND_REACHED             =           nint(VAR(68))       
NR_CONVERGENCE           =           nint(VAR(69))       
TOO                                 =           VAR(70)       

ELSEIF (INDEX==2)THEN


VAR(1)           =      x 
VAR(2)           =      xo
VAR(3)           =      xd
VAR(4)           =      xdo
VAR(5)           =      zt
VAR(6)           =      zto
VAR(7)           =      ztd            
VAR(8)           =      ztdo               
VAR(9)           =      dxd                            
VAR(10)         =      real(flag_fwd)             
VAR(11)         =      real(flag_rev)             
VAR(12)         =      real(transformation)    
VAR(13:18)    =      et(1:6) 
VAR(19:24)    =      eto(1:6)                       
VAR(25:30)    =      ep(1:6)
VAR(31:36)    =      epo(1:6)                        
VAR(37:42)    =      et_tr(1:6)                  
VAR(43:48)    =      Lamdat_r(1:6)            
VAR(49:54)    =      Backstress(1:6)          
VAR(55:60)    =      e(1:6)
VAR(61:66)    =      eo(1:6)
VAR(67)         =      RPLC
VAR(68)         =      real(BOUND_REACHED)
VAR(69)         =      real(NR_CONVERGENCE)
VAR(70)         =      TOO

ENDIF

END SUBROUTINE   

SUBROUTINE ELASTIC_PREDICTOR(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,PARAM,NPARAM,VAR,NVAR)
IMPLICIT NONE

!DEFINITION OF ABAQUS PROVIDED VARIABLES
REAL*8 STRAN(6),DSTRAN(6),COORDS(3),DROT(3,3),TIME(2),DFGRD0(3,3),DFGRD1(3,3)
REAL*8 TEMP,DTEMP,DTIME
REAL*8 STATEV(NSTATV),PROPS(NPROPS)
INTEGER JSTEP(4)
INTEGER NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,NSTATV,NPROPS


!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
INTEGER NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


REAL*8 Tn, Tn1,VARIABLE

! variables using in ELASTIC_PREDICTOR_subroutine
REAL*8 PHI_Fwd_Cur,PHI_Fwd_Pre,PHI_Rev_Cur,PHI_Rev_Pre,DelPHI_Fwd,DelPHI_Rev,stress(6)
INTEGER load_DIR, checknan(10), checkinf(10), index1, index2

checknan=0
checkinf=0

CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)
!CALL STRESS_UPDATE(PARAM,NPARAM,VAR,NVAR,temp+dtemp,stress)

if (temp<=0) then
    temp=1.0e-8
endif
Tn=temp                                !previous temp

! if (dtemp==0) then
! dtemp=1.0e-8
! endif
Tn1=temp+dtemp                  !current temp
! eto=et
! epo=ep

!CALL smart_props(propsUR,propsUI,Sa,Sm,Ms,Mf,As,Af,aA,aM,v,rcA,rcM,kt,Hmax,rdso,mos,Yo,d1b,d2b,d3b,d4b,d5b,d1d,d2d,c1,c2,m1,too,x_initial)
CALL PHI_Cal(PARAM,NPARAM,x,xd,zt,ztd,e,et,ep,et_tr,backstress,Lamdat_r,Tn1,too,PHI_Fwd_Cur,1)
CALL check_nan1(PHI_Fwd_Cur,1,checknan(1))
CALL check_inf1(PHI_Fwd_Cur,1,checkinf(1))
      if(checknan(1)==1)then
      write(*,*)'NAN  in PHI_Fwd_Cur'
      endif
      if(checkinf(1)==1)then
      write(*,*)'INF  in PHI_Fwd_Cur'
      endif
 ! if(checknan(1)==1 .or. checkinf(1)) then
 	! Read(*,*),variable
! endif
	

CALL PHI_Cal(PARAM,NPARAM,xo,xdo,zto,ztdo,eo,eto,epo,et_tr,backstress,Lamdat_r,Tn,too,PHI_Fwd_Pre,1)
CALL check_nan1(PHI_Fwd_Pre,1,checknan(2))
CALL check_inf1(PHI_Fwd_Pre,1,checkinf(2))
      if(checknan(2)==1)then
      write(*,*)'NAN  in PHI_Fwd_Pre'
      endif
      if(checkinf(2)==1)then
      write(*,*)'INF  in PHI_Fwd_Pre'
      endif


CALL PHI_Cal(PARAM,NPARAM,x,xd,zt,ztd,e,et,ep,et_tr,backstress,Lamdat_r,Tn1,too,PHI_Rev_Cur,-1)
CALL check_nan1(PHI_Rev_Cur,1,checknan(3))
CALL check_inf1(PHI_Rev_Cur,1,checkinf(3))      
      if(checknan(3)==1)then
      write(*,*)'NAN  in PHI_Rev_Cur'
      endif
      if(checkinf(3)==1)then
      write(*,*)'INF  in PHI_Rev_Cur'
      endif      


CALL PHI_Cal(PARAM,NPARAM,xo,xdo,zto,ztdo,eo,eto,epo,et_tr,backstress,Lamdat_r,Tn,too,PHI_Rev_Pre,-1)
CALL check_nan1(PHI_Rev_Pre,1,checknan(4))
CALL check_inf1(PHI_Rev_Pre,1,checkinf(4))      
      if(checknan(4)==1)then
      write(*,*)'NAN  in PHI_Rev_Pre'
      endif
      if(checkinf(4)==1)then
      write(*,*)'INF  in PHI_Rev_Pre'
      endif

	  
if (maxval(checknan)==1 .or. maxval(checkinf)==1 )then


	! write(*,*)'Current time is***************', propsUR(29)
	! write(*,*)"kinc",kinc
	! write(*,*)"npt",npt
    ! write(*,*)"noel",noel
	write(*,*)"Volume Fraction",x
	write(*,*)"transformation",transformation
	write(*,*)"flag_fwd",flag_fwd
	write(*,*)"flag_rev",flag_rev	
    write(*,*)"ep",ep
    write(*,*)"et",et
	! write(*,*)"stress",stress
    write(*,*)"lamdat_r",lamdat_r
	! write(*,*)'JACOBIAN_MATRIX'
	! write(*,*),DDSDDE
	!read(*,*),variable

endif



DelPHI_Fwd = PHI_Fwd_Cur-PHI_Fwd_Pre
DelPHI_Rev = PHI_Rev_Cur-PHI_Rev_Pre



! Specify loading direction
if (DelPHI_Fwd>0 .and. DelPHI_Rev<=0 ) then           ! Loading
	Load_DIR=1

elseif (DelPHI_Fwd<=0 .and. DelPHI_Rev>0 ) then   ! Unloading
	Load_DIR=-1

elseif (DelPHI_Fwd==0 .and. DelPHI_Rev==0 ) then ! stay the same
	Load_DIR=0

else                                                                        ! wierd uncertain loading condition
	Load_DIR=2
	
endif

! Specify whether transformation happens
IF(Load_DIR==0)then
		transformation=0
		
ELSEIF (Load_DIR==1)then

	if(PHI_Fwd_Cur<=0) then
		transformation=0

	elseif(PHI_Fwd_Cur>0)then
		IF(flag_fwd==0)then
		    transformation=1
		else    ! Forward finished
		    transformation=0
		   
		endif
	endif


ELSEIF(Load_DIR==-1)then
	if(PHI_Rev_Cur<=0) then
		transformation=0

	elseif(PHI_Rev_Cur>0)then
		if(flag_rev==0)then
		   transformation=-1
			
	    else   ! Reverse finished
		   transformation=0
		
		endif
	endif

ELSEIF(Load_DIR==2)then


	if(PHI_Fwd_Cur<0 .and. PHI_Rev_Cur<0)then
       transformation=0
	   

	elseif(PHI_Fwd_Cur>0 .and. PHI_Rev_Cur<0)then
	
	      if(flag_fwd==0 .and. DelPHI_Fwd>0)THEN
             transformation=1    
          else
             transformation=0
          endif
		  
	   
	elseif(PHI_Fwd_Cur<0 .and. PHI_Rev_Cur>0)then
	      
	      if(flag_rev==0 .and. DelPHI_Rev>0)THEN
             transformation=-1    
          else
             transformation=0
          endif
		  
	   
	elseif(PHI_Fwd_Cur>0 .and. PHI_Rev_Cur>0)then
	    
	     if(nint(STATEV(8))==1) then		 
	   		
	          if(flag_fwd==0) then
                 transformation=1    
              else
                 transformation=0
              endif			
			
	    elseif (nint(STATEV(8))==-1) then
		
	          if(flag_rev==0 ) then
                 transformation=-1    
              else
                 transformation=0
              endif			
        
		elseif (nint(STATEV(8))==0) then
				
			 if(DelPHI_Fwd>=DelPHI_Rev)then
	            if(flag_fwd==0) then
                   transformation=1    
                else
                   transformation=0
                endif	                  

             elseif(DelPHI_Fwd<DelPHI_Rev)then
				if(flag_rev==0 ) then
                      transformation=-1    
                else
                      transformation=0
                endif
             endif	  
	    endif 
	
 
	endif
ENDIF

IF(TRANSFORMATION==-1.AND.BOUND_REACHED==1) THEN
BOUND_REACHED=0;
!READ(*,*),VARIABLE
ELSEIF(TRANSFORMATION==1.AND.BOUND_REACHED==-1)THEN
BOUND_REACHED=0;
ENDIF

CALL VAR_ASSIGNMENTS(2,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

! write(*,*)'***************************************time(1)',time(1)
! write(*,*)'transformation',transformation
! write(*,*)"flag_fwd",flag_fwd
! write(*,*)"flag_rev",flag_rev	
! write(*,*)'PHI_Rev_Cur',PHI_Rev_Cur
! write(*,*)'PHI_Rev_Pre',PHI_Rev_Pre
! write(*,*)'PHI_FWD_Cur',PHI_FWD_Cur
! write(*,*)'PHI_FWD_Pre',PHI_FWD_Pre	
! write(*,*)"lamdat_r",lamdat_r


end SUBROUTINE

SUBROUTINE TRANS_CORRECTOR(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,PARAM,NPARAM,VAR,NVAR)
IMPLICIT NONE

!DEFINITION OF ABAQUS PROVIDED VARIABLES
REAL*8 STRAN(6),DSTRAN(6),COORDS(3),DROT(3,3),TIME(2)
REAL*8 TEMP,DTEMP,DTIME
REAL*8 STATEV(NSTATV),PROPS(NPROPS)
INTEGER JSTEP(4)
INTEGER NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,NSTATV,NPROPS


!DEFINITION OF  INPUT MODEL PARAMETERS
REAL*8  SA,SM,V,ALPHA,CA,CM,MS,MF,AS,AF,HMAX,KT,N1,N2,N3,N4,A1,A2,A3,RDSO,RDUO,DC,YO,EB,C1,C2,L1,TOO,X_UP_BOUND,X_LOW_BOUND,RATIO,KP
INTEGER NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 X,XO,XD,XDO,DXD,DXDO,ZT,ZTO,ZTD,ZTDO
REAL*8 RPLC
REAL*8 E(6),EO(6),ET(6),ETO(6),EP(6),EPO(6),LAMDAT_R(6),ET_TR(6),BACKSTRESS(6)
INTEGER FLAG_FWD,FLAG_REV,TRANSFORMATION,BOUND_REACHED,NR_CONVERGENCE
REAL*8 VAR(NVAR)
INTEGER NVAR


! Variables used specific for this subroutine
REAL*8 TN, TN1,VARIABLE,VARIABLE2
REAL*8 NR_RE(13), NR_RE_R(12), NR_JAC(13,13),NR_JAC_INV(13,13), NR_JAC_R(12,12),INV_NR_JAC_R(12,12),DU(13),DU_R(12),NR_JAC1(13,13),RESIDUAL(2)
INTEGER ITER, II, JJ,OUTPUT(4)

TN=TEMP
TN1=TEMP+DTEMP

CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,SA,SM,V,ALPHA,CA,CM,MS,MF,AS,AF,HMAX,KT,N1,N2,N3,N4,A1,A2,A3,RDSO,RDUO,DC,YO,EB,C1,C2,L1,TOO,X_UP_BOUND,X_LOW_BOUND,RATIO,NLGEOM,ELASTIC,COUPLING,MODEL,KP) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,X,XO,XD,XDO,DXD,DXDO,ZT,ZTO,ZTD,ZTDO,FLAG_FWD,FLAG_REV,TRANSFORMATION,ET,ETO,EP,EPO,ET_TR,LAMDAT_R,BACKSTRESS,E,EO,RPLC,BOUND_REACHED,NR_CONVERGENCE,TOO)



NR_CONVERGENCE=0
ITER=0
NR_RE=0
NR_RE_R=0
NR_JAC=0
NR_JAC_R=0
NR_JAC_INV=0
INV_NR_JAC_R=0
DU=0
DU_R=0
RESIDUAL(1)=10
RESIDUAL(2)=10



IF (TRANSFORMATION/=0) THEN

! IF (TIME(2)>=1.53) THEN

DO while (NR_CONVERGENCE==0)  
      
! IF(ITER>=10)THEN
! if (JSTEP(1)==2.AND.time(1)>=0.270) then
  ! READ(*,*),VARIABLE
! Read(*,*),variable2
! ENDIF
 
 
CALL N_R_Residual(PARAM,NPARAM,VAR,NVAR,TN1,TRANSFORMATION,NR_RE)
CALL N_R_JAC(PARAM,NPARAM,VAR,NVAR,TN1,TRANSFORMATION,NR_JAC)
    
	!WRITE(*,*)"************NR_JAC=",NR_JAC
	
      NR_RE_R(1:12)=NR_RE(1:12)                           ! Bound reached NR_RE
      NR_JAC_R(1:12,1:12)=NR_JAC(1:12,1:12)       ! Bound reached NR_JAC
      
      IF (BOUND_REACHED==0)THEN
      RESIDUAL(1)=maxval(abs(NR_RE(1:12))) 
      RESIDUAL(2)=abs(NR_RE(12+1))
	  
      ELSE
      RESIDUAL(1)=maxval(abs(NR_RE_R(1:12)))      ! Bound reached case
      RESIDUAL(2)=0
      ENDIF


!***************Debug Part***************
OUTPUT=0

call check_inf1(NR_RE,13,OUTPUT(1))
call check_inf2(NR_JAC,13,OUTPUT(2))
call check_nan1(NR_RE,13,OUTPUT(3))
call check_nan2(NR_JAC,13,OUTPUT(4))
if(OUTPUT(1)==1.or.OUTPUT(2)==1.or.OUTPUT(3)==1.or.OUTPUT(4)==1)then
                   write(*,*)' Total time at beginning of current step',time(2)  
                   write(*,*)"output",output 
                   write(*,*)"transformation",transformation				   
                   write(*,*)"kinc",kinc
                   write(*,*)"iteration",Iter
                   write(*,*)"Element number",noel
                   write(*,*)"Integration number",npt
                   write(*,*)"Volume Frac.",x
                   write(*,*)"Current temperature temp+dtemp",Tn1			   
                   write(*,*)"total strain E********************************************"
				   write(*,*),e
                   write(*,*)"ET********************************************"
				   write(*,*),ET            
                   write(*,*)"Ep********************************************"
				   write(*,*),Ep
                   write(*,*)"du********************************************"
				   write(*,*),du
                   write(*,*)"ET_tr********************************************"
				   write(*,*),ET_tr
                   write(*,*)"lamdat_r********************************************"
				   write(*,*),LAMDAT_R		
                   write(*,*)"TRANSFORMATION",TRANSFORMATION
                   ! write(*,*)"sa",Sa
                   ! write(*,*)"sm",Sm
                   ! write(*,*)"poisson",v
                   ! write(*,*)"alpha",alpha
                   ! write(*,*)"hmax",hmax
                   ! write(*,*)"kt",kt
                   ! WRITE(*,*)"n1",n1
                   ! WRITE(*,*)"n2",n2
                   ! WRITE(*,*)"n3",n3
                   ! WRITE(*,*)"n4",n4
                   ! write(*,*)"DC",DC
                   ! write(*,*)"rdso",rdso
                   ! write(*,*)"a1",a1
                   ! write(*,*)"a2",a2
                   ! write(*,*)"a3",a3
                   ! write(*,*)"rduo",rduo
                   ! write(*,*)"yo",yo  
                   ! write(*,*)"as",as
                   ! write(*,*)"af",af 
                   ! write(*,*)"ms",ms 
                   ! write(*,*)"mf",mf 
                   ! write(*,*)"ca",ca
                   ! write(*,*)"cm",cm                     
                   do ii=1,13
                       do jj=1,13
                            write(*,*)"NR_JAC",NR_JAC(ii,jj),ii,jj
                       enddo
                   enddo                 
stop
endif
      
       ! IF(KINC/=0)THEN 
       ! CALL DEBUG(5,2,0,NR_RE,NTENS,STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,
        ! &NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PARAM,NPARAM,VAR,NVAR)     
        ! CALL DEBUG(6,1,NR_JAC,0,NTENS,STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,
         ! &NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PARAM,NPARAM,VAR,NVAR)    
      ! ENDIF
      
IF(RESIDUAL(1)<1E-6_8 .and. RESIDUAL(2)<1E-6_8)then
          
      NR_CONVERGENCE=1
      

	  ELSE
          
      IF (BOUND_REACHED==0) THEN
      
           CALL INVERSE(NR_JAC,NR_JAC_INV,13)
           DU=-matmul(NR_JAC_INV,NR_RE)
           X=X+DU(13)    
           
           IF(TRANSFORMATION==1 .and.X>X_UP_BOUND)then
               X=X_UP_BOUND
               BOUND_REACHED=1  
			   FLAG_FWD=1                                 ! Forward transformation finished
			   FLAG_REV= 0 			   
           ELSEif(TRANSFORMATION==1 .and.X<X_LOW_BOUND)then 
               X=X_LOW_BOUND
               BOUND_REACHED=-1
			   	FLAG_FWD=0
	            FLAG_REV=1 
           ELSEif(TRANSFORMATION==-1.and. X<X_LOW_BOUND)then
                X=X_LOW_BOUND
                BOUND_REACHED=-1
				FLAG_FWD=0
	            FLAG_REV=1                                 ! Reverse transformation finished     
           ELSEif(TRANSFORMATION==-1.and. X>X_UP_BOUND)then 
                X=X_UP_BOUND
                BOUND_REACHED=1
				FLAG_FWD=1
	            FLAG_REV=0	
           ELSE
               BOUND_REACHED=0 
               ET(1:6)=ET(1:6)+DU(1:6)
			   EP(1:6)=EP(1:6)+DU(7:12)		   
           ENDIF 
			IF(ITER>=10.and. maxval(abs(DU(1:13))) <1e-6 )THEN				! In case get stuck when x close to 0 and 1.
							NR_CONVERGENCE=1
			ENDIF
     ELSE
         CALL INVERSE(NR_JAC_R,INV_NR_JAC_R,12)
         DU_R=-matmul(INV_NR_JAC_R,NR_RE_R)
         ET(1:6)=ET(1:6)+DU_R(1:6)
		 EP(1:6)=EP(1:6)+DU_R(7:12)	
         ENDIF
         
         ITER=ITER+1 !ITERATIONS COUNTER
      
     IF(ITER>15)THEN
		NR_CONVERGENCE=-1
        WRITE(*,*)' Total time at beginning of current step',time(2)
		WRITE(*,*)"x =",x		
		WRITE(*,*)"************delta x =",DU(13)		
		WRITE(*,*)"************max dET and dEP =",maxval(abs(NR_RE_R(1:12))) 			
		WRITE(*,*)"Hard to Converge in Newton_Raphson,  Iter=",ITER

     ENDIF
      
ENDIF

if(TRANSFORMATION==1)then
LAMDAT_R=ET/X ! THOSE VALUES ARE UPDATED ONLY DURING FORWARD TRANSFORMATION
ET_TR=ET
endif      
      
CALL VAR_ASSIGNMENTS(2,VAR,NVAR,X,XO,XD,XDO,DXD,DXDO,ZT,ZTO,ZTD,ZTDO,FLAG_FWD,FLAG_REV,TRANSFORMATION,ET,ETO,EP,EPO,ET_TR,LAMDAT_R,BACKSTRESS,E,EO,RPLC,BOUND_REACHED,NR_CONVERGENCE,TOO)
      
END DO  

ENDIF

end

Subroutine PHI_Cal(PARAM,NPARAM,x,xd,zt,ztd,e,et,ep,et_tr,backstress,lamdat_r,Tn,too,phi,index1)

implicit real*8 (t)

!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
integer NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

real*8   Tn

! Changing state variables
real*8 e(6),et(6),ep(6),et_tr(6),Lamdat_r(6),backstress(6),stress(6)
real*8 x,xo,xd,xdo,zt,zto,ztd,ztdo

real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6
real*8 lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 phi

CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 


if (too<=0)then
too=1.0d-10
endif


if (x<1.0e-8)then
x=1.0e-8
endif

if (x>=1.0d0)then
x=1.0d0-1.0e-8
endif

   phi=0

   

   e1=e(1)
   e2=e(2)
   e3=e(3)
   e4=e(4)
   e5=e(5)
   e6=e(6)
    if(maxval(abs(e))<=1.0e-8)then
       e1=1.0e-8_8
    end if

   
   et1=et(1)
   et2=et(2)
   et3=et(3)
   et4=et(4)
   et5=et(5)
   et6=et(6)
   if(maxval(abs(et))<=1.0e-9)then
       et1=1.0e-9_8
    end if
   

   et_tr1=et_tr(1)
   et_tr2=et_tr(2)
   et_tr3=et_tr(3)
   et_tr4=et_tr(4)
   et_tr5=et_tr(5)
   et_tr6=et_tr(6)
   if(maxval(abs(et_tr))<=1.0e-10)then
       et_tr1=1.0e-10_8
    end if

   ep1=ep(1)
   ep2=ep(2)
   ep3=ep(3)
   ep4=ep(4)
   ep5=ep(5)
   ep6=ep(6)

   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)
   
   
!Backstress Calculation   
bs1=backstress(1)
bs2=backstress(2)
bs3=backstress(3)
bs4=backstress(4)
bs5=backstress(5)
bs6=backstress(6)
!DEC$ NOFREEFORM   


      if(index1==1) then    !Forward
      t2 = Sm*x
      t11 = Sa*x
      t3 = Sa+t2-t11
      t4 = 1.0D0/t3
      t5 = v**2
      t6 = t5*2.0D0
      t7 = t6+v-1.0D0
      t8 = 1.0D0/t7
      t9 = Tn-too
      t10 = alpha*t9
      t12 = v-1.0D0
      t13 = -e1+ep1+et1+t10
      t14 = -e2+ep2+et2+t10
      t15 = -e3+ep3+et3+t10
      t16 = Sa-Sm
      t17 = t4*t8*t15*v
      t18 = t4*t8*t13*v
      t19 = t4*t8*t14*v
      t20 = t4*t8*t15*v*(1.0D0/2.0D0)
      t30 = t4*t8*t12*t14
      t21 = t17+t18-t30
      t28 = t4*t8*t12*t13
      t22 = t17+t19-t28
      t27 = t4*t8*t12*t15
      t23 = t18+t19-t27
      t24 = t16*t23*v
      t25 = t4*t8*t13*v*(1.0D0/2.0D0)
      t26 = t4*t8*t14*v*(1.0D0/2.0D0)
      t29 = t16*t22*v
      t31 = t16*t21*v
      t33 = v+1.0D0
      t34 = 1.0D0/t33
      t61 = -e4+ep4+et4
      t62 = t4*t34*t61*(1.0D0/2.0D0)
      t32 = bs4-t62
      t65 = -e5+ep5+et5
      t66 = t4*t34*t65*(1.0D0/2.0D0)
      t35 = bs5-t66
      t69 = -e6+ep6+et6
      t70 = t4*t34*t69*(1.0D0/2.0D0)
      t36 = bs6-t70
      t40 = Sa*bs1
      t41 = Sa*bs1*v
      t42 = Sa*bs1*x
      t43 = Sm*bs1*x
      t44 = Sa*bs1*v*x
      t45 = Sm*bs1*v*x
      t47 = Sa*bs2
      t49 = Sa*bs2*v
      t51 = Sa*bs2*x
      t53 = Sm*bs2*x
      t55 = Sa*bs2*v*x
      t57 = Sm*bs2*v*x
      t37 = e1-e2-ep1+ep2-et1+et2+t40+t41-t42+t43-t44+t45-t47-t49+t51-t5
     &3+t55-t57
      t38 = 1.0D0/t33**2
      t39 = 1.0D0/t3**2
      t48 = Sa*bs3
      t50 = Sa*bs3*v
      t52 = Sa*bs3*x
      t54 = Sm*bs3*x
      t56 = Sa*bs3*v*x
      t58 = Sm*bs3*v*x
      t46 = e1-e3-ep1+ep3-et1+et3+t40+t41-t42+t43-t44+t45-t48-t50+t52-t5
     &4+t56-t58
      t59 = e2-e3-ep2+ep3-et2+et3+t47-t48+t49-t50-t51+t52+t53-t54-t55+t5
     &6+t57-t58
      t60 = sqrt(2.0D0)
      t63 = t32**2
      t64 = t63*6.0D0
      t67 = t35**2
      t68 = t67*6.0D0
      t71 = t36**2
      t72 = t71*6.0D0
      t73 = t37**2
      t74 = t38*t39*t73
      t75 = t46**2
      t76 = t38*t39*t75
      t77 = t59**2
      t78 = t38*t39*t77
      t79 = t64+t68+t72+t74+t76+t78
      t80 = abs(t79)
      t81 = 1.0D0/sqrt(t80)
      t82 = sqrt(t80)
      t85 = kt*t60*t82*(1.0D0/2.0D0)
      t83 = exp(-t85)
      t84 = t83-1.0D0
      t86 = v*2.0D0
      t87 = t86+2.0D0
      t88 = ep1*2.0D0
      t89 = et1*2.0D0
      t90 = Sa*bs1*x*2.0D0
      t91 = Sa*bs1*v*x*2.0D0
      t107 = e1*2.0D0
      t108 = Sa*bs1*2.0D0
      t109 = Sa*bs1*v*2.0D0
      t110 = Sm*bs1*x*2.0D0
      t111 = Sm*bs1*v*x*2.0D0
      t92 = e2+e3-ep2-ep3-et2-et3+t47+t48+t49+t50-t51-t52+t53+t54-t55-t5
     &6+t57+t58+t88+t89+t90+t91-t107-t108-t109-t110-t111
      t93 = ep2*2.0D0
      t94 = et2*2.0D0
      t95 = Sa*bs2*x*2.0D0
      t96 = Sa*bs2*v*x*2.0D0
      t112 = e2*2.0D0
      t113 = Sa*bs2*2.0D0
      t114 = Sa*bs2*v*2.0D0
      t115 = Sm*bs2*x*2.0D0
      t116 = Sm*bs2*v*x*2.0D0
      t97 = e1+e3-ep1-ep3-et1-et3+t40+t41-t42+t43-t44+t45+t48+t50-t52+t5
     &4-t56+t58+t93+t94+t95+t96-t112-t113-t114-t115-t116
      t98 = ep3*2.0D0
      t99 = et3*2.0D0
      t100 = Sa*bs3*x*2.0D0
      t101 = Sa*bs3*v*x*2.0D0
      t117 = e3*2.0D0
      t118 = Sa*bs3*2.0D0
      t119 = Sa*bs3*v*2.0D0
      t120 = Sm*bs3*x*2.0D0
      t121 = Sm*bs3*v*x*2.0D0
      t102 = e1+e2-ep1-ep2-et1-et2+t40+t41-t42+t43-t44+t45+t47+t49-t51+t
     &53-t55+t57+t98+t99+t100+t101-t117-t118-t119-t120-t121
      t103 = t84**2
      t104 = c2*ztd
      t105 = t104+1.0D0
      t106 = 1.0D0/t105
      t0 = -Yo-a3-rduo+Tn*rdso+(t20+t25-t4*t8*t12*t14*(1.0D0/2.0D0))*(t2
     &4+t29-t16*t21)+(t20+t26-t4*t8*t12*t13*(1.0D0/2.0D0))*(t24+t31-t16*
     &t22)+(t25+t26-t4*t8*t12*t15*(1.0D0/2.0D0))*(t29+t31-t16*t23)-a1*(-
     &(-x+1.0D0)**n2+x**n1+1.0D0)*(1.0D0/2.0D0)-t16*t38*t39*t61**2*t87*(
     &1.0D0/8.0D0)-t16*t38*t39*t65**2*t87*(1.0D0/8.0D0)-t16*t38*t39*t69*
     &*2*t87*(1.0D0/8.0D0)-Hmax*t60*t63*t81*t84*3.0D0-Hmax*t60*t67*t81*t
     &84*3.0D0-Hmax*t60*t71*t81*t84*3.0D0+Hmax*t4*t34*t60*t81*t84*t92*(b
     &s1+t17+t19-t28)*(1.0D0/2.0D0)+Hmax*t4*t34*t60*t81*t84*t97*(bs2+t17
     &+t18-t30)*(1.0D0/2.0D0)+Hmax*t4*t34*t60*t81*t84*t102*(bs3+t18+t19-
     &t27)*(1.0D0/2.0D0)-Dc*Hmax*t4*t32*t34*t60*t61*t81*t84*(3.0D0/2.0D0
     &)-Dc*Hmax*t4*t34*t35*t60*t65*t81*t84*(3.0D0/2.0D0)-Dc*Hmax*t4*t34*
     &t36*t60*t69*t81*t84*(3.0D0/2.0D0)-Dc*Hmax*t4*t22*t34*t60*t81*t84*t
     &92*(1.0D0/2.0D0)-Dc*Hmax*t4*t21*t34*t60*t81*t84*t97*(1.0D0/2.0D0)-
     &Dc*Hmax*t4*t23*t34*t60*t81*t84*t102*(1.0D0/2.0D0)-c1*c2*t4*t32*t34
     &*t60*t61*t81*t103*t106*(3.0D0/2.0D0)-c1*c2*t4*t34*t35*t60*t65*t81*
     &t103*t106*(3.0D0/2.0D0)-c1*c2*t4*t34*t36*t60*t69*t81*t103*t106*(3.
     &0D0/2.0D0)-c1*c2*t4*t22*t34*t60*t81*t92*t103*t106*(1.0D0/2.0D0)-c1
     &*c2*t4*t21*t34*t60*t81*t97*t103*t106*(1.0D0/2.0D0)-c1*c2*t4*t23*t3
     &4*t60*t81*t102*t103*t106*(1.0D0/2.0D0)



	 
	    phi=t0
 
       elseif (index1==-1)then ! Reverse
      t2 = v+1.0D0
      t3 = 1.0D0/t2
      t4 = Sm*x
      t7 = Sa*x
      t5 = Sa+t4-t7
      t6 = 1.0D0/t5
      t8 = v**2
      t9 = t8*2.0D0
      t10 = t9+v-1.0D0
      t11 = 1.0D0/t10
      t12 = Tn*alpha
      t13 = v-1.0D0
      t15 = alpha*too
      t14 = -e2+ep2+et2+t12-t15
      t16 = -e3+ep3+et3+t12-t15
      t17 = -e1+ep1+et1+t12-t15
      t18 = t6*t11*t17*v
      t19 = t6*t11*t14*v
      t20 = Tn-too
      t21 = alpha*t20
      t22 = -e2+ep2+et2+t21
      t23 = -e1+ep1+et1+t21
      t24 = -e3+ep3+et3+t21
      t25 = t6*t11*t24*v
      t26 = t6*t11*t23*v
      t27 = t6*t11*t22*v
      t28 = Sa-Sm
      t29 = 1.0D0/t5**2
      t30 = e1*v
      t31 = ep2*v
      t32 = ep3*v
      t33 = et2*v
      t34 = et3*v
      t35 = Tn*alpha*v
      t36 = e2*v
      t37 = ep1*v
      t38 = et1*v
      t39 = -e4+ep4+et4
      t40 = -e5+ep5+et5
      t41 = -e6+ep6+et6
      t42 = v*2.0D0
      t43 = t42+2.0D0
      t44 = 1.0D0/t2**2
      t70 = t3*t6*t39*(1.0D0/2.0D0)
      t45 = bs4-t70
      t73 = t3*t6*t40*(1.0D0/2.0D0)
      t46 = bs5-t73
      t76 = t3*t6*t41*(1.0D0/2.0D0)
      t47 = bs6-t76
      t49 = Sa*bs1
      t50 = Sa*bs1*v
      t51 = Sa*bs1*x
      t52 = Sm*bs1*x
      t53 = Sa*bs1*v*x
      t54 = Sm*bs1*v*x
      t56 = Sa*bs2
      t58 = Sa*bs2*v
      t60 = Sa*bs2*x
      t62 = Sm*bs2*x
      t64 = Sa*bs2*v*x
      t66 = Sm*bs2*v*x
      t48 = e1-e2-ep1+ep2-et1+et2+t49+t50-t51+t52-t53+t54-t56-t58+t60-t6
     &2+t64-t66
      t57 = Sa*bs3
      t59 = Sa*bs3*v
      t61 = Sa*bs3*x
      t63 = Sm*bs3*x
      t65 = Sa*bs3*v*x
      t67 = Sm*bs3*v*x
      t55 = e1-e3-ep1+ep3-et1+et3+t49+t50-t51+t52-t53+t54-t57-t59+t61-t6
     &3+t65-t67
      t68 = e2-e3-ep2+ep3-et2+et3+t56-t57+t58-t59-t60+t61+t62-t63-t64+t6
     &5+t66-t67
      t69 = sqrt(2.0D0)
      t71 = t45**2
      t72 = t71*6.0D0
      t74 = t46**2
      t75 = t74*6.0D0
      t77 = t47**2
      t78 = t77*6.0D0
      t79 = t48**2
      t80 = t29*t44*t79
      t81 = t55**2
      t82 = t29*t44*t81
      t83 = t68**2
      t84 = t29*t44*t83
      t85 = t72+t75+t78+t80+t82+t84
      t86 = abs(t85)
      t89 = sqrt(t86)
      t90 = kt*t69*t89*(1.0D0/2.0D0)
      t91 = exp(-t90)
      t87 = t91-1.0D0
      t88 = 1.0D0/sqrt(t86)
      t92 = t87**2
      t93 = c2*ztd
      t94 = t93+1.0D0
      t95 = 1.0D0/t94
      t96 = e3*v
      t97 = -e3+ep3+et3+t12-t15-t30+t31-t32+t33-t34+t35-t36+t37+t38+t96-
     &alpha*too*v
      t0 = -Yo-a3+rduo-Tn*rdso-lamdat_r4*t45-lamdat_r5*t46-lamdat_r6*t47
     &-lamdat_r2*(bs2+t18-t6*t11*t13*t14+t6*t11*t16*v)-lamdat_r3*(bs3+t1
     &8+t19-t6*t11*t13*t16)-lamdat_r1*(bs1+t19-t6*t11*t13*(-e1+ep1+et1+t
     &12-alpha*too)+t6*t11*v*(-e3+ep3+et3+t12-alpha*too))+a2*(-(-x+1.0D0
     &)**n4+x**n3+1.0D0)*(1.0D0/2.0D0)-Dc*lamdat_r2*(t25+t26-t6*t11*t13*
     &t22)-Dc*lamdat_r1*(t25+t27-t6*t11*t13*t23)-Dc*lamdat_r3*(t26+t27-t
     &6*t11*t13*t24)+t28*t29*t39**2*t43*t44*(1.0D0/8.0D0)+t28*t29*t40**2
     &*t43*t44*(1.0D0/8.0D0)+t28*t29*t41**2*t43*t44*(1.0D0/8.0D0)+Dc*lam
     &dat_r4*t3*t6*t39*(1.0D0/2.0D0)+Dc*lamdat_r5*t3*t6*t40*(1.0D0/2.0D0
     &)+Dc*lamdat_r6*t3*t6*t41*(1.0D0/2.0D0)-t11*t14*t28*t29*(-e2+ep2+et
     &2+t12-t15-t30-t31+t32-t33+t34+t35+t36+t37+t38-e3*v-alpha*too*v)*(1
     &.0D0/2.0D0)-t11*t16*t28*t29*t97*(1.0D0/2.0D0)-t11*t17*t28*t29*(-e1
     &+ep1+et1+t12-t15+t30+t31+t32+t33+t34+t35-e2*v-e3*v-ep1*v-et1*v-alp
     &ha*too*v)*(1.0D0/2.0D0)-c1*c2*t3*t6*t39*t45*t69*t88*t92*t95*(3.0D0
     &/2.0D0)-c1*c2*t3*t6*t40*t46*t69*t88*t92*t95*(3.0D0/2.0D0)-c1*c2*t3
     &*t6*t41*t47*t69*t88*t92*t95*(3.0D0/2.0D0)-c1*c2*t3*t11*t29*t69*t88
     &*t92*t95*(-e2+ep2+et2+t12-t15-t30-t31+t32-t33+t34+t35+t36+t37+t38-
     &t96-alpha*too*v)*(e1-e2*2.0D0+e3-ep1+ep2*2.0D0-ep3-et1+et2*2.0D0-e
     &t3+t49+t50-t51+t52-t53+t54+t57+t59-t61+t63-t65+t67-Sa*bs2*2.0D0-Sa
     &*bs2*v*2.0D0+Sa*bs2*x*2.0D0-Sm*bs2*x*2.0D0+Sa*bs2*v*x*2.0D0-Sm*bs2
     &*v*x*2.0D0)*(1.0D0/2.0D0)-c1*c2*t3*t11*t29*t69*t88*t92*t95*(-e1+ep
     &1+et1+t12-t15+t30+t31+t32+t33+t34+t35-t36-t37-t38-t96-alpha*too*v)
     &*(e1*(-2.0D0)+e2+e3+ep1*2.0D0-ep2-ep3+et1*2.0D0-et2-et3+t56+t57+t5
     &8+t59-t60-t61+t62+t63-t64-t65+t66+t67-Sa*bs1*2.0D0-Sa*bs1*v*2.0D0+
     &Sa*bs1*x*2.0D0-Sm*bs1*x*2.0D0+Sa*bs1*v*x*2.0D0-Sm*bs1*v*x*2.0D0)*(
     &1.0D0/2.0D0)-c1*c2*t3*t11*t29*t69*t88*t92*t95*t97*(e1+e2-e3*2.0D0-
     &ep1-ep2+ep3*2.0D0-et1-et2+et3*2.0D0+t49+t50-t51+t52-t53+t54+t56+t5
     &8-t60+t62-t64+t66-Sa*bs3*2.0D0-Sa*bs3*v*2.0D0+Sa*bs3*x*2.0D0-Sm*bs
     &3*x*2.0D0+Sa*bs3*v*x*2.0D0-Sm*bs3*v*x*2.0D0)*(1.0D0/2.0D0)





       phi=t0

!DEC$ FREEFORM


  

	ENDif
	
	
! write(*,*)'Mart. Vol x',x
! write(*,*)'Accum Mart. Vol ztd',ztd
! write(*,*)'total   strain e',e
! write(*,*)'trans. strain et',et
! write(*,*)'plastic. strain et',ep
! write(*,*)'lamdat_r',lamdat_r


end

Subroutine N_R_JAC(PARAM,NPARAM,VAR,NVAR,Tn1,index1,NR_JAC)

implicit real*8 (t)

!DEFINITION OF  INPUT MODEL PARAMETERS
REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
integer NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


! FEA integer information
integer index1
real*8   Tn


! Changing state variables
real*8 NR_JAC(13,13),A0(13,13)


real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6
real*8 eto1,eto2,eto3,eto4,eto5,eto6,epo1,epo2,epo3,epo4,epo5,epo6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

   NR_JAC=0
   A0=0


   if(too<=0)then
   too=1e-10
   endif

   
   e1=e(1)
   e2=e(2)
   e3=e(3)
   e4=e(4)
   e5=e(5)
   e6=e(6)
    if(maxval(abs(e))<=1.0e-8)then
       e1=10.0e-8_8
    end if

	
   et1=et(1)
   et2=et(2)
   et3=et(3)
   et4=et(4)
   et5=et(5)
   et6=et(6)
   if(maxval(abs(et))<=1.0e-8)then
       et1=1.0e-8_8
    end if

   et_tr1=et_tr(1)
   et_tr2=et_tr(2)
   et_tr3=et_tr(3)
   et_tr4=et_tr(4)
   et_tr5=et_tr(5)
   et_tr6=et_tr(6)
   if(maxval(abs(et_tr))<=1.0e-8)then
       et_tr1=1.0e-8_8
    end if
	
   eto1=eto(1)
   eto2=eto(2)
   eto3=eto(3)
   eto4=eto(4)
   eto5=eto(5)
   eto6=eto(6)
    if(maxval(abs(eto))<=1.0e-8)then
       eto1=1.0e-8_8
    end if
	
   ep1=ep(1)
   ep2=ep(2)
   ep3=ep(3)
   ep4=ep(4)
   ep5=ep(5)
   ep6=ep(6)


   epo1=epo(1)
   epo2=epo(2)
   epo3=epo(3)
   epo4=epo(4)
   epo5=epo(5)
   epo6=epo(6)

   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)  
   
   !Backstress Calculation   
   bs1=backstress(1)
   bs2=backstress(2)
   bs3=backstress(3)
   bs4=backstress(4)
   bs5=backstress(5)
   bs6=backstress(6)
   

   if(abs(x-xo)<=1e-15)then
       x=xo+1.0e-15_8
   end if

   

!DEC$ NOFREEFORM
 
       if(index1==1)then    !Forward
       t3 = v+1.0D0
      t4 = 1.0D0/t3
      t5 = Sa*x
      t6 = Sm*x
      t7 = Sa-t5+t6
      t8 = 1.0D0/t7
      t35 = -e4+ep4+et4
      t36 = t4*t8*t35*(1.0D0/2.0D0)
      t2 = bs4-t36
      t39 = -e5+ep5+et5
      t40 = t4*t8*t39*(1.0D0/2.0D0)
      t9 = bs5-t40
      t43 = -e6+ep6+et6
      t44 = t4*t8*t43*(1.0D0/2.0D0)
      t10 = bs6-t44
      t14 = Sa*bs1
      t15 = Sa*bs1*v
      t16 = Sa*bs1*x
      t17 = Sm*bs1*x
      t18 = Sa*bs1*v*x
      t19 = Sm*bs1*v*x
      t21 = Sa*bs2
      t23 = Sa*bs2*v
      t25 = Sa*bs2*x
      t27 = Sm*bs2*x
      t29 = Sa*bs2*v*x
      t31 = Sm*bs2*v*x
      t11 = e1-e2-ep1+ep2-et1+et2+t14+t15-t16+t17-t18+t19-t21-t23+t25-t2
     &7+t29-t31
      t12 = 1.0D0/t3**2
      t13 = 1.0D0/t7**2
      t22 = Sa*bs3
      t24 = Sa*bs3*v
      t26 = Sa*bs3*x
      t28 = Sm*bs3*x
      t30 = Sa*bs3*v*x
      t32 = Sm*bs3*v*x
      t20 = e1-e3-ep1+ep3-et1+et3+t14+t15-t16+t17-t18+t19-t22-t24+t26-t2
     &8+t30-t32
      t33 = e2-e3-ep2+ep3-et2+et3+t21-t22+t23-t24-t25+t26+t27-t28-t29+t3
     &0+t31-t32
      t34 = sqrt(2.0D0)
      t37 = t2**2
      t38 = t37*6.0D0
      t41 = t9**2
      t42 = t41*6.0D0
      t45 = t10**2
      t46 = t45*6.0D0
      t47 = t11**2
      t48 = t12*t13*t47
      t49 = t20**2
      t50 = t12*t13*t49
      t51 = t33**2
      t52 = t12*t13*t51
      t53 = t38+t42+t46+t48+t50+t52
      t54 = abs(t53)
      t55 = sqrt(t54)
      t79 = kt*t34*t55*(1.0D0/2.0D0)
      t56 = exp(-t79)
      t57 = e1*2.0D0
      t58 = ep1*2.0D0
      t59 = et1*2.0D0
      t60 = Sa*bs1*x*2.0D0
      t61 = Sa*bs1*v*x*2.0D0
      t62 = x-xo
      t63 = (t53/abs(t53))
      t64 = e2*2.0D0
      t65 = ep2*2.0D0
      t66 = et2*2.0D0
      t67 = Sa*bs1*2.0D0
      t68 = Sa*bs1*v*2.0D0
      t69 = Sa*bs2*x*2.0D0
      t70 = Sm*bs1*x*2.0D0
      t71 = Sa*bs2*v*x*2.0D0
      t72 = Sm*bs1*v*x*2.0D0
      t83 = Sa*bs2*2.0D0
      t84 = Sa*bs2*v*2.0D0
      t85 = Sm*bs2*x*2.0D0
      t86 = Sm*bs2*v*x*2.0D0
      t73 = t57-t58-t59-t60-t61-t64+t65+t66+t67+t68+t69+t70+t71+t72-t83-
     &t84-t85-t86
      t74 = e3*2.0D0
      t75 = ep3*2.0D0
      t76 = et3*2.0D0
      t77 = Sa*bs3*x*2.0D0
      t78 = Sa*bs3*v*x*2.0D0
      t80 = t56-1.0D0
      t81 = 1.0D0/sqrt(t54)
      t82 = 1.0D0/t54
      t87 = t12*t13*t73
      t88 = e2+e3-ep2-ep3-et2-et3+t21+t22+t23+t24-t25-t26+t27+t28-t29-t3
     &0+t31+t32-t57+t58+t59+t60+t61-t67-t68-t70-t72
      t89 = 1.0D0/t54**(3.0D0/2.0D0)
      t92 = Sa*bs3*2.0D0
      t93 = Sa*bs3*v*2.0D0
      t94 = Sm*bs3*x*2.0D0
      t95 = Sm*bs3*v*x*2.0D0
      t90 = t64-t65-t66-t69-t71-t74+t75+t76+t77+t78+t83+t84+t85+t86-t92-
     &t93-t94-t95
      t96 = t12*t13*t90
      t91 = t87-t96
      t97 = t57-t58-t59-t60-t61+t67+t68+t70+t72-t74+t75+t76+t77+t78-t92-
     &t93-t94-t95
      t98 = t12*t13*t97
      t99 = t96+t98
      t100 = Hmax*t4*t8*t34*t62*t80*t81
      t101 = t87+t98
      t103 = Hmax*t4*t8*t34*t62*t80*t81*(1.0D0/2.0D0)
      t102 = -t103-Hmax*kt*t4*t8*t56*t62*t63*t82*t88*t91*(1.0D0/4.0D0)-H
     &max*t4*t8*t34*t62*t63*t80*t88*t89*t91*(1.0D0/4.0D0)
      t104 = Hmax*t2*t12*t13*t34*t62*t63*t80*t88*t89*(3.0D0/2.0D0)
      t105 = Hmax*kt*t2*t12*t13*t56*t62*t63*t82*t88*(3.0D0/2.0D0)
      t106 = t104+t105
      t107 = Hmax*t9*t12*t13*t34*t62*t63*t80*t88*t89*(3.0D0/2.0D0)
      t108 = Hmax*kt*t9*t12*t13*t56*t62*t63*t82*t88*(3.0D0/2.0D0)
      t109 = t107+t108
      t110 = Hmax*t10*t12*t13*t34*t62*t63*t80*t88*t89*(3.0D0/2.0D0)
      t111 = Hmax*kt*t10*t12*t13*t56*t62*t63*t82*t88*(3.0D0/2.0D0)
      t112 = t110+t111
      t113 = Sa-Sm
      t114 = 1.0D0/t7**3
      t115 = Sm*bs2
      t116 = Sm*bs3
      t117 = Sm*bs2*v
      t118 = Sm*bs3*v
      t121 = Sm*bs1
      t122 = Sm*bs1*v
      t119 = t14+t15-t21-t23+t115+t117-t121-t122
      t120 = t11*t12*t13*t119*2.0D0
      t123 = t21-t22+t23-t24-t115+t116-t117+t118
      t124 = t12*t13*t33*t123*2.0D0
      t125 = t2*t4*t13*t35*t113*6.0D0
      t126 = t4*t9*t13*t39*t113*6.0D0
      t127 = t4*t10*t13*t43*t113*6.0D0
      t128 = e1+e3-ep1-ep3-et1-et3+t14+t15-t16+t17-t18+t19+t22+t24-t26+t
     &28-t30+t32-t64+t65+t66+t69+t71-t83-t84-t85-t86
      t129 = Hmax*kt*t4*t8*t56*t62*t63*t82*t101*t128*(1.0D0/4.0D0)
      t130 = Hmax*t4*t8*t34*t62*t63*t80*t89*t101*t128*(1.0D0/4.0D0)
      t131 = -t103+t129+t130
      t132 = -t103-Hmax*kt*t4*t8*t56*t62*t63*t82*t99*t128*(1.0D0/4.0D0)-
     &Hmax*t4*t8*t34*t62*t63*t80*t89*t99*t128*(1.0D0/4.0D0)
      t133 = Hmax*t2*t12*t13*t34*t62*t63*t80*t89*t128*(3.0D0/2.0D0)
      t134 = Hmax*kt*t2*t12*t13*t56*t62*t63*t82*t128*(3.0D0/2.0D0)
      t135 = t133+t134
      t136 = Hmax*t9*t12*t13*t34*t62*t63*t80*t89*t128*(3.0D0/2.0D0)
      t137 = Hmax*kt*t9*t12*t13*t56*t62*t63*t82*t128*(3.0D0/2.0D0)
      t138 = t136+t137
      t139 = Hmax*t10*t12*t13*t34*t62*t63*t80*t89*t128*(3.0D0/2.0D0)
      t140 = Hmax*kt*t10*t12*t13*t56*t62*t63*t82*t128*(3.0D0/2.0D0)
      t141 = t139+t140
      t142 = t14+t15-t22-t24+t116+t118-t121-t122
      t143 = t12*t13*t20*t142*2.0D0
      t145 = t12*t47*t113*t114*2.0D0
      t146 = t12*t49*t113*t114*2.0D0
      t147 = t12*t51*t113*t114*2.0D0
      t144 = t120+t124+t125+t126+t127+t143-t145-t146-t147
      t148 = e1+e2-ep1-ep2-et1-et2+t14+t15-t16+t17-t18+t19+t21+t23-t25+t
     &27-t29+t31-t74+t75+t76+t77+t78-t92-t93-t94-t95
      t149 = Hmax*kt*t4*t8*t56*t62*t63*t82*t101*t148*(1.0D0/4.0D0)
      t150 = Hmax*t4*t8*t34*t62*t63*t80*t89*t101*t148*(1.0D0/4.0D0)
      t151 = -t103+t149+t150
      t152 = -t103-Hmax*kt*t4*t8*t56*t62*t63*t82*t91*t148*(1.0D0/4.0D0)-
     &Hmax*t4*t8*t34*t62*t63*t80*t89*t91*t148*(1.0D0/4.0D0)
      t153 = Hmax*t2*t12*t13*t34*t62*t63*t80*t89*t148*(3.0D0/2.0D0)
      t154 = Hmax*kt*t2*t12*t13*t56*t62*t63*t82*t148*(3.0D0/2.0D0)
      t155 = t153+t154
      t156 = Hmax*t9*t12*t13*t34*t62*t63*t80*t89*t148*(3.0D0/2.0D0)
      t157 = Hmax*kt*t9*t12*t13*t56*t62*t63*t82*t148*(3.0D0/2.0D0)
      t158 = t156+t157
      t159 = Hmax*t10*t12*t13*t34*t62*t63*t80*t89*t148*(3.0D0/2.0D0)
      t160 = Hmax*kt*t10*t12*t13*t56*t62*t63*t82*t148*(3.0D0/2.0D0)
      t161 = t159+t160
      t162 = Hmax*kt*t2*t56*t62*t63*t82*t101*(-3.0D0/2.0D0)-Hmax*t2*t34*
     &t62*t63*t80*t89*t101*(3.0D0/2.0D0)
      t163 = Hmax*t2*t34*t62*t63*t80*t89*t91*(3.0D0/2.0D0)
      t164 = Hmax*kt*t2*t56*t62*t63*t82*t91*(3.0D0/2.0D0)
      t165 = t163+t164
      t166 = Hmax*t2*t34*t62*t63*t80*t89*t99*(3.0D0/2.0D0)
      t167 = Hmax*kt*t2*t56*t62*t63*t82*t99*(3.0D0/2.0D0)
      t168 = t166+t167
      t169 = Hmax*t4*t8*t34*t62*t80*t81*(3.0D0/2.0D0)
      t172 = Hmax*t2*t4*t8*t9*t34*t62*t63*t80*t89*9.0D0
      t173 = Hmax*kt*t2*t4*t8*t9*t56*t62*t63*t82*9.0D0
      t170 = -t172-t173
      t182 = Hmax*t2*t4*t8*t10*t34*t62*t63*t80*t89*9.0D0
      t183 = Hmax*kt*t2*t4*t8*t10*t56*t62*t63*t82*9.0D0
      t171 = -t182-t183
      t174 = Hmax*kt*t9*t56*t62*t63*t82*t101*(-3.0D0/2.0D0)-Hmax*t9*t34*
     &t62*t63*t80*t89*t101*(3.0D0/2.0D0)
      t175 = Hmax*t9*t34*t62*t63*t80*t89*t91*(3.0D0/2.0D0)
      t176 = Hmax*kt*t9*t56*t62*t63*t82*t91*(3.0D0/2.0D0)
      t177 = t175+t176
      t178 = Hmax*t9*t34*t62*t63*t80*t89*t99*(3.0D0/2.0D0)
      t179 = Hmax*kt*t9*t56*t62*t63*t82*t99*(3.0D0/2.0D0)
      t180 = t178+t179
      t184 = Hmax*t4*t8*t9*t10*t34*t62*t63*t80*t89*9.0D0
      t185 = Hmax*kt*t4*t8*t9*t10*t56*t62*t63*t82*9.0D0
      t181 = -t184-t185
      t186 = Hmax*kt*t10*t56*t62*t63*t82*t101*(-3.0D0/2.0D0)-Hmax*t10*t3
     &4*t62*t63*t80*t89*t101*(3.0D0/2.0D0)
      t187 = Hmax*t10*t34*t62*t63*t80*t89*t91*(3.0D0/2.0D0)
      t188 = Hmax*kt*t10*t56*t62*t63*t82*t91*(3.0D0/2.0D0)
      t189 = t187+t188
      t190 = Hmax*t10*t34*t62*t63*t80*t89*t99*(3.0D0/2.0D0)
      t191 = Hmax*kt*t10*t56*t62*t63*t82*t99*(3.0D0/2.0D0)
      t192 = t190+t191
      t193 = t80**2
      t194 = c2*ztd
      t195 = t194+1.0D0
      t196 = 1.0D0/t195
      t197 = c1*c2*t4*t8*t34*t62*t81*t193*t196*(1.0D0/2.0D0)
      t198 = c1*c2*t4*t8*t34*t62*t63*t88*t89*t91*t193*t196*(1.0D0/4.0D0)
      t199 = c1*c2*kt*t4*t8*t56*t62*t63*t80*t82*t88*t91*t196*(1.0D0/2.0D
     &0)
      t200 = t197+t198+t199
      t201 = c1*c2*t4*t8*t34*t62*t63*t88*t89*t99*t193*t196*(1.0D0/4.0D0)
      t202 = c1*c2*kt*t4*t8*t56*t62*t63*t80*t82*t88*t99*t196*(1.0D0/2.0D
     &0)
      t203 = t197+t201+t202
      t204 = c1*c2*t2*t12*t13*t34*t62*t63*t88*t89*t193*t196*(-3.0D0/2.0D
     &0)-c1*c2*kt*t2*t12*t13*t56*t62*t63*t80*t82*t88*t196*3.0D0
      t205 = c1*c2*t9*t12*t13*t34*t62*t63*t88*t89*t193*t196*(-3.0D0/2.0D
     &0)-c1*c2*kt*t9*t12*t13*t56*t62*t63*t80*t82*t88*t196*3.0D0
      t206 = c1*c2*t10*t12*t13*t34*t62*t63*t88*t89*t193*t196*(-3.0D0/2.0
     &D0)-c1*c2*kt*t10*t12*t13*t56*t62*t63*t80*t82*t88*t196*3.0D0
      t207 = Sm*bs1*2.0D0
      t208 = Sm*bs1*v*2.0D0
      t209 = t197-c1*c2*t4*t8*t34*t62*t63*t89*t101*t128*t193*t196*(1.0D0
     &/4.0D0)-c1*c2*kt*t4*t8*t56*t62*t63*t80*t82*t101*t128*t196*(1.0D0/2
     &.0D0)
      t210 = c1*c2*t4*t8*t34*t62*t63*t89*t91*t128*t193*t196*(1.0D0/4.0D0
     &)
      t211 = c1*c2*kt*t4*t8*t56*t62*t63*t80*t82*t91*t128*t196*(1.0D0/2.0
     &D0)
      t212 = c1*c2*t4*t8*t34*t62*t63*t89*t99*t128*t193*t196*(1.0D0/4.0D0
     &)
      t213 = c1*c2*kt*t4*t8*t56*t62*t63*t80*t82*t99*t128*t196*(1.0D0/2.0
     &D0)
      t214 = t197+t212+t213
      t215 = c1*c2*t2*t12*t13*t34*t62*t63*t89*t128*t193*t196*(-3.0D0/2.0
     &D0)-c1*c2*kt*t2*t12*t13*t56*t62*t63*t80*t82*t128*t196*3.0D0
      t216 = c1*c2*t9*t12*t13*t34*t62*t63*t89*t128*t193*t196*(-3.0D0/2.0
     &D0)-c1*c2*kt*t9*t12*t13*t56*t62*t63*t80*t82*t128*t196*3.0D0
      t217 = c1*c2*t10*t12*t13*t34*t62*t63*t89*t128*t193*t196*(-3.0D0/2.
     &0D0)-c1*c2*kt*t10*t12*t13*t56*t62*t63*t80*t82*t128*t196*3.0D0
      t218 = Sm*bs2*2.0D0
      t219 = Sm*bs2*v*2.0D0
      t220 = t14+t15+t22+t24-t83-t84-t116-t118-t121-t122+t218+t219
      t221 = t197-c1*c2*t4*t8*t34*t62*t63*t89*t101*t148*t193*t196*(1.0D0
     &/4.0D0)-c1*c2*kt*t4*t8*t56*t62*t63*t80*t82*t101*t148*t196*(1.0D0/2
     &.0D0)
      t222 = c1*c2*t4*t8*t34*t62*t63*t89*t91*t148*t193*t196*(1.0D0/4.0D0
     &)
      t223 = c1*c2*kt*t4*t8*t56*t62*t63*t80*t82*t91*t148*t196*(1.0D0/2.0
     &D0)
      t224 = t197+t222+t223
      t225 = c1*c2*t4*t8*t34*t62*t63*t89*t99*t148*t193*t196*(1.0D0/4.0D0
     &)
      t226 = c1*c2*kt*t4*t8*t56*t62*t63*t80*t82*t99*t148*t196*(1.0D0/2.0
     &D0)
      t227 = c1*c2*t2*t12*t13*t34*t62*t63*t89*t148*t193*t196*(-3.0D0/2.0
     &D0)-c1*c2*kt*t2*t12*t13*t56*t62*t63*t80*t82*t148*t196*3.0D0
      t228 = c1*c2*t9*t12*t13*t34*t62*t63*t89*t148*t193*t196*(-3.0D0/2.0
     &D0)-c1*c2*kt*t9*t12*t13*t56*t62*t63*t80*t82*t148*t196*3.0D0
      t229 = c1*c2*t10*t12*t13*t34*t62*t63*t89*t148*t193*t196*(-3.0D0/2.
     &0D0)-c1*c2*kt*t10*t12*t13*t56*t62*t63*t80*t82*t148*t196*3.0D0
      t230 = Sm*bs3*2.0D0
      t231 = Sm*bs3*v*2.0D0
      t232 = t14+t15+t21+t23-t92-t93-t115-t117-t121-t122+t230+t231
      t233 = c1*c2*t2*t34*t62*t63*t89*t101*t193*t196*(3.0D0/2.0D0)
      t234 = c1*c2*kt*t2*t56*t62*t63*t80*t82*t101*t196*3.0D0
      t235 = t233+t234
      t236 = c1*c2*t2*t34*t62*t63*t89*t91*t193*t196*(-3.0D0/2.0D0)-c1*c2
     &*kt*t2*t56*t62*t63*t80*t82*t91*t196*3.0D0
      t237 = c1*c2*t2*t34*t62*t63*t89*t99*t193*t196*(-3.0D0/2.0D0)-c1*c2
     &*kt*t2*t56*t62*t63*t80*t82*t99*t196*3.0D0
      t238 = c1*c2*t4*t8*t34*t37*t62*t63*t89*t193*t196*9.0D0
      t239 = c1*c2*kt*t4*t8*t37*t56*t62*t63*t80*t82*t196*1.8D1
      t240 = c1*c2*t2*t4*t8*t9*t34*t62*t63*t89*t193*t196*9.0D0
      t241 = c1*c2*kt*t2*t4*t8*t9*t56*t62*t63*t80*t82*t196*1.8D1
      t242 = t240+t241
      t243 = c1*c2*t2*t4*t8*t10*t34*t62*t63*t89*t193*t196*9.0D0
      t244 = c1*c2*kt*t2*t4*t8*t10*t56*t62*t63*t80*t82*t196*1.8D1
      t245 = t243+t244
      t246 = c1*c2*t9*t34*t62*t63*t89*t101*t193*t196*(3.0D0/2.0D0)
      t247 = c1*c2*kt*t9*t56*t62*t63*t80*t82*t101*t196*3.0D0
      t248 = t246+t247
      t249 = c1*c2*t9*t34*t62*t63*t89*t91*t193*t196*(-3.0D0/2.0D0)-c1*c2
     &*kt*t9*t56*t62*t63*t80*t82*t91*t196*3.0D0
      t250 = c1*c2*t9*t34*t62*t63*t89*t99*t193*t196*(-3.0D0/2.0D0)-c1*c2
     &*kt*t9*t56*t62*t63*t80*t82*t99*t196*3.0D0
      t251 = c1*c2*t4*t8*t34*t41*t62*t63*t89*t193*t196*9.0D0
      t252 = c1*c2*kt*t4*t8*t41*t56*t62*t63*t80*t82*t196*1.8D1
      t253 = c1*c2*t4*t8*t9*t10*t34*t62*t63*t89*t193*t196*9.0D0
      t254 = c1*c2*kt*t4*t8*t9*t10*t56*t62*t63*t80*t82*t196*1.8D1
      t255 = t253+t254
      t256 = c1*c2*t10*t34*t62*t63*t89*t101*t193*t196*(3.0D0/2.0D0)
      t257 = c1*c2*kt*t10*t56*t62*t63*t80*t82*t101*t196*3.0D0
      t258 = t256+t257
      t259 = c1*c2*t10*t34*t62*t63*t89*t91*t193*t196*(-3.0D0/2.0D0)-c1*c
     &2*kt*t10*t56*t62*t63*t80*t82*t91*t196*3.0D0
      t260 = c1*c2*t10*t34*t62*t63*t89*t99*t193*t196*(-3.0D0/2.0D0)-c1*c
     &2*kt*t10*t56*t62*t63*t80*t82*t99*t196*3.0D0
      t261 = c1*c2*t4*t8*t34*t45*t62*t63*t89*t193*t196*9.0D0
      t262 = c1*c2*kt*t4*t8*t45*t56*t62*t63*t80*t82*t196*1.8D1
      t263 = v**2
      t264 = t263*2.0D0
      t265 = t264+v-1.0D0
      t266 = 1.0D0/t265
      t267 = Tn-too
      t268 = alpha*t267
      t269 = v-1.0D0
      t270 = -e3+ep3+et3+t268
      t271 = -e1+ep1+et1+t268
      t272 = t8*t266*t271*v*(1.0D0/2.0D0)
      t273 = -e2+ep2+et2+t268
      t274 = t8*t113*t266*v
      t275 = t8*t113*t266*t269*v
      t295 = t8*t113*t263*t266
      t276 = t274+t275-t295
      t277 = t8*t266*t273*v*(1.0D0/2.0D0)
      t278 = t8*t266*t270*v*(1.0D0/2.0D0)
      t279 = t8*t266*t270*v
      t280 = t8*t266*t271*v
      t281 = t8*t266*t273*v
      t289 = t8*t266*t269*t270
      t282 = t280+t281-t289
      t286 = t8*t266*t269*t271
      t283 = t279+t281-t286
      t284 = t113*t283*v
      t287 = t8*t266*t269*t273
      t285 = t279+t280-t287
      t288 = t113*t285*v
      t290 = t113*t282*v
      t291 = bs1+t279+t281-t286
      t292 = bs2+t279+t280-t287
      t293 = bs3+t280+t281-t289
      t308 = t8*t266*t269*t271*(1.0D0/2.0D0)
      t294 = t277+t278-t308
      t310 = t8*t266*t269*t270*(1.0D0/2.0D0)
      t296 = t272+t277-t310
      t297 = t8*t113*t266*t269
      t298 = t8*t113*t263*t266*2.0D0
      t299 = t297+t298
      t309 = t8*t266*t269*t273*(1.0D0/2.0D0)
      t300 = t272+t278-t309
      t311 = t113*t283
      t301 = t288+t290-t311
      t315 = t113*t282
      t302 = t284+t288-t315
      t303 = t8*t266*t302*v*(1.0D0/2.0D0)
      t313 = t113*t285
      t304 = t284+t290-t313
      t305 = Dc*Hmax*t4*t8*t34*t80*t81*t282*(1.0D0/2.0D0)
      t306 = Hmax*t4*t13*t34*t80*t81*t148*t266*v*(1.0D0/2.0D0)
      t307 = c1*c2*t4*t8*t34*t81*t193*t196*t282*(1.0D0/2.0D0)
      t312 = t8*t266*t301*v*(1.0D0/2.0D0)
      t314 = t8*t266*t304*v*(1.0D0/2.0D0)
      t316 = Dc*Hmax*t4*t8*t34*t80*t81*t283*(1.0D0/2.0D0)
      t317 = Dc*Hmax*t4*t8*t34*t80*t81*t285*(1.0D0/2.0D0)
      t318 = Hmax*t4*t13*t34*t80*t81*t128*t266*v*(1.0D0/2.0D0)
      t319 = Hmax*t4*t13*t34*t80*t81*t88*t266*v*(1.0D0/2.0D0)
      t320 = c1*c2*t4*t8*t34*t81*t193*t196*t283*(1.0D0/2.0D0)
      t321 = c1*c2*t4*t8*t34*t81*t193*t196*t285*(1.0D0/2.0D0)
      t322 = v*2.0D0
      t323 = t322+2.0D0
      t324 = t294*t299
      t325 = Hmax*t4*t8*t34*t80*t81*t291
      t326 = Dc*Hmax*t4*t13*t34*t80*t81*t88*t266*t269*(1.0D0/2.0D0)
      t327 = Hmax*kt*t4*t8*t56*t63*t82*t88*t101*t291*(1.0D0/4.0D0)
      t328 = Hmax*kt*t4*t8*t56*t63*t82*t101*t128*t292*(1.0D0/4.0D0)
      t329 = Hmax*kt*t4*t8*t56*t63*t82*t101*t148*t293*(1.0D0/4.0D0)
      t330 = Hmax*t4*t8*t34*t63*t80*t88*t89*t101*t291*(1.0D0/4.0D0)
      t331 = Hmax*t4*t8*t34*t63*t80*t89*t101*t128*t292*(1.0D0/4.0D0)
      t332 = Hmax*t4*t8*t34*t63*t80*t89*t101*t148*t293*(1.0D0/4.0D0)
      t333 = c1*c2*t4*t13*t34*t81*t88*t193*t196*t266*t269*(1.0D0/2.0D0)
      t335 = t276*t296
      t338 = Hmax*t4*t8*t34*t80*t81*t293*(1.0D0/2.0D0)
      t346 = Dc*Hmax*t4*t13*t34*t80*t81*t148*t266*v*(1.0D0/2.0D0)
      t350 = c1*c2*t4*t13*t34*t81*t148*t193*t196*t266*v*(1.0D0/2.0D0)
      t373 = t276*t300
      t375 = Hmax*t4*t8*t34*t80*t81*t292*(1.0D0/2.0D0)
      t384 = Dc*Hmax*t4*t13*t34*t80*t81*t128*t266*v*(1.0D0/2.0D0)
      t388 = c1*c2*t4*t13*t34*t81*t128*t193*t196*t266*v*(1.0D0/2.0D0)
      t334 = t303+t305+t306+t307+t314+t317+t318+t321+t324+t325+t326+t327
     &+t328+t329+t330+t331+t332+t333-t335-t338-t346-t350-t373-t375-t384-
     &t388-t8*t266*t269*t301*(1.0D0/2.0D0)-Hmax*kt*t37*t56*t63*t82*t101*
     &(3.0D0/2.0D0)-Hmax*kt*t41*t56*t63*t82*t101*(3.0D0/2.0D0)-Hmax*kt*t
     &45*t56*t63*t82*t101*(3.0D0/2.0D0)-Hmax*t34*t37*t63*t80*t89*t101*(3
     &.0D0/2.0D0)-Hmax*t34*t41*t63*t80*t89*t101*(3.0D0/2.0D0)-Hmax*t34*t
     &45*t63*t80*t89*t101*(3.0D0/2.0D0)-Dc*Hmax*t4*t8*t34*t80*t81*t283-H
     &max*t4*t13*t34*t80*t81*t88*t266*t269*(1.0D0/2.0D0)-c1*c2*t4*t8*t34
     &*t81*t193*t196*t283-Dc*Hmax*t2*t4*t8*t34*t35*t63*t80*t89*t101*(3.0
     &D0/4.0D0)-Dc*Hmax*t4*t8*t9*t34*t39*t63*t80*t89*t101*(3.0D0/4.0D0)-
     &Dc*Hmax*t4*t8*t10*t34*t43*t63*t80*t89*t101*(3.0D0/4.0D0)-Dc*Hmax*t
     &4*t8*t34*t63*t80*t88*t89*t101*t283*(1.0D0/4.0D0)-Dc*Hmax*t4*t8*t34
     &*t63*t80*t89*t101*t128*t285*(1.0D0/4.0D0)-Dc*Hmax*t4*t8*t34*t63*t8
     &0*t89*t101*t148*t282*(1.0D0/4.0D0)-Dc*Hmax*kt*t2*t4*t8*t35*t56*t63
     &*t82*t101*(3.0D0/4.0D0)-Dc*Hmax*kt*t4*t8*t9*t39*t56*t63*t82*t101*(
     &3.0D0/4.0D0)-Dc*Hmax*kt*t4*t8*t10*t43*t56*t63*t82*t101*(3.0D0/4.0D
     &0)-Dc*Hmax*kt*t4*t8*t56*t63*t82*t88*t101*t283*(1.0D0/4.0D0)-Dc*Hma
     &x*kt*t4*t8*t56*t63*t82*t101*t128*t285*(1.0D0/4.0D0)-Dc*Hmax*kt*t4*
     &t8*t56*t63*t82*t101*t148*t282*(1.0D0/4.0D0)-c1*c2*t2*t4*t8*t34*t35
     &*t63*t89*t101*t193*t196*(3.0D0/4.0D0)-c1*c2*t4*t8*t9*t34*t39*t63*t
     &89*t101*t193*t196*(3.0D0/4.0D0)-c1*c2*t4*t8*t10*t34*t43*t63*t89*t1
     &01*t193*t196*(3.0D0/4.0D0)-c1*c2*t4*t8*t34*t63*t88*t89*t101*t193*t
     &196*t283*(1.0D0/4.0D0)-c1*c2*t4*t8*t34*t63*t89*t101*t128*t193*t196
     &*t285*(1.0D0/4.0D0)-c1*c2*t4*t8*t34*t63*t89*t101*t148*t193*t196*t2
     &82*(1.0D0/4.0D0)-c1*c2*kt*t2*t4*t8*t35*t56*t63*t80*t82*t101*t196*(
     &3.0D0/2.0D0)-c1*c2*kt*t4*t8*t9*t39*t56*t63*t80*t82*t101*t196*(3.0D
     &0/2.0D0)-c1*c2*kt*t4*t8*t10*t43*t56*t63*t80*t82*t101*t196*(3.0D0/2
     &.0D0)-c1*c2*kt*t4*t8*t56*t63*t80*t82*t88*t101*t196*t283*(1.0D0/2.0
     &D0)-c1*c2*kt*t4*t8*t56*t63*t80*t82*t101*t128*t196*t285*(1.0D0/2.0D
     &0)-c1*c2*kt*t4*t8*t56*t63*t80*t82*t101*t148*t196*t282*(1.0D0/2.0D0
     &)
      t336 = t299*t300
      t337 = Hmax*t4*t8*t34*t80*t81*t292
      t339 = Hmax*kt*t37*t56*t63*t82*t91*(3.0D0/2.0D0)
      t340 = Hmax*kt*t41*t56*t63*t82*t91*(3.0D0/2.0D0)
      t341 = Hmax*kt*t45*t56*t63*t82*t91*(3.0D0/2.0D0)
      t342 = Hmax*t34*t37*t63*t80*t89*t91*(3.0D0/2.0D0)
      t343 = Hmax*t34*t41*t63*t80*t89*t91*(3.0D0/2.0D0)
      t344 = Hmax*t34*t45*t63*t80*t89*t91*(3.0D0/2.0D0)
      t345 = Dc*Hmax*t4*t13*t34*t80*t81*t128*t266*t269*(1.0D0/2.0D0)
      t347 = Dc*Hmax*t2*t4*t8*t34*t35*t63*t80*t89*t91*(3.0D0/4.0D0)
      t348 = Dc*Hmax*t4*t8*t9*t34*t39*t63*t80*t89*t91*(3.0D0/4.0D0)
      t349 = Dc*Hmax*t4*t8*t10*t34*t43*t63*t80*t89*t91*(3.0D0/4.0D0)
      t351 = Dc*Hmax*kt*t4*t8*t56*t63*t82*t88*t91*t283*(1.0D0/4.0D0)
      t352 = Dc*Hmax*kt*t4*t8*t56*t63*t82*t91*t128*t285*(1.0D0/4.0D0)
      t353 = Dc*Hmax*kt*t4*t8*t56*t63*t82*t91*t148*t282*(1.0D0/4.0D0)
      t354 = c1*c2*t4*t13*t34*t81*t128*t193*t196*t266*t269*(1.0D0/2.0D0)
      t355 = Dc*Hmax*t4*t8*t34*t63*t80*t88*t89*t91*t283*(1.0D0/4.0D0)
      t356 = Dc*Hmax*t4*t8*t34*t63*t80*t89*t91*t128*t285*(1.0D0/4.0D0)
      t357 = Dc*Hmax*t4*t8*t34*t63*t80*t89*t91*t148*t282*(1.0D0/4.0D0)
      t358 = Dc*Hmax*kt*t2*t4*t8*t35*t56*t63*t82*t91*(3.0D0/4.0D0)
      t359 = Dc*Hmax*kt*t4*t8*t9*t39*t56*t63*t82*t91*(3.0D0/4.0D0)
      t360 = Dc*Hmax*kt*t4*t8*t10*t43*t56*t63*t82*t91*(3.0D0/4.0D0)
      t361 = c1*c2*t4*t8*t34*t63*t88*t89*t91*t193*t196*t283*(1.0D0/4.0D0
     &)
      t362 = c1*c2*t4*t8*t34*t63*t89*t91*t128*t193*t196*t285*(1.0D0/4.0D
     &0)
      t363 = c1*c2*t4*t8*t34*t63*t89*t91*t148*t193*t196*t282*(1.0D0/4.0D
     &0)
      t364 = c1*c2*t2*t4*t8*t34*t35*t63*t89*t91*t193*t196*(3.0D0/4.0D0)
      t365 = c1*c2*t4*t8*t9*t34*t39*t63*t89*t91*t193*t196*(3.0D0/4.0D0)
      t366 = c1*c2*t4*t8*t10*t34*t43*t63*t89*t91*t193*t196*(3.0D0/4.0D0)
      t367 = c1*c2*kt*t4*t8*t56*t63*t80*t82*t88*t91*t196*t283*(1.0D0/2.0
     &D0)
      t368 = c1*c2*kt*t4*t8*t56*t63*t80*t82*t91*t128*t196*t285*(1.0D0/2.
     &0D0)
      t369 = c1*c2*kt*t4*t8*t56*t63*t80*t82*t91*t148*t196*t282*(1.0D0/2.
     &0D0)
      t370 = c1*c2*kt*t2*t4*t8*t35*t56*t63*t80*t82*t91*t196*(3.0D0/2.0D0
     &)
      t371 = c1*c2*kt*t4*t8*t9*t39*t56*t63*t80*t82*t91*t196*(3.0D0/2.0D0
     &)
      t372 = c1*c2*kt*t4*t8*t10*t43*t56*t63*t80*t82*t91*t196*(3.0D0/2.0D
     &0)
      t374 = t296*t299
      t376 = Hmax*t4*t8*t34*t80*t81*t293
      t377 = Hmax*kt*t37*t56*t63*t82*t99*(3.0D0/2.0D0)
      t378 = Hmax*kt*t41*t56*t63*t82*t99*(3.0D0/2.0D0)
      t379 = Hmax*kt*t45*t56*t63*t82*t99*(3.0D0/2.0D0)
      t380 = Hmax*t34*t37*t63*t80*t89*t99*(3.0D0/2.0D0)
      t381 = Hmax*t34*t41*t63*t80*t89*t99*(3.0D0/2.0D0)
      t382 = Hmax*t34*t45*t63*t80*t89*t99*(3.0D0/2.0D0)
      t383 = Dc*Hmax*t4*t13*t34*t80*t81*t148*t266*t269*(1.0D0/2.0D0)
      t385 = Dc*Hmax*t2*t4*t8*t34*t35*t63*t80*t89*t99*(3.0D0/4.0D0)
      t386 = Dc*Hmax*t4*t8*t9*t34*t39*t63*t80*t89*t99*(3.0D0/4.0D0)
      t387 = Dc*Hmax*t4*t8*t10*t34*t43*t63*t80*t89*t99*(3.0D0/4.0D0)
      t389 = Dc*Hmax*kt*t4*t8*t56*t63*t82*t88*t99*t283*(1.0D0/4.0D0)
      t390 = Dc*Hmax*kt*t4*t8*t56*t63*t82*t99*t128*t285*(1.0D0/4.0D0)
      t391 = Dc*Hmax*kt*t4*t8*t56*t63*t82*t99*t148*t282*(1.0D0/4.0D0)
      t392 = c1*c2*t4*t13*t34*t81*t148*t193*t196*t266*t269*(1.0D0/2.0D0)
      t393 = Dc*Hmax*t4*t8*t34*t63*t80*t88*t89*t99*t283*(1.0D0/4.0D0)
      t394 = Dc*Hmax*t4*t8*t34*t63*t80*t89*t99*t128*t285*(1.0D0/4.0D0)
      t395 = Dc*Hmax*t4*t8*t34*t63*t80*t89*t99*t148*t282*(1.0D0/4.0D0)
      t396 = Dc*Hmax*kt*t2*t4*t8*t35*t56*t63*t82*t99*(3.0D0/4.0D0)
      t397 = Dc*Hmax*kt*t4*t8*t9*t39*t56*t63*t82*t99*(3.0D0/4.0D0)
      t398 = Dc*Hmax*kt*t4*t8*t10*t43*t56*t63*t82*t99*(3.0D0/4.0D0)
      t399 = c1*c2*t4*t8*t34*t63*t88*t89*t99*t193*t196*t283*(1.0D0/4.0D0
     &)
      t400 = c1*c2*t4*t8*t34*t63*t89*t99*t128*t193*t196*t285*(1.0D0/4.0D
     &0)
      t401 = c1*c2*t4*t8*t34*t63*t89*t99*t148*t193*t196*t282*(1.0D0/4.0D
     &0)
      t402 = c1*c2*t2*t4*t8*t34*t35*t63*t89*t99*t193*t196*(3.0D0/4.0D0)
      t403 = c1*c2*t4*t8*t9*t34*t39*t63*t89*t99*t193*t196*(3.0D0/4.0D0)
      t404 = c1*c2*t4*t8*t10*t34*t43*t63*t89*t99*t193*t196*(3.0D0/4.0D0)
      t405 = c1*c2*kt*t4*t8*t56*t63*t80*t82*t88*t99*t196*t283*(1.0D0/2.0
     &D0)
      t406 = c1*c2*kt*t4*t8*t56*t63*t80*t82*t99*t128*t196*t285*(1.0D0/2.
     &0D0)
      t407 = c1*c2*kt*t4*t8*t56*t63*t80*t82*t99*t148*t196*t282*(1.0D0/2.
     &0D0)
      t408 = c1*c2*kt*t2*t4*t8*t35*t56*t63*t80*t82*t99*t196*(3.0D0/2.0D0
     &)
      t409 = c1*c2*kt*t4*t8*t9*t39*t56*t63*t80*t82*t99*t196*(3.0D0/2.0D0
     &)
      t410 = c1*c2*kt*t4*t8*t10*t43*t56*t63*t80*t82*t99*t196*(3.0D0/2.0D
     &0)
      t411 = ep4*2.0D0
      t412 = et4*2.0D0
      t413 = e4*(-2.0D0)+t411+t412
      t414 = Hmax*t2*t4*t8*t34*t80*t81*3.0D0
      t415 = Dc*Hmax*t12*t13*t34*t35*t80*t81*(3.0D0/4.0D0)
      t416 = c1*c2*t12*t13*t34*t35*t81*t193*t196*(3.0D0/4.0D0)
      t417 = Hmax*t2*t12*t13*t34*t63*t80*t88*t89*t291*(3.0D0/2.0D0)
      t418 = Hmax*t2*t12*t13*t34*t63*t80*t89*t128*t292*(3.0D0/2.0D0)
      t419 = Hmax*t2*t12*t13*t34*t63*t80*t89*t148*t293*(3.0D0/2.0D0)
      t420 = Hmax*kt*t2*t12*t13*t56*t63*t82*t88*t291*(3.0D0/2.0D0)
      t421 = Hmax*kt*t2*t12*t13*t56*t63*t82*t128*t292*(3.0D0/2.0D0)
      t422 = Hmax*kt*t2*t12*t13*t56*t63*t82*t148*t293*(3.0D0/2.0D0)
      t423 = t414+t415+t416+t417+t418+t419+t420+t421+t422-t12*t13*t113*t
     &323*t413*(1.0D0/8.0D0)-Dc*Hmax*t2*t4*t8*t34*t80*t81*(3.0D0/2.0D0)-
     &Hmax*kt*t2*t4*t8*t37*t56*t63*t82*9.0D0-Hmax*kt*t2*t4*t8*t41*t56*t6
     &3*t82*9.0D0-Hmax*kt*t2*t4*t8*t45*t56*t63*t82*9.0D0-Hmax*t2*t4*t8*t
     &34*t37*t63*t80*t89*9.0D0-Hmax*t2*t4*t8*t34*t41*t63*t80*t89*9.0D0-H
     &max*t2*t4*t8*t34*t45*t63*t80*t89*9.0D0-c1*c2*t2*t4*t8*t34*t81*t193
     &*t196*(3.0D0/2.0D0)-Dc*Hmax*kt*t12*t13*t35*t37*t56*t63*t82*(9.0D0/
     &2.0D0)-Dc*Hmax*t12*t13*t34*t35*t37*t63*t80*t89*(9.0D0/2.0D0)-Dc*Hm
     &ax*t2*t9*t12*t13*t34*t39*t63*t80*t89*(9.0D0/2.0D0)-Dc*Hmax*t2*t10*
     &t12*t13*t34*t43*t63*t80*t89*(9.0D0/2.0D0)-Dc*Hmax*t2*t12*t13*t34*t
     &63*t80*t88*t89*t283*(3.0D0/2.0D0)-Dc*Hmax*t2*t12*t13*t34*t63*t80*t
     &89*t128*t285*(3.0D0/2.0D0)-Dc*Hmax*t2*t12*t13*t34*t63*t80*t89*t148
     &*t282*(3.0D0/2.0D0)-c1*c2*t12*t13*t34*t35*t37*t63*t89*t193*t196*(9
     &.0D0/2.0D0)-Dc*Hmax*kt*t2*t9*t12*t13*t39*t56*t63*t82*(9.0D0/2.0D0)
     &-Dc*Hmax*kt*t2*t10*t12*t13*t43*t56*t63*t82*(9.0D0/2.0D0)-Dc*Hmax*k
     &t*t2*t12*t13*t56*t63*t82*t88*t283*(3.0D0/2.0D0)-Dc*Hmax*kt*t2*t12*
     &t13*t56*t63*t82*t128*t285*(3.0D0/2.0D0)-Dc*Hmax*kt*t2*t12*t13*t56*
     &t63*t82*t148*t282*(3.0D0/2.0D0)-c1*c2*kt*t12*t13*t35*t37*t56*t63*t
     &80*t82*t196*9.0D0-c1*c2*t2*t9*t12*t13*t34*t39*t63*t89*t193*t196*(9
     &.0D0/2.0D0)-c1*c2*t2*t10*t12*t13*t34*t43*t63*t89*t193*t196*(9.0D0/
     &2.0D0)-c1*c2*t2*t12*t13*t34*t63*t88*t89*t193*t196*t283*(3.0D0/2.0D
     &0)-c1*c2*t2*t12*t13*t34*t63*t89*t128*t193*t196*t285*(3.0D0/2.0D0)-
     &c1*c2*t2*t12*t13*t34*t63*t89*t148*t193*t196*t282*(3.0D0/2.0D0)-c1*
     &c2*kt*t2*t9*t12*t13*t39*t56*t63*t80*t82*t196*9.0D0-c1*c2*kt*t2*t10
     &*t12*t13*t43*t56*t63*t80*t82*t196*9.0D0-c1*c2*kt*t2*t12*t13*t56*t6
     &3*t80*t82*t88*t196*t283*3.0D0-c1*c2*kt*t2*t12*t13*t56*t63*t80*t82*
     &t128*t196*t285*3.0D0-c1*c2*kt*t2*t12*t13*t56*t63*t80*t82*t148*t196
     &*t282*3.0D0
      t424 = ep5*2.0D0
      t425 = et5*2.0D0
      t426 = e5*(-2.0D0)+t424+t425
      t427 = Hmax*t4*t8*t9*t34*t80*t81*3.0D0
      t428 = Dc*Hmax*t12*t13*t34*t39*t80*t81*(3.0D0/4.0D0)
      t429 = c1*c2*t12*t13*t34*t39*t81*t193*t196*(3.0D0/4.0D0)
      t430 = Hmax*t9*t12*t13*t34*t63*t80*t88*t89*t291*(3.0D0/2.0D0)
      t431 = Hmax*t9*t12*t13*t34*t63*t80*t89*t128*t292*(3.0D0/2.0D0)
      t432 = Hmax*t9*t12*t13*t34*t63*t80*t89*t148*t293*(3.0D0/2.0D0)
      t433 = Hmax*kt*t9*t12*t13*t56*t63*t82*t88*t291*(3.0D0/2.0D0)
      t434 = Hmax*kt*t9*t12*t13*t56*t63*t82*t128*t292*(3.0D0/2.0D0)
      t435 = Hmax*kt*t9*t12*t13*t56*t63*t82*t148*t293*(3.0D0/2.0D0)
      t436 = t427+t428+t429+t430+t431+t432+t433+t434+t435-t12*t13*t113*t
     &323*t426*(1.0D0/8.0D0)-Dc*Hmax*t4*t8*t9*t34*t80*t81*(3.0D0/2.0D0)-
     &Hmax*kt*t4*t8*t9*t37*t56*t63*t82*9.0D0-Hmax*kt*t4*t8*t9*t41*t56*t6
     &3*t82*9.0D0-Hmax*kt*t4*t8*t9*t45*t56*t63*t82*9.0D0-Hmax*t4*t8*t9*t
     &34*t37*t63*t80*t89*9.0D0-Hmax*t4*t8*t9*t34*t41*t63*t80*t89*9.0D0-H
     &max*t4*t8*t9*t34*t45*t63*t80*t89*9.0D0-c1*c2*t4*t8*t9*t34*t81*t193
     &*t196*(3.0D0/2.0D0)-Dc*Hmax*kt*t12*t13*t39*t41*t56*t63*t82*(9.0D0/
     &2.0D0)-Dc*Hmax*t12*t13*t34*t39*t41*t63*t80*t89*(9.0D0/2.0D0)-Dc*Hm
     &ax*t2*t9*t12*t13*t34*t35*t63*t80*t89*(9.0D0/2.0D0)-Dc*Hmax*t9*t10*
     &t12*t13*t34*t43*t63*t80*t89*(9.0D0/2.0D0)-Dc*Hmax*t9*t12*t13*t34*t
     &63*t80*t88*t89*t283*(3.0D0/2.0D0)-Dc*Hmax*t9*t12*t13*t34*t63*t80*t
     &89*t128*t285*(3.0D0/2.0D0)-Dc*Hmax*t9*t12*t13*t34*t63*t80*t89*t148
     &*t282*(3.0D0/2.0D0)-c1*c2*t12*t13*t34*t39*t41*t63*t89*t193*t196*(9
     &.0D0/2.0D0)-Dc*Hmax*kt*t2*t9*t12*t13*t35*t56*t63*t82*(9.0D0/2.0D0)
     &-Dc*Hmax*kt*t9*t10*t12*t13*t43*t56*t63*t82*(9.0D0/2.0D0)-Dc*Hmax*k
     &t*t9*t12*t13*t56*t63*t82*t88*t283*(3.0D0/2.0D0)-Dc*Hmax*kt*t9*t12*
     &t13*t56*t63*t82*t128*t285*(3.0D0/2.0D0)-Dc*Hmax*kt*t9*t12*t13*t56*
     &t63*t82*t148*t282*(3.0D0/2.0D0)-c1*c2*kt*t12*t13*t39*t41*t56*t63*t
     &80*t82*t196*9.0D0-c1*c2*t2*t9*t12*t13*t34*t35*t63*t89*t193*t196*(9
     &.0D0/2.0D0)-c1*c2*t9*t10*t12*t13*t34*t43*t63*t89*t193*t196*(9.0D0/
     &2.0D0)-c1*c2*t9*t12*t13*t34*t63*t88*t89*t193*t196*t283*(3.0D0/2.0D
     &0)-c1*c2*t9*t12*t13*t34*t63*t89*t128*t193*t196*t285*(3.0D0/2.0D0)-
     &c1*c2*t9*t12*t13*t34*t63*t89*t148*t193*t196*t282*(3.0D0/2.0D0)-c1*
     &c2*kt*t2*t9*t12*t13*t35*t56*t63*t80*t82*t196*9.0D0-c1*c2*kt*t9*t10
     &*t12*t13*t43*t56*t63*t80*t82*t196*9.0D0-c1*c2*kt*t9*t12*t13*t56*t6
     &3*t80*t82*t88*t196*t283*3.0D0-c1*c2*kt*t9*t12*t13*t56*t63*t80*t82*
     &t128*t196*t285*3.0D0-c1*c2*kt*t9*t12*t13*t56*t63*t80*t82*t148*t196
     &*t282*3.0D0
      t437 = ep6*2.0D0
      t438 = et6*2.0D0
      t439 = e6*(-2.0D0)+t437+t438
      t440 = Hmax*t4*t8*t10*t34*t80*t81*3.0D0
      t441 = Dc*Hmax*t12*t13*t34*t43*t80*t81*(3.0D0/4.0D0)
      t442 = c1*c2*t12*t13*t34*t43*t81*t193*t196*(3.0D0/4.0D0)
      t443 = Hmax*t10*t12*t13*t34*t63*t80*t88*t89*t291*(3.0D0/2.0D0)
      t444 = Hmax*t10*t12*t13*t34*t63*t80*t89*t128*t292*(3.0D0/2.0D0)
      t445 = Hmax*t10*t12*t13*t34*t63*t80*t89*t148*t293*(3.0D0/2.0D0)
      t446 = Hmax*kt*t10*t12*t13*t56*t63*t82*t88*t291*(3.0D0/2.0D0)
      t447 = Hmax*kt*t10*t12*t13*t56*t63*t82*t128*t292*(3.0D0/2.0D0)
      t448 = Hmax*kt*t10*t12*t13*t56*t63*t82*t148*t293*(3.0D0/2.0D0)
      t449 = t440+t441+t442+t443+t444+t445+t446+t447+t448-t12*t13*t113*t
     &323*t439*(1.0D0/8.0D0)-Dc*Hmax*t4*t8*t10*t34*t80*t81*(3.0D0/2.0D0)
     &-Hmax*kt*t4*t8*t10*t37*t56*t63*t82*9.0D0-Hmax*kt*t4*t8*t10*t41*t56
     &*t63*t82*9.0D0-Hmax*kt*t4*t8*t10*t45*t56*t63*t82*9.0D0-Hmax*t4*t8*
     &t10*t34*t37*t63*t80*t89*9.0D0-Hmax*t4*t8*t10*t34*t41*t63*t80*t89*9
     &.0D0-Hmax*t4*t8*t10*t34*t45*t63*t80*t89*9.0D0-c1*c2*t4*t8*t10*t34*
     &t81*t193*t196*(3.0D0/2.0D0)-Dc*Hmax*kt*t12*t13*t43*t45*t56*t63*t82
     &*(9.0D0/2.0D0)-Dc*Hmax*t12*t13*t34*t43*t45*t63*t80*t89*(9.0D0/2.0D
     &0)-Dc*Hmax*t2*t10*t12*t13*t34*t35*t63*t80*t89*(9.0D0/2.0D0)-Dc*Hma
     &x*t9*t10*t12*t13*t34*t39*t63*t80*t89*(9.0D0/2.0D0)-Dc*Hmax*t10*t12
     &*t13*t34*t63*t80*t88*t89*t283*(3.0D0/2.0D0)-Dc*Hmax*t10*t12*t13*t3
     &4*t63*t80*t89*t128*t285*(3.0D0/2.0D0)-Dc*Hmax*t10*t12*t13*t34*t63*
     &t80*t89*t148*t282*(3.0D0/2.0D0)-c1*c2*t12*t13*t34*t43*t45*t63*t89*
     &t193*t196*(9.0D0/2.0D0)-Dc*Hmax*kt*t2*t10*t12*t13*t35*t56*t63*t82*
     &(9.0D0/2.0D0)-Dc*Hmax*kt*t9*t10*t12*t13*t39*t56*t63*t82*(9.0D0/2.0
     &D0)-Dc*Hmax*kt*t10*t12*t13*t56*t63*t82*t88*t283*(3.0D0/2.0D0)-Dc*H
     &max*kt*t10*t12*t13*t56*t63*t82*t128*t285*(3.0D0/2.0D0)-Dc*Hmax*kt*
     &t10*t12*t13*t56*t63*t82*t148*t282*(3.0D0/2.0D0)-c1*c2*kt*t12*t13*t
     &43*t45*t56*t63*t80*t82*t196*9.0D0-c1*c2*t2*t10*t12*t13*t34*t35*t63
     &*t89*t193*t196*(9.0D0/2.0D0)-c1*c2*t9*t10*t12*t13*t34*t39*t63*t89*
     &t193*t196*(9.0D0/2.0D0)-c1*c2*t10*t12*t13*t34*t63*t88*t89*t193*t19
     &6*t283*(3.0D0/2.0D0)-c1*c2*t10*t12*t13*t34*t63*t89*t128*t193*t196*
     &t285*(3.0D0/2.0D0)-c1*c2*t10*t12*t13*t34*t63*t89*t148*t193*t196*t2
     &82*(3.0D0/2.0D0)-c1*c2*kt*t2*t10*t12*t13*t35*t56*t63*t80*t82*t196*
     &9.0D0-c1*c2*kt*t9*t10*t12*t13*t39*t56*t63*t80*t82*t196*9.0D0-c1*c2
     &*kt*t10*t12*t13*t56*t63*t80*t82*t88*t196*t283*3.0D0-c1*c2*kt*t10*t
     &12*t13*t56*t63*t80*t82*t128*t196*t285*3.0D0-c1*c2*kt*t10*t12*t13*t
     &56*t63*t80*t82*t148*t196*t282*3.0D0
      t450 = t13*t113*t266*t271*v*(1.0D0/2.0D0)
      t451 = t13*t113*t266*t273*v*(1.0D0/2.0D0)
      t452 = t13*t113*t266*t270*v*(1.0D0/2.0D0)
      t453 = t13*t113*t266*t271*v
      t454 = t13*t113*t266*t273*v
      t455 = t13*t113*t266*t270*v
      t463 = t13*t113*t266*t269*t273
      t456 = t453+t455-t463
      t461 = t13*t113*t266*t269*t270
      t457 = t453+t454-t461
      t460 = t13*t113*t266*t269*t271
      t458 = t454+t455-t460
      t459 = t113*t458*v
      t462 = t113*t457*v
      t464 = t113*t456*v
      t465 = t113**2
      t466 = t21+t22+t23+t24-t67-t68-t115-t116-t117-t118+t207+t208
      t467 = t35**2
      t468 = t39**2
      t469 = t43**2
      A0(1,1) = t100+Hmax*kt*t4*t8*t56*t62*t63*t82*(t87+t12*t13*(t57-t58
     &-t59-t60-t61-t74+t75+t76+t77+t78+Sa*bs1*2.0D0-Sa*bs3*2.0D0+Sa*bs1*
     &v*2.0D0-Sa*bs3*v*2.0D0+Sm*bs1*x*2.0D0-Sm*bs3*x*2.0D0+Sm*bs1*v*x*2.
     &0D0-Sm*bs3*v*x*2.0D0))*(e2+e3-ep2-ep3-et2-et3+t21+t22+t23+t24-t25-
     &t26+t27+t28-t29-t30+t31+t32-t57+t58+t59+t60+t61-Sa*bs1*2.0D0-Sa*bs
     &1*v*2.0D0-Sm*bs1*x*2.0D0-Sm*bs1*v*x*2.0D0)*(1.0D0/4.0D0)+Hmax*t4*t
     &8*t34*t62*t63*t80*t88*t89*(t87+t12*t13*(t57-t58-t59-t60-t61+t67+t6
     &8+t70+t72-t74+t75+t76+t77+t78-Sa*bs3*2.0D0-Sa*bs3*v*2.0D0-Sm*bs3*x
     &*2.0D0-Sm*bs3*v*x*2.0D0))*(1.0D0/4.0D0)-1.0D0
      A0(1,2) = t102
      A0(1,3) = Hmax*t4*t8*t34*t62*t80*t81*(-1.0D0/2.0D0)-Hmax*kt*t4*t8*
     &t56*t62*t63*t82*t88*t99*(1.0D0/4.0D0)-Hmax*t4*t8*t34*t62*t63*t80*t
     &88*t89*t99*(1.0D0/4.0D0)
      A0(1,4) = t106
      A0(1,5) = t109
      A0(1,6) = t112
      A0(1,7) = t100+Hmax*kt*t4*t8*t56*t62*t63*t82*t88*t101*(1.0D0/4.0D0
     &)+Hmax*t4*t8*t34*t62*t63*t80*t88*t89*t101*(1.0D0/4.0D0)
      A0(1,8) = t102
      A0(1,9) = -t103-Hmax*kt*t4*t8*t56*t62*t63*t82*t88*t99*(1.0D0/4.0D0
     &)-Hmax*t4*t8*t34*t62*t63*t80*t88*t89*t99*(1.0D0/4.0D0)
      A0(1,10) = t106
      A0(1,11) = t109
      A0(1,12) = t112
      A0(1,13) = Hmax*t4*t8*t34*t80*t81*t88*(1.0D0/2.0D0)-Hmax*t4*t8*t34
     &*t62*t80*t81*(t21+t22+t23+t24-t67-t68+t207+t208-Sm*bs2-Sm*bs3-Sm*b
     &s2*v-Sm*bs3*v)*(1.0D0/2.0D0)+Hmax*t4*t13*t34*t62*t80*t81*t88*t113*
     &(1.0D0/2.0D0)+Hmax*t4*t8*t34*t62*t63*t80*t88*t89*(t120+t124+t125+t
     &126+t127-t12*t47*t113*t114*2.0D0-t12*t49*t113*t114*2.0D0-t12*t51*t
     &113*t114*2.0D0+t12*t13*t20*(t14+t15-t22-t24+t116+t118-Sm*bs1-Sm*bs
     &1*v)*2.0D0)*(1.0D0/4.0D0)+Hmax*kt*t4*t8*t56*t62*t63*t82*t88*t144*(
     &1.0D0/4.0D0)
      A0(2,1) = t131
      A0(2,2) = t100-Hmax*kt*t4*t8*t56*t62*t63*t82*t91*t128*(1.0D0/4.0D0
     &)-Hmax*t4*t8*t34*t62*t63*t80*t89*t91*t128*(1.0D0/4.0D0)-1.0D0
      A0(2,3) = t132
      A0(2,4) = t135
      A0(2,5) = t138
      A0(2,6) = t141
      A0(2,7) = t131
      A0(2,8) = t100-Hmax*kt*t4*t8*t56*t62*t63*t82*t91*t128*(1.0D0/4.0D0
     &)-Hmax*t4*t8*t34*t62*t63*t80*t89*t91*t128*(1.0D0/4.0D0)
      A0(2,9) = t132
      A0(2,10) = t135
      A0(2,11) = t138
      A0(2,12) = t141
      A0(2,13) = Hmax*t4*t8*t34*t80*t81*t128*(1.0D0/2.0D0)-Hmax*t4*t8*t3
     &4*t62*t80*t81*t220*(1.0D0/2.0D0)+Hmax*t4*t13*t34*t62*t80*t81*t113*
     &t128*(1.0D0/2.0D0)+Hmax*kt*t4*t8*t56*t62*t63*t82*t128*t144*(1.0D0/
     &4.0D0)+Hmax*t4*t8*t34*t62*t63*t80*t89*t128*t144*(1.0D0/4.0D0)
      A0(3,1) = t151
      A0(3,2) = t152
      A0(3,3) = t100-Hmax*kt*t4*t8*t56*t62*t63*t82*t99*t148*(1.0D0/4.0D0
     &)-Hmax*t4*t8*t34*t62*t63*t80*t89*t99*t148*(1.0D0/4.0D0)-1.0D0
      A0(3,4) = t155
      A0(3,5) = t158
      A0(3,6) = t161
      A0(3,7) = t151
      A0(3,8) = t152
      A0(3,9) = t100-Hmax*kt*t4*t8*t56*t62*t63*t82*t99*t148*(1.0D0/4.0D0
     &)-Hmax*t4*t8*t34*t62*t63*t80*t89*t99*t148*(1.0D0/4.0D0)
      A0(3,10) = t155
      A0(3,11) = t158
      A0(3,12) = t161
      A0(3,13) = Hmax*t4*t8*t34*t80*t81*t148*(1.0D0/2.0D0)-Hmax*t4*t8*t3
     &4*t62*t80*t81*t232*(1.0D0/2.0D0)+Hmax*t4*t13*t34*t62*t80*t81*t113*
     &t148*(1.0D0/2.0D0)+Hmax*kt*t4*t8*t56*t62*t63*t82*t144*t148*(1.0D0/
     &4.0D0)+Hmax*t4*t8*t34*t62*t63*t80*t89*t144*t148*(1.0D0/4.0D0)
      A0(4,1) = t162
      A0(4,2) = t165
      A0(4,3) = t168
      A0(4,4) = t169-Hmax*kt*t4*t8*t37*t56*t62*t63*t82*9.0D0-Hmax*t4*t8*
     &t34*t37*t62*t63*t80*t89*9.0D0-1.0D0
      A0(4,5) = t170
      A0(4,6) = t171
      A0(4,7) = t162
      A0(4,8) = t165
      A0(4,9) = t168
      A0(4,10) = t169-Hmax*kt*t4*t8*t37*t56*t62*t63*t82*9.0D0-Hmax*t4*t8
     &*t34*t37*t62*t63*t80*t89*9.0D0
      A0(4,11) = t170
      A0(4,12) = t171
      A0(4,13) = Hmax*t2*t34*t80*t81*(-3.0D0)-Hmax*kt*t2*t56*t62*t63*t82
     &*t144*(3.0D0/2.0D0)-Hmax*t2*t34*t62*t63*t80*t89*t144*(3.0D0/2.0D0)
     &+Hmax*t4*t13*t34*t35*t62*t80*t81*t113*(3.0D0/2.0D0)
      A0(5,1) = t174
      A0(5,2) = t177
      A0(5,3) = t180
      A0(5,4) = t170
      A0(5,5) = t169-Hmax*kt*t4*t8*t41*t56*t62*t63*t82*9.0D0-Hmax*t4*t8*
     &t34*t41*t62*t63*t80*t89*9.0D0-1.0D0
      A0(5,6) = t181
      A0(5,7) = t174
      A0(5,8) = t177
      A0(5,9) = t180
      A0(5,10) = t170
      A0(5,11) = t169-Hmax*kt*t4*t8*t41*t56*t62*t63*t82*9.0D0-Hmax*t4*t8
     &*t34*t41*t62*t63*t80*t89*9.0D0
      A0(5,12) = t181
      A0(5,13) = Hmax*t9*t34*t80*t81*(-3.0D0)-Hmax*kt*t9*t56*t62*t63*t82
     &*t144*(3.0D0/2.0D0)-Hmax*t9*t34*t62*t63*t80*t89*t144*(3.0D0/2.0D0)
     &+Hmax*t4*t13*t34*t39*t62*t80*t81*t113*(3.0D0/2.0D0)
      A0(6,1) = t186
      A0(6,2) = t189
      A0(6,3) = t192
      A0(6,4) = t171
      A0(6,5) = t181
      A0(6,6) = t169-Hmax*kt*t4*t8*t45*t56*t62*t63*t82*9.0D0-Hmax*t4*t8*
     &t34*t45*t62*t63*t80*t89*9.0D0-1.0D0
      A0(6,7) = t186
      A0(6,8) = t189
      A0(6,9) = t192
      A0(6,10) = t171
      A0(6,11) = t181
      A0(6,12) = t169-Hmax*kt*t4*t8*t45*t56*t62*t63*t82*9.0D0-Hmax*t4*t8
     &*t34*t45*t62*t63*t80*t89*9.0D0
      A0(6,13) = Hmax*t10*t34*t80*t81*(-3.0D0)-Hmax*kt*t10*t56*t62*t63*t
     &82*t144*(3.0D0/2.0D0)-Hmax*t10*t34*t62*t63*t80*t89*t144*(3.0D0/2.0
     &D0)+Hmax*t4*t13*t34*t43*t62*t80*t81*t113*(3.0D0/2.0D0)
      A0(7,1) = -c1*c2*t4*t8*t34*t62*t81*t193*t196-c1*c2*t4*t8*t34*t62*t
     &63*t88*t89*t101*t193*t196*(1.0D0/4.0D0)-c1*c2*kt*t4*t8*t56*t62*t63
     &*t80*t82*t88*t101*t196*(1.0D0/2.0D0)
      A0(7,2) = t200
      A0(7,3) = t203
      A0(7,4) = t204
      A0(7,5) = t205
      A0(7,6) = t206
      A0(7,7) = -c1*c2*t4*t8*t34*t62*t81*t193*t196-c1*c2*t4*t8*t34*t62*t
     &63*t88*t89*t101*t193*t196*(1.0D0/4.0D0)-c1*c2*kt*t4*t8*t56*t62*t63
     &*t80*t82*t88*t101*t196*(1.0D0/2.0D0)-1.0D0
      A0(7,8) = t200
      A0(7,9) = t203
      A0(7,10) = t204
      A0(7,11) = t205
      A0(7,12) = t206
      A0(7,13) = c1*c2*t4*t8*t34*t81*t88*t193*t196*(-1.0D0/2.0D0)+c1*c2*
     &t4*t8*t34*t62*t81*t193*t196*t466*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t6
     &2*t81*t88*t113*t193*t196*(1.0D0/2.0D0)-c1*c2*t4*t8*t34*t62*t63*t88
     &*t89*t144*t193*t196*(1.0D0/4.0D0)-c1*c2*kt*t4*t8*t56*t62*t63*t80*t
     &82*t88*t144*t196*(1.0D0/2.0D0)
      A0(8,1) = t209
      A0(8,2) = t210+t211-c1*c2*t4*t8*t34*t62*t81*t193*t196
      A0(8,3) = t214
      A0(8,4) = t215
      A0(8,5) = t216
      A0(8,6) = t217
      A0(8,7) = t209
      A0(8,8) = t210+t211-c1*c2*t4*t8*t34*t62*t81*t193*t196-1.0D0
      A0(8,9) = t214
      A0(8,10) = t215
      A0(8,11) = t216
      A0(8,12) = t217
      A0(8,13) = c1*c2*t4*t8*t34*t81*t128*t193*t196*(-1.0D0/2.0D0)+c1*c2
     &*t4*t8*t34*t62*t81*t193*t196*t220*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t
     &62*t81*t113*t128*t193*t196*(1.0D0/2.0D0)-c1*c2*t4*t8*t34*t62*t63*t
     &89*t128*t144*t193*t196*(1.0D0/4.0D0)-c1*c2*kt*t4*t8*t56*t62*t63*t8
     &0*t82*t128*t144*t196*(1.0D0/2.0D0)
      A0(9,1) = t221
      A0(9,2) = t224
      A0(9,3) = t225+t226-c1*c2*t4*t8*t34*t62*t81*t193*t196
      A0(9,4) = t227
      A0(9,5) = t228
      A0(9,6) = t229
      A0(9,7) = t221
      A0(9,8) = t224
      A0(9,9) = t225+t226-c1*c2*t4*t8*t34*t62*t81*t193*t196-1.0D0
      A0(9,10) = t227
      A0(9,11) = t228
      A0(9,12) = t229
      A0(9,13) = c1*c2*t4*t8*t34*t81*t148*t193*t196*(-1.0D0/2.0D0)+c1*c2
     &*t4*t8*t34*t62*t81*t193*t196*t232*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t
     &62*t81*t113*t148*t193*t196*(1.0D0/2.0D0)-c1*c2*t4*t8*t34*t62*t63*t
     &89*t144*t148*t193*t196*(1.0D0/4.0D0)-c1*c2*kt*t4*t8*t56*t62*t63*t8
     &0*t82*t144*t148*t196*(1.0D0/2.0D0)
      A0(10,1) = t235
      A0(10,2) = t236
      A0(10,3) = t237
      A0(10,4) = t238+t239-c1*c2*t4*t8*t34*t62*t81*t193*t196*(3.0D0/2.0D
     &0)
      A0(10,5) = t242
      A0(10,6) = t245
      A0(10,7) = t235
      A0(10,8) = t236
      A0(10,9) = t237
      A0(10,10) = t238+t239-c1*c2*t4*t8*t34*t62*t81*t193*t196*(3.0D0/2.0
     &D0)-1.0D0
      A0(10,11) = t242
      A0(10,12) = t245
      A0(10,13) = c1*c2*t2*t34*t81*t193*t196*3.0D0+c1*c2*t2*t34*t62*t63*
     &t89*t144*t193*t196*(3.0D0/2.0D0)+c1*c2*kt*t2*t56*t62*t63*t80*t82*t
     &144*t196*3.0D0-c1*c2*t4*t13*t34*t35*t62*t81*t113*t193*t196*(3.0D0/
     &2.0D0)
      A0(11,1) = t248
      A0(11,2) = t249
      A0(11,3) = t250
      A0(11,4) = t242
      A0(11,5) = t251+t252-c1*c2*t4*t8*t34*t62*t81*t193*t196*(3.0D0/2.0D
     &0)
      A0(11,6) = t255
      A0(11,7) = t248
      A0(11,8) = t249
      A0(11,9) = t250
      A0(11,10) = t242
      A0(11,11) = t251+t252-c1*c2*t4*t8*t34*t62*t81*t193*t196*(3.0D0/2.0
     &D0)-1.0D0
      A0(11,12) = t255
      A0(11,13) = c1*c2*t9*t34*t81*t193*t196*3.0D0+c1*c2*t9*t34*t62*t63*
     &t89*t144*t193*t196*(3.0D0/2.0D0)+c1*c2*kt*t9*t56*t62*t63*t80*t82*t
     &144*t196*3.0D0-c1*c2*t4*t13*t34*t39*t62*t81*t113*t193*t196*(3.0D0/
     &2.0D0)
      A0(12,1) = t258
      A0(12,2) = t259
      A0(12,3) = t260
      A0(12,4) = t245
      A0(12,5) = t255
      A0(12,6) = t261+t262-c1*c2*t4*t8*t34*t62*t81*t193*t196*(3.0D0/2.0D
     &0)
      A0(12,7) = t258
      A0(12,8) = t259
      A0(12,9) = t260
      A0(12,10) = t245
      A0(12,11) = t255
      A0(12,12) = t261+t262-c1*c2*t4*t8*t34*t62*t81*t193*t196*(3.0D0/2.0
     &D0)-1.0D0
      A0(12,13) = c1*c2*t10*t34*t81*t193*t196*3.0D0+c1*c2*t10*t34*t62*t6
     &3*t89*t144*t193*t196*(3.0D0/2.0D0)+c1*c2*kt*t10*t56*t62*t63*t80*t8
     &2*t144*t196*3.0D0-c1*c2*t4*t13*t34*t43*t62*t81*t113*t193*t196*(3.0
     &D0/2.0D0)
      A0(13,1) = t334
      A0(13,2) = t303+t305+t306+t307+t312+t316+t319+t320+t336+t337+t339+
     &t340+t341+t342+t343+t344+t345+t347+t348+t349+t351+t352+t353+t354+t
     &355+t356+t357+t358+t359+t360+t361+t362+t363+t364+t365+t366+t367+t3
     &68+t369+t370+t371+t372-t276*t294-t276*t296-t8*t266*t269*t304*(1.0D
     &0/2.0D0)-Hmax*t4*t8*t34*t80*t81*t291*(1.0D0/2.0D0)-Hmax*t4*t8*t34*
     &t80*t81*t293*(1.0D0/2.0D0)-Dc*Hmax*t4*t8*t34*t80*t81*t285-Hmax*t4*
     &t13*t34*t80*t81*t128*t266*t269*(1.0D0/2.0D0)-c1*c2*t4*t8*t34*t81*t
     &193*t196*t285-Dc*Hmax*t4*t13*t34*t80*t81*t88*t266*v*(1.0D0/2.0D0)-
     &Dc*Hmax*t4*t13*t34*t80*t81*t148*t266*v*(1.0D0/2.0D0)-Hmax*kt*t4*t8
     &*t56*t63*t82*t88*t91*t291*(1.0D0/4.0D0)-Hmax*kt*t4*t8*t56*t63*t82*
     &t91*t128*t292*(1.0D0/4.0D0)-Hmax*kt*t4*t8*t56*t63*t82*t91*t148*t29
     &3*(1.0D0/4.0D0)-Hmax*t4*t8*t34*t63*t80*t88*t89*t91*t291*(1.0D0/4.0
     &D0)-Hmax*t4*t8*t34*t63*t80*t89*t91*t128*t292*(1.0D0/4.0D0)-Hmax*t4
     &*t8*t34*t63*t80*t89*t91*t148*t293*(1.0D0/4.0D0)-c1*c2*t4*t13*t34*t
     &81*t88*t193*t196*t266*v*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t81*t148*t1
     &93*t196*t266*v*(1.0D0/2.0D0)
      A0(13,3) = t312+t314+t316+t317+t318+t319+t320+t321+t374+t376+t377+
     &t378+t379+t380+t381+t382+t383+t385+t386+t387+t389+t390+t391+t392+t
     &393+t394+t395+t396+t397+t398+t399+t400+t401+t402+t403+t404+t405+t4
     &06+t407+t408+t409+t410-t276*t294-t276*t300-t8*t266*t269*t302*(1.0D
     &0/2.0D0)-Hmax*t4*t8*t34*t80*t81*t291*(1.0D0/2.0D0)-Hmax*t4*t8*t34*
     &t80*t81*t292*(1.0D0/2.0D0)-Dc*Hmax*t4*t8*t34*t80*t81*t282-Hmax*t4*
     &t13*t34*t80*t81*t148*t266*t269*(1.0D0/2.0D0)-c1*c2*t4*t8*t34*t81*t
     &193*t196*t282-Dc*Hmax*t4*t13*t34*t80*t81*t88*t266*v*(1.0D0/2.0D0)-
     &Dc*Hmax*t4*t13*t34*t80*t81*t128*t266*v*(1.0D0/2.0D0)-Hmax*kt*t4*t8
     &*t56*t63*t82*t88*t99*t291*(1.0D0/4.0D0)-Hmax*kt*t4*t8*t56*t63*t82*
     &t99*t128*t292*(1.0D0/4.0D0)-Hmax*kt*t4*t8*t56*t63*t82*t99*t148*t29
     &3*(1.0D0/4.0D0)-Hmax*t4*t8*t34*t63*t80*t88*t89*t99*t291*(1.0D0/4.0
     &D0)-Hmax*t4*t8*t34*t63*t80*t89*t99*t128*t292*(1.0D0/4.0D0)-Hmax*t4
     &*t8*t34*t63*t80*t89*t99*t148*t293*(1.0D0/4.0D0)-c1*c2*t4*t13*t34*t
     &81*t88*t193*t196*t266*v*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t81*t128*t1
     &93*t196*t266*v*(1.0D0/2.0D0)
      A0(13,4) = t423
      A0(13,5) = t436
      A0(13,6) = t449
      A0(13,7) = t334
      A0(13,8) = t303+t305+t306+t307+t312+t316+t319+t320-t335+t336+t337-
     &t338+t339+t340+t341+t342+t343+t344+t345-t346+t347+t348+t349-t350+t
     &351+t352+t353+t354+t355+t356+t357+t358+t359+t360+t361+t362+t363+t3
     &64+t365+t366+t367+t368+t369+t370+t371+t372-t276*t294-t8*t266*t269*
     &t304*(1.0D0/2.0D0)-Hmax*t4*t8*t34*t80*t81*t291*(1.0D0/2.0D0)-Dc*Hm
     &ax*t4*t8*t34*t80*t81*t285-Hmax*t4*t13*t34*t80*t81*t128*t266*t269*(
     &1.0D0/2.0D0)-c1*c2*t4*t8*t34*t81*t193*t196*t285-Dc*Hmax*t4*t13*t34
     &*t80*t81*t88*t266*v*(1.0D0/2.0D0)-Hmax*kt*t4*t8*t56*t63*t82*t88*t9
     &1*t291*(1.0D0/4.0D0)-Hmax*kt*t4*t8*t56*t63*t82*t91*t128*t292*(1.0D
     &0/4.0D0)-Hmax*kt*t4*t8*t56*t63*t82*t91*t148*t293*(1.0D0/4.0D0)-Hma
     &x*t4*t8*t34*t63*t80*t88*t89*t91*t291*(1.0D0/4.0D0)-Hmax*t4*t8*t34*
     &t63*t80*t89*t91*t128*t292*(1.0D0/4.0D0)-Hmax*t4*t8*t34*t63*t80*t89
     &*t91*t148*t293*(1.0D0/4.0D0)-c1*c2*t4*t13*t34*t81*t88*t193*t196*t2
     &66*v*(1.0D0/2.0D0)
      A0(13,9) = t312+t314+t316+t317+t318+t319+t320+t321-t373+t374-t375+
     &t376+t377+t378+t379+t380+t381+t382+t383-t384+t385+t386+t387-t388+t
     &389+t390+t391+t392+t393+t394+t395+t396+t397+t398+t399+t400+t401+t4
     &02+t403+t404+t405+t406+t407+t408+t409+t410-t276*t294-t8*t266*t269*
     &t302*(1.0D0/2.0D0)-Hmax*t4*t8*t34*t80*t81*t291*(1.0D0/2.0D0)-Dc*Hm
     &ax*t4*t8*t34*t80*t81*t282-Hmax*t4*t13*t34*t80*t81*t148*t266*t269*(
     &1.0D0/2.0D0)-c1*c2*t4*t8*t34*t81*t193*t196*t282-Dc*Hmax*t4*t13*t34
     &*t80*t81*t88*t266*v*(1.0D0/2.0D0)-Hmax*kt*t4*t8*t56*t63*t82*t88*t9
     &9*t291*(1.0D0/4.0D0)-Hmax*kt*t4*t8*t56*t63*t82*t99*t128*t292*(1.0D
     &0/4.0D0)-Hmax*kt*t4*t8*t56*t63*t82*t99*t148*t293*(1.0D0/4.0D0)-Hma
     &x*t4*t8*t34*t63*t80*t88*t89*t99*t291*(1.0D0/4.0D0)-Hmax*t4*t8*t34*
     &t63*t80*t89*t99*t128*t292*(1.0D0/4.0D0)-Hmax*t4*t8*t34*t63*t80*t89
     &*t99*t148*t293*(1.0D0/4.0D0)-c1*c2*t4*t13*t34*t81*t88*t193*t196*t2
     &66*v*(1.0D0/2.0D0)
      A0(13,10) = t423
      A0(13,11) = t436
      A0(13,12) = t449
      A0(13,13) = t296*(t459+t464-t113*t457)+t300*(t459+t462-t113*t456)+
     &t294*(t462+t464-t113*t458)+t302*(t450+t451-t13*t113*t266*t269*t270
     &*(1.0D0/2.0D0))+t301*(t451+t452-t13*t113*t266*t269*t271*(1.0D0/2.0
     &D0))+t304*(t450+t452-t13*t113*t266*t269*t273*(1.0D0/2.0D0))-a1*(n2
     &*(-x+1.0D0)**(n2-1.0D0)+n1*x**(n1-1.0D0))*(1.0D0/2.0D0)-t12*t114*t
     &323*t465*t467*(1.0D0/4.0D0)-t12*t114*t323*t465*t468*(1.0D0/4.0D0)-
     &t12*t114*t323*t465*t469*(1.0D0/4.0D0)-Hmax*kt*t37*t56*t63*t82*t144
     &*(3.0D0/2.0D0)-Hmax*kt*t41*t56*t63*t82*t144*(3.0D0/2.0D0)-Hmax*kt*
     &t45*t56*t63*t82*t144*(3.0D0/2.0D0)-Hmax*t34*t37*t63*t80*t89*t144*(
     &3.0D0/2.0D0)-Hmax*t34*t41*t63*t80*t89*t144*(3.0D0/2.0D0)-Hmax*t34*
     &t45*t63*t80*t89*t144*(3.0D0/2.0D0)-Hmax*t4*t8*t34*t80*t81*t220*t29
     &2*(1.0D0/2.0D0)-Hmax*t4*t8*t34*t80*t81*t232*t293*(1.0D0/2.0D0)+Hma
     &x*t4*t8*t34*t80*t81*t88*t458*(1.0D0/2.0D0)+Hmax*t4*t8*t34*t80*t81*
     &t128*t456*(1.0D0/2.0D0)+Hmax*t4*t8*t34*t80*t81*t148*t457*(1.0D0/2.
     &0D0)-Hmax*t4*t8*t34*t80*t81*t291*t466*(1.0D0/2.0D0)+Dc*Hmax*t4*t8*
     &t34*t80*t81*t220*t285*(1.0D0/2.0D0)+Dc*Hmax*t4*t8*t34*t80*t81*t232
     &*t282*(1.0D0/2.0D0)-Dc*Hmax*t4*t8*t34*t80*t81*t88*t458*(1.0D0/2.0D
     &0)-Dc*Hmax*t4*t8*t34*t80*t81*t128*t456*(1.0D0/2.0D0)-Dc*Hmax*t4*t8
     &*t34*t80*t81*t148*t457*(1.0D0/2.0D0)+Dc*Hmax*t12*t34*t80*t81*t113*
     &t114*t467*(3.0D0/4.0D0)+Dc*Hmax*t12*t34*t80*t81*t113*t114*t468*(3.
     &0D0/4.0D0)+Dc*Hmax*t12*t34*t80*t81*t113*t114*t469*(3.0D0/4.0D0)+Hm
     &ax*t2*t4*t13*t34*t35*t80*t81*t113*3.0D0+Hmax*t4*t9*t13*t34*t39*t80
     &*t81*t113*3.0D0+Hmax*t4*t10*t13*t34*t43*t80*t81*t113*3.0D0+Hmax*t4
     &*t13*t34*t80*t81*t88*t113*t291*(1.0D0/2.0D0)+Hmax*t4*t13*t34*t80*t
     &81*t113*t128*t292*(1.0D0/2.0D0)+Hmax*t4*t13*t34*t80*t81*t113*t148*
     &t293*(1.0D0/2.0D0)+Dc*Hmax*t4*t8*t34*t80*t81*t283*(t21+t22+t23+t24
     &-t67-t68-t115-t116-t117-t118+t207+t208)*(1.0D0/2.0D0)+c1*c2*t4*t8*
     &t34*t81*t193*t196*t283*(t21+t22+t23+t24-t67-t68-t115-t116-t117-t11
     &8+t207+t208)*(1.0D0/2.0D0)-Dc*Hmax*t2*t4*t13*t34*t35*t80*t81*t113*
     &(3.0D0/2.0D0)-Dc*Hmax*t4*t9*t13*t34*t39*t80*t81*t113*(3.0D0/2.0D0)
     &-Dc*Hmax*t4*t10*t13*t34*t43*t80*t81*t113*(3.0D0/2.0D0)-Dc*Hmax*t4*
     &t13*t34*t80*t81*t88*t113*t283*(1.0D0/2.0D0)-Dc*Hmax*t4*t13*t34*t80
     &*t81*t113*t128*t285*(1.0D0/2.0D0)-Dc*Hmax*t4*t13*t34*t80*t81*t113*
     &t148*t282*(1.0D0/2.0D0)+Hmax*kt*t4*t8*t56*t63*t82*t88*t144*t291*(1
     &.0D0/4.0D0)+Hmax*kt*t4*t8*t56*t63*t82*t128*t144*t292*(1.0D0/4.0D0)
     &+Hmax*kt*t4*t8*t56*t63*t82*t144*t148*t293*(1.0D0/4.0D0)+Hmax*t4*t8
     &*t34*t63*t80*t88*t89*t144*t291*(1.0D0/4.0D0)+Hmax*t4*t8*t34*t63*t8
     &0*t89*t128*t144*t292*(1.0D0/4.0D0)+Hmax*t4*t8*t34*t63*t80*t89*t144
     &*t148*t293*(1.0D0/4.0D0)+c1*c2*t4*t8*t34*t81*t193*t196*t220*t285*(
     &1.0D0/2.0D0)+c1*c2*t4*t8*t34*t81*t193*t196*t232*t282*(1.0D0/2.0D0)
     &-c1*c2*t4*t8*t34*t81*t88*t193*t196*t458*(1.0D0/2.0D0)-c1*c2*t4*t8*
     &t34*t81*t128*t193*t196*t456*(1.0D0/2.0D0)-c1*c2*t4*t8*t34*t81*t148
     &*t193*t196*t457*(1.0D0/2.0D0)+c1*c2*t12*t34*t81*t113*t114*t193*t19
     &6*t467*(3.0D0/4.0D0)+c1*c2*t12*t34*t81*t113*t114*t193*t196*t468*(3
     &.0D0/4.0D0)+c1*c2*t12*t34*t81*t113*t114*t193*t196*t469*(3.0D0/4.0D
     &0)-Dc*Hmax*t2*t4*t8*t34*t35*t63*t80*t89*t144*(3.0D0/4.0D0)-Dc*Hmax
     &*t4*t8*t9*t34*t39*t63*t80*t89*t144*(3.0D0/4.0D0)-Dc*Hmax*t4*t8*t10
     &*t34*t43*t63*t80*t89*t144*(3.0D0/4.0D0)-Dc*Hmax*t4*t8*t34*t63*t80*
     &t88*t89*t144*t283*(1.0D0/4.0D0)-Dc*Hmax*t4*t8*t34*t63*t80*t89*t128
     &*t144*t285*(1.0D0/4.0D0)-Dc*Hmax*t4*t8*t34*t63*t80*t89*t144*t148*t
     &282*(1.0D0/4.0D0)-c1*c2*t2*t4*t13*t34*t35*t81*t113*t193*t196*(3.0D
     &0/2.0D0)-c1*c2*t4*t9*t13*t34*t39*t81*t113*t193*t196*(3.0D0/2.0D0)-
     &c1*c2*t4*t10*t13*t34*t43*t81*t113*t193*t196*(3.0D0/2.0D0)-c1*c2*t4
     &*t13*t34*t81*t88*t113*t193*t196*t283*(1.0D0/2.0D0)-c1*c2*t4*t13*t3
     &4*t81*t113*t128*t193*t196*t285*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t81*
     &t113*t148*t193*t196*t282*(1.0D0/2.0D0)-Dc*Hmax*kt*t2*t4*t8*t35*t56
     &*t63*t82*t144*(3.0D0/4.0D0)-Dc*Hmax*kt*t4*t8*t9*t39*t56*t63*t82*t1
     &44*(3.0D0/4.0D0)-Dc*Hmax*kt*t4*t8*t10*t43*t56*t63*t82*t144*(3.0D0/
     &4.0D0)-Dc*Hmax*kt*t4*t8*t56*t63*t82*t88*t144*t283*(1.0D0/4.0D0)-Dc
     &*Hmax*kt*t4*t8*t56*t63*t82*t128*t144*t285*(1.0D0/4.0D0)-Dc*Hmax*kt
     &*t4*t8*t56*t63*t82*t144*t148*t282*(1.0D0/4.0D0)-c1*c2*t2*t4*t8*t34
     &*t35*t63*t89*t144*t193*t196*(3.0D0/4.0D0)-c1*c2*t4*t8*t9*t34*t39*t
     &63*t89*t144*t193*t196*(3.0D0/4.0D0)-c1*c2*t4*t8*t10*t34*t43*t63*t8
     &9*t144*t193*t196*(3.0D0/4.0D0)-c1*c2*t4*t8*t34*t63*t88*t89*t144*t1
     &93*t196*t283*(1.0D0/4.0D0)-c1*c2*t4*t8*t34*t63*t89*t128*t144*t193*
     &t196*t285*(1.0D0/4.0D0)-c1*c2*t4*t8*t34*t63*t89*t144*t148*t193*t19
     &6*t282*(1.0D0/4.0D0)-c1*c2*kt*t2*t4*t8*t35*t56*t63*t80*t82*t144*t1
     &96*(3.0D0/2.0D0)-c1*c2*kt*t4*t8*t9*t39*t56*t63*t80*t82*t144*t196*(
     &3.0D0/2.0D0)-c1*c2*kt*t4*t8*t10*t43*t56*t63*t80*t82*t144*t196*(3.0
     &D0/2.0D0)-c1*c2*kt*t4*t8*t56*t63*t80*t82*t88*t144*t196*t283*(1.0D0
     &/2.0D0)-c1*c2*kt*t4*t8*t56*t63*t80*t82*t128*t144*t196*t285*(1.0D0/
     &2.0D0)-c1*c2*kt*t4*t8*t56*t63*t80*t82*t144*t148*t196*t282*(1.0D0/2
     &.0D0)





       NR_JAC=A0
       !NR_JAC(13,)=t0
!   **************** New Elemnts *********!**********************************************************************
!   **********************************************************************
   

       elseif (index1==-1)then ! Reverse
      t3 = v+1.0D0
      t4 = 1.0D0/t3
      t5 = Sa*x
      t6 = Sm*x
      t7 = Sa-t5+t6
      t8 = 1.0D0/t7
      t35 = -e4+ep4+et4
      t36 = t4*t8*t35*(1.0D0/2.0D0)
      t2 = bs4-t36
      t39 = -e5+ep5+et5
      t40 = t4*t8*t39*(1.0D0/2.0D0)
      t9 = bs5-t40
      t43 = -e6+ep6+et6
      t44 = t4*t8*t43*(1.0D0/2.0D0)
      t10 = bs6-t44
      t14 = Sa*bs1
      t15 = Sa*bs1*v
      t16 = Sa*bs1*x
      t17 = Sm*bs1*x
      t18 = Sa*bs1*v*x
      t19 = Sm*bs1*v*x
      t21 = Sa*bs2
      t23 = Sa*bs2*v
      t25 = Sa*bs2*x
      t27 = Sm*bs2*x
      t29 = Sa*bs2*v*x
      t31 = Sm*bs2*v*x
      t11 = e1-e2-ep1+ep2-et1+et2+t14+t15-t16+t17-t18+t19-t21-t23+t25-t2
     &7+t29-t31
      t12 = 1.0D0/t3**2
      t13 = 1.0D0/t7**2
      t22 = Sa*bs3
      t24 = Sa*bs3*v
      t26 = Sa*bs3*x
      t28 = Sm*bs3*x
      t30 = Sa*bs3*v*x
      t32 = Sm*bs3*v*x
      t20 = e1-e3-ep1+ep3-et1+et3+t14+t15-t16+t17-t18+t19-t22-t24+t26-t2
     &8+t30-t32
      t33 = e2-e3-ep2+ep3-et2+et3+t21-t22+t23-t24-t25+t26+t27-t28-t29+t3
     &0+t31-t32
      t34 = sqrt(2.0D0)
      t37 = t2**2
      t38 = t37*6.0D0
      t41 = t9**2
      t42 = t41*6.0D0
      t45 = t10**2
      t46 = t45*6.0D0
      t47 = t11**2
      t48 = t12*t13*t47
      t49 = t20**2
      t50 = t12*t13*t49
      t51 = t33**2
      t52 = t12*t13*t51
      t53 = t38+t42+t46+t48+t50+t52
      t54 = abs(t53)
      t61 = sqrt(t54)
      t62 = kt*t34*t61*(1.0D0/2.0D0)
      t63 = exp(-t62)
      t55 = t63-1.0D0
      t56 = e1*2.0D0
      t57 = ep1*2.0D0
      t58 = et1*2.0D0
      t59 = Sa*bs1*x*2.0D0
      t60 = Sa*bs1*v*x*2.0D0
      t64 = t55**2
      t65 = c2*ztd
      t66 = t65+1.0D0
      t67 = 1.0D0/t66
      t68 = x-xo
      t69 = (t53/abs(t53))
      t70 = e2*2.0D0
      t71 = ep2*2.0D0
      t72 = et2*2.0D0
      t73 = Sa*bs1*2.0D0
      t74 = Sa*bs1*v*2.0D0
      t75 = Sa*bs2*x*2.0D0
      t76 = Sm*bs1*x*2.0D0
      t77 = Sa*bs2*v*x*2.0D0
      t78 = Sm*bs1*v*x*2.0D0
      t87 = Sa*bs2*2.0D0
      t88 = Sa*bs2*v*2.0D0
      t89 = Sm*bs2*x*2.0D0
      t90 = Sm*bs2*v*x*2.0D0
      t79 = t56-t57-t58-t59-t60-t70+t71+t72+t73+t74+t75+t76+t77+t78-t87-
     &t88-t89-t90
      t80 = e3*2.0D0
      t81 = ep3*2.0D0
      t82 = et3*2.0D0
      t83 = Sa*bs3*x*2.0D0
      t84 = Sa*bs3*v*x*2.0D0
      t85 = 1.0D0/sqrt(t54)
      t86 = 1.0D0/t54**(3.0D0/2.0D0)
      t91 = t12*t13*t79
      t92 = e2+e3-ep2-ep3-et2-et3+t21+t22+t23+t24-t25-t26+t27+t28-t29-t3
     &0+t31+t32-t56+t57+t58+t59+t60-t73-t74-t76-t78
      t93 = 1.0D0/t54
      t96 = Sa*bs3*2.0D0
      t97 = Sa*bs3*v*2.0D0
      t98 = Sm*bs3*x*2.0D0
      t99 = Sm*bs3*v*x*2.0D0
      t94 = t70-t71-t72-t75-t77-t80+t81+t82+t83+t84+t87+t88+t89+t90-t96-
     &t97-t98-t99
      t100 = t12*t13*t94
      t95 = t91-t100
      t101 = t56-t57-t58-t59-t60+t73+t74+t76+t78-t80+t81+t82+t83+t84-t96
     &-t97-t98-t99
      t102 = t12*t13*t101
      t103 = t100+t102
      t104 = c1*c2*t4*t8*t34*t64*t67*t68*t85
      t105 = t91+t102
      t107 = c1*c2*t4*t8*t34*t64*t67*t68*t85*(1.0D0/2.0D0)
      t106 = -t107-c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t92*t95*(1.0D0/4.
     &0D0)-c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t92*t93*t95*(1.0D0/2.0D0)
      t108 = c1*c2*t2*t12*t13*t34*t64*t67*t68*t69*t86*t92*(3.0D0/2.0D0)
      t109 = c1*c2*kt*t2*t12*t13*t55*t63*t67*t68*t69*t92*t93*3.0D0
      t110 = t108+t109
      t111 = c1*c2*t9*t12*t13*t34*t64*t67*t68*t69*t86*t92*(3.0D0/2.0D0)
      t112 = c1*c2*kt*t9*t12*t13*t55*t63*t67*t68*t69*t92*t93*3.0D0
      t113 = t111+t112
      t114 = c1*c2*t10*t12*t13*t34*t64*t67*t68*t69*t86*t92*(3.0D0/2.0D0)
      t115 = c1*c2*kt*t10*t12*t13*t55*t63*t67*t68*t69*t92*t93*3.0D0
      t116 = t114+t115
      t117 = Sa-Sm
      t118 = 1.0D0/t7**3
      t119 = Sm*bs2
      t120 = Sm*bs3
      t121 = Sm*bs2*v
      t122 = Sm*bs3*v
      t125 = Sm*bs1
      t126 = Sm*bs1*v
      t123 = t14+t15-t21-t23+t119+t121-t125-t126
      t124 = t11*t12*t13*t123*2.0D0
      t127 = t21-t22+t23-t24-t119+t120-t121+t122
      t128 = t12*t13*t33*t127*2.0D0
      t129 = t2*t4*t13*t35*t117*6.0D0
      t130 = t4*t9*t13*t39*t117*6.0D0
      t131 = t4*t10*t13*t43*t117*6.0D0
      t132 = e1+e3-ep1-ep3-et1-et3+t14+t15-t16+t17-t18+t19+t22+t24-t26+t
     &28-t30+t32-t70+t71+t72+t75+t77-t87-t88-t89-t90
      t133 = c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t105*t132*(1.0D0/4.0D0)
      t134 = c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93*t105*t132*(1.0D0/2.0
     &D0)
      t135 = -t107+t133+t134
      t136 = -t107-c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t103*t132*(1.0D0/
     &4.0D0)-c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93*t103*t132*(1.0D0/2.0
     &D0)
      t137 = c1*c2*t2*t12*t13*t34*t64*t67*t68*t69*t86*t132*(3.0D0/2.0D0)
      t138 = c1*c2*kt*t2*t12*t13*t55*t63*t67*t68*t69*t93*t132*3.0D0
      t139 = t137+t138
      t140 = c1*c2*t9*t12*t13*t34*t64*t67*t68*t69*t86*t132*(3.0D0/2.0D0)
      t141 = c1*c2*kt*t9*t12*t13*t55*t63*t67*t68*t69*t93*t132*3.0D0
      t142 = t140+t141
      t143 = c1*c2*t10*t12*t13*t34*t64*t67*t68*t69*t86*t132*(3.0D0/2.0D0
     &)
      t144 = c1*c2*kt*t10*t12*t13*t55*t63*t67*t68*t69*t93*t132*3.0D0
      t145 = t143+t144
      t146 = t14+t15-t22-t24+t120+t122-t125-t126
      t147 = t12*t13*t20*t146*2.0D0
      t149 = t12*t47*t117*t118*2.0D0
      t150 = t12*t49*t117*t118*2.0D0
      t151 = t12*t51*t117*t118*2.0D0
      t148 = t124+t128+t129+t130+t131+t147-t149-t150-t151
      t152 = e1+e2-ep1-ep2-et1-et2+t14+t15-t16+t17-t18+t19+t21+t23-t25+t
     &27-t29+t31-t80+t81+t82+t83+t84-t96-t97-t98-t99
      t153 = c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t105*t152*(1.0D0/4.0D0)
      t154 = c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93*t105*t152*(1.0D0/2.0
     &D0)
      t155 = -t107+t153+t154
      t156 = -t107-c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t95*t152*(1.0D0/4
     &.0D0)-c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93*t95*t152*(1.0D0/2.0D0
     &)
      t157 = c1*c2*t2*t12*t13*t34*t64*t67*t68*t69*t86*t152*(3.0D0/2.0D0)
      t158 = c1*c2*kt*t2*t12*t13*t55*t63*t67*t68*t69*t93*t152*3.0D0
      t159 = t157+t158
      t160 = c1*c2*t9*t12*t13*t34*t64*t67*t68*t69*t86*t152*(3.0D0/2.0D0)
      t161 = c1*c2*kt*t9*t12*t13*t55*t63*t67*t68*t69*t93*t152*3.0D0
      t162 = t160+t161
      t163 = c1*c2*t10*t12*t13*t34*t64*t67*t68*t69*t86*t152*(3.0D0/2.0D0
     &)
      t164 = c1*c2*kt*t10*t12*t13*t55*t63*t67*t68*t69*t93*t152*3.0D0
      t165 = t163+t164
      t166 = c1*c2*t2*t34*t64*t67*t68*t69*t86*t105*(-3.0D0/2.0D0)-c1*c2*
     &kt*t2*t55*t63*t67*t68*t69*t93*t105*3.0D0
      t167 = c1*c2*t2*t34*t64*t67*t68*t69*t86*t95*(3.0D0/2.0D0)
      t168 = c1*c2*kt*t2*t55*t63*t67*t68*t69*t93*t95*3.0D0
      t169 = t167+t168
      t170 = c1*c2*t2*t34*t64*t67*t68*t69*t86*t103*(3.0D0/2.0D0)
      t171 = c1*c2*kt*t2*t55*t63*t67*t68*t69*t93*t103*3.0D0
      t172 = t170+t171
      t173 = c1*c2*t4*t8*t34*t64*t67*t68*t85*(3.0D0/2.0D0)
      t176 = c1*c2*t2*t4*t8*t9*t34*t64*t67*t68*t69*t86*9.0D0
      t177 = c1*c2*kt*t2*t4*t8*t9*t55*t63*t67*t68*t69*t93*1.8D1
      t174 = -t176-t177
      t186 = c1*c2*t2*t4*t8*t10*t34*t64*t67*t68*t69*t86*9.0D0
      t187 = c1*c2*kt*t2*t4*t8*t10*t55*t63*t67*t68*t69*t93*1.8D1
      t175 = -t186-t187
      t178 = c1*c2*t9*t34*t64*t67*t68*t69*t86*t105*(-3.0D0/2.0D0)-c1*c2*
     &kt*t9*t55*t63*t67*t68*t69*t93*t105*3.0D0
      t179 = c1*c2*t9*t34*t64*t67*t68*t69*t86*t95*(3.0D0/2.0D0)
      t180 = c1*c2*kt*t9*t55*t63*t67*t68*t69*t93*t95*3.0D0
      t181 = t179+t180
      t182 = c1*c2*t9*t34*t64*t67*t68*t69*t86*t103*(3.0D0/2.0D0)
      t183 = c1*c2*kt*t9*t55*t63*t67*t68*t69*t93*t103*3.0D0
      t184 = t182+t183
      t188 = c1*c2*t4*t8*t9*t10*t34*t64*t67*t68*t69*t86*9.0D0
      t189 = c1*c2*kt*t4*t8*t9*t10*t55*t63*t67*t68*t69*t93*1.8D1
      t185 = -t188-t189
      t190 = c1*c2*t10*t34*t64*t67*t68*t69*t86*t105*(-3.0D0/2.0D0)-c1*c2
     &*kt*t10*t55*t63*t67*t68*t69*t93*t105*3.0D0
      t191 = c1*c2*t10*t34*t64*t67*t68*t69*t86*t95*(3.0D0/2.0D0)
      t192 = c1*c2*kt*t10*t55*t63*t67*t68*t69*t93*t95*3.0D0
      t193 = t191+t192
      t194 = c1*c2*t10*t34*t64*t67*t68*t69*t86*t103*(3.0D0/2.0D0)
      t195 = c1*c2*kt*t10*t55*t63*t67*t68*t69*t93*t103*3.0D0
      t196 = t194+t195
      t197 = v**2
      t198 = t197*2.0D0
      t199 = t198+v-1.0D0
      t200 = 1.0D0/t199
      t201 = Tn*alpha
      t202 = v-1.0D0
      t203 = e1*v
      t204 = ep2*v
      t205 = ep3*v
      t206 = et2*v
      t207 = et3*v
      t208 = Tn*alpha*v
      t210 = alpha*too
      t211 = e2*v
      t212 = e3*v
      t213 = ep1*v
      t214 = et1*v
      t215 = alpha*too*v
      t209 = -e1+ep1+et1+t201+t203+t204+t205+t206+t207+t208-t210-t211-t2
     &12-t213-t214-t215
      t216 = -e2+ep2+et2+t201-t203-t204+t205-t206+t207+t208-t210+t211-t2
     &12+t213+t214-t215
      t217 = -e3+ep3+et3+t201-t203+t204-t205+t206-t207+t208-t210-t211+t2
     &12+t213+t214-t215
      t218 = c1*c2*t4*t13*t34*t64*t67*t85*t200*t217*(1.0D0/2.0D0)
      t219 = -e1+ep1+et1+t201-t210
      t220 = -e2+ep2+et2+t201-t210
      t221 = -e3+ep3+et3+t201-t210
      t222 = c1*c2*t4*t13*t34*t64*t67*t85*t200*t209*(1.0D0/2.0D0)
      t223 = c1*c2*t4*t13*t34*t64*t67*t85*t200*t216*(1.0D0/2.0D0)
      t224 = v*2.0D0
      t225 = t224+2.0D0
      t226 = lamdat_r1*t8*t200*t202
      t227 = Dc*lamdat_r1*t8*t200*t202
      t228 = c1*c2*t4*t13*t34*t64*t67*t85*t92*t200*t202*(1.0D0/2.0D0)
      t229 = lamdat_r2*t8*t200*t202
      t230 = t13*t117*t200*t202*t220*(1.0D0/2.0D0)
      t231 = Dc*lamdat_r2*t8*t200*t202
      t232 = c1*c2*t4*t13*t34*t64*t67*t85*t132*t200*t202*(1.0D0/2.0D0)
      t233 = c1*c2*t2*t4*t8*t34*t35*t64*t67*t69*t86*t95*(3.0D0/4.0D0)
      t234 = c1*c2*t4*t8*t9*t34*t39*t64*t67*t69*t86*t95*(3.0D0/4.0D0)
      t235 = c1*c2*t4*t8*t10*t34*t43*t64*t67*t69*t86*t95*(3.0D0/4.0D0)
      t236 = c1*c2*kt*t2*t4*t8*t35*t55*t63*t67*t69*t93*t95*(3.0D0/2.0D0)
      t237 = c1*c2*kt*t4*t8*t9*t39*t55*t63*t67*t69*t93*t95*(3.0D0/2.0D0)
      t238 = c1*c2*kt*t4*t8*t10*t43*t55*t63*t67*t69*t93*t95*(3.0D0/2.0D0
     &)
      t239 = c1*c2*t4*t13*t34*t64*t67*t69*t86*t92*t95*t200*t209*(1.0D0/4
     &.0D0)
      t240 = c1*c2*t4*t13*t34*t64*t67*t69*t86*t95*t132*t200*t216*(1.0D0/
     &4.0D0)
      t241 = c1*c2*t4*t13*t34*t64*t67*t69*t86*t95*t152*t200*t217*(1.0D0/
     &4.0D0)
      t242 = c1*c2*kt*t4*t13*t55*t63*t67*t69*t92*t93*t95*t200*t209*(1.0D
     &0/2.0D0)
      t243 = c1*c2*kt*t4*t13*t55*t63*t67*t69*t93*t95*t132*t200*t216*(1.0
     &D0/2.0D0)
      t244 = c1*c2*kt*t4*t13*t55*t63*t67*t69*t93*t95*t152*t200*t217*(1.0
     &D0/2.0D0)
      t246 = lamdat_r1*t8*t200*v
      t248 = t13*t117*t200*t219*v*(1.0D0/2.0D0)
      t249 = Dc*lamdat_r1*t8*t200*v
      t252 = c1*c2*t4*t13*t34*t64*t67*t85*t92*t200*v*(1.0D0/2.0D0)
      t245 = t218+t222+t229+t230+t231+t232+t233+t234+t235+t236+t237+t238
     &+t239+t240+t241+t242+t243+t244-t246-t248-t249-t252-lamdat_r3*t8*t2
     &00*v-t13*t117*t200*t216*(1.0D0/2.0D0)-Dc*lamdat_r3*t8*t200*v-t13*t
     &117*t200*t221*v*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t64*t67*t85*t200*t2
     &16-c1*c2*t4*t13*t34*t64*t67*t85*t152*t200*v*(1.0D0/2.0D0)
      t247 = lamdat_r3*t8*t200*t202
      t250 = t13*t117*t200*t202*t221*(1.0D0/2.0D0)
      t251 = Dc*lamdat_r3*t8*t200*t202
      t253 = c1*c2*t4*t13*t34*t64*t67*t85*t152*t200*t202*(1.0D0/2.0D0)
      t254 = c1*c2*t2*t4*t8*t34*t35*t64*t67*t69*t86*t103*(3.0D0/4.0D0)
      t255 = c1*c2*t4*t8*t9*t34*t39*t64*t67*t69*t86*t103*(3.0D0/4.0D0)
      t256 = c1*c2*t4*t8*t10*t34*t43*t64*t67*t69*t86*t103*(3.0D0/4.0D0)
      t257 = c1*c2*kt*t2*t4*t8*t35*t55*t63*t67*t69*t93*t103*(3.0D0/2.0D0
     &)
      t258 = c1*c2*kt*t4*t8*t9*t39*t55*t63*t67*t69*t93*t103*(3.0D0/2.0D0
     &)
      t259 = c1*c2*kt*t4*t8*t10*t43*t55*t63*t67*t69*t93*t103*(3.0D0/2.0D
     &0)
      t260 = c1*c2*t4*t13*t34*t64*t67*t69*t86*t92*t103*t200*t209*(1.0D0/
     &4.0D0)
      t261 = c1*c2*t4*t13*t34*t64*t67*t69*t86*t103*t132*t200*t216*(1.0D0
     &/4.0D0)
      t262 = c1*c2*t4*t13*t34*t64*t67*t69*t86*t103*t152*t200*t217*(1.0D0
     &/4.0D0)
      t263 = c1*c2*kt*t4*t13*t55*t63*t67*t69*t92*t93*t103*t200*t209*(1.0
     &D0/2.0D0)
      t264 = c1*c2*kt*t4*t13*t55*t63*t67*t69*t93*t103*t132*t200*t216*(1.
     &0D0/2.0D0)
      t265 = c1*c2*kt*t4*t13*t55*t63*t67*t69*t93*t103*t152*t200*t217*(1.
     &0D0/2.0D0)
      t266 = lamdat_r4*t4*t8*(1.0D0/2.0D0)
      t267 = Dc*lamdat_r4*t4*t8*(1.0D0/2.0D0)
      t268 = ep4*2.0D0
      t269 = et4*2.0D0
      t270 = e4*(-2.0D0)+t268+t269
      t271 = t12*t13*t117*t225*t270*(1.0D0/8.0D0)
      t272 = c1*c2*t12*t13*t34*t35*t64*t67*t85*(3.0D0/4.0D0)
      t273 = t266+t267+t271+t272-c1*c2*t2*t4*t8*t34*t64*t67*t85*(3.0D0/2
     &.0D0)-c1*c2*t12*t13*t34*t35*t37*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*c
     &2*kt*t12*t13*t35*t37*t55*t63*t67*t69*t93*9.0D0-c1*c2*t2*t9*t12*t13
     &*t34*t39*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*c2*t2*t10*t12*t13*t34*t4
     &3*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*c2*kt*t2*t9*t12*t13*t39*t55*t63
     &*t67*t69*t93*9.0D0-c1*c2*kt*t2*t10*t12*t13*t43*t55*t63*t67*t69*t93
     &*9.0D0-c1*c2*t2*t12*t34*t64*t67*t69*t86*t92*t118*t200*t209*(3.0D0/
     &2.0D0)-c1*c2*t2*t12*t34*t64*t67*t69*t86*t118*t132*t200*t216*(3.0D0
     &/2.0D0)-c1*c2*t2*t12*t34*t64*t67*t69*t86*t118*t152*t200*t217*(3.0D
     &0/2.0D0)-c1*c2*kt*t2*t12*t55*t63*t67*t69*t92*t93*t118*t200*t209*3.
     &0D0-c1*c2*kt*t2*t12*t55*t63*t67*t69*t93*t118*t132*t200*t216*3.0D0-
     &c1*c2*kt*t2*t12*t55*t63*t67*t69*t93*t118*t152*t200*t217*3.0D0
      t274 = lamdat_r5*t4*t8*(1.0D0/2.0D0)
      t275 = Dc*lamdat_r5*t4*t8*(1.0D0/2.0D0)
      t276 = ep5*2.0D0
      t277 = et5*2.0D0
      t278 = e5*(-2.0D0)+t276+t277
      t279 = t12*t13*t117*t225*t278*(1.0D0/8.0D0)
      t280 = c1*c2*t12*t13*t34*t39*t64*t67*t85*(3.0D0/4.0D0)
      t281 = t274+t275+t279+t280-c1*c2*t4*t8*t9*t34*t64*t67*t85*(3.0D0/2
     &.0D0)-c1*c2*t12*t13*t34*t39*t41*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*c
     &2*kt*t12*t13*t39*t41*t55*t63*t67*t69*t93*9.0D0-c1*c2*t2*t9*t12*t13
     &*t34*t35*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*c2*t9*t10*t12*t13*t34*t4
     &3*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*c2*kt*t2*t9*t12*t13*t35*t55*t63
     &*t67*t69*t93*9.0D0-c1*c2*kt*t9*t10*t12*t13*t43*t55*t63*t67*t69*t93
     &*9.0D0-c1*c2*t9*t12*t34*t64*t67*t69*t86*t92*t118*t200*t209*(3.0D0/
     &2.0D0)-c1*c2*t9*t12*t34*t64*t67*t69*t86*t118*t132*t200*t216*(3.0D0
     &/2.0D0)-c1*c2*t9*t12*t34*t64*t67*t69*t86*t118*t152*t200*t217*(3.0D
     &0/2.0D0)-c1*c2*kt*t9*t12*t55*t63*t67*t69*t92*t93*t118*t200*t209*3.
     &0D0-c1*c2*kt*t9*t12*t55*t63*t67*t69*t93*t118*t132*t200*t216*3.0D0-
     &c1*c2*kt*t9*t12*t55*t63*t67*t69*t93*t118*t152*t200*t217*3.0D0
      t282 = lamdat_r6*t4*t8*(1.0D0/2.0D0)
      t283 = Dc*lamdat_r6*t4*t8*(1.0D0/2.0D0)
      t284 = ep6*2.0D0
      t285 = et6*2.0D0
      t286 = e6*(-2.0D0)+t284+t285
      t287 = t12*t13*t117*t225*t286*(1.0D0/8.0D0)
      t288 = c1*c2*t12*t13*t34*t43*t64*t67*t85*(3.0D0/4.0D0)
      t289 = t282+t283+t287+t288-c1*c2*t4*t8*t10*t34*t64*t67*t85*(3.0D0/
     &2.0D0)-c1*c2*t12*t13*t34*t43*t45*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*
     &c2*kt*t12*t13*t43*t45*t55*t63*t67*t69*t93*9.0D0-c1*c2*t2*t10*t12*t
     &13*t34*t35*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*c2*t9*t10*t12*t13*t34*
     &t39*t64*t67*t69*t86*(9.0D0/2.0D0)-c1*c2*kt*t2*t10*t12*t13*t35*t55*
     &t63*t67*t69*t93*9.0D0-c1*c2*kt*t9*t10*t12*t13*t39*t55*t63*t67*t69*
     &t93*9.0D0-c1*c2*t10*t12*t34*t64*t67*t69*t86*t92*t118*t200*t209*(3.
     &0D0/2.0D0)-c1*c2*t10*t12*t34*t64*t67*t69*t86*t118*t132*t200*t216*(
     &3.0D0/2.0D0)-c1*c2*t10*t12*t34*t64*t67*t69*t86*t118*t152*t200*t217
     &*(3.0D0/2.0D0)-c1*c2*kt*t10*t12*t55*t63*t67*t69*t92*t93*t118*t200*
     &t209*3.0D0-c1*c2*kt*t10*t12*t55*t63*t67*t69*t93*t118*t132*t200*t21
     &6*3.0D0-c1*c2*kt*t10*t12*t55*t63*t67*t69*t93*t118*t152*t200*t217*3
     &.0D0
      t290 = t13*t117*t200*t221*v
      t291 = t13*t117*t200*t219*v
      t292 = t13*t117*t200*t220*v
      t293 = Tn-too
      t294 = alpha*t293
      t295 = -e1+ep1+et1+t294
      t296 = -e3+ep3+et3+t294
      t297 = t13*t117*t200*t296*v
      t298 = -e2+ep2+et2+t294
      t299 = t13*t117*t200*t295*v
      t300 = t13*t117*t200*t298*v
      t301 = t117**2
      t302 = t35**2
      t303 = t39**2
      t304 = t43**2
      t305 = Sm*bs1*2.0D0
      t306 = Sm*bs1*v*2.0D0
      t307 = Sm*bs2*2.0D0
      t308 = Sm*bs2*v*2.0D0
      t309 = t14+t15+t22+t24-t87-t88-t120-t122-t125-t126+t307+t308
      t310 = Sm*bs3*2.0D0
      t311 = Sm*bs3*v*2.0D0
      t312 = t14+t15+t21+t23-t96-t97-t119-t121-t125-t126+t310+t311
      A0(1,1) = -1.0D0
      A0(1,13) = lamdat_r1
      A0(2,2) = -1.0D0
      A0(2,13) = lamdat_r2
      A0(3,3) = -1.0D0
      A0(3,13) = lamdat_r3
      A0(4,4) = -1.0D0
      A0(4,13) = lamdat_r4
      A0(5,5) = -1.0D0
      A0(5,13) = lamdat_r5
      A0(6,6) = -1.0D0
      A0(6,13) = lamdat_r6
      A0(7,1) = t104+c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*(t91+t12*t13*(t
     &56-t57-t58-t59-t60-t80+t81+t82+t83+t84+Sa*bs1*2.0D0-Sa*bs3*2.0D0+S
     &a*bs1*v*2.0D0-Sa*bs3*v*2.0D0+Sm*bs1*x*2.0D0-Sm*bs3*x*2.0D0+Sm*bs1*
     &v*x*2.0D0-Sm*bs3*v*x*2.0D0))*(e2+e3-ep2-ep3-et2-et3+t21+t22+t23+t2
     &4-t25-t26+t27+t28-t29-t30+t31+t32-t56+t57+t58+t59+t60-Sa*bs1*2.0D0
     &-Sa*bs1*v*2.0D0-Sm*bs1*x*2.0D0-Sm*bs1*v*x*2.0D0)*(1.0D0/4.0D0)+c1*
     &c2*kt*t4*t8*t55*t63*t67*t68*t69*t92*t93*(t91+t12*t13*(t56-t57-t58-
     &t59-t60+t73+t74+t76+t78-t80+t81+t82+t83+t84-Sa*bs3*2.0D0-Sa*bs3*v*
     &2.0D0-Sm*bs3*x*2.0D0-Sm*bs3*v*x*2.0D0))*(1.0D0/2.0D0)
      A0(7,2) = t106
      A0(7,3) = c1*c2*t4*t8*t34*t64*t67*t68*t85*(-1.0D0/2.0D0)-c1*c2*t4*
     &t8*t34*t64*t67*t68*t69*t86*t92*t103*(1.0D0/4.0D0)-c1*c2*kt*t4*t8*t
     &55*t63*t67*t68*t69*t92*t93*t103*(1.0D0/2.0D0)
      A0(7,4) = t110
      A0(7,5) = t113
      A0(7,6) = t116
      A0(7,7) = t104+c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t92*t105*(1.0D0
     &/4.0D0)+c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t92*t93*t105*(1.0D0/2.0
     &D0)-1.0D0
      A0(7,8) = t106
      A0(7,9) = -t107-c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t92*t103*(1.0D
     &0/4.0D0)-c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t92*t93*t103*(1.0D0/2.
     &0D0)
      A0(7,10) = t110
      A0(7,11) = t113
      A0(7,12) = t116
      A0(7,13) = c1*c2*t4*t8*t34*t64*t67*t85*t92*(1.0D0/2.0D0)-c1*c2*t4*
     &t8*t34*t64*t67*t68*t85*(t21+t22+t23+t24-t73-t74+t305+t306-Sm*bs2-S
     &m*bs3-Sm*bs2*v-Sm*bs3*v)*(1.0D0/2.0D0)+c1*c2*t4*t13*t34*t64*t67*t6
     &8*t85*t92*t117*(1.0D0/2.0D0)+c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t
     &92*(t124+t128+t129+t130+t131-t12*t47*t117*t118*2.0D0-t12*t49*t117*
     &t118*2.0D0-t12*t51*t117*t118*2.0D0+t12*t13*t20*(t14+t15-t22-t24+t1
     &20+t122-Sm*bs1-Sm*bs1*v)*2.0D0)*(1.0D0/4.0D0)+c1*c2*kt*t4*t8*t55*t
     &63*t67*t68*t69*t92*t93*t148*(1.0D0/2.0D0)
      A0(8,1) = t135
      A0(8,2) = t104-c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t95*t132*(1.0D0
     &/4.0D0)-c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93*t95*t132*(1.0D0/2.0
     &D0)
      A0(8,3) = t136
      A0(8,4) = t139
      A0(8,5) = t142
      A0(8,6) = t145
      A0(8,7) = t135
      A0(8,8) = t104-c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t95*t132*(1.0D0
     &/4.0D0)-c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93*t95*t132*(1.0D0/2.0
     &D0)-1.0D0
      A0(8,9) = t136
      A0(8,10) = t139
      A0(8,11) = t142
      A0(8,12) = t145
      A0(8,13) = c1*c2*t4*t8*t34*t64*t67*t85*t132*(1.0D0/2.0D0)-c1*c2*t4
     &*t8*t34*t64*t67*t68*t85*t309*(1.0D0/2.0D0)+c1*c2*t4*t13*t34*t64*t6
     &7*t68*t85*t117*t132*(1.0D0/2.0D0)+c1*c2*t4*t8*t34*t64*t67*t68*t69*
     &t86*t132*t148*(1.0D0/4.0D0)+c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93
     &*t132*t148*(1.0D0/2.0D0)
      A0(9,1) = t155
      A0(9,2) = t156
      A0(9,3) = t104-c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t103*t152*(1.0D
     &0/4.0D0)-c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93*t103*t152*(1.0D0/2
     &.0D0)
      A0(9,4) = t159
      A0(9,5) = t162
      A0(9,6) = t165
      A0(9,7) = t155
      A0(9,8) = t156
      A0(9,9) = t104-c1*c2*t4*t8*t34*t64*t67*t68*t69*t86*t103*t152*(1.0D
     &0/4.0D0)-c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93*t103*t152*(1.0D0/2
     &.0D0)-1.0D0
      A0(9,10) = t159
      A0(9,11) = t162
      A0(9,12) = t165
      A0(9,13) = c1*c2*t4*t8*t34*t64*t67*t85*t152*(1.0D0/2.0D0)-c1*c2*t4
     &*t8*t34*t64*t67*t68*t85*t312*(1.0D0/2.0D0)+c1*c2*t4*t13*t34*t64*t6
     &7*t68*t85*t117*t152*(1.0D0/2.0D0)+c1*c2*t4*t8*t34*t64*t67*t68*t69*
     &t86*t148*t152*(1.0D0/4.0D0)+c1*c2*kt*t4*t8*t55*t63*t67*t68*t69*t93
     &*t148*t152*(1.0D0/2.0D0)
      A0(10,1) = t166
      A0(10,2) = t169
      A0(10,3) = t172
      A0(10,4) = t173-c1*c2*t4*t8*t34*t37*t64*t67*t68*t69*t86*9.0D0-c1*c
     &2*kt*t4*t8*t37*t55*t63*t67*t68*t69*t93*1.8D1
      A0(10,5) = t174
      A0(10,6) = t175
      A0(10,7) = t166
      A0(10,8) = t169
      A0(10,9) = t172
      A0(10,10) = t173-c1*c2*t4*t8*t34*t37*t64*t67*t68*t69*t86*9.0D0-c1*
     &c2*kt*t4*t8*t37*t55*t63*t67*t68*t69*t93*1.8D1-1.0D0
      A0(10,11) = t174
      A0(10,12) = t175
      A0(10,13) = c1*c2*t2*t34*t64*t67*t85*(-3.0D0)-c1*c2*t2*t34*t64*t67
     &*t68*t69*t86*t148*(3.0D0/2.0D0)-c1*c2*kt*t2*t55*t63*t67*t68*t69*t9
     &3*t148*3.0D0+c1*c2*t4*t13*t34*t35*t64*t67*t68*t85*t117*(3.0D0/2.0D
     &0)
      A0(11,1) = t178
      A0(11,2) = t181
      A0(11,3) = t184
      A0(11,4) = t174
      A0(11,5) = t173-c1*c2*t4*t8*t34*t41*t64*t67*t68*t69*t86*9.0D0-c1*c
     &2*kt*t4*t8*t41*t55*t63*t67*t68*t69*t93*1.8D1
      A0(11,6) = t185
      A0(11,7) = t178
      A0(11,8) = t181
      A0(11,9) = t184
      A0(11,10) = t174
      A0(11,11) = t173-c1*c2*t4*t8*t34*t41*t64*t67*t68*t69*t86*9.0D0-c1*
     &c2*kt*t4*t8*t41*t55*t63*t67*t68*t69*t93*1.8D1-1.0D0
      A0(11,12) = t185
      A0(11,13) = c1*c2*t9*t34*t64*t67*t85*(-3.0D0)-c1*c2*t9*t34*t64*t67
     &*t68*t69*t86*t148*(3.0D0/2.0D0)-c1*c2*kt*t9*t55*t63*t67*t68*t69*t9
     &3*t148*3.0D0+c1*c2*t4*t13*t34*t39*t64*t67*t68*t85*t117*(3.0D0/2.0D
     &0)
      A0(12,1) = t190
      A0(12,2) = t193
      A0(12,3) = t196
      A0(12,4) = t175
      A0(12,5) = t185
      A0(12,6) = t173-c1*c2*t4*t8*t34*t45*t64*t67*t68*t69*t86*9.0D0-c1*c
     &2*kt*t4*t8*t45*t55*t63*t67*t68*t69*t93*1.8D1
      A0(12,7) = t190
      A0(12,8) = t193
      A0(12,9) = t196
      A0(12,10) = t175
      A0(12,11) = t185
      A0(12,12) = t173-c1*c2*t4*t8*t34*t45*t64*t67*t68*t69*t86*9.0D0-c1*
     &c2*kt*t4*t8*t45*t55*t63*t67*t68*t69*t93*1.8D1-1.0D0
      A0(12,13) = c1*c2*t10*t34*t64*t67*t85*(-3.0D0)-c1*c2*t10*t34*t64*t
     &67*t68*t69*t86*t148*(3.0D0/2.0D0)-c1*c2*kt*t10*t55*t63*t67*t68*t69
     &*t93*t148*3.0D0+c1*c2*t4*t13*t34*t43*t64*t67*t68*t85*t117*(3.0D0/2
     &.0D0)
      A0(13,1) = t218+t223+t226+t227+t228-lamdat_r2*t8*t200*v-lamdat_r3*
     &t8*t200*v-t13*t117*t200*t209*(1.0D0/2.0D0)-Dc*lamdat_r2*t8*t200*v-
     &Dc*lamdat_r3*t8*t200*v+t13*t117*t200*t202*(-e1+ep1+et1+t201-alpha*
     &too)*(1.0D0/2.0D0)-t13*t117*t200*v*(-e2+ep2+et2+t201-alpha*too)*(1
     &.0D0/2.0D0)-t13*t117*t200*v*(-e3+ep3+et3+t201-alpha*too)*(1.0D0/2.
     &0D0)-c1*c2*t4*t13*t34*t64*t67*t85*t200*t209-c1*c2*t4*t13*t34*t64*t
     &67*t85*t132*t200*v*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t64*t67*t85*t152
     &*t200*v*(1.0D0/2.0D0)-c1*c2*t2*t4*t8*t34*t35*t64*t67*t69*t86*t105*
     &(3.0D0/4.0D0)-c1*c2*t4*t8*t9*t34*t39*t64*t67*t69*t86*t105*(3.0D0/4
     &.0D0)-c1*c2*t4*t8*t10*t34*t43*t64*t67*t69*t86*t105*(3.0D0/4.0D0)-c
     &1*c2*kt*t2*t4*t8*t35*t55*t63*t67*t69*t93*t105*(3.0D0/2.0D0)-c1*c2*
     &kt*t4*t8*t9*t39*t55*t63*t67*t69*t93*t105*(3.0D0/2.0D0)-c1*c2*kt*t4
     &*t8*t10*t43*t55*t63*t67*t69*t93*t105*(3.0D0/2.0D0)-c1*c2*t4*t13*t3
     &4*t64*t67*t69*t86*t92*t105*t200*t209*(1.0D0/4.0D0)-c1*c2*t4*t13*t3
     &4*t64*t67*t69*t86*t105*t132*t200*t216*(1.0D0/4.0D0)-c1*c2*t4*t13*t
     &34*t64*t67*t69*t86*t105*t152*t200*t217*(1.0D0/4.0D0)-c1*c2*kt*t4*t
     &13*t55*t63*t67*t69*t92*t93*t105*t200*t209*(1.0D0/2.0D0)-c1*c2*kt*t
     &4*t13*t55*t63*t67*t69*t93*t105*t132*t200*t216*(1.0D0/2.0D0)-c1*c2*
     &kt*t4*t13*t55*t63*t67*t69*t93*t105*t152*t200*t217*(1.0D0/2.0D0)
      A0(13,2) = t245
      A0(13,3) = t222+t223+t247+t250+t251+t253+t254+t255+t256+t257+t258+
     &t259+t260+t261+t262+t263+t264+t265-lamdat_r1*t8*t200*v-lamdat_r2*t
     &8*t200*v-t13*t117*t200*t217*(1.0D0/2.0D0)-Dc*lamdat_r1*t8*t200*v-D
     &c*lamdat_r2*t8*t200*v-t13*t117*t200*t219*v*(1.0D0/2.0D0)-t13*t117*
     &t200*t220*v*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t64*t67*t85*t200*t217-c
     &1*c2*t4*t13*t34*t64*t67*t85*t92*t200*v*(1.0D0/2.0D0)-c1*c2*t4*t13*
     &t34*t64*t67*t85*t132*t200*v*(1.0D0/2.0D0)
      A0(13,4) = t273
      A0(13,5) = t281
      A0(13,6) = t289
      A0(13,7) = t218+t223+t226+t227+t228-lamdat_r2*t8*t200*v-lamdat_r3*
     &t8*t200*v-t13*t117*t200*t209*(1.0D0/2.0D0)-Dc*lamdat_r2*t8*t200*v-
     &Dc*lamdat_r3*t8*t200*v+t13*t117*t200*t202*t219*(1.0D0/2.0D0)-t13*t
     &117*t200*t220*v*(1.0D0/2.0D0)-t13*t117*t200*t221*v*(1.0D0/2.0D0)-c
     &1*c2*t4*t13*t34*t64*t67*t85*t200*t209-c1*c2*t4*t13*t34*t64*t67*t85
     &*t132*t200*v*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t64*t67*t85*t152*t200*
     &v*(1.0D0/2.0D0)-c1*c2*t2*t4*t8*t34*t35*t64*t67*t69*t86*t105*(3.0D0
     &/4.0D0)-c1*c2*t4*t8*t9*t34*t39*t64*t67*t69*t86*t105*(3.0D0/4.0D0)-
     &c1*c2*t4*t8*t10*t34*t43*t64*t67*t69*t86*t105*(3.0D0/4.0D0)-c1*c2*k
     &t*t2*t4*t8*t35*t55*t63*t67*t69*t93*t105*(3.0D0/2.0D0)-c1*c2*kt*t4*
     &t8*t9*t39*t55*t63*t67*t69*t93*t105*(3.0D0/2.0D0)-c1*c2*kt*t4*t8*t1
     &0*t43*t55*t63*t67*t69*t93*t105*(3.0D0/2.0D0)-c1*c2*t4*t13*t34*t64*
     &t67*t69*t86*t92*t105*t200*t209*(1.0D0/4.0D0)-c1*c2*t4*t13*t34*t64*
     &t67*t69*t86*t105*t132*t200*t216*(1.0D0/4.0D0)-c1*c2*t4*t13*t34*t64
     &*t67*t69*t86*t105*t152*t200*t217*(1.0D0/4.0D0)-c1*c2*kt*t4*t13*t55
     &*t63*t67*t69*t92*t93*t105*t200*t209*(1.0D0/2.0D0)-c1*c2*kt*t4*t13*
     &t55*t63*t67*t69*t93*t105*t132*t200*t216*(1.0D0/2.0D0)-c1*c2*kt*t4*
     &t13*t55*t63*t67*t69*t93*t105*t152*t200*t217*(1.0D0/2.0D0)
      A0(13,8) = t245
      A0(13,9) = t222+t223-t246+t247-t248-t249+t250+t251-t252+t253+t254+
     &t255+t256+t257+t258+t259+t260+t261+t262+t263+t264+t265-lamdat_r2*t
     &8*t200*v-t13*t117*t200*t217*(1.0D0/2.0D0)-Dc*lamdat_r2*t8*t200*v-t
     &13*t117*t200*t220*v*(1.0D0/2.0D0)-c1*c2*t4*t13*t34*t64*t67*t85*t20
     &0*t217-c1*c2*t4*t13*t34*t64*t67*t85*t132*t200*v*(1.0D0/2.0D0)
      A0(13,10) = t273
      A0(13,11) = t281
      A0(13,12) = t289
      A0(13,13) = -lamdat_r1*(t290+t292-t13*t117*t200*t202*t219)-lamdat_
     &r2*(t290+t291-t13*t117*t200*t202*t220)-lamdat_r3*(t291+t292-t13*t1
     &17*t200*t202*t221)+a2*(n4*(-x+1.0D0)**(n4-1.0D0)+n3*x**(n3-1.0D0))
     &*(1.0D0/2.0D0)-Dc*lamdat_r1*(t297+t300-t13*t117*t200*t202*t295)-Dc
     &*lamdat_r2*(t297+t299-t13*t117*t200*t202*t298)-Dc*lamdat_r3*(t299+
     &t300-t13*t117*t200*t202*t296)+lamdat_r4*t4*t13*t35*t117*(1.0D0/2.0
     &D0)+lamdat_r5*t4*t13*t39*t117*(1.0D0/2.0D0)+lamdat_r6*t4*t13*t43*t
     &117*(1.0D0/2.0D0)+t12*t118*t225*t301*t302*(1.0D0/4.0D0)+t12*t118*t
     &225*t301*t303*(1.0D0/4.0D0)+t12*t118*t225*t301*t304*(1.0D0/4.0D0)-
     &t118*t200*t209*t219*t301-t118*t200*t216*t220*t301-t118*t200*t217*t
     &221*t301+Dc*lamdat_r4*t4*t13*t35*t117*(1.0D0/2.0D0)+Dc*lamdat_r5*t
     &4*t13*t39*t117*(1.0D0/2.0D0)+Dc*lamdat_r6*t4*t13*t43*t117*(1.0D0/2
     &.0D0)+c1*c2*t12*t34*t64*t67*t85*t117*t118*t302*(3.0D0/4.0D0)+c1*c2
     &*t12*t34*t64*t67*t85*t117*t118*t303*(3.0D0/4.0D0)+c1*c2*t12*t34*t6
     &4*t67*t85*t117*t118*t304*(3.0D0/4.0D0)-c1*c2*t2*t4*t13*t34*t35*t64
     &*t67*t85*t117*(3.0D0/2.0D0)-c1*c2*t4*t9*t13*t34*t39*t64*t67*t85*t1
     &17*(3.0D0/2.0D0)-c1*c2*t4*t10*t13*t34*t43*t64*t67*t85*t117*(3.0D0/
     &2.0D0)+c1*c2*t4*t13*t34*t64*t67*t85*t200*t216*t309*(1.0D0/2.0D0)+c
     &1*c2*t4*t13*t34*t64*t67*t85*t200*t217*t312*(1.0D0/2.0D0)+c1*c2*t4*
     &t13*t34*t64*t67*t85*t200*t209*(t21+t22+t23+t24-t73-t74-t119-t120-t
     &121-t122+t305+t306)*(1.0D0/2.0D0)-c1*c2*t2*t4*t8*t34*t35*t64*t67*t
     &69*t86*t148*(3.0D0/4.0D0)-c1*c2*t4*t8*t9*t34*t39*t64*t67*t69*t86*t
     &148*(3.0D0/4.0D0)-c1*c2*t4*t8*t10*t34*t43*t64*t67*t69*t86*t148*(3.
     &0D0/4.0D0)-c1*c2*t4*t34*t64*t67*t85*t92*t117*t118*t200*t209-c1*c2*
     &t4*t34*t64*t67*t85*t117*t118*t132*t200*t216-c1*c2*t4*t34*t64*t67*t
     &85*t117*t118*t152*t200*t217-c1*c2*kt*t2*t4*t8*t35*t55*t63*t67*t69*
     &t93*t148*(3.0D0/2.0D0)-c1*c2*kt*t4*t8*t9*t39*t55*t63*t67*t69*t93*t
     &148*(3.0D0/2.0D0)-c1*c2*kt*t4*t8*t10*t43*t55*t63*t67*t69*t93*t148*
     &(3.0D0/2.0D0)-c1*c2*t4*t13*t34*t64*t67*t69*t86*t92*t148*t200*t209*
     &(1.0D0/4.0D0)-c1*c2*t4*t13*t34*t64*t67*t69*t86*t132*t148*t200*t216
     &*(1.0D0/4.0D0)-c1*c2*t4*t13*t34*t64*t67*t69*t86*t148*t152*t200*t21
     &7*(1.0D0/4.0D0)-c1*c2*kt*t4*t13*t55*t63*t67*t69*t92*t93*t148*t200*
     &t209*(1.0D0/2.0D0)-c1*c2*kt*t4*t13*t55*t63*t67*t69*t93*t132*t148*t
     &200*t216*(1.0D0/2.0D0)-c1*c2*kt*t4*t13*t55*t63*t67*t69*t93*t148*t1
     &52*t200*t217*(1.0D0/2.0D0)




      NR_JAC=A0


      ENDIF
!DEC$ FREEFORM

end subroutine

subroutine N_R_Residual(PARAM,NPARAM,VAR,NVAR,Tn,index1,NR_RE)

implicit real*8 (t)



!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
integer NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


! NR Variables
real*8 NR_RE(13,1),A0(13,1)


real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6
real*8 eto1,eto2,eto3,eto4,eto5,eto6,epo1,epo2,epo3,epo4,epo5,epo6
real*8 lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 phi, Tn

integer index1


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

   phi=0
   
   if(too<=0)then
   too=1e-10
   endif

   
   e1=e(1)
   e2=e(2)
   e3=e(3)
   e4=e(4)
   e5=e(5)
   e6=e(6)
    if(maxval(abs(e))<=1.0e-8)then
       e1=1.0e-8_8
    end if

	
   et1=et(1)
   et2=et(2)
   et3=et(3)
   et4=et(4)
   et5=et(5)
   et6=et(6)
   if(maxval(abs(et))<=1.0e-9)then
       et1=1.0e-9_8
    end if

   et_tr1=et_tr(1)
   et_tr2=et_tr(2)
   et_tr3=et_tr(3)
   et_tr4=et_tr(4)
   et_tr5=et_tr(5)
   et_tr6=et_tr(6)
   if(maxval(abs(et_tr))<=1.0e-10)then
       et_tr1=1.0e-10_8
    end if
	
   eto1=eto(1)
   eto2=eto(2)
   eto3=eto(3)
   eto4=eto(4)
   eto5=eto(5)
   eto6=eto(6)
    if(maxval(abs(eto))<=1.0e-11)then
       eto1=1.0e-11_8
    end if
	
   ep1=ep(1)
   ep2=ep(2)
   ep3=ep(3)
   ep4=ep(4)
   ep5=ep(5)
   ep6=ep(6)


   epo1=epo(1)
   epo2=epo(2)
   epo3=epo(3)
   epo4=epo(4)
   epo5=epo(5)
   epo6=epo(6)

   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)  
   
   
   bs1=backstress(1)
   bs2=backstress(2)
   bs3=backstress(3)
   bs4=backstress(4)
   bs5=backstress(5)
   bs6=backstress(6)
   

   if(abs(x-xo)<=0)then
       x=xo+1.0e-16_8
   end if

!DEC$ NOFREEFORM   


      if(index1==1)then    !Forward
      t3 = v+1.0D0
      t4 = 1.0D0/t3
      t5 = Sa*x
      t6 = Sm*x
      t7 = Sa-t5+t6
      t8 = 1.0D0/t7
      t35 = -e4+ep4+et4
      t36 = t4*t8*t35*(1.0D0/2.0D0)
      t2 = bs4-t36
      t39 = -e5+ep5+et5
      t40 = t4*t8*t39*(1.0D0/2.0D0)
      t9 = bs5-t40
      t43 = -e6+ep6+et6
      t44 = t4*t8*t43*(1.0D0/2.0D0)
      t10 = bs6-t44
      t14 = Sa*bs1
      t15 = Sa*bs1*v
      t16 = Sa*bs1*x
      t17 = Sm*bs1*x
      t18 = Sa*bs1*v*x
      t19 = Sm*bs1*v*x
      t21 = Sa*bs2
      t23 = Sa*bs2*v
      t25 = Sa*bs2*x
      t27 = Sm*bs2*x
      t29 = Sa*bs2*v*x
      t31 = Sm*bs2*v*x
      t11 = e1-e2-ep1+ep2-et1+et2+t14+t15-t16+t17-t18+t19-t21-t23+t25-t2
     &7+t29-t31
      t12 = 1.0D0/t3**2
      t13 = 1.0D0/t7**2
      t22 = Sa*bs3
      t24 = Sa*bs3*v
      t26 = Sa*bs3*x
      t28 = Sm*bs3*x
      t30 = Sa*bs3*v*x
      t32 = Sm*bs3*v*x
      t20 = e1-e3-ep1+ep3-et1+et3+t14+t15-t16+t17-t18+t19-t22-t24+t26-t2
     &8+t30-t32
      t33 = e2-e3-ep2+ep3-et2+et3+t21-t22+t23-t24-t25+t26+t27-t28-t29+t3
     &0+t31-t32
      t34 = sqrt(2.0D0)
      t37 = t2**2
      t38 = t37*6.0D0
      t41 = t9**2
      t42 = t41*6.0D0
      t45 = t10**2
      t46 = t45*6.0D0
      t47 = t11**2
      t48 = t12*t13*t47
      t49 = t20**2
      t50 = t12*t13*t49
      t51 = t33**2
      t52 = t12*t13*t51
      t53 = t38+t42+t46+t48+t50+t52
      t54 = abs(t53)
      t55 = 1.0D0/sqrt(t54)
      t56 = sqrt(t54)
      t60 = kt*t34*t56*(1.0D0/2.0D0)
      t57 = exp(-t60)
      t58 = t57-1.0D0
      t59 = x-xo
      t61 = ep1*2.0D0
      t62 = et1*2.0D0
      t63 = Sa*bs1*x*2.0D0
      t64 = Sa*bs1*v*x*2.0D0
      t108 = e1*2.0D0
      t109 = Sa*bs1*2.0D0
      t110 = Sa*bs1*v*2.0D0
      t111 = Sm*bs1*x*2.0D0
      t112 = Sm*bs1*v*x*2.0D0
      t65 = e2+e3-ep2-ep3-et2-et3+t21+t22+t23+t24-t25-t26+t27+t28-t29-t3
     &0+t31+t32+t61+t62+t63+t64-t108-t109-t110-t111-t112
      t66 = t58**2
      t67 = c2*ztd
      t68 = t67+1.0D0
      t69 = 1.0D0/t68
      t70 = ep2*2.0D0
      t71 = et2*2.0D0
      t72 = Sa*bs2*x*2.0D0
      t73 = Sa*bs2*v*x*2.0D0
      t113 = e2*2.0D0
      t114 = Sa*bs2*2.0D0
      t115 = Sa*bs2*v*2.0D0
      t116 = Sm*bs2*x*2.0D0
      t117 = Sm*bs2*v*x*2.0D0
      t74 = e1+e3-ep1-ep3-et1-et3+t14+t15-t16+t17-t18+t19+t22+t24-t26+t2
     &8-t30+t32+t70+t71+t72+t73-t113-t114-t115-t116-t117
      t75 = ep3*2.0D0
      t76 = et3*2.0D0
      t77 = Sa*bs3*x*2.0D0
      t78 = Sa*bs3*v*x*2.0D0
      t118 = e3*2.0D0
      t119 = Sa*bs3*2.0D0
      t120 = Sa*bs3*v*2.0D0
      t121 = Sm*bs3*x*2.0D0
      t122 = Sm*bs3*v*x*2.0D0
      t79 = e1+e2-ep1-ep2-et1-et2+t14+t15-t16+t17-t18+t19+t21+t23-t25+t2
     &7-t29+t31+t75+t76+t77+t78-t118-t119-t120-t121-t122
      t80 = v**2
      t81 = t80*2.0D0
      t82 = t81+v-1.0D0
      t83 = 1.0D0/t82
      t84 = Tn-too
      t85 = alpha*t84
      t86 = v-1.0D0
      t87 = -e1+ep1+et1+t85
      t88 = -e2+ep2+et2+t85
      t89 = -e3+ep3+et3+t85
      t90 = Sa-Sm
      t91 = t8*t83*t89*v
      t92 = t8*t83*t87*v
      t93 = t8*t83*t88*v
      t94 = t8*t83*t89*v*(1.0D0/2.0D0)
      t104 = t8*t83*t86*t88
      t95 = t91+t92-t104
      t102 = t8*t83*t86*t87
      t96 = t91+t93-t102
      t101 = t8*t83*t86*t89
      t97 = t92+t93-t101
      t98 = t90*t97*v
      t99 = t8*t83*t87*v*(1.0D0/2.0D0)
      t100 = t8*t83*t88*v*(1.0D0/2.0D0)
      t103 = t90*t96*v
      t105 = t90*t95*v
      t106 = v*2.0D0
      t107 = t106+2.0D0
      A0(1,1) = -et1+eto1+Hmax*t4*t8*t34*t55*t58*t59*t65*(1.0D0/2.0D0)
      A0(2,1) = -et2+eto2+Hmax*t4*t8*t34*t55*t58*t59*t74*(1.0D0/2.0D0)
      A0(3,1) = -et3+eto3+Hmax*t4*t8*t34*t55*t58*t59*t79*(1.0D0/2.0D0)
      A0(4,1) = -et4+eto4-Hmax*t2*t34*t55*t58*t59*3.0D0
      A0(5,1) = -et5+eto5-Hmax*t9*t34*t55*t58*t59*3.0D0
      A0(6,1) = -et6+eto6-Hmax*t10*t34*t55*t58*t59*3.0D0
      A0(7,1) = -ep1+epo1-c1*c2*t4*t8*t34*t55*t59*t65*t66*t69*(1.0D0/2.0
     &D0)
      A0(8,1) = -ep2+epo2-c1*c2*t4*t8*t34*t55*t59*t66*t69*t74*(1.0D0/2.0
     &D0)
      A0(9,1) = -ep3+epo3-c1*c2*t4*t8*t34*t55*t59*t66*t69*t79*(1.0D0/2.0
     &D0)
      A0(10,1) = -ep4+epo4+c1*c2*t2*t34*t55*t59*t66*t69*3.0D0
      A0(11,1) = -ep5+epo5+c1*c2*t9*t34*t55*t59*t66*t69*3.0D0
      A0(12,1) = -ep6+epo6+c1*c2*t10*t34*t55*t59*t66*t69*3.0D0
      A0(13,1) = -Yo-a3-rduo+Tn*rdso+(t94+t99-t8*t83*t86*t88*(1.0D0/2.0D
     &0))*(t98+t103-t90*t95)+(t94+t100-t8*t83*t86*t87*(1.0D0/2.0D0))*(t9
     &8+t105-t90*t96)+(t99+t100-t8*t83*t86*t89*(1.0D0/2.0D0))*(t103+t105
     &-t90*t97)-a1*(-(-x+1.0D0)**n2+x**n1+1.0D0)*(1.0D0/2.0D0)-t12*t13*t
     &35**2*t90*t107*(1.0D0/8.0D0)-t12*t13*t39**2*t90*t107*(1.0D0/8.0D0)
     &-t12*t13*t43**2*t90*t107*(1.0D0/8.0D0)-Hmax*t34*t37*t55*t58*3.0D0-
     &Hmax*t34*t41*t55*t58*3.0D0-Hmax*t34*t45*t55*t58*3.0D0+Hmax*t4*t8*t
     &34*t55*t58*t65*(bs1+t91+t93-t102)*(1.0D0/2.0D0)+Hmax*t4*t8*t34*t55
     &*t58*t74*(bs2+t91+t92-t104)*(1.0D0/2.0D0)+Hmax*t4*t8*t34*t55*t58*t
     &79*(bs3+t92+t93-t101)*(1.0D0/2.0D0)-Dc*Hmax*t2*t4*t8*t34*t35*t55*t
     &58*(3.0D0/2.0D0)-Dc*Hmax*t4*t8*t9*t34*t39*t55*t58*(3.0D0/2.0D0)-Dc
     &*Hmax*t4*t8*t10*t34*t43*t55*t58*(3.0D0/2.0D0)-Dc*Hmax*t4*t8*t34*t5
     &5*t58*t65*t96*(1.0D0/2.0D0)-Dc*Hmax*t4*t8*t34*t55*t58*t74*t95*(1.0
     &D0/2.0D0)-Dc*Hmax*t4*t8*t34*t55*t58*t79*t97*(1.0D0/2.0D0)-c1*c2*t2
     &*t4*t8*t34*t35*t55*t66*t69*(3.0D0/2.0D0)-c1*c2*t4*t8*t9*t34*t39*t5
     &5*t66*t69*(3.0D0/2.0D0)-c1*c2*t4*t8*t10*t34*t43*t55*t66*t69*(3.0D0
     &/2.0D0)-c1*c2*t4*t8*t34*t55*t65*t66*t69*t96*(1.0D0/2.0D0)-c1*c2*t4
     &*t8*t34*t55*t66*t69*t74*t95*(1.0D0/2.0D0)-c1*c2*t4*t8*t34*t55*t66*
     &t69*t79*t97*(1.0D0/2.0D0)



 
        NR_RE=A0

       elseif(index1==-1)then ! Reverse
       t2 = x-xo
      t4 = v+1.0D0
      t5 = 1.0D0/t4
      t6 = Sa*x
      t7 = Sm*x
      t8 = Sa-t6+t7
      t9 = 1.0D0/t8
      t36 = -e4+ep4+et4
      t37 = t5*t9*t36*(1.0D0/2.0D0)
      t3 = bs4-t37
      t40 = -e5+ep5+et5
      t41 = t5*t9*t40*(1.0D0/2.0D0)
      t10 = bs5-t41
      t44 = -e6+ep6+et6
      t45 = t5*t9*t44*(1.0D0/2.0D0)
      t11 = bs6-t45
      t15 = Sa*bs1
      t16 = Sa*bs1*v
      t17 = Sa*bs1*x
      t18 = Sm*bs1*x
      t19 = Sa*bs1*v*x
      t20 = Sm*bs1*v*x
      t22 = Sa*bs2
      t24 = Sa*bs2*v
      t26 = Sa*bs2*x
      t28 = Sm*bs2*x
      t30 = Sa*bs2*v*x
      t32 = Sm*bs2*v*x
      t12 = e1-e2-ep1+ep2-et1+et2+t15+t16-t17+t18-t19+t20-t22-t24+t26-t2
     &8+t30-t32
      t13 = 1.0D0/t4**2
      t14 = 1.0D0/t8**2
      t23 = Sa*bs3
      t25 = Sa*bs3*v
      t27 = Sa*bs3*x
      t29 = Sm*bs3*x
      t31 = Sa*bs3*v*x
      t33 = Sm*bs3*v*x
      t21 = e1-e3-ep1+ep3-et1+et3+t15+t16-t17+t18-t19+t20-t23-t25+t27-t2
     &9+t31-t33
      t34 = e2-e3-ep2+ep3-et2+et3+t22-t23+t24-t25-t26+t27+t28-t29-t30+t3
     &1+t32-t33
      t35 = sqrt(2.0D0)
      t38 = t3**2
      t39 = t38*6.0D0
      t42 = t10**2
      t43 = t42*6.0D0
      t46 = t11**2
      t47 = t46*6.0D0
      t48 = t12**2
      t49 = t13*t14*t48
      t50 = t21**2
      t51 = t13*t14*t50
      t52 = t34**2
      t53 = t13*t14*t52
      t54 = t39+t43+t47+t49+t51+t53
      t55 = abs(t54)
      t58 = sqrt(t55)
      t59 = kt*t35*t58*(1.0D0/2.0D0)
      t60 = exp(-t59)
      t56 = t60-1.0D0
      t57 = 1.0D0/sqrt(t55)
      t61 = t56**2
      t62 = c2*ztd
      t63 = t62+1.0D0
      t64 = 1.0D0/t63
      t65 = v**2
      t66 = t65*2.0D0
      t67 = t66+v-1.0D0
      t68 = 1.0D0/t67
      t69 = Tn*alpha
      t70 = v-1.0D0
      t72 = alpha*too
      t71 = -e2+ep2+et2+t69-t72
      t73 = -e3+ep3+et3+t69-t72
      t74 = -e1+ep1+et1+t69-t72
      t75 = t9*t68*t74*v
      t76 = t9*t68*t71*v
      t77 = Tn-too
      t78 = alpha*t77
      t79 = -e2+ep2+et2+t78
      t80 = -e1+ep1+et1+t78
      t81 = -e3+ep3+et3+t78
      t82 = t9*t68*t81*v
      t83 = t9*t68*t80*v
      t84 = t9*t68*t79*v
      t85 = Sa-Sm
      t86 = e1*v
      t87 = ep2*v
      t88 = ep3*v
      t89 = et2*v
      t90 = et3*v
      t91 = Tn*alpha*v
      t92 = e2*v
      t93 = ep1*v
      t94 = et1*v
      t95 = v*2.0D0
      t96 = t95+2.0D0
      t97 = e3*v
      t98 = ep1*2.0D0
      t99 = et1*2.0D0
      t100 = Sa*bs1*x*2.0D0
      t101 = Sa*bs1*v*x*2.0D0
      t102 = e1*(-2.0D0)+e2+e3-ep2-ep3-et2-et3+t22+t23+t24+t25-t26-t27+t
     &28+t29-t30-t31+t32+t33+t98+t99+t100+t101-Sa*bs1*2.0D0-Sa*bs1*v*2.0
     &D0-Sm*bs1*x*2.0D0-Sm*bs1*v*x*2.0D0
      t103 = ep2*2.0D0
      t104 = et2*2.0D0
      t105 = Sa*bs2*x*2.0D0
      t106 = Sa*bs2*v*x*2.0D0
      t107 = e1-e2*2.0D0+e3-ep1-ep3-et1-et3+t15+t16-t17+t18-t19+t20+t23+
     &t25-t27+t29-t31+t33+t103+t104+t105+t106-Sa*bs2*2.0D0-Sa*bs2*v*2.0D
     &0-Sm*bs2*x*2.0D0-Sm*bs2*v*x*2.0D0
      t108 = -e3+ep3+et3+t69-t72-t86+t87-t88+t89-t90+t91-t92+t93+t94+t97
     &-alpha*too*v
      t109 = ep3*2.0D0
      t110 = et3*2.0D0
      t111 = Sa*bs3*x*2.0D0
      t112 = Sa*bs3*v*x*2.0D0
      t113 = e1+e2-e3*2.0D0-ep1-ep2-et1-et2+t15+t16-t17+t18-t19+t20+t22+
     &t24-t26+t28-t30+t32+t109+t110+t111+t112-Sa*bs3*2.0D0-Sa*bs3*v*2.0D
     &0-Sm*bs3*x*2.0D0-Sm*bs3*v*x*2.0D0
      A0(1,1) = -et1+eto1+lamdat_r1*t2
      A0(2,1) = -et2+eto2+lamdat_r2*t2
      A0(3,1) = -et3+eto3+lamdat_r3*t2
      A0(4,1) = -et4+eto4+lamdat_r4*t2
      A0(5,1) = -et5+eto5+lamdat_r5*t2
      A0(6,1) = -et6+eto6+lamdat_r6*t2
      A0(7,1) = -ep1+epo1+c1*c2*t2*t5*t9*t35*t57*t61*t64*t102*(1.0D0/2.0
     &D0)
      A0(8,1) = -ep2+epo2+c1*c2*t2*t5*t9*t35*t57*t61*t64*t107*(1.0D0/2.0
     &D0)
      A0(9,1) = -ep3+epo3+c1*c2*t2*t5*t9*t35*t57*t61*t64*t113*(1.0D0/2.0
     &D0)
      A0(10,1) = -ep4+epo4-c1*c2*t2*t3*t35*t57*t61*t64*3.0D0
      A0(11,1) = -ep5+epo5-c1*c2*t2*t10*t35*t57*t61*t64*3.0D0
      A0(12,1) = -ep6+epo6-c1*c2*t2*t11*t35*t57*t61*t64*3.0D0
      A0(13,1) = -Yo-a3+rduo-Tn*rdso-lamdat_r4*t3-lamdat_r5*t10-lamdat_r
     &6*t11-lamdat_r2*(bs2+t75-t9*t68*t70*t71+t9*t68*t73*v)-lamdat_r3*(b
     &s3+t75+t76-t9*t68*t70*t73)-lamdat_r1*(bs1+t76-t9*t68*t70*(-e1+ep1+
     &et1+t69-alpha*too)+t9*t68*v*(-e3+ep3+et3+t69-alpha*too))+a2*(-(-x+
     &1.0D0)**n4+x**n3+1.0D0)*(1.0D0/2.0D0)-Dc*lamdat_r2*(t82+t83-t9*t68
     &*t70*t79)-Dc*lamdat_r1*(t82+t84-t9*t68*t70*t80)-Dc*lamdat_r3*(t83+
     &t84-t9*t68*t70*t81)+t13*t14*t36**2*t85*t96*(1.0D0/8.0D0)+t13*t14*t
     &40**2*t85*t96*(1.0D0/8.0D0)+t13*t14*t44**2*t85*t96*(1.0D0/8.0D0)+D
     &c*lamdat_r4*t5*t9*t36*(1.0D0/2.0D0)+Dc*lamdat_r5*t5*t9*t40*(1.0D0/
     &2.0D0)+Dc*lamdat_r6*t5*t9*t44*(1.0D0/2.0D0)-t14*t68*t71*t85*(-e2+e
     &p2+et2+t69-t72-t86-t87+t88-t89+t90+t91+t92+t93+t94-e3*v-alpha*too*
     &v)*(1.0D0/2.0D0)-t14*t68*t73*t85*t108*(1.0D0/2.0D0)-t14*t68*t74*t8
     &5*(-e1+ep1+et1+t69-t72+t86+t87+t88+t89+t90+t91-e2*v-e3*v-ep1*v-et1
     &*v-alpha*too*v)*(1.0D0/2.0D0)-c1*c2*t3*t5*t9*t35*t36*t57*t61*t64*(
     &3.0D0/2.0D0)-c1*c2*t5*t9*t10*t35*t40*t57*t61*t64*(3.0D0/2.0D0)-c1*
     &c2*t5*t9*t11*t35*t44*t57*t61*t64*(3.0D0/2.0D0)-c1*c2*t5*t14*t35*t5
     &7*t61*t64*t68*t108*t113*(1.0D0/2.0D0)-c1*c2*t5*t14*t35*t57*t61*t64
     &*t68*t102*(-e1+ep1+et1+t69-t72+t86+t87+t88+t89+t90+t91-t92-t93-t94
     &-t97-alpha*too*v)*(1.0D0/2.0D0)-c1*c2*t5*t14*t35*t57*t61*t64*t68*t
     &107*(-e2+ep2+et2+t69-t72-t86-t87+t88-t89+t90+t91+t92+t93+t94-t97-a
     &lpha*too*v)*(1.0D0/2.0D0)


 
 	 NR_RE=A0

!DEC$ FREEFORM
      endif

END subroutine

subroutine STRESS_UPDATE(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,PARAM,NPARAM,VAR,NVAR,STRESS)

implicit real*8(t)


!DEFINITION OF ABAQUS PROVIDED VARIABLES
REAL*8 STRAN(6),DSTRAN(6),COORDS(3),DROT(3,3),TIME(2)
REAL*8 TEMP,DTEMP,DTIME
REAL*8 STATEV(NSTATV),PROPS(NPROPS)
INTEGER JSTEP(4)
INTEGER NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,NSTATV,NPROPS


!DEFINITION OF  INPUT MODEL PARAMETERS
REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
integer NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


real*8   Tn,stress(6)
! Changing state variables
real*8 A0(6,1)


real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6


Tn=temp+dtemp;

CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)


   e1=e(1)
   e2=e(2)
   e3=e(3)
   e4=e(4)
   e5=e(5)
   e6=e(6)


   et1=et(1)
   et2=et(2)
   et3=et(3)
   et4=et(4)
   et5=et(5)
   et6=et(6)

   ep1=ep(1)
   ep2=ep(2)
   ep3=ep(3)
   ep4=ep(4)
   ep5=ep(5)
   ep6=ep(6)

      t2 = Sm*x
      t16 = Sa*x
      t3 = Sa+t2-t16
      t4 = 1.0D0/t3
      t5 = v**2
      t6 = t5*2.0D0
      t7 = t6+v-1.0D0
      t8 = 1.0D0/t7
      t9 = Tn*alpha
      t10 = e1*v
      t11 = ep2*v
      t12 = ep3*v
      t13 = et2*v
      t14 = et3*v
      t15 = Tn*alpha*v
      t17 = e2*v
      t18 = ep1*v
      t19 = et1*v
      t20 = v+1.0D0
      t21 = 1.0D0/t20
      A0(1,1) = t4*t8*(-e1+ep1+et1+t9+t10+t11+t12+t13+t14+t15-alpha*too-e2*v-e3*v-ep1*v-et1*v-alpha*too*v)
      A0(2,1) = t4*t8*(-e2+ep2+et2+t9-t10-t11+t12-t13+t14+t15+t17+t18+t19-alpha*too-e3*v-alpha*too*v)
      A0(3,1) = t4*t8*(-e3+ep3+et3+t9-t10+t11-t12+t13-t14+t15-t17+t18+t19-alpha*too+e3*v-alpha*too*v)
      A0(4,1) = t4*t21*(-e4+ep4+et4)*(-1.0D0/2.0D0)
      A0(5,1) = t4*t21*(-e5+ep5+et5)*(-1.0D0/2.0D0)
      A0(6,1) = t4*t21*(-e6+ep6+et6)*(-1.0D0/2.0D0)

	  stress(1:6)=A0(1:6,1)

end

subroutine JACOBIAN_MATRIX(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,PARAM,NPARAM,VAR,NVAR,STRESS,ddsdde)

implicit real*8(t)

!DEFINITION OF  INPUT MODEL PARAMETERS
REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
integer NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR

real*8   Tn           ! temp+dtemp

! Changing state variables
real*8 stress(6)
real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,stress1,stress2,stress3,stress4,stress5,stress6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6,lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6

real*8 AAA(6,1),BBB(1,6),LLL(6,6),CCC,A0(6,6),M(6,6),M_inv(6,6),N(6,1),dpds(1,6),dpdx,ddsdde(6,6),S_inv(6,6),AAA_rev(6,1)

Tn=temp+dtemp

CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

A0=0
M=0
M_INV=0
N=0
dpds=0
dpdx=0





   
   e1=e(1)
   e2=e(2)
   e3=e(3)
   e4=e(4)
   e5=e(5)
   e6=e(6)
    if(maxval(abs(e))<=1.0e-8)then
       e1=1.0e-8_8
    end if

	
   et1=et(1)
   et2=et(2)
   et3=et(3)
   et4=et(4)
   et5=et(5)
   et6=et(6)
   if(maxval(abs(et))<=1.0e-9)then
       et1=1.0e-9_8
    end if

   et_tr1=et_tr(1)
   et_tr2=et_tr(2)
   et_tr3=et_tr(3)
   et_tr4=et_tr(4)
   et_tr5=et_tr(5)
   et_tr6=et_tr(6)
   if(maxval(abs(et_tr))<=1.0e-10)then
       et_tr1=1.0e-10_8
    end if
	
   eto1=eto(1)
   eto2=eto(2)
   eto3=eto(3)
   eto4=eto(4)
   eto5=eto(5)
   eto6=eto(6)
    if(maxval(abs(eto))<=1.0e-11)then
       eto1=1.0e-11_8
    end if
	
   ep1=ep(1)
   ep2=ep(2)
   ep3=ep(3)
   ep4=ep(4)
   ep5=ep(5)
   ep6=ep(6)


   epo1=epo(1)
   epo2=epo(2)
   epo3=epo(3)
   epo4=epo(4)
   epo5=epo(5)
   epo6=epo(6)

   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)  
 
   !Backstress Calculation    
   bs1=backstress(1)
   bs2=backstress(2)
   bs3=backstress(3)
   bs4=backstress(4)
   bs5=backstress(5)
   bs6=backstress(6)
   
   stress1=stress(1)
   stress2=stress(2)
   stress3=stress(3)
   stress4=stress(4)
   stress5=stress(5)
   stress6=stress(6)



   if(abs(x-xo)<=0)then
       x=xo+1.0e-12_8
   end if







if (transformation == 0 .OR. (TRANSFORMATION/=0 .AND. BOUND_REACHED/=0)) then
!DEC$ NOFREEFORM
      t2 = Sm*x
      t9 = Sa*x
      t3 = Sa+t2-t9
      t4 = 1.0D0/t3
      t5 = v**2
      t6 = t5*2.0D0
      t7 = t6+v-1.0D0
      t8 = 1.0D0/t7
      t10 = v-1.0D0
      t11 = t4*t8*t10
      t12 = v+1.0D0
      t13 = 1.0D0/t12
      t14 = t4*t13*(1.0D0/2.0D0)
      A0(1,1) = t11
      A0(1,2) = -t4*t8*v
      A0(1,3) = -t4*t8*v
      A0(2,1) = -t4*t8*v
      A0(2,2) = t11
      A0(2,3) = -t4*t8*v
      A0(3,1) = -t4*t8*v
      A0(3,2) = -t4*t8*v
      A0(3,3) = t11
      A0(4,4) = t14
      A0(5,5) = t14
      A0(6,6) = t14
 
      ddsdde=A0

!DEC$ FREEFORM



elseif (transformation == 1) then   ! Forward Transformaion

!DEC$ NOFREEFORM
      t2 = stress1-stress2
      t3 = stress1-stress3
      t4 = stress2-stress3
      t5 = bs1-bs2+stress1-stress2
      t6 = bs1-bs3+stress1-stress3
      t7 = bs2-bs3+stress2-stress3
      t8 = bs4+stress4
      t9 = bs5+stress5
      t10 = bs6+stress6
      t11 = t2**2
      t12 = t11*(1.0D0/2.0D0)
      t13 = t3**2
      t14 = t13*(1.0D0/2.0D0)
      t15 = t4**2
      t16 = t15*(1.0D0/2.0D0)
      t17 = stress4**2
      t18 = t17*3.0D0
      t19 = stress5**2
      t20 = t19*3.0D0
      t21 = stress6**2
      t22 = t21*3.0D0
      t23 = t12+t14+t16+t18+t20+t22
      t24 = sqrt(t23)
      t40 = kt*t24
      t25 = exp(-t40)
      t26 = t25-1.0D0
      t27 = t5**2
      t28 = t27*(1.0D0/2.0D0)
      t29 = t6**2
      t30 = t29*(1.0D0/2.0D0)
      t31 = t7**2
      t32 = t31*(1.0D0/2.0D0)
      t33 = t8**2
      t34 = t33*3.0D0
      t35 = t9**2
      t36 = t35*3.0D0
      t37 = t10**2
      t38 = t37*3.0D0
      t39 = t28+t30+t32+t34+t36+t38
      t41 = bs2*(1.0D0/3.0D0)
      t42 = bs3*(1.0D0/3.0D0)
      t43 = stress2*(1.0D0/3.0D0)
      t44 = stress3*(1.0D0/3.0D0)
      t64 = bs1*(2.0D0/3.0D0)
      t65 = stress1*(2.0D0/3.0D0)
      t45 = t41+t42+t43+t44-t64-t65
      t46 = 1.0D0/sqrt(t39)
      t47 = kp*t24
      t48 = exp(t47)
      t49 = c2*t48*zt
      t50 = t49+1.0D0
      t51 = 1.0D0/t50
      t54 = stress1*2.0D0
      t88 = bs1*2.0D0
      t52 = bs2+bs3+stress2+stress3-t54-t88
      t53 = 1.0D0/t39**(3.0D0/2.0D0)
      t55 = bs2*(1.0D0/2.0D0)
      t56 = bs3*(1.0D0/2.0D0)
      t57 = stress2*(1.0D0/2.0D0)
      t58 = stress3*(1.0D0/2.0D0)
      t59 = -bs1-stress1+t55+t56+t57+t58
      t60 = 1.0D0/sqrt(t23)
      t61 = stress2+stress3-t54
      t62 = Sa-Sm
      t63 = x-xo
      t67 = stress2*2.0D0
      t97 = bs2*2.0D0
      t66 = bs1+bs3+stress1+stress3-t67-t97
      t68 = c2**2
      t69 = kp*t24*2.0D0
      t70 = exp(t69)
      t71 = 1.0D0/t50**2
      t72 = stress1+stress3-t67
      t87 = t62*x
      t73 = Sa-t87
      t74 = Hmax*t26*t46*(1.0D0/2.0D0)
      t76 = stress3*2.0D0
      t99 = bs3*2.0D0
      t75 = bs1+bs2+stress1+stress2-t76-t99
      t77 = stress1+stress2-t76
      t78 = bs4*6.0D0
      t79 = stress4*6.0D0
      t80 = t78+t79
      t81 = bs5*6.0D0
      t82 = stress5*6.0D0
      t83 = t81+t82
      t84 = bs6*6.0D0
      t85 = stress6*6.0D0
      t86 = t84+t85
      t89 = bs1*(1.0D0/3.0D0)
      t90 = stress1*(1.0D0/3.0D0)
      t95 = bs2*(2.0D0/3.0D0)
      t96 = stress2*(2.0D0/3.0D0)
      t91 = t42+t44+t89+t90-t95-t96
      t92 = bs1*(1.0D0/2.0D0)
      t93 = stress1*(1.0D0/2.0D0)
      t94 = -bs2-stress2+t56+t58+t92+t93
      t98 = c1*c2*t46*t48*t51
      t102 = bs3*(2.0D0/3.0D0)
      t103 = stress3*(2.0D0/3.0D0)
      t100 = t41+t43+t89+t90-t102-t103
      t101 = -bs3-stress3+t55+t57+t92+t93
      t104 = bs4*(3.0D0/2.0D0)
      t105 = stress4*(3.0D0/2.0D0)
      t106 = t104+t105
      t107 = bs5*(3.0D0/2.0D0)
      t108 = stress5*(3.0D0/2.0D0)
      t109 = t107+t108
      t110 = v*2.0D0
      t111 = t110+2.0D0
      t112 = t73*t111
      t113 = c1*c2*t46*t48*t51*3.0D0
      t114 = bs6*(3.0D0/2.0D0)
      t115 = stress6*(3.0D0/2.0D0)
      t116 = t114+t115
      A0(1,1) = Sa+t98+t63*(-Hmax*t26*t46+Hmax*t26*t45*t52*t53*(3.0D0/4.
     &0D0)+Hmax*kt*t25*t45*t46*t60*(stress1*(-2.0D0)+stress2+stress3)*(3
     &.0D0/4.0D0))-t62*x-c1*c2*t48*t51*t52*t53*t59*(1.0D0/2.0D0)+c1*c2*k
     &p*t46*t48*t51*t59*t60*t61*(1.0D0/2.0D0)-c1*kp*t46*t59*t60*t61*t68*
     &t70*t71*zt*(1.0D0/2.0D0)
      A0(1,2) = t63*(t74+Hmax*t26*t45*t53*t66*(3.0D0/4.0D0)+Hmax*kt*t25*
     &t45*t46*t60*(stress1-stress2*2.0D0+stress3)*(3.0D0/4.0D0))-t73*v-c
     &1*c2*t46*t48*t51*(1.0D0/2.0D0)-c1*c2*t48*t51*t53*t59*t66*(1.0D0/2.
     &0D0)+c1*c2*kp*t46*t48*t51*t59*t60*t72*(1.0D0/2.0D0)-c1*kp*t46*t59*
     &t60*t68*t70*t71*t72*zt*(1.0D0/2.0D0)
      A0(1,3) = t63*(t74+Hmax*t26*t45*t53*t75*(3.0D0/4.0D0)+Hmax*kt*t25*
     &t45*t46*t60*(stress1+stress2-stress3*2.0D0)*(3.0D0/4.0D0))-t73*v-c
     &1*c2*t46*t48*t51*(1.0D0/2.0D0)-c1*c2*t48*t51*t53*t59*t75*(1.0D0/2.
     &0D0)+c1*c2*kp*t46*t48*t51*t59*t60*t77*(1.0D0/2.0D0)-c1*kp*t46*t59*
     &t60*t68*t70*t71*t77*zt*(1.0D0/2.0D0)
      A0(1,4) = -t63*(Hmax*t26*t45*t53*t80*(3.0D0/4.0D0)+Hmax*kt*stress4
     &*t25*t45*t46*t60*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t59*t80*(1.0D0/2
     &.0D0)-c1*c2*kp*stress4*t46*t48*t51*t59*t60*3.0D0+c1*kp*stress4*t46
     &*t59*t60*t68*t70*t71*zt*3.0D0
      A0(1,5) = -t63*(Hmax*t26*t45*t53*t83*(3.0D0/4.0D0)+Hmax*kt*stress5
     &*t25*t45*t46*t60*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t59*t83*(1.0D0/2
     &.0D0)-c1*c2*kp*stress5*t46*t48*t51*t59*t60*3.0D0+c1*kp*stress5*t46
     &*t59*t60*t68*t70*t71*zt*3.0D0
      A0(1,6) = -t63*(Hmax*t26*t45*t53*t86*(3.0D0/4.0D0)+Hmax*kt*stress6
     &*t25*t45*t46*t60*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t59*t86*(1.0D0/2
     &.0D0)-c1*c2*kp*stress6*t46*t48*t51*t59*t60*3.0D0+c1*kp*stress6*t46
     &*t59*t60*t68*t70*t71*zt*3.0D0
      A0(2,1) = -t73*v+t63*(t74+Hmax*t26*t52*t53*t91*(3.0D0/4.0D0)+Hmax*
     &kt*t25*t46*t60*t61*t91*(3.0D0/4.0D0))-c1*c2*t46*t48*t51*(1.0D0/2.0
     &D0)-c1*c2*t48*t51*t52*t53*t94*(1.0D0/2.0D0)+c1*c2*kp*t46*t48*t51*t
     &60*t61*t94*(1.0D0/2.0D0)-c1*kp*t46*t60*t61*t68*t70*t71*t94*zt*(1.0
     &D0/2.0D0)
      A0(2,2) = Sa-t87+t98+t63*(-Hmax*t26*t46+Hmax*t26*t53*t66*t91*(3.0D
     &0/4.0D0)+Hmax*kt*t25*t46*t60*t72*t91*(3.0D0/4.0D0))-c1*c2*t48*t51*
     &t53*t66*t94*(1.0D0/2.0D0)+c1*c2*kp*t46*t48*t51*t60*t72*t94*(1.0D0/
     &2.0D0)-c1*kp*t46*t60*t68*t70*t71*t72*t94*zt*(1.0D0/2.0D0)
      A0(2,3) = -t73*v+t63*(t74+Hmax*t26*t53*t75*t91*(3.0D0/4.0D0)+Hmax*
     &kt*t25*t46*t60*t77*t91*(3.0D0/4.0D0))-c1*c2*t46*t48*t51*(1.0D0/2.0
     &D0)-c1*c2*t48*t51*t53*t75*t94*(1.0D0/2.0D0)+c1*c2*kp*t46*t48*t51*t
     &60*t77*t94*(1.0D0/2.0D0)-c1*kp*t46*t60*t68*t70*t71*t77*t94*zt*(1.0
     &D0/2.0D0)
      A0(2,4) = -t63*(Hmax*t26*t53*t80*t91*(3.0D0/4.0D0)+Hmax*kt*stress4
     &*t25*t46*t60*t91*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t80*t94*(1.0D0/2
     &.0D0)-c1*c2*kp*stress4*t46*t48*t51*t60*t94*3.0D0+c1*kp*stress4*t46
     &*t60*t68*t70*t71*t94*zt*3.0D0
      A0(2,5) = -t63*(Hmax*t26*t53*t83*t91*(3.0D0/4.0D0)+Hmax*kt*stress5
     &*t25*t46*t60*t91*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t83*t94*(1.0D0/2
     &.0D0)-c1*c2*kp*stress5*t46*t48*t51*t60*t94*3.0D0+c1*kp*stress5*t46
     &*t60*t68*t70*t71*t94*zt*3.0D0
      A0(2,6) = -t63*(Hmax*t26*t53*t86*t91*(3.0D0/4.0D0)+Hmax*kt*stress6
     &*t25*t46*t60*t91*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t86*t94*(1.0D0/2
     &.0D0)-c1*c2*kp*stress6*t46*t48*t51*t60*t94*3.0D0+c1*kp*stress6*t46
     &*t60*t68*t70*t71*t94*zt*3.0D0
      A0(3,1) = -t73*v+t63*(t74+Hmax*t26*t52*t53*t100*(3.0D0/4.0D0)+Hmax
     &*kt*t25*t46*t60*t61*t100*(3.0D0/4.0D0))-c1*c2*t46*t48*t51*(1.0D0/2
     &.0D0)-c1*c2*t48*t51*t52*t53*t101*(1.0D0/2.0D0)+c1*c2*kp*t46*t48*t5
     &1*t60*t61*t101*(1.0D0/2.0D0)-c1*kp*t46*t60*t61*t68*t70*t71*t101*zt
     &*(1.0D0/2.0D0)
      A0(3,2) = -t73*v+t63*(t74+Hmax*t26*t53*t66*t100*(3.0D0/4.0D0)+Hmax
     &*kt*t25*t46*t60*t72*t100*(3.0D0/4.0D0))-c1*c2*t46*t48*t51*(1.0D0/2
     &.0D0)-c1*c2*t48*t51*t53*t66*t101*(1.0D0/2.0D0)+c1*c2*kp*t46*t48*t5
     &1*t60*t72*t101*(1.0D0/2.0D0)-c1*kp*t46*t60*t68*t70*t71*t72*t101*zt
     &*(1.0D0/2.0D0)
      A0(3,3) = Sa-t87+t98+t63*(-Hmax*t26*t46+Hmax*t26*t53*t75*t100*(3.0
     &D0/4.0D0)+Hmax*kt*t25*t46*t60*t77*t100*(3.0D0/4.0D0))-c1*c2*t48*t5
     &1*t53*t75*t101*(1.0D0/2.0D0)+c1*c2*kp*t46*t48*t51*t60*t77*t101*(1.
     &0D0/2.0D0)-c1*kp*t46*t60*t68*t70*t71*t77*t101*zt*(1.0D0/2.0D0)
      A0(3,4) = -t63*(Hmax*t26*t53*t80*t100*(3.0D0/4.0D0)+Hmax*kt*stress
     &4*t25*t46*t60*t100*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t80*t101*(1.0D
     &0/2.0D0)-c1*c2*kp*stress4*t46*t48*t51*t60*t101*3.0D0+c1*kp*stress4
     &*t46*t60*t68*t70*t71*t101*zt*3.0D0
      A0(3,5) = -t63*(Hmax*t26*t53*t83*t100*(3.0D0/4.0D0)+Hmax*kt*stress
     &5*t25*t46*t60*t100*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t83*t101*(1.0D
     &0/2.0D0)-c1*c2*kp*stress5*t46*t48*t51*t60*t101*3.0D0+c1*kp*stress5
     &*t46*t60*t68*t70*t71*t101*zt*3.0D0
      A0(3,6) = -t63*(Hmax*t26*t53*t86*t100*(3.0D0/4.0D0)+Hmax*kt*stress
     &6*t25*t46*t60*t100*(9.0D0/2.0D0))+c1*c2*t48*t51*t53*t86*t101*(1.0D
     &0/2.0D0)-c1*c2*kp*stress6*t46*t48*t51*t60*t101*3.0D0+c1*kp*stress6
     &*t46*t60*t68*t70*t71*t101*zt*3.0D0
      A0(4,1) = -t63*(Hmax*t8*t26*t52*t53*(3.0D0/2.0D0)+Hmax*kt*t8*t25*t
     &46*t60*t61*(3.0D0/2.0D0))+c1*c2*t48*t51*t52*t53*t106-c1*c2*kp*t46*
     &t48*t51*t60*t61*t106+c1*kp*t46*t60*t61*t68*t70*t71*t106*zt
      A0(4,2) = -t63*(Hmax*t8*t26*t53*t66*(3.0D0/2.0D0)+Hmax*kt*t8*t25*t
     &46*t60*t72*(3.0D0/2.0D0))+c1*c2*t48*t51*t53*t66*t106-c1*c2*kp*t46*
     &t48*t51*t60*t72*t106+c1*kp*t46*t60*t68*t70*t71*t72*t106*zt
      A0(4,3) = -t63*(Hmax*t8*t26*t53*t75*(3.0D0/2.0D0)+Hmax*kt*t8*t25*t
     &46*t60*t77*(3.0D0/2.0D0))+c1*c2*t48*t51*t53*t75*t106-c1*c2*kp*t46*
     &t48*t51*t60*t77*t106+c1*kp*t46*t60*t68*t70*t71*t77*t106*zt
      A0(4,4) = t112+t113+t63*(Hmax*t26*t46*(-3.0D0)+Hmax*t8*t26*t53*t80
     &*(3.0D0/2.0D0)+Hmax*kt*stress4*t8*t25*t46*t60*9.0D0)-c1*c2*t48*t51
     &*t53*t80*t106+c1*c2*kp*stress4*t46*t48*t51*t60*t106*6.0D0-c1*kp*st
     &ress4*t46*t60*t68*t70*t71*t106*zt*6.0D0
      A0(4,5) = t63*(Hmax*t8*t26*t53*t83*(3.0D0/2.0D0)+Hmax*kt*stress5*t
     &8*t25*t46*t60*9.0D0)-c1*c2*t48*t51*t53*t83*t106+c1*c2*kp*stress5*t
     &46*t48*t51*t60*t106*6.0D0-c1*kp*stress5*t46*t60*t68*t70*t71*t106*z
     &t*6.0D0
      A0(4,6) = t63*(Hmax*t8*t26*t53*t86*(3.0D0/2.0D0)+Hmax*kt*stress6*t
     &8*t25*t46*t60*9.0D0)-c1*c2*t48*t51*t53*t86*t106+c1*c2*kp*stress6*t
     &46*t48*t51*t60*t106*6.0D0-c1*kp*stress6*t46*t60*t68*t70*t71*t106*z
     &t*6.0D0
      A0(5,1) = -t63*(Hmax*t9*t26*t52*t53*(3.0D0/2.0D0)+Hmax*kt*t9*t25*t
     &46*t60*t61*(3.0D0/2.0D0))+c1*c2*t48*t51*t52*t53*t109-c1*c2*kp*t46*
     &t48*t51*t60*t61*t109+c1*kp*t46*t60*t61*t68*t70*t71*t109*zt
      A0(5,2) = -t63*(Hmax*t9*t26*t53*t66*(3.0D0/2.0D0)+Hmax*kt*t9*t25*t
     &46*t60*t72*(3.0D0/2.0D0))+c1*c2*t48*t51*t53*t66*t109-c1*c2*kp*t46*
     &t48*t51*t60*t72*t109+c1*kp*t46*t60*t68*t70*t71*t72*t109*zt
      A0(5,3) = -t63*(Hmax*t9*t26*t53*t75*(3.0D0/2.0D0)+Hmax*kt*t9*t25*t
     &46*t60*t77*(3.0D0/2.0D0))+c1*c2*t48*t51*t53*t75*t109-c1*c2*kp*t46*
     &t48*t51*t60*t77*t109+c1*kp*t46*t60*t68*t70*t71*t77*t109*zt
      A0(5,4) = t63*(Hmax*t9*t26*t53*t80*(3.0D0/2.0D0)+Hmax*kt*stress4*t
     &9*t25*t46*t60*9.0D0)-c1*c2*t48*t51*t53*t80*t109+c1*c2*kp*stress4*t
     &46*t48*t51*t60*t109*6.0D0-c1*kp*stress4*t46*t60*t68*t70*t71*t109*z
     &t*6.0D0
      A0(5,5) = t112+t113+t63*(Hmax*t26*t46*(-3.0D0)+Hmax*t9*t26*t53*t83
     &*(3.0D0/2.0D0)+Hmax*kt*stress5*t9*t25*t46*t60*9.0D0)-c1*c2*t48*t51
     &*t53*t83*t109+c1*c2*kp*stress5*t46*t48*t51*t60*t109*6.0D0-c1*kp*st
     &ress5*t46*t60*t68*t70*t71*t109*zt*6.0D0
      A0(5,6) = t63*(Hmax*t9*t26*t53*t86*(3.0D0/2.0D0)+Hmax*kt*stress6*t
     &9*t25*t46*t60*9.0D0)-c1*c2*t48*t51*t53*t86*t109+c1*c2*kp*stress6*t
     &46*t48*t51*t60*t109*6.0D0-c1*kp*stress6*t46*t60*t68*t70*t71*t109*z
     &t*6.0D0
      A0(6,1) = -t63*(Hmax*t10*t26*t52*t53*(3.0D0/2.0D0)+Hmax*kt*t10*t25
     &*t46*t60*t61*(3.0D0/2.0D0))+c1*c2*t48*t51*t52*t53*t116-c1*c2*kp*t4
     &6*t48*t51*t60*t61*t116+c1*kp*t46*t60*t61*t68*t70*t71*t116*zt
      A0(6,2) = -t63*(Hmax*t10*t26*t53*t66*(3.0D0/2.0D0)+Hmax*kt*t10*t25
     &*t46*t60*t72*(3.0D0/2.0D0))+c1*c2*t48*t51*t53*t66*t116-c1*c2*kp*t4
     &6*t48*t51*t60*t72*t116+c1*kp*t46*t60*t68*t70*t71*t72*t116*zt
      A0(6,3) = -t63*(Hmax*t10*t26*t53*t75*(3.0D0/2.0D0)+Hmax*kt*t10*t25
     &*t46*t60*t77*(3.0D0/2.0D0))+c1*c2*t48*t51*t53*t75*t116-c1*c2*kp*t4
     &6*t48*t51*t60*t77*t116+c1*kp*t46*t60*t68*t70*t71*t77*t116*zt
      A0(6,4) = t63*(Hmax*t10*t26*t53*t80*(3.0D0/2.0D0)+Hmax*kt*stress4*
     &t10*t25*t46*t60*9.0D0)-c1*c2*t48*t51*t53*t80*t116+c1*c2*kp*stress4
     &*t46*t48*t51*t60*t116*6.0D0-c1*kp*stress4*t46*t60*t68*t70*t71*t116
     &*zt*6.0D0
      A0(6,5) = t63*(Hmax*t10*t26*t53*t83*(3.0D0/2.0D0)+Hmax*kt*stress5*
     &t10*t25*t46*t60*9.0D0)-c1*c2*t48*t51*t53*t83*t116+c1*c2*kp*stress5
     &*t46*t48*t51*t60*t116*6.0D0-c1*kp*stress5*t46*t60*t68*t70*t71*t116
     &*zt*6.0D0
      A0(6,6) = t112+t113+t63*(Hmax*t26*t46*(-3.0D0)+Hmax*t10*t26*t53*t8
     &6*(3.0D0/2.0D0)+Hmax*kt*stress6*t10*t25*t46*t60*9.0D0)-c1*c2*t48*t
     &51*t53*t86*t116+c1*c2*kp*stress6*t46*t48*t51*t60*t116*6.0D0-c1*kp*
     &stress6*t46*t60*t68*t70*t71*t116*zt*6.0D0

 
 	   M_INV=A0
 	   A0=0
 	   CALL inverse(M_INV,M,6)
 !********************************************************************************************************
 !********************************************************************************************************
      t2 = Sa-Sm
      t3 = stress1-stress2
      t4 = stress1-stress3
      t5 = stress2-stress3
      t6 = bs1-bs2+stress1-stress2
      t7 = bs1-bs3+stress1-stress3
      t8 = bs2-bs3+stress2-stress3
      t9 = bs4+stress4
      t10 = bs5+stress5
      t11 = bs6+stress6
      t12 = t3**2
      t13 = t12*(1.0D0/2.0D0)
      t14 = t4**2
      t15 = t14*(1.0D0/2.0D0)
      t16 = t5**2
      t17 = t16*(1.0D0/2.0D0)
      t18 = stress4**2
      t19 = t18*3.0D0
      t20 = stress5**2
      t21 = t20*3.0D0
      t22 = stress6**2
      t23 = t22*3.0D0
      t24 = t13+t15+t17+t19+t21+t23
      t25 = sqrt(t24)
      t26 = kp*t25
      t27 = exp(t26)
      t28 = t6**2
      t29 = t28*(1.0D0/2.0D0)
      t30 = t7**2
      t31 = t30*(1.0D0/2.0D0)
      t32 = t8**2
      t33 = t32*(1.0D0/2.0D0)
      t34 = t9**2
      t35 = t34*3.0D0
      t36 = t10**2
      t37 = t36*3.0D0
      t38 = t11**2
      t39 = t38*3.0D0
      t40 = t29+t31+t33+t35+t37+t39
      t41 = 1.0D0/sqrt(t40)
      t42 = stress3*t2*v
      t54 = kt*t25
      t43 = exp(-t54)
      t44 = t43-1.0D0
      t45 = bs3*(1.0D0/3.0D0)
      t46 = stress3*(1.0D0/3.0D0)
      t47 = c2*t27*zt
      t48 = t47+1.0D0
      t49 = 1.0D0/t48
      t50 = bs3*(1.0D0/2.0D0)
      t51 = stress3*(1.0D0/2.0D0)
      t52 = stress1*t2*v
      t53 = stress2*t2*v
      t55 = bs1*(1.0D0/3.0D0)
      t56 = bs2*(1.0D0/3.0D0)
      t57 = stress1*(1.0D0/3.0D0)
      t58 = stress2*(1.0D0/3.0D0)
      t59 = bs1*(1.0D0/2.0D0)
      t60 = bs2*(1.0D0/2.0D0)
      t61 = stress1*(1.0D0/2.0D0)
      t62 = stress2*(1.0D0/2.0D0)
      t63 = v*2.0D0
      t64 = t63+2.0D0
      A0(1,1) = t42+t53-stress1*t2+Hmax*t41*t44*(bs1*(-2.0D0/3.0D0)-stre
     &ss1*(2.0D0/3.0D0)+t45+t46+t56+t58)*(3.0D0/2.0D0)-c1*c2*t27*t41*t49
     &*(-bs1-stress1+t50+t51+t60+t62)
      A0(2,1) = t42+t52-stress2*t2+Hmax*t41*t44*(bs2*(-2.0D0/3.0D0)-stre
     &ss2*(2.0D0/3.0D0)+t45+t46+t55+t57)*(3.0D0/2.0D0)-c1*c2*t27*t41*t49
     &*(-bs2-stress2+t50+t51+t59+t61)
      A0(3,1) = t52+t53-stress3*t2+Hmax*t41*t44*(bs3*(-2.0D0/3.0D0)-stre
     &ss3*(2.0D0/3.0D0)+t55+t56+t57+t58)*(3.0D0/2.0D0)-c1*c2*t27*t41*t49
     &*(-bs3-stress3+t59+t60+t61+t62)
      A0(4,1) = -stress4*t2*t64-Hmax*t9*t41*t44*3.0D0+c1*c2*t27*t41*t49*
     &(bs4*(3.0D0/2.0D0)+stress4*(3.0D0/2.0D0))*2.0D0
      A0(5,1) = -stress5*t2*t64-Hmax*t10*t41*t44*3.0D0+c1*c2*t27*t41*t49
     &*(bs5*(3.0D0/2.0D0)+stress5*(3.0D0/2.0D0))*2.0D0
      A0(6,1) = -stress6*t2*t64-Hmax*t11*t41*t44*3.0D0+c1*c2*t27*t41*t49
     &*(bs6*(3.0D0/2.0D0)+stress6*(3.0D0/2.0D0))*2.0D0

 
       N(1:6,1)=A0(1:6,1)
 	   A0=0
 !********************************************************************************************************
 !********************************************************************************************************
       t2 = Sa-Sm
      t3 = stress1-stress2
      t4 = stress1-stress3
      t5 = stress2-stress3
      t6 = bs1-bs2+stress1-stress2
      t7 = bs1-bs3+stress1-stress3
      t8 = bs2-bs3+stress2-stress3
      t9 = bs4+stress4
      t10 = bs5+stress5
      t11 = bs6+stress6
      t12 = t3**2
      t13 = t12*(1.0D0/2.0D0)
      t14 = t4**2
      t15 = t14*(1.0D0/2.0D0)
      t16 = t5**2
      t17 = t16*(1.0D0/2.0D0)
      t18 = stress4**2
      t19 = t18*3.0D0
      t20 = stress5**2
      t21 = t20*3.0D0
      t22 = stress6**2
      t23 = t22*3.0D0
      t24 = t13+t15+t17+t19+t21+t23
      t25 = sqrt(t24)
      t42 = kt*t25
      t26 = exp(-t42)
      t27 = t26-1.0D0
      t28 = t6**2
      t29 = t28*(1.0D0/2.0D0)
      t30 = t7**2
      t31 = t30*(1.0D0/2.0D0)
      t32 = t8**2
      t33 = t32*(1.0D0/2.0D0)
      t34 = t9**2
      t35 = t34*3.0D0
      t36 = t10**2
      t37 = t36*3.0D0
      t38 = t11**2
      t39 = t38*3.0D0
      t40 = t29+t31+t33+t35+t37+t39
      t41 = 1.0D0/sqrt(t40)
      t45 = bs1*2.0D0
      t46 = stress1*2.0D0
      t43 = bs2+bs3+stress2+stress3-t45-t46
      t44 = 1.0D0/t40**(3.0D0/2.0D0)
      t47 = bs2*(1.0D0/3.0D0)
      t48 = bs3*(1.0D0/3.0D0)
      t49 = stress2*(1.0D0/3.0D0)
      t50 = stress3*(1.0D0/3.0D0)
      t58 = bs1*(2.0D0/3.0D0)
      t59 = stress1*(2.0D0/3.0D0)
      t51 = t47+t48+t49+t50-t58-t59
      t52 = kp*t25
      t53 = exp(t52)
      t54 = c2*t53*zt
      t55 = t54+1.0D0
      t56 = 1.0D0/t55
      t57 = bs1+stress1
      t60 = bs2+stress2
      t61 = bs3+stress3
      t62 = bs1*(1.0D0/3.0D0)
      t63 = stress1*(1.0D0/3.0D0)
      t64 = stress2+stress3-t46
      t65 = 1.0D0/sqrt(t24)
      t68 = bs2*(2.0D0/3.0D0)
      t69 = stress2*(2.0D0/3.0D0)
      t66 = t48+t50+t62+t63-t68-t69
      t70 = bs3*(2.0D0/3.0D0)
      t71 = stress3*(2.0D0/3.0D0)
      t67 = t47+t49+t62+t63-t70-t71
      t72 = bs2*(1.0D0/2.0D0)
      t73 = bs3*(1.0D0/2.0D0)
      t74 = stress2*(1.0D0/2.0D0)
      t75 = stress3*(1.0D0/2.0D0)
      t76 = -bs1-stress1+t72+t73+t74+t75
      t77 = bs1*(1.0D0/2.0D0)
      t78 = stress1*(1.0D0/2.0D0)
      t79 = bs4*(3.0D0/2.0D0)
      t80 = stress4*(3.0D0/2.0D0)
      t81 = t79+t80
      t82 = bs5*(3.0D0/2.0D0)
      t83 = stress5*(3.0D0/2.0D0)
      t84 = t82+t83
      t85 = bs6*(3.0D0/2.0D0)
      t86 = stress6*(3.0D0/2.0D0)
      t87 = t85+t86
      t88 = -bs2-stress2+t73+t75+t77+t78
      t89 = -bs3-stress3+t72+t74+t77+t78
      t90 = c2**2
      t91 = kp*t25*2.0D0
      t92 = exp(t91)
      t93 = 1.0D0/t55**2
      t94 = stress3*t2*v
      t95 = Hmax*t27*t41*t61*(1.0D0/2.0D0)
      t97 = bs2*2.0D0
      t98 = stress2*2.0D0
      t96 = bs1+bs3+stress1+stress3-t97-t98
      t99 = stress1+stress3-t98
      t100 = stress1*t2*v
      t101 = stress2*t2*v
      t102 = Hmax*t27*t41*t57*(1.0D0/2.0D0)
      t103 = Hmax*t27*t41*t60*(1.0D0/2.0D0)
      t105 = bs3*2.0D0
      t106 = stress3*2.0D0
      t104 = bs1+bs2+stress1+stress2-t105-t106
      t107 = stress1+stress2-t106
      t108 = bs4*6.0D0
      t109 = stress4*6.0D0
      t110 = t108+t109
      t111 = v*2.0D0
      t112 = t111+2.0D0
      t113 = bs5*6.0D0
      t114 = stress5*6.0D0
      t115 = t113+t114
      t116 = bs6*6.0D0
      t117 = stress6*6.0D0
      t118 = t116+t117
      A0(1,1) = t94+t95+t101+t103-stress1*t2+Hmax*t27*t41*t51*(3.0D0/2.0
     &D0)-Hmax*t27*t41*t57+Dc*Hmax*stress1*t27*t41-Dc*Hmax*stress2*t27*t
     &41*(1.0D0/2.0D0)-Dc*Hmax*stress3*t27*t41*(1.0D0/2.0D0)-Dc*Hmax*t27
     &*t41*t51*(3.0D0/2.0D0)-Hmax*t27*t34*t43*t44*(3.0D0/2.0D0)-Hmax*t27
     &*t36*t43*t44*(3.0D0/2.0D0)-Hmax*t27*t38*t43*t44*(3.0D0/2.0D0)+Hmax
     &*t27*t43*t44*t51*t57*(3.0D0/4.0D0)+Hmax*t27*t43*t44*t60*t66*(3.0D0
     &/4.0D0)+Hmax*t27*t43*t44*t61*t67*(3.0D0/4.0D0)+c1*c2*stress1*t41*t
     &53*t56-c1*c2*stress2*t41*t53*t56*(1.0D0/2.0D0)-c1*c2*stress3*t41*t
     &53*t56*(1.0D0/2.0D0)-c1*c2*t41*t53*t56*t76+Dc*Hmax*stress4*t9*t27*
     &t43*t44*(3.0D0/2.0D0)+Dc*Hmax*stress5*t10*t27*t43*t44*(3.0D0/2.0D0
     &)+Dc*Hmax*stress6*t11*t27*t43*t44*(3.0D0/2.0D0)-Dc*Hmax*stress1*t2
     &7*t43*t44*t51*(3.0D0/4.0D0)-Dc*Hmax*stress2*t27*t43*t44*t66*(3.0D0
     &/4.0D0)-Dc*Hmax*stress3*t27*t43*t44*t67*(3.0D0/4.0D0)-Hmax*kt*t26*
     &t34*t41*t64*t65*(3.0D0/2.0D0)-Hmax*kt*t26*t36*t41*t64*t65*(3.0D0/2
     &.0D0)-Hmax*kt*t26*t38*t41*t64*t65*(3.0D0/2.0D0)+Hmax*kt*t26*t41*t5
     &1*t57*t64*t65*(3.0D0/4.0D0)+Hmax*kt*t26*t41*t60*t64*t65*t66*(3.0D0
     &/4.0D0)+Hmax*kt*t26*t41*t61*t64*t65*t67*(3.0D0/4.0D0)-c1*c2*stress
     &1*t43*t44*t53*t56*t76*(1.0D0/2.0D0)+c1*c2*stress4*t43*t44*t53*t56*
     &t81+c1*c2*stress5*t43*t44*t53*t56*t84-c1*c2*stress2*t43*t44*t53*t5
     &6*t88*(1.0D0/2.0D0)-c1*c2*stress3*t43*t44*t53*t56*t89*(1.0D0/2.0D0
     &)+c1*c2*stress6*t43*t44*t53*t56*t87+Dc*Hmax*kt*stress4*t9*t26*t41*
     &t64*t65*(3.0D0/2.0D0)+Dc*Hmax*kt*stress5*t10*t26*t41*t64*t65*(3.0D
     &0/2.0D0)+Dc*Hmax*kt*stress6*t11*t26*t41*t64*t65*(3.0D0/2.0D0)-Dc*H
     &max*kt*stress1*t26*t41*t51*t64*t65*(3.0D0/4.0D0)-Dc*Hmax*kt*stress
     &2*t26*t41*t64*t65*t66*(3.0D0/4.0D0)-Dc*Hmax*kt*stress3*t26*t41*t64
     &*t65*t67*(3.0D0/4.0D0)+c1*c2*kp*stress1*t41*t53*t56*t64*t65*t76*(1
     &.0D0/2.0D0)-c1*c2*kp*stress4*t41*t53*t56*t64*t65*t81-c1*c2*kp*stre
     &ss5*t41*t53*t56*t64*t65*t84+c1*c2*kp*stress2*t41*t53*t56*t64*t65*t
     &88*(1.0D0/2.0D0)+c1*c2*kp*stress3*t41*t53*t56*t64*t65*t89*(1.0D0/2
     &.0D0)-c1*c2*kp*stress6*t41*t53*t56*t64*t65*t87-c1*kp*stress1*t41*t
     &64*t65*t76*t90*t92*t93*zt*(1.0D0/2.0D0)+c1*kp*stress4*t41*t64*t65*
     &t81*t90*t92*t93*zt+c1*kp*stress5*t41*t64*t65*t84*t90*t92*t93*zt-c1
     &*kp*stress2*t41*t64*t65*t88*t90*t92*t93*zt*(1.0D0/2.0D0)-c1*kp*str
     &ess3*t41*t64*t65*t89*t90*t92*t93*zt*(1.0D0/2.0D0)+c1*kp*stress6*t4
     &1*t64*t65*t87*t90*t92*t93*zt
      A0(1,2) = t94+t95+t100+t102-stress2*t2-Hmax*t27*t41*t60+Hmax*t27*t
     &41*t66*(3.0D0/2.0D0)-Dc*Hmax*stress1*t27*t41*(1.0D0/2.0D0)+Dc*Hmax
     &*stress2*t27*t41-Dc*Hmax*stress3*t27*t41*(1.0D0/2.0D0)-Dc*Hmax*t27
     &*t41*t66*(3.0D0/2.0D0)-Hmax*t27*t34*t44*t96*(3.0D0/2.0D0)-Hmax*t27
     &*t36*t44*t96*(3.0D0/2.0D0)-Hmax*t27*t38*t44*t96*(3.0D0/2.0D0)+Hmax
     &*t27*t44*t51*t57*t96*(3.0D0/4.0D0)+Hmax*t27*t44*t60*t66*t96*(3.0D0
     &/4.0D0)+Hmax*t27*t44*t61*t67*t96*(3.0D0/4.0D0)-c1*c2*stress1*t41*t
     &53*t56*(1.0D0/2.0D0)+c1*c2*stress2*t41*t53*t56-c1*c2*stress3*t41*t
     &53*t56*(1.0D0/2.0D0)-c1*c2*t41*t53*t56*t88+Dc*Hmax*stress4*t9*t27*
     &t44*t96*(3.0D0/2.0D0)+Dc*Hmax*stress5*t10*t27*t44*t96*(3.0D0/2.0D0
     &)+Dc*Hmax*stress6*t11*t27*t44*t96*(3.0D0/2.0D0)-Dc*Hmax*stress1*t2
     &7*t44*t51*t96*(3.0D0/4.0D0)-Dc*Hmax*stress2*t27*t44*t66*t96*(3.0D0
     &/4.0D0)-Dc*Hmax*stress3*t27*t44*t67*t96*(3.0D0/4.0D0)-Hmax*kt*t26*
     &t34*t41*t65*t99*(3.0D0/2.0D0)-Hmax*kt*t26*t36*t41*t65*t99*(3.0D0/2
     &.0D0)-Hmax*kt*t26*t38*t41*t65*t99*(3.0D0/2.0D0)+Hmax*kt*t26*t41*t5
     &1*t57*t65*t99*(3.0D0/4.0D0)+Hmax*kt*t26*t41*t60*t65*t66*t99*(3.0D0
     &/4.0D0)+Hmax*kt*t26*t41*t61*t65*t67*t99*(3.0D0/4.0D0)-c1*c2*stress
     &1*t44*t53*t56*t76*t96*(1.0D0/2.0D0)+c1*c2*stress4*t44*t53*t56*t81*
     &t96+c1*c2*stress5*t44*t53*t56*t84*t96-c1*c2*stress2*t44*t53*t56*t8
     &8*t96*(1.0D0/2.0D0)-c1*c2*stress3*t44*t53*t56*t89*t96*(1.0D0/2.0D0
     &)+c1*c2*stress6*t44*t53*t56*t87*t96+Dc*Hmax*kt*stress4*t9*t26*t41*
     &t65*t99*(3.0D0/2.0D0)+Dc*Hmax*kt*stress5*t10*t26*t41*t65*t99*(3.0D
     &0/2.0D0)+Dc*Hmax*kt*stress6*t11*t26*t41*t65*t99*(3.0D0/2.0D0)-Dc*H
     &max*kt*stress1*t26*t41*t51*t65*t99*(3.0D0/4.0D0)-Dc*Hmax*kt*stress
     &2*t26*t41*t65*t66*t99*(3.0D0/4.0D0)-Dc*Hmax*kt*stress3*t26*t41*t65
     &*t67*t99*(3.0D0/4.0D0)+c1*c2*kp*stress1*t41*t53*t56*t65*t76*t99*(1
     &.0D0/2.0D0)-c1*c2*kp*stress4*t41*t53*t56*t65*t81*t99-c1*c2*kp*stre
     &ss5*t41*t53*t56*t65*t84*t99+c1*c2*kp*stress2*t41*t53*t56*t65*t88*t
     &99*(1.0D0/2.0D0)+c1*c2*kp*stress3*t41*t53*t56*t65*t89*t99*(1.0D0/2
     &.0D0)-c1*c2*kp*stress6*t41*t53*t56*t65*t87*t99-c1*kp*stress1*t41*t
     &65*t76*t90*t92*t93*t99*zt*(1.0D0/2.0D0)+c1*kp*stress4*t41*t65*t81*
     &t90*t92*t93*t99*zt+c1*kp*stress5*t41*t65*t84*t90*t92*t93*t99*zt-c1
     &*kp*stress2*t41*t65*t88*t90*t92*t93*t99*zt*(1.0D0/2.0D0)-c1*kp*str
     &ess3*t41*t65*t89*t90*t92*t93*t99*zt*(1.0D0/2.0D0)+c1*kp*stress6*t4
     &1*t65*t87*t90*t92*t93*t99*zt
      A0(1,3) = t100+t101+t102+t103-stress3*t2-Hmax*t27*t41*t61+Hmax*t27
     &*t41*t67*(3.0D0/2.0D0)-Dc*Hmax*stress1*t27*t41*(1.0D0/2.0D0)-Dc*Hm
     &ax*stress2*t27*t41*(1.0D0/2.0D0)+Dc*Hmax*stress3*t27*t41-Dc*Hmax*t
     &27*t41*t67*(3.0D0/2.0D0)-Hmax*t27*t34*t44*t104*(3.0D0/2.0D0)-Hmax*
     &t27*t36*t44*t104*(3.0D0/2.0D0)-Hmax*t27*t38*t44*t104*(3.0D0/2.0D0)
     &+Hmax*t27*t44*t51*t57*t104*(3.0D0/4.0D0)+Hmax*t27*t44*t60*t66*t104
     &*(3.0D0/4.0D0)+Hmax*t27*t44*t61*t67*t104*(3.0D0/4.0D0)-c1*c2*stres
     &s1*t41*t53*t56*(1.0D0/2.0D0)-c1*c2*stress2*t41*t53*t56*(1.0D0/2.0D
     &0)+c1*c2*stress3*t41*t53*t56-c1*c2*t41*t53*t56*t89+Dc*Hmax*stress4
     &*t9*t27*t44*t104*(3.0D0/2.0D0)+Dc*Hmax*stress5*t10*t27*t44*t104*(3
     &.0D0/2.0D0)+Dc*Hmax*stress6*t11*t27*t44*t104*(3.0D0/2.0D0)-Dc*Hmax
     &*stress1*t27*t44*t51*t104*(3.0D0/4.0D0)-Dc*Hmax*stress2*t27*t44*t6
     &6*t104*(3.0D0/4.0D0)-Dc*Hmax*stress3*t27*t44*t67*t104*(3.0D0/4.0D0
     &)-Hmax*kt*t26*t34*t41*t65*t107*(3.0D0/2.0D0)-Hmax*kt*t26*t36*t41*t
     &65*t107*(3.0D0/2.0D0)-Hmax*kt*t26*t38*t41*t65*t107*(3.0D0/2.0D0)+H
     &max*kt*t26*t41*t51*t57*t65*t107*(3.0D0/4.0D0)+Hmax*kt*t26*t41*t60*
     &t65*t66*t107*(3.0D0/4.0D0)+Hmax*kt*t26*t41*t61*t65*t67*t107*(3.0D0
     &/4.0D0)-c1*c2*stress1*t44*t53*t56*t76*t104*(1.0D0/2.0D0)+c1*c2*str
     &ess4*t44*t53*t56*t81*t104+c1*c2*stress5*t44*t53*t56*t84*t104-c1*c2
     &*stress2*t44*t53*t56*t88*t104*(1.0D0/2.0D0)-c1*c2*stress3*t44*t53*
     &t56*t89*t104*(1.0D0/2.0D0)+c1*c2*stress6*t44*t53*t56*t87*t104+Dc*H
     &max*kt*stress4*t9*t26*t41*t65*t107*(3.0D0/2.0D0)+Dc*Hmax*kt*stress
     &5*t10*t26*t41*t65*t107*(3.0D0/2.0D0)+Dc*Hmax*kt*stress6*t11*t26*t4
     &1*t65*t107*(3.0D0/2.0D0)-Dc*Hmax*kt*stress1*t26*t41*t51*t65*t107*(
     &3.0D0/4.0D0)-Dc*Hmax*kt*stress2*t26*t41*t65*t66*t107*(3.0D0/4.0D0)
     &-Dc*Hmax*kt*stress3*t26*t41*t65*t67*t107*(3.0D0/4.0D0)+c1*c2*kp*st
     &ress1*t41*t53*t56*t65*t76*t107*(1.0D0/2.0D0)-c1*c2*kp*stress4*t41*
     &t53*t56*t65*t81*t107-c1*c2*kp*stress5*t41*t53*t56*t65*t84*t107+c1*
     &c2*kp*stress2*t41*t53*t56*t65*t88*t107*(1.0D0/2.0D0)+c1*c2*kp*stre
     &ss3*t41*t53*t56*t65*t89*t107*(1.0D0/2.0D0)-c1*c2*kp*stress6*t41*t5
     &3*t56*t65*t87*t107-c1*kp*stress1*t41*t65*t76*t90*t92*t93*t107*zt*(
     &1.0D0/2.0D0)+c1*kp*stress4*t41*t65*t81*t90*t92*t93*t107*zt+c1*kp*s
     &tress5*t41*t65*t84*t90*t92*t93*t107*zt-c1*kp*stress2*t41*t65*t88*t
     &90*t92*t93*t107*zt*(1.0D0/2.0D0)-c1*kp*stress3*t41*t65*t89*t90*t92
     &*t93*t107*zt*(1.0D0/2.0D0)+c1*kp*stress6*t41*t65*t87*t90*t92*t93*t
     &107*zt
      A0(1,4) = -stress4*t2*t112-Hmax*t27*t41*(bs4*2.0D0+stress4*2.0D0)*
     &3.0D0+Dc*Hmax*stress4*t27*t41*3.0D0+Dc*Hmax*t9*t27*t41*3.0D0+Hmax*
     &t27*t34*t44*t110*(3.0D0/2.0D0)+Hmax*t27*t36*t44*t110*(3.0D0/2.0D0)
     &+Hmax*t27*t38*t44*t110*(3.0D0/2.0D0)-Hmax*t27*t44*t51*t57*t110*(3.
     &0D0/4.0D0)-Hmax*t27*t44*t60*t66*t110*(3.0D0/4.0D0)-Hmax*t27*t44*t6
     &1*t67*t110*(3.0D0/4.0D0)+c1*c2*stress4*t41*t53*t56*3.0D0+c1*c2*t41
     &*t53*t56*t81*2.0D0-Dc*Hmax*stress4*t9*t27*t44*t110*(3.0D0/2.0D0)-D
     &c*Hmax*stress5*t10*t27*t44*t110*(3.0D0/2.0D0)-Dc*Hmax*stress6*t11*
     &t27*t44*t110*(3.0D0/2.0D0)+Dc*Hmax*stress1*t27*t44*t51*t110*(3.0D0
     &/4.0D0)+Dc*Hmax*stress2*t27*t44*t66*t110*(3.0D0/4.0D0)+Dc*Hmax*str
     &ess3*t27*t44*t67*t110*(3.0D0/4.0D0)+Hmax*kt*stress4*t26*t34*t41*t6
     &5*9.0D0+Hmax*kt*stress4*t26*t36*t41*t65*9.0D0+Hmax*kt*stress4*t26*
     &t38*t41*t65*9.0D0-Dc*Hmax*kt*t9*t18*t26*t41*t65*9.0D0-Hmax*kt*stre
     &ss4*t26*t41*t51*t57*t65*(9.0D0/2.0D0)-Hmax*kt*stress4*t26*t41*t60*
     &t65*t66*(9.0D0/2.0D0)-Hmax*kt*stress4*t26*t41*t61*t65*t67*(9.0D0/2
     &.0D0)+c1*c2*stress1*t44*t53*t56*t76*t110*(1.0D0/2.0D0)-c1*c2*stres
     &s4*t44*t53*t56*t81*t110-c1*c2*stress5*t44*t53*t56*t84*t110+c1*c2*s
     &tress2*t44*t53*t56*t88*t110*(1.0D0/2.0D0)+c1*c2*stress3*t44*t53*t5
     &6*t89*t110*(1.0D0/2.0D0)-c1*c2*stress6*t44*t53*t56*t87*t110-Dc*Hma
     &x*kt*stress4*stress5*t10*t26*t41*t65*9.0D0-Dc*Hmax*kt*stress4*stre
     &ss6*t11*t26*t41*t65*9.0D0+Dc*Hmax*kt*stress1*stress4*t26*t41*t51*t
     &65*(9.0D0/2.0D0)+Dc*Hmax*kt*stress2*stress4*t26*t41*t65*t66*(9.0D0
     &/2.0D0)+Dc*Hmax*kt*stress3*stress4*t26*t41*t65*t67*(9.0D0/2.0D0)+c
     &1*c2*kp*t18*t41*t53*t56*t65*t81*6.0D0-c1*c2*kp*stress1*stress4*t41
     &*t53*t56*t65*t76*3.0D0+c1*c2*kp*stress4*stress5*t41*t53*t56*t65*t8
     &4*6.0D0-c1*c2*kp*stress2*stress4*t41*t53*t56*t65*t88*3.0D0-c1*c2*k
     &p*stress3*stress4*t41*t53*t56*t65*t89*3.0D0+c1*c2*kp*stress4*stres
     &s6*t41*t53*t56*t65*t87*6.0D0-c1*kp*t18*t41*t65*t81*t90*t92*t93*zt*
     &6.0D0+c1*kp*stress1*stress4*t41*t65*t76*t90*t92*t93*zt*3.0D0-c1*kp
     &*stress4*stress5*t41*t65*t84*t90*t92*t93*zt*6.0D0+c1*kp*stress2*st
     &ress4*t41*t65*t88*t90*t92*t93*zt*3.0D0+c1*kp*stress3*stress4*t41*t
     &65*t89*t90*t92*t93*zt*3.0D0-c1*kp*stress4*stress6*t41*t65*t87*t90*
     &t92*t93*zt*6.0D0
      A0(1,5) = -stress5*t2*t112-Hmax*t27*t41*(bs5*2.0D0+stress5*2.0D0)*
     &3.0D0+Dc*Hmax*stress5*t27*t41*3.0D0+Dc*Hmax*t10*t27*t41*3.0D0+Hmax
     &*t27*t34*t44*t115*(3.0D0/2.0D0)+Hmax*t27*t36*t44*t115*(3.0D0/2.0D0
     &)+Hmax*t27*t38*t44*t115*(3.0D0/2.0D0)-Hmax*t27*t44*t51*t57*t115*(3
     &.0D0/4.0D0)-Hmax*t27*t44*t60*t66*t115*(3.0D0/4.0D0)-Hmax*t27*t44*t
     &61*t67*t115*(3.0D0/4.0D0)+c1*c2*stress5*t41*t53*t56*3.0D0+c1*c2*t4
     &1*t53*t56*t84*2.0D0-Dc*Hmax*stress4*t9*t27*t44*t115*(3.0D0/2.0D0)-
     &Dc*Hmax*stress5*t10*t27*t44*t115*(3.0D0/2.0D0)-Dc*Hmax*stress6*t11
     &*t27*t44*t115*(3.0D0/2.0D0)+Dc*Hmax*stress1*t27*t44*t51*t115*(3.0D
     &0/4.0D0)+Dc*Hmax*stress2*t27*t44*t66*t115*(3.0D0/4.0D0)+Dc*Hmax*st
     &ress3*t27*t44*t67*t115*(3.0D0/4.0D0)+Hmax*kt*stress5*t26*t34*t41*t
     &65*9.0D0+Hmax*kt*stress5*t26*t36*t41*t65*9.0D0+Hmax*kt*stress5*t26
     &*t38*t41*t65*9.0D0-Dc*Hmax*kt*t10*t20*t26*t41*t65*9.0D0-Hmax*kt*st
     &ress5*t26*t41*t51*t57*t65*(9.0D0/2.0D0)-Hmax*kt*stress5*t26*t41*t6
     &0*t65*t66*(9.0D0/2.0D0)-Hmax*kt*stress5*t26*t41*t61*t65*t67*(9.0D0
     &/2.0D0)+c1*c2*stress1*t44*t53*t56*t76*t115*(1.0D0/2.0D0)-c1*c2*str
     &ess4*t44*t53*t56*t81*t115-c1*c2*stress5*t44*t53*t56*t84*t115+c1*c2
     &*stress2*t44*t53*t56*t88*t115*(1.0D0/2.0D0)+c1*c2*stress3*t44*t53*
     &t56*t89*t115*(1.0D0/2.0D0)-c1*c2*stress6*t44*t53*t56*t87*t115-Dc*H
     &max*kt*stress4*stress5*t9*t26*t41*t65*9.0D0-Dc*Hmax*kt*stress5*str
     &ess6*t11*t26*t41*t65*9.0D0+Dc*Hmax*kt*stress1*stress5*t26*t41*t51*
     &t65*(9.0D0/2.0D0)+Dc*Hmax*kt*stress2*stress5*t26*t41*t65*t66*(9.0D
     &0/2.0D0)+Dc*Hmax*kt*stress3*stress5*t26*t41*t65*t67*(9.0D0/2.0D0)+
     &c1*c2*kp*t20*t41*t53*t56*t65*t84*6.0D0-c1*c2*kp*stress1*stress5*t4
     &1*t53*t56*t65*t76*3.0D0+c1*c2*kp*stress4*stress5*t41*t53*t56*t65*t
     &81*6.0D0-c1*c2*kp*stress2*stress5*t41*t53*t56*t65*t88*3.0D0-c1*c2*
     &kp*stress3*stress5*t41*t53*t56*t65*t89*3.0D0+c1*c2*kp*stress5*stre
     &ss6*t41*t53*t56*t65*t87*6.0D0-c1*kp*t20*t41*t65*t84*t90*t92*t93*zt
     &*6.0D0+c1*kp*stress1*stress5*t41*t65*t76*t90*t92*t93*zt*3.0D0-c1*k
     &p*stress4*stress5*t41*t65*t81*t90*t92*t93*zt*6.0D0+c1*kp*stress2*s
     &tress5*t41*t65*t88*t90*t92*t93*zt*3.0D0+c1*kp*stress3*stress5*t41*
     &t65*t89*t90*t92*t93*zt*3.0D0-c1*kp*stress5*stress6*t41*t65*t87*t90
     &*t92*t93*zt*6.0D0
      A0(1,6) = -stress6*t2*t112-Hmax*t27*t41*(bs6*2.0D0+stress6*2.0D0)*
     &3.0D0+Dc*Hmax*stress6*t27*t41*3.0D0+Dc*Hmax*t11*t27*t41*3.0D0+Hmax
     &*t27*t34*t44*t118*(3.0D0/2.0D0)+Hmax*t27*t36*t44*t118*(3.0D0/2.0D0
     &)+Hmax*t27*t38*t44*t118*(3.0D0/2.0D0)-Hmax*t27*t44*t51*t57*t118*(3
     &.0D0/4.0D0)-Hmax*t27*t44*t60*t66*t118*(3.0D0/4.0D0)-Hmax*t27*t44*t
     &61*t67*t118*(3.0D0/4.0D0)+c1*c2*stress6*t41*t53*t56*3.0D0+c1*c2*t4
     &1*t53*t56*t87*2.0D0-Dc*Hmax*stress4*t9*t27*t44*t118*(3.0D0/2.0D0)-
     &Dc*Hmax*stress5*t10*t27*t44*t118*(3.0D0/2.0D0)-Dc*Hmax*stress6*t11
     &*t27*t44*t118*(3.0D0/2.0D0)+Dc*Hmax*stress1*t27*t44*t51*t118*(3.0D
     &0/4.0D0)+Dc*Hmax*stress2*t27*t44*t66*t118*(3.0D0/4.0D0)+Dc*Hmax*st
     &ress3*t27*t44*t67*t118*(3.0D0/4.0D0)+Hmax*kt*stress6*t26*t34*t41*t
     &65*9.0D0+Hmax*kt*stress6*t26*t36*t41*t65*9.0D0+Hmax*kt*stress6*t26
     &*t38*t41*t65*9.0D0-Dc*Hmax*kt*t11*t22*t26*t41*t65*9.0D0-Hmax*kt*st
     &ress6*t26*t41*t51*t57*t65*(9.0D0/2.0D0)-Hmax*kt*stress6*t26*t41*t6
     &0*t65*t66*(9.0D0/2.0D0)-Hmax*kt*stress6*t26*t41*t61*t65*t67*(9.0D0
     &/2.0D0)+c1*c2*stress1*t44*t53*t56*t76*t118*(1.0D0/2.0D0)-c1*c2*str
     &ess4*t44*t53*t56*t81*t118-c1*c2*stress5*t44*t53*t56*t84*t118+c1*c2
     &*stress2*t44*t53*t56*t88*t118*(1.0D0/2.0D0)+c1*c2*stress3*t44*t53*
     &t56*t89*t118*(1.0D0/2.0D0)-c1*c2*stress6*t44*t53*t56*t87*t118-Dc*H
     &max*kt*stress4*stress6*t9*t26*t41*t65*9.0D0-Dc*Hmax*kt*stress5*str
     &ess6*t10*t26*t41*t65*9.0D0+Dc*Hmax*kt*stress1*stress6*t26*t41*t51*
     &t65*(9.0D0/2.0D0)+Dc*Hmax*kt*stress2*stress6*t26*t41*t65*t66*(9.0D
     &0/2.0D0)+Dc*Hmax*kt*stress3*stress6*t26*t41*t65*t67*(9.0D0/2.0D0)+
     &c1*c2*kp*t22*t41*t53*t56*t65*t87*6.0D0-c1*c2*kp*stress1*stress6*t4
     &1*t53*t56*t65*t76*3.0D0+c1*c2*kp*stress4*stress6*t41*t53*t56*t65*t
     &81*6.0D0+c1*c2*kp*stress5*stress6*t41*t53*t56*t65*t84*6.0D0-c1*c2*
     &kp*stress2*stress6*t41*t53*t56*t65*t88*3.0D0-c1*c2*kp*stress3*stre
     &ss6*t41*t53*t56*t65*t89*3.0D0-c1*kp*t22*t41*t65*t87*t90*t92*t93*zt
     &*6.0D0+c1*kp*stress1*stress6*t41*t65*t76*t90*t92*t93*zt*3.0D0-c1*k
     &p*stress4*stress6*t41*t65*t81*t90*t92*t93*zt*6.0D0-c1*kp*stress5*s
     &tress6*t41*t65*t84*t90*t92*t93*zt*6.0D0+c1*kp*stress2*stress6*t41*
     &t65*t88*t90*t92*t93*zt*3.0D0+c1*kp*stress3*stress6*t41*t65*t89*t90
     &*t92*t93*zt*3.0D0

 
       dpds(1,1:6)=A0(1,1:6)
       A0=0
 !********************************************************************************************************
 !********************************************************************************************************
      t0 = a1*(n2*(-x+1.0D0)**(n2-1.0D0)+n1*x**(n1-1.0D0))*(-1.0D0/2.0D0
     &)


 
 		dpdx=t0

!********************************************************************************************************
!********************************************************************************************************
!DEC$ FREEFORM


AAA=matmul(M,N)
BBB=matmul(dpds,M)
ccc=0
do ii=1,6
    ccc=ccc+dpds(1,ii)*AAA(ii,1)
enddo
!CCC=matmul(dpds,AAA)
!real*8 AAA(6,1),BBB(1,6),CCC,A0(6,6),M(6,6),N(6,1),dpds(1,6),dpdx


do ii=1,6
	do jj=1,6
		LLL(ii,jj)=AAA(ii,1)*BBB(1,jj)

		!write(*,*)"invjac",invjac(ii,jj),ii,jj
	enddo
enddo

ddsdde=M-LLL/(CCC-dpdx)



elseif (transformation == -1) then                  ! Reverse Transformaion
! Read(*,*),variable
!DEC$ NOFREEFORM
       t2 = stress1-stress2
      t3 = stress1-stress3
      t4 = stress2-stress3
      t5 = t2**2
      t6 = t5*(1.0D0/2.0D0)
      t7 = t3**2
      t8 = t7*(1.0D0/2.0D0)
      t9 = t4**2
      t10 = t9*(1.0D0/2.0D0)
      t11 = stress4**2
      t12 = t11*3.0D0
      t13 = stress5**2
      t14 = t13*3.0D0
      t15 = stress6**2
      t16 = t15*3.0D0
      t17 = t6+t8+t10+t12+t14+t16
      t18 = sqrt(t17)
      t19 = kp*t18
      t20 = exp(t19)
      t21 = c2*t20*zt
      t22 = t21+1.0D0
      t35 = stress1*2.0D0
      t23 = stress2+stress3-t35
      t24 = 1.0D0/sqrt(t17)
      t25 = Sa-Sm
      t26 = 1.0D0/t22
      t27 = c2**2
      t28 = kp*t18*2.0D0
      t29 = exp(t28)
      t30 = 1.0D0/t22**2
      t36 = stress2*2.0D0
      t31 = stress1+stress3-t36
      t34 = t25*x
      t32 = Sa-t34
      t37 = stress3*2.0D0
      t33 = stress1+stress2-t37
      t38 = v*2.0D0
      t39 = t38+2.0D0
      t40 = t32*t39
      A0(1,1) = Sa-t25*x+c1*c2*kp*lamdat_r1*t20*t23*t24*t26*(1.0D0/2.0D0
     &)-c1*kp*lamdat_r1*t23*t24*t27*t29*t30*zt*(1.0D0/2.0D0)
      A0(1,2) = -t32*v+c1*c2*kp*lamdat_r1*t20*t24*t26*t31*(1.0D0/2.0D0)-
     &c1*kp*lamdat_r1*t24*t27*t29*t30*t31*zt*(1.0D0/2.0D0)
      A0(1,3) = -t32*v+c1*c2*kp*lamdat_r1*t20*t24*t26*t33*(1.0D0/2.0D0)-
     &c1*kp*lamdat_r1*t24*t27*t29*t30*t33*zt*(1.0D0/2.0D0)
      A0(1,4) = c1*c2*kp*lamdat_r1*stress4*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r1*stress4*t24*t27*t29*t30*zt*3.0D0
      A0(1,5) = c1*c2*kp*lamdat_r1*stress5*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r1*stress5*t24*t27*t29*t30*zt*3.0D0
      A0(1,6) = c1*c2*kp*lamdat_r1*stress6*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r1*stress6*t24*t27*t29*t30*zt*3.0D0
      A0(2,1) = -t32*v+c1*c2*kp*lamdat_r2*t20*t23*t24*t26*(1.0D0/2.0D0)-
     &c1*kp*lamdat_r2*t23*t24*t27*t29*t30*zt*(1.0D0/2.0D0)
      A0(2,2) = Sa-t34+c1*c2*kp*lamdat_r2*t20*t24*t26*t31*(1.0D0/2.0D0)-
     &c1*kp*lamdat_r2*t24*t27*t29*t30*t31*zt*(1.0D0/2.0D0)
      A0(2,3) = -t32*v+c1*c2*kp*lamdat_r2*t20*t24*t26*t33*(1.0D0/2.0D0)-
     &c1*kp*lamdat_r2*t24*t27*t29*t30*t33*zt*(1.0D0/2.0D0)
      A0(2,4) = c1*c2*kp*lamdat_r2*stress4*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r2*stress4*t24*t27*t29*t30*zt*3.0D0
      A0(2,5) = c1*c2*kp*lamdat_r2*stress5*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r2*stress5*t24*t27*t29*t30*zt*3.0D0
      A0(2,6) = c1*c2*kp*lamdat_r2*stress6*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r2*stress6*t24*t27*t29*t30*zt*3.0D0
      A0(3,1) = -t32*v+c1*c2*kp*lamdat_r3*t20*t23*t24*t26*(1.0D0/2.0D0)-
     &c1*kp*lamdat_r3*t23*t24*t27*t29*t30*zt*(1.0D0/2.0D0)
      A0(3,2) = -t32*v+c1*c2*kp*lamdat_r3*t20*t24*t26*t31*(1.0D0/2.0D0)-
     &c1*kp*lamdat_r3*t24*t27*t29*t30*t31*zt*(1.0D0/2.0D0)
      A0(3,3) = Sa-t34+c1*c2*kp*lamdat_r3*t20*t24*t26*t33*(1.0D0/2.0D0)-
     &c1*kp*lamdat_r3*t24*t27*t29*t30*t33*zt*(1.0D0/2.0D0)
      A0(3,4) = c1*c2*kp*lamdat_r3*stress4*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r3*stress4*t24*t27*t29*t30*zt*3.0D0
      A0(3,5) = c1*c2*kp*lamdat_r3*stress5*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r3*stress5*t24*t27*t29*t30*zt*3.0D0
      A0(3,6) = c1*c2*kp*lamdat_r3*stress6*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r3*stress6*t24*t27*t29*t30*zt*3.0D0
      A0(4,1) = c1*c2*kp*lamdat_r4*t20*t23*t24*t26*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r4*t23*t24*t27*t29*t30*zt*(1.0D0/2.0D0)
      A0(4,2) = c1*c2*kp*lamdat_r4*t20*t24*t26*t31*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r4*t24*t27*t29*t30*t31*zt*(1.0D0/2.0D0)
      A0(4,3) = c1*c2*kp*lamdat_r4*t20*t24*t26*t33*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r4*t24*t27*t29*t30*t33*zt*(1.0D0/2.0D0)
      A0(4,4) = t40-c1*c2*kp*lamdat_r4*stress4*t20*t24*t26*3.0D0+c1*kp*l
     &amdat_r4*stress4*t24*t27*t29*t30*zt*3.0D0
      A0(4,5) = c1*c2*kp*lamdat_r4*stress5*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r4*stress5*t24*t27*t29*t30*zt*3.0D0
      A0(4,6) = c1*c2*kp*lamdat_r4*stress6*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r4*stress6*t24*t27*t29*t30*zt*3.0D0
      A0(5,1) = c1*c2*kp*lamdat_r5*t20*t23*t24*t26*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r5*t23*t24*t27*t29*t30*zt*(1.0D0/2.0D0)
      A0(5,2) = c1*c2*kp*lamdat_r5*t20*t24*t26*t31*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r5*t24*t27*t29*t30*t31*zt*(1.0D0/2.0D0)
      A0(5,3) = c1*c2*kp*lamdat_r5*t20*t24*t26*t33*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r5*t24*t27*t29*t30*t33*zt*(1.0D0/2.0D0)
      A0(5,4) = c1*c2*kp*lamdat_r5*stress4*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r5*stress4*t24*t27*t29*t30*zt*3.0D0
      A0(5,5) = t40-c1*c2*kp*lamdat_r5*stress5*t20*t24*t26*3.0D0+c1*kp*l
     &amdat_r5*stress5*t24*t27*t29*t30*zt*3.0D0
      A0(5,6) = c1*c2*kp*lamdat_r5*stress6*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r5*stress6*t24*t27*t29*t30*zt*3.0D0
      A0(6,1) = c1*c2*kp*lamdat_r6*t20*t23*t24*t26*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r6*t23*t24*t27*t29*t30*zt*(1.0D0/2.0D0)
      A0(6,2) = c1*c2*kp*lamdat_r6*t20*t24*t26*t31*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r6*t24*t27*t29*t30*t31*zt*(1.0D0/2.0D0)
      A0(6,3) = c1*c2*kp*lamdat_r6*t20*t24*t26*t33*(1.0D0/2.0D0)-c1*kp*l
     &amdat_r6*t24*t27*t29*t30*t33*zt*(1.0D0/2.0D0)
      A0(6,4) = c1*c2*kp*lamdat_r6*stress4*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r6*stress4*t24*t27*t29*t30*zt*3.0D0
      A0(6,5) = c1*c2*kp*lamdat_r6*stress5*t20*t24*t26*(-3.0D0)+c1*kp*la
     &mdat_r6*stress5*t24*t27*t29*t30*zt*3.0D0
      A0(6,6) = t40-c1*c2*kp*lamdat_r6*stress6*t20*t24*t26*3.0D0+c1*kp*l
     &amdat_r6*stress6*t24*t27*t29*t30*zt*3.0D0



 
 	   M_INV=A0
 	   A0=0
 	   CALL inverse(M_INV,M,6)
 !********************************************************************************************************
 !********************************************************************************************************
      t2 = Sa-Sm
      t3 = stress1-stress2
      t4 = stress1-stress3
      t5 = stress2-stress3
      t6 = t3**2
      t7 = t6*(1.0D0/2.0D0)
      t8 = t4**2
      t9 = t8*(1.0D0/2.0D0)
      t10 = t5**2
      t11 = t10*(1.0D0/2.0D0)
      t12 = stress4**2
      t13 = t12*3.0D0
      t14 = stress5**2
      t15 = t14*3.0D0
      t16 = stress6**2
      t17 = t16*3.0D0
      t18 = t7+t9+t11+t13+t15+t17
      t19 = sqrt(t18)
      t20 = kp*t19
      t21 = exp(t20)
      t22 = stress3*t2*v
      t23 = c2*t21*zt
      t24 = t23+1.0D0
      t25 = 1.0D0/t24
      t26 = stress1*t2*v
      t27 = stress2*t2*v
      t28 = v*2.0D0
      t29 = t28+2.0D0
      A0(1,1) = lamdat_r1+t22+t27-stress1*t2-c1*c2*lamdat_r1*t21*t25
      A0(2,1) = lamdat_r2+t22+t26-stress2*t2-c1*c2*lamdat_r2*t21*t25
      A0(3,1) = lamdat_r3+t26+t27-stress3*t2-c1*c2*lamdat_r3*t21*t25
      A0(4,1) = lamdat_r4-stress4*t2*t29-c1*c2*lamdat_r4*t21*t25
      A0(5,1) = lamdat_r5-stress5*t2*t29-c1*c2*lamdat_r5*t21*t25
      A0(6,1) = lamdat_r6-stress6*t2*t29-c1*c2*lamdat_r6*t21*t25

	  
       N(1:6,1)=A0(1:6,1)
 	   A0=0
 !********************************************************************************************************
 !********************************************************************************************************
      t2 = Sa-Sm
      t3 = stress1-stress2
      t4 = stress1-stress3
      t5 = stress2-stress3
      t6 = t3**2
      t7 = t6*(1.0D0/2.0D0)
      t8 = t4**2
      t9 = t8*(1.0D0/2.0D0)
      t10 = t5**2
      t11 = t10*(1.0D0/2.0D0)
      t12 = stress4**2
      t13 = t12*3.0D0
      t14 = stress5**2
      t15 = t14*3.0D0
      t16 = stress6**2
      t17 = t16*3.0D0
      t18 = t7+t9+t11+t13+t15+t17
      t19 = sqrt(t18)
      t20 = kp*t19
      t21 = exp(t20)
      t22 = c2*t21*zt
      t23 = t22+1.0D0
      t24 = 1.0D0/t23
      t27 = stress1*2.0D0
      t25 = stress2+stress3-t27
      t26 = 1.0D0/sqrt(t18)
      t28 = c2**2
      t29 = kp*t19*2.0D0
      t30 = exp(t29)
      t31 = 1.0D0/t23**2
      t33 = stress2*2.0D0
      t32 = stress1+stress3-t33
      t35 = stress3*2.0D0
      t34 = stress1+stress2-t35
      t36 = v*2.0D0
      t37 = t36+2.0D0
      A0(1,1) = -lamdat_r1+stress1*t2-stress2*t2*v-stress3*t2*v+c1*c2*la
     &mdat_r1*t21*t24-c1*c2*kp*lamdat_r1*stress1*t21*t24*t25*t26*(1.0D0/
     &2.0D0)-c1*c2*kp*lamdat_r2*stress2*t21*t24*t25*t26*(1.0D0/2.0D0)-c1
     &*c2*kp*lamdat_r3*stress3*t21*t24*t25*t26*(1.0D0/2.0D0)-c1*c2*kp*la
     &mdat_r4*stress4*t21*t24*t25*t26*(1.0D0/2.0D0)-c1*c2*kp*lamdat_r5*s
     &tress5*t21*t24*t25*t26*(1.0D0/2.0D0)-c1*c2*kp*lamdat_r6*stress6*t2
     &1*t24*t25*t26*(1.0D0/2.0D0)+c1*kp*lamdat_r1*stress1*t25*t26*t28*t3
     &0*t31*zt*(1.0D0/2.0D0)+c1*kp*lamdat_r2*stress2*t25*t26*t28*t30*t31
     &*zt*(1.0D0/2.0D0)+c1*kp*lamdat_r3*stress3*t25*t26*t28*t30*t31*zt*(
     &1.0D0/2.0D0)+c1*kp*lamdat_r4*stress4*t25*t26*t28*t30*t31*zt*(1.0D0
     &/2.0D0)+c1*kp*lamdat_r5*stress5*t25*t26*t28*t30*t31*zt*(1.0D0/2.0D
     &0)+c1*kp*lamdat_r6*stress6*t25*t26*t28*t30*t31*zt*(1.0D0/2.0D0)
      A0(1,2) = -lamdat_r2+stress2*t2-stress1*t2*v-stress3*t2*v+c1*c2*la
     &mdat_r2*t21*t24-c1*c2*kp*lamdat_r1*stress1*t21*t24*t26*t32*(1.0D0/
     &2.0D0)-c1*c2*kp*lamdat_r2*stress2*t21*t24*t26*t32*(1.0D0/2.0D0)-c1
     &*c2*kp*lamdat_r3*stress3*t21*t24*t26*t32*(1.0D0/2.0D0)-c1*c2*kp*la
     &mdat_r4*stress4*t21*t24*t26*t32*(1.0D0/2.0D0)-c1*c2*kp*lamdat_r5*s
     &tress5*t21*t24*t26*t32*(1.0D0/2.0D0)-c1*c2*kp*lamdat_r6*stress6*t2
     &1*t24*t26*t32*(1.0D0/2.0D0)+c1*kp*lamdat_r1*stress1*t26*t28*t30*t3
     &1*t32*zt*(1.0D0/2.0D0)+c1*kp*lamdat_r2*stress2*t26*t28*t30*t31*t32
     &*zt*(1.0D0/2.0D0)+c1*kp*lamdat_r3*stress3*t26*t28*t30*t31*t32*zt*(
     &1.0D0/2.0D0)+c1*kp*lamdat_r4*stress4*t26*t28*t30*t31*t32*zt*(1.0D0
     &/2.0D0)+c1*kp*lamdat_r5*stress5*t26*t28*t30*t31*t32*zt*(1.0D0/2.0D
     &0)+c1*kp*lamdat_r6*stress6*t26*t28*t30*t31*t32*zt*(1.0D0/2.0D0)
      A0(1,3) = -lamdat_r3+stress3*t2-stress1*t2*v-stress2*t2*v+c1*c2*la
     &mdat_r3*t21*t24-c1*c2*kp*lamdat_r1*stress1*t21*t24*t26*t34*(1.0D0/
     &2.0D0)-c1*c2*kp*lamdat_r2*stress2*t21*t24*t26*t34*(1.0D0/2.0D0)-c1
     &*c2*kp*lamdat_r3*stress3*t21*t24*t26*t34*(1.0D0/2.0D0)-c1*c2*kp*la
     &mdat_r4*stress4*t21*t24*t26*t34*(1.0D0/2.0D0)-c1*c2*kp*lamdat_r5*s
     &tress5*t21*t24*t26*t34*(1.0D0/2.0D0)-c1*c2*kp*lamdat_r6*stress6*t2
     &1*t24*t26*t34*(1.0D0/2.0D0)+c1*kp*lamdat_r1*stress1*t26*t28*t30*t3
     &1*t34*zt*(1.0D0/2.0D0)+c1*kp*lamdat_r2*stress2*t26*t28*t30*t31*t34
     &*zt*(1.0D0/2.0D0)+c1*kp*lamdat_r3*stress3*t26*t28*t30*t31*t34*zt*(
     &1.0D0/2.0D0)+c1*kp*lamdat_r4*stress4*t26*t28*t30*t31*t34*zt*(1.0D0
     &/2.0D0)+c1*kp*lamdat_r5*stress5*t26*t28*t30*t31*t34*zt*(1.0D0/2.0D
     &0)+c1*kp*lamdat_r6*stress6*t26*t28*t30*t31*t34*zt*(1.0D0/2.0D0)
      A0(1,4) = -lamdat_r4+stress4*t2*t37+c1*c2*lamdat_r4*t21*t24+c1*c2*
     &kp*lamdat_r4*t12*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r1*stress1*stre
     &ss4*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r2*stress2*stress4*t21*t24*t
     &26*3.0D0+c1*c2*kp*lamdat_r3*stress3*stress4*t21*t24*t26*3.0D0+c1*c
     &2*kp*lamdat_r5*stress4*stress5*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r
     &6*stress4*stress6*t21*t24*t26*3.0D0-c1*kp*lamdat_r4*t12*t26*t28*t3
     &0*t31*zt*3.0D0-c1*kp*lamdat_r1*stress1*stress4*t26*t28*t30*t31*zt*
     &3.0D0-c1*kp*lamdat_r2*stress2*stress4*t26*t28*t30*t31*zt*3.0D0-c1*
     &kp*lamdat_r3*stress3*stress4*t26*t28*t30*t31*zt*3.0D0-c1*kp*lamdat
     &_r5*stress4*stress5*t26*t28*t30*t31*zt*3.0D0-c1*kp*lamdat_r6*stres
     &s4*stress6*t26*t28*t30*t31*zt*3.0D0
      A0(1,5) = -lamdat_r5+stress5*t2*t37+c1*c2*lamdat_r5*t21*t24+c1*c2*
     &kp*lamdat_r5*t14*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r1*stress1*stre
     &ss5*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r2*stress2*stress5*t21*t24*t
     &26*3.0D0+c1*c2*kp*lamdat_r3*stress3*stress5*t21*t24*t26*3.0D0+c1*c
     &2*kp*lamdat_r4*stress4*stress5*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r
     &6*stress5*stress6*t21*t24*t26*3.0D0-c1*kp*lamdat_r5*t14*t26*t28*t3
     &0*t31*zt*3.0D0-c1*kp*lamdat_r1*stress1*stress5*t26*t28*t30*t31*zt*
     &3.0D0-c1*kp*lamdat_r2*stress2*stress5*t26*t28*t30*t31*zt*3.0D0-c1*
     &kp*lamdat_r3*stress3*stress5*t26*t28*t30*t31*zt*3.0D0-c1*kp*lamdat
     &_r4*stress4*stress5*t26*t28*t30*t31*zt*3.0D0-c1*kp*lamdat_r6*stres
     &s5*stress6*t26*t28*t30*t31*zt*3.0D0
      A0(1,6) = -lamdat_r6+stress6*t2*t37+c1*c2*lamdat_r6*t21*t24+c1*c2*
     &kp*lamdat_r6*t16*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r1*stress1*stre
     &ss6*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r2*stress2*stress6*t21*t24*t
     &26*3.0D0+c1*c2*kp*lamdat_r3*stress3*stress6*t21*t24*t26*3.0D0+c1*c
     &2*kp*lamdat_r4*stress4*stress6*t21*t24*t26*3.0D0+c1*c2*kp*lamdat_r
     &5*stress5*stress6*t21*t24*t26*3.0D0-c1*kp*lamdat_r6*t16*t26*t28*t3
     &0*t31*zt*3.0D0-c1*kp*lamdat_r1*stress1*stress6*t26*t28*t30*t31*zt*
     &3.0D0-c1*kp*lamdat_r2*stress2*stress6*t26*t28*t30*t31*zt*3.0D0-c1*
     &kp*lamdat_r3*stress3*stress6*t26*t28*t30*t31*zt*3.0D0-c1*kp*lamdat
     &_r4*stress4*stress6*t26*t28*t30*t31*zt*3.0D0-c1*kp*lamdat_r5*stres
     &s5*stress6*t26*t28*t30*t31*zt*3.0D0

 
       dpds(1,1:6)=A0(1,1:6)
       A0=0
 !********************************************************************************************************
 !********************************************************************************************************
      t0 = a2*(n4*(-x+1.0D0)**(n4-1.0D0)+n3*x**(n3-1.0D0))*(1.0D0/2.0D0)
 
 		dpdx=t0

!********************************************************************************************************
!********************************************************************************************************






!DEC$ FREEFORM


AAA=matmul(M,N)
BBB=matmul(dpds,M)
ccc=0
do ii=1,6
    ccc=ccc+dpds(1,ii)*AAA(ii,1)
enddo
!CCC=matmul(dpds,AAA)
!real*8 AAA(6,1),BBB(1,6),CCC,A0(6,6),M(6,6),N(6,1),dpds(1,6),dpdx


do ii=1,6
	do jj=1,6
		LLL(ii,jj)=AAA(ii,1)*BBB(1,jj)

		!write(*,*)"invjac",invjac(ii,jj),ii,jj
	enddo
enddo

ddsdde=M-LLL/(CCC-dpdx)




endif




end

subroutine  Accumulated_zt(PARAM,NPARAM,VAR,NVAR)
implicit real*8(t)


!DEFINITION OF  INPUT MODEL PARAMETERS
REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
integer NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


! Changing state variables
real*8 A0(4,1)
real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 et_initial1,et_initial2,et_initial3,et_initial4,et_initial5,et_initial6


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)



       et_initial1=0
       et_initial2=0
       et_initial3=0
       et_initial4=0
       et_initial5=0
       et_initial6=0

       e1=e(1)
       e2=e(2)
       e3=e(3)
       e4=e(4)
       e5=e(5)
       e6=e(6)
	   if (maxval(e)<=1.0e-8)then
          e1 = 1.0e-8_8
      end if

      
      et1=et(1)
      et2=et(2)
      et3=et(3)
      et4=et(4)
      et5=et(5)
      et6=et(6)	  
      if (maxval(et)<=1.0e-8)then
          et1 = 1.0e-8_8
      end if
	  
	  
   ep1=ep(1)
   ep2=ep(2)
   ep3=ep(3)
   ep4=ep(4)
   ep5=ep(5)
   if(maxval(abs(ep))<=1.0e-10)then
       ep1=1.0e-10_8
    end if

	  !Backstress Calculation   
   bs1=backstress(1)
   bs2=backstress(2)
   bs3=backstress(3)
   bs4=backstress(4)
   bs5=backstress(5)
   bs6=backstress(6)


!DEC$ NOFREEFORM
 !******************************************************************************
      t3 = v+1.0D0
      t4 = 1.0D0/t3
      t5 = Sa*x
      t6 = Sm*x
      t7 = Sa-t5+t6
      t8 = 1.0D0/t7
      t30 = -e4+ep4+et4
      t31 = t4*t8*t30*(1.0D0/2.0D0)
      t2 = bs4-t31
      t34 = -e5+ep5+et5
      t35 = t4*t8*t34*(1.0D0/2.0D0)
      t9 = bs5-t35
      t38 = -e6+ep6+et6
      t39 = t4*t8*t38*(1.0D0/2.0D0)
      t10 = bs6-t39
      t11 = v-1.0D0
      t12 = v**2
      t13 = t12*2.0D0
      t14 = t13+v-1.0D0
      t15 = 1.0D0/t14
      t16 = Tn-too
      t17 = alpha*t16
      t18 = -e1+ep1+et1+t17
      t19 = -e2+ep2+et2+t17
      t21 = t8*t11*t15*t18
      t22 = t8*t15*t18*v
      t25 = t8*t11*t15*t19
      t27 = t8*t15*t19*v
      t20 = bs1-bs2-t21-t22+t25+t27
      t23 = -e3+ep3+et3+t17
      t26 = t8*t11*t15*t23
      t28 = t8*t15*t23*v
      t24 = bs1-bs3-t21-t22+t26+t28
      t29 = bs2-bs3-t25+t26-t27+t28
      t32 = t2**2
      t33 = t32*3.0D0
      t36 = t9**2
      t37 = t36*3.0D0
      t40 = t10**2
      t41 = t40*3.0D0
      t42 = t20**2
      t43 = t42*(1.0D0/2.0D0)
      t44 = t24**2
      t45 = t44*(1.0D0/2.0D0)
      t46 = t29**2
      t47 = t46*(1.0D0/2.0D0)
      t48 = t33+t37+t41+t43+t45+t47
      t49 = sqrt(t48)
      t53 = kt*t49
      t50 = exp(-t53)
      t51 = t50-1.0D0
      t52 = x-xo
      A0(1,1) = xdo-t51*t52
      A0(2,1) = -t51*t52
      A0(3,1) = zto+abs(t52)
      A0(4,1) = ztdo+abs(t51*t52)


       xd = A0(1,1)
       dxd=A0(2,1) 
       zt  = A0(3,1) 
       ztd=A0(4,1) 
!DEC$ FREEFORM	   
       	  
CALL VAR_ASSIGNMENTS(2,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)


!
!
 end

 subroutine Backstress_Update(PARAM,NPARAM,VAR,NVAR)
 
implicit real*8 (t)

!DEFINITION OF  INPUT MODEL PARAMETERS
REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
integer NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


! Changing state variables
real*8 A0(6,1)


real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6
real*8 eto1,eto2,eto3,eto4,eto5,eto6,epo1,epo2,epo3,epo4,epo5,epo6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)



   e1=e(1)
   e2=e(2)
   e3=e(3)
   e4=e(4)
   e5=e(5)
   e6=e(6)
    if(maxval(abs(e))<=1.0e-8)then
       e1=1.0e-8_8
    end if

   
   et1=et(1)
   et2=et(2)
   et3=et(3)
   et4=et(4)
   et5=et(5)
   et6=et(6)
   if(maxval(abs(et))<=1.0e-9)then
       et1=1.0e-9_8
    end if
   

   et_tr1=et_tr(1)
   et_tr2=et_tr(2)
   et_tr3=et_tr(3)
   et_tr4=et_tr(4)
   et_tr5=et_tr(5)
   et_tr6=et_tr(6)
   if(maxval(abs(et_tr))<=1.0e-10)then
       et_tr1=1.0e-10_8
    end if

   ep1=ep(1)
   ep2=ep(2)
   ep3=ep(3)
   ep4=ep(4)
   ep5=ep(5)
   if(maxval(abs(ep))<=1.0e-10)then
       ep1=1.0e-10_8
    end if

   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)
   
   
   !Backstress Calculation   
   bs1=backstress(1)
   bs2=backstress(2)
   bs3=backstress(3)
   bs4=backstress(4)
   bs5=backstress(5)
   bs6=backstress(6)

      t3 = v+1.0D0
      t4 = 1.0D0/t3
      t5 = Sa*x
      t6 = Sm*x
      t7 = Sa-t5+t6
      t8 = 1.0D0/t7
      t32 = -e4+ep4+et4
      t33 = t4*t8*t32*(1.0D0/2.0D0)
      t2 = bs4-t33
      t36 = -e5+ep5+et5
      t37 = t4*t8*t36*(1.0D0/2.0D0)
      t9 = bs5-t37
      t40 = -e6+ep6+et6
      t41 = t4*t8*t40*(1.0D0/2.0D0)
      t10 = bs6-t41
      t11 = v-1.0D0
      t12 = v**2
      t13 = t12*2.0D0
      t14 = t13+v-1.0D0
      t15 = 1.0D0/t14
      t16 = Tn-too
      t17 = alpha*t16
      t18 = -e1+ep1+et1+t17
      t19 = -e2+ep2+et2+t17
      t21 = t8*t11*t15*t18
      t22 = t8*t15*t18*v
      t25 = t8*t11*t15*t19
      t27 = t8*t15*t19*v
      t20 = bs1-bs2-t21-t22+t25+t27
      t23 = -e3+ep3+et3+t17
      t26 = t8*t11*t15*t23
      t28 = t8*t15*t23*v
      t24 = bs1-bs3-t21-t22+t26+t28
      t29 = bs2-bs3-t25+t26-t27+t28
      t51 = l1*ztd
      t30 = exp(-t51)
      t31 = t30-1.0D0
      t34 = t2**2
      t35 = t34*6.0D0
      t38 = t9**2
      t39 = t38*6.0D0
      t42 = t10**2
      t43 = t42*6.0D0
      t44 = t20**2
      t45 = t24**2
      t46 = t29**2
      t47 = t35+t39+t43+t44+t45+t46
      t48 = abs(t47)
      t49 = t48*(1.0D0/2.0D0)
      t50 = 1.0D0/sqrt(t49)
      A0(1,1) = -Eb*t31*t50*(bs1-t21+t27+t28)
      A0(2,1) = -Eb*t31*t50*(bs2+t22-t25+t28)
      A0(3,1) = -Eb*t31*t50*(bs3+t22-t26+t27)
      A0(4,1) = -Eb*t2*t31*t50
      A0(5,1) = -Eb*t9*t31*t50
      A0(6,1) = -Eb*t10*t31*t50



      backstress(1)=A0(1,1)
      backstress(2)=A0(2,1)
      backstress(3)=A0(3,1)
      backstress(4)=A0(4,1)
      backstress(5)=A0(5,1)
      backstress(6)=A0(6,1)


CALL VAR_ASSIGNMENTS(2,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

 end subroutine
 
 
subroutine SDV_UPDATE(PARAM,NPARAM,VAR,NVAR,STATEV,NSTATV)
implicit none


!DEFINITION OF  INPUT MODEL PARAMETERS
REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,kp
INTEGER NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),ep(6),epo(6),Lamdat_r(6),et_tr(6),backstress(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


INTEGER NSTATV
REAL*8 STATEV(NSTATV)


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,Eb,C1,C2,L1,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL,kp) 
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,xd,xdo,dxd,dxdo,zt,zto,ztd,ztdo,flag_fwd,flag_rev,Transformation,et,eto,ep,epo,et_tr,lamdat_r,backstress,e,eo,RPLC,Bound_Reached,NR_Convergence,too)



	 statev(1)           =      x
     statev(2)           =      xd
     statev(3)           =      zt
     statev(4)           =      ztd
     statev(5)           =      dxd
     statev(6)           =      real(flag_fwd)
     statev(7)           =      real(flag_rev)
     statev(8)           =      real(transformation)
     statev(9:14)      =      et(1:6)
     statev(15:20)    =      ep(1:6)
     statev(21:26)    =      et_tr(1:6)
	 statev(27:32)    =      Lamdat_r(1:6)
	 statev(33:38)	  =      backstress(1:6)
     statev(39)         =      RPLC
	 statev(40)         =      BOUND_REACHED
	 statev(41)         =      NR_CONVERGENCE
	 statev(42)         =      TOO 
	 
end

!
!! !######################## Utlity Subroutines ##############################
!! !##################################################################
        
subroutine check_nan1(A,K,OUTPUT)  ! A IS THE INPUT VECTOR , K IS THE DIMENSION OF THE VECTOR,OUTPUT TAKES 0 VALUE IF THERE IS NO NAN AND 1 IF THERE IS A NAN
IMPLICIT NONE
    real*8 A(K)
    integer K,i,OUTPUT
    logical c(K),d
    
OUTPUT=0
i=1
do while(i<=K.and.d/=.true.)
c(i)=isnan(a(i))
d=c(i)
IF(D==.TRUE.)THEN
    OUTPUT=1
ENDIF

i=i+1
enddo
    
    
    ENDSUBROUTINE
    
    
subroutine check_nan2(A,k,OUTPUT)  ! A IS THE INPUT MATRIX , K IS THE DIMENSION OF THE SQUARE MATRIX,OUTPUT TAKES 0 VALUE IF THERE IS NO NAN AND 1 IF THERE IS A NAN
IMPLICIT NONE
    real*8 A(K,K)
    integer K,i,J,OUTPUT
    logical c(K,K),d
    

    
OUTPUT=0
i=1
J=1
do while(i<=K.and.d/=.true.)
c(I,J)=isnan(a(i,J))
d=c(i,J)

IF(D==.TRUE.)THEN
    OUTPUT=1
ENDIF

J=J+1
IF(J==K+1)THEN
    J=1
    I=I+1
    ENDIF
enddo
    

    
    ENDSUBROUTINE
    
    
    subroutine check_inf1(A,K,OUTPUT) ! A IS THE INPUT ARRAY , K IS THE DIMENSION OF THE ARRAY,OUTPUT TAKES 0 VALUE IF THERE IS NO NAN AND 1 IF THERE IS A NAN
IMPLICIT NONE
    real*8 A(K),inf_no,b
    integer K,i,OUTPUT

b=1
inf_no=b/0
OUTPUT=0
i=1
do while(i<=K.and.output/=1)

      if(abs(A(i))<abs(inf_no))then
          output=0
      else
          output=1
          endif

i=i+1
enddo
    
    
    ENDSUBROUTINE
    
    
subroutine check_inf2(A,k,OUTPUT) ! A IS THE INPUT MATRIX , K IS THE DIMENSION OF THE SQUARE MATRIX,OUTPUT TAKES 0 VALUE IF THERE IS NO INF AND 1 IF THERE IS AN INF
IMPLICIT NONE
    real*8 A(K,K),b,inf_no
    integer K,i,J,OUTPUT

    b=1
inf_no=b/0
OUTPUT=0
i=1
J=1
do while(i<=K.and.output/=1)

      if(abs(A(i,j))<abs(inf_no))then
          output=0
      else
          output=1
          endif

J=J+1
IF(J==K+1)THEN
    J=1
    I=I+1
    ENDIF
enddo
    
    
    ENDSUBROUTINE
!
!DEC$ NOFREEFORM
!! !######################################################
!! !######################################################
!! !######################################################
!! !######################################################
!! !######################################################
!! !######################################################
!
!****************************************************************************
!   This is the key subroutine is to find the incremental rotation matrix RLOG and log spin tensor
!   SPINLOG defined in large deformation.
      subroutine FIND_LOG(RLOG,FN,FN1)
      !
      ! Definition of Variables
      !
      ! RLOG(3.3)         :  incremental log roation matrix_  output
      ! SPINLOG(3,3)    :  incremental log spin tensor, used to find RLOG(3,3)
      ! FN (3,3)            : Deformation Gradient correspond to X(n) configuration, input as DFGRD0
      ! FN_INV (3,3)     : Inverse of FN correspond to X(n) configuration
      ! FN1(3,3)           : Deformation Gradient correspond to X(n+1) configuration input as DFGRD1
      ! Delta_F(3,3)             :  the incremental deformation gradient can transfer X(n) to X(n+1)
      ! Delta_F_INV(3,3)      : inverse of Delta_F(3,3)
      !	det_Del_F                : the determinant of Del_F(3,3); det_Del_F > 0
      !
      !	B(3,3):              the incremental left Cauchy-Green tensor;
      !	VL(3,3):            Velocity gradient
      !	D(3,3):              Streching tensor
      ! W(3,3):             Spin tensor
      ! B_P(3)              :the principal stretches of B(3,3)
      !	B_sqrtval(3)      :the squares roots of the B_eigval(3)
      !	B_E(3,3)           :matrix of eigenvectors of B(3,3)
      !	E(3,3)               :the logarithmic strain tensor; output
      !R_N(3,3)           :Rotation matrix at N configuration
      !R_N1(3,3)         :Rotation matrix at N+1 configuration
      !F_INV(3,3)        :the inverse of 0.5D0*(FN+FN1)
      !q(3),q0, q_star, norm                  :Variable  associated with exponetional MAP
      ! istat:        success flag, istat=0 for a failed attempt; output
 
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8  det_Del_F , det_temp, det_FN, q0, q_star, norm, N_LOG
 
 !     DIMENSION STATEV(NSTATEV)                               ! input R_N
      DIMENSION RLOG(3,3), SPINLOG(3,3)                   ! input as initial value, output after call
      DIMENSION FN(3,3), FN1(3,3), FN_INV(3,3),FN_AVG(3,3),Delta_F(3,3) 
      DIMENSION B_VEC(3,3), B_PV(3), lgtt1(3,3),N_LOG(3,3), EXP_M(3,3)
      DIMENSION q(3),FN_AVG_INV(3,3)
 
 	  DIMENSION rlog1(3,3), rlog2(3,3)
      DIMENSION VL(3,3), D(3,3), D1(3,3),W(3,3),B(3,3),B1(3,3),B2(3,3)
 
 
      TOL=0.1e-30
 
      N_LOG=0.d0
 	  FN_AVG=FN1
 !      FN_AVG=0.5d0*(FN+FN1)
      Delta_F=FN1-FN
      call inverse(FN_AVG,FN_AVG_INV,3)
      VL=matmul(Delta_F,FN_AVG_INV)
      D=0.5d0*(VL+transpose(VL))
      D1=D
      W=0.5d0*(VL-transpose(VL))
      B=matmul(FN_AVG,transpose(FN_AVG))
      B1=B
      B2=B
      call MAT_P(B1,3,3,B_PV,B_VEC)
      call FIND_N_LOG(B_PV,B2,D1,N_LOG)
      SPINLOG = W + N_LOG                       !logarithmic spin
 
 
 ! ****** Exponential Map to calculate ΔR********
      EXP_M=0.0D0
      q=0.0D0
      q(1)=-SPINLOG(2,3)
      q(2)= SPINLOG(1,3)
      q(3)=-SPINLOG(1,2)
      norm = sqrt ( q(1)**2 + q(2)**2 + q(3)**2)
 !         WRITE(*,*)'This is norm'
 !         WRITE(*,*),norm
       q0      = cos(0.5D0*norm)
       q_star  = sin (0.5D0*norm)
 !
      IF (abs(q_star) .GT. TOL) THEN
           q_star=0.5D0*sin(0.5D0*norm)/(0.5D0*norm)
      ELSE
          q_star=0.5D0*(1.0D0-norm**2/24.0D0+norm**4/1920.0D0)
      ENDIF
      q=q_star*q
      EXP_M(1,1)=q0**2+q(1)**2-0.5D0
      EXP_M(2,2)=q0**2+q(2)**2-0.5D0
      EXP_M(3,3)=q0**2+q(3)**2-0.5D0
      EXP_M(1,2)=q(1)*q(2)-q(3)*q0
      EXP_M(1,3)=q(1)*q(3)+q(2)*q0
      EXP_M(2,1)=q(2)*q(1)+q(3)*q0
      EXP_M(2,3)=q(2)*q(3)-q(1)*q0
      EXP_M(3,1)=q(3)*q(1)-q(2)*q0
      EXP_M(3,2)=q(3)*q(2)+q(1)*q0
 
      EXP_M=2.0D0*EXP_M
 
      RLOG=EXP_M
 
      return
      end
 
 
 !**************************************************************************
 !  this subroutine is to fin the N_LOG used to calculate SPIN_LOG
      SUBROUTINE FIND_N_LOG(B_P,B,D,N_LOG)
 
      IMPLICIT REAL*8 (A-H,M-N,O-Z)
      real*8 one,two,three, tempv, delta
      PARAMETER(one=1.d0,two=2.d0,three=3.d0)
      DIMENSION B_P(3), B(3,3), D(3,3), N_LOG(3,3)
      DIMENSION temp_v(3), eps(3)
      integer i,j
 
      TOL=0.00000000001D0                                                   ! Tolerance to determine
 
      temp_v= 0.0D0
      N_LOG = 0.0D0
 
 
      IF ((ABS(B_P(1)-B_P(2)).LE.TOL).AND.(ABS(B_P(2)-B_P(3))
     *     .LE.TOL).AND.(ABS(B_P(1)-B_P(3)).LE.TOL)) THEN                   !Case 1: three eigenvalues equals each other
 !             write(*,*)'CASE 1'
           N_LOG = 0.0D0
 
      ELSEIF((ABS(B_P(1)-B_P(2)).GE.TOL).AND.(ABS(B_P(2)-B_P(3))
     *      .GE.TOL).AND.(ABS(B_P(1)-B_P(3)).GE.TOL)) THEN
 !        write(*,*)'CASE 3'                                             !Case 3: None of eigenvalues equals each other
      eps(1) = B_P(2)/B_P(3)
      eps(2) = B_P(3)/B_P(1)
      eps(3) = B_P(1)/B_P(2)
      delta   =  (B_P(1)-B_P(2))*(B_P(2)-B_P(3))*(B_P(3)-B_P(1))
      DO j=1,3
        DO i=1,3
             temp_v(j) = temp_v(j)- one/delta*(-B_P(i))**(3-j)*((1       !!!******注意的地方，符号的正负
     *     +eps(i))/(1-eps(i))+two/log(eps(i)))
        end do
 !         WRITE(*,*)'This is temp_v(j)'
 !         WRITE(*,*),temp_v(j)
      end do
      N_LOG=temp_v(1)*(MATMUL(B,D)-MATMUL(D,B)) +
     *        temp_v(2)*(MATMUL(MATMUL(B,B),D)-MATMUL(D,MATMUL(B,B)))+
     *        temp_v(3)*(MATMUL(MATMUL(MATMUL(B,B),D),B)-
     *        MATMUL(B,MATMUL(D,MATMUL(B,B))))
 !          WRITE(*,*)'This is N_LOG'
 !          WRITE(*,*),N_LOG
 
      ELSEIF((ABS(B_P(1)-B_P(2)).GT.TOL).AND.(ABS(B_P(2)-B_P(3))        ! Two eigenvalues equals each other
     * .LE.TOL)) THEN                                                   ! Case 2-1: b1 != b2; b2=b3
         tempv =one/(B_P(1)-B_P(2))*((one+B_P(1)/B_P(2))/
     *   (one-B_P(1)/B_P(2))+two/(log(B_P(1)/B_P(2))))
 !            write(*,*)'CASE 2-1'
          N_LOG=tempv*(MATMUL(B,D)-MATMUL(D,B))
 
      ELSEIF((ABS(B_P(1)-B_P(2)).LE.TOL).AND.(ABS(B_P(2)-B_P(3))
     * .GT.TOL)) THEN
 !          write(*,*)'CASE 2-2'                                         ! Case 2-2: !Case 2-1: b1= b2; b2!=b3
        tempv =one/(B_P(3)-B_P(2))*((one+B_P(3)/B_P(2))/
     *  (one-B_P(3)/B_P(2))+two/(log(B_P(3)/B_P(2))))
 
          N_LOG=tempv*(MATMUL(B,D)-MATMUL(D,B))
 
      ELSEIF((ABS(B_P(1)-B_P(3)).LE.TOL).AND.(ABS(B_P(1)-B_P(2))
     *      .GE.TOL)) THEN
 !          write(*,*)'CASE 2-3'                                         ! Case 2-3: !Case 2-1: b1= b3; b1!=b2
       tempv =one/(B_P(2)-B_P(3))*((one+B_P(2)/B_P(3))/
     *       (one-B_P(2)/B_P(3))+two/(log(B_P(2)/B_P(3))))
 
          N_LOG=tempv*(MATMUL(B,D)-MATMUL(D,B))
 
 
      ELSE
        WRITE(*,*)'ERROR in SUNROUTINE FIND_N_LOG'
        WRITE(*,*)'This is B_P'
        WRITE(*,*),B_P
 
      ENDIF
 
      RETURN
      END
 
 
 !**************************************************************************
 !  This subroutine is to transfer Cauchy stress into Kirchhoff stress
       SUBROUTINE STRESS_K_C(ST,DFGR)
      !
 	  ! ST:   		passed in as cauchy stress, output as Kirchhoff stress
 	  ! DFG  : 		the current deformation gradient which in UMAT is DFGRD1
 	  ! detF :       determinant of F
 
      !
      IMPLICIT REAL*8 (A-H,O-Z)
      !
      DIMENSION ST(6) ,DFGR(3,3)
 !           detF=1d0
         TOL=0.00000001D0
      CALL mt_det(DFGR,detF)     ! determinant of DFGR
 
      IF (detF.LE.TOL) THEN
      WRITE(*,*)'ERR in Sub_STRESS_K_C, detF of DFGRD1,',detF
      ENDIF
 
      ST=ST/detF
 !        WRITE(*,*),'After K_C stress component',ST
 !
 
       RETURN
 
       END
 
 !**************************************************************************
 !  This subroutine is to transfer Cauchy stress into Kirchhoff stress
      SUBROUTINE STRESS_C_K(ST,DFGR)
      !
 	  ! ST:   		passed in as cauchy stress, output as Kirchhoff stress
 	  ! DFGR  : 		the current deformation gradient which in UMAT is DFGRD1
 	  ! detF :       determinant of F
 
      !
        IMPLICIT REAL*8 (A-H,O-Z)
 !      !
        DIMENSION ST(6), DFGR(3,3)
 !          detF=1d0
        TOL=0.00000001D0
 
      CALL mt_det(DFGR,detF)
 !       WRITE(*,*)'DetF of DFGRD1',detF
      IF (detF.LE.TOL) THEN
      WRITE(*,*)'ERR in Sub_STRESS_C_K, detF of DFGRD1,',detF
      ENDIF
 
       ST=ST*detF
 
       RETURN
       END
 
 !     ******************************************************************
        subroutine FIND_E(F,e_total)
 
        implicit none
        real*8  F(3,3), e_total(6), HenckyStr(3,3)
 
        CALL skinem(F,HenckyStr)
 
        e_total(1)  = HenckyStr(1,1)
        e_total(2)  = HenckyStr(2,2)
        e_total(3)  = HenckyStr(3,3)
        e_total(4)  = 2.d0*HenckyStr(1,2)
        e_total(5)  = 2.d0*HenckyStr(1,3)
        e_total(6)  = 2.d0*HenckyStr(2,3)
 
        end subroutine
 
 !****************************************************************************
 !     The next subroutine calculates various kinematical quantities
 !      associated with the deformation gradient - used in large
 !      deformation analyses.
 !****************************************************************************
 
      subroutine skinem(F,E)
      !
      ! This subroutine performs the left polar decomposition F = VR of the deformation gradient F into a rotation R
      ! and the left stretch tensor V.  The logarithmic strain E = ln(V) is also returned.
 
      !	F(3,3):       the deformation gradient; input
      !	detF:         the determinant of F; detF > 0
      !	R(3,3):       the rotation matrix; output
      !	V(3,3):       the left stretch tensor; output
      !	Vinv(3,3):    the inverse of V
      !	b(3,3):       the right Cauchy-Green tensor
      !	omega(3):     the squares of the principal stretches
      ! Veigval(3):   the principal stretches
      !	eigvec(3,3):  matrix of eigenvectors of V
      !	E(3,3):       the logarithmic strain tensor; output
      ! istat:        success flag, istat=0 for a failed attempt; output
      !
      implicit none
      !
      integer istat
      !
      real*8 F(3,3),b(3,3),omega(3),Veigval(3),eigvec(3,3),
     +  V(3,3),E(3,3),Vinv(3,3),R(3,3),detF
 
 
      !	Store the identity matrix in R, V, and Vinv
      !
      call onem(R)
      call onem(V)
      call onem(Vinv)
 
      ! Store the zero matrix in E
      !
      E = 0.d0
 
      ! Check if the determinant of F is greater than zero.
      !  If not, then print a diagnostic and cut back the
      !  time increment.
      !
      call mt_det(F,detF)
      if (detF.le.0.d0) then
        write(*,'(/5X,A/)') '--problem in kinematics-- the',
     +       ' determinant of F is not greater than 0'
        istat = 0
        return
      end if
 
      ! Calculate the left Cauchy-Green tensor b
      b = matmul(F,transpose(F))
      !b = matmul(transpose(F),F)
 
      ! Calculate the eigenvalues and eigenvectors of b
      call spectral(b,omega,eigvec,istat)
 
 
      ! Calculate the principal values of V and E
      !
      Veigval(1) = dsqrt(omega(1))
      Veigval(2) = dsqrt(omega(2))
      Veigval(3) = dsqrt(omega(3))
      !
      V(1,1) = Veigval(1)
      V(2,2) = Veigval(2)
      V(3,3) = Veigval(3)
      !
      E(1,1) = dlog(Veigval(1))
      E(2,2) = dlog(Veigval(2))
      E(3,3) = dlog(Veigval(3))
 
      V = matmul(matmul(eigvec,V),transpose(eigvec))
      E = matmul(matmul(eigvec,E),transpose(eigvec))
 
 
      return
      end subroutine skinem
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !****************************************************************************%%%%%%%%%%%
 ! %%%%%%    									  Utility subroutines                                                   %%%%%%%%%%%
 !****************************************************************************%%%%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !****************************************************************************%%%%%%%%%%%
 
      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)
 
 
      do i=1,3
         do J=1,3
 	    if (i .eq. j) then
              A(i,j) = 1.d0
            else
              A(i,j) = 0.d0
            end if
         end do
      end do
 
 
      return
      end subroutine onem
 
 !****************************************************************************
 
 
 !****************************************************************************
 
      subroutine jacobi(A,n,np,D,V,nrot,istat)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     +  sm,tresh,G,T,H,theta,S,C,tau
 
 
      ! Initialize V to the identity matrix
      !
      call onem(V)
 
 
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
    	B(ip) = A(ip,ip)
    	D(ip) = B(ip)
    	Z(ip) = 0.d0
      end do
 
 
      ! Begin iteration
      !
      nrot = 0
      do i=1,50
          !
          ! Sum off-diagonal elements
          !
          sm = 0.d0
          do ip=1,n-1
            do iq=ip+1,n
 	      sm = sm + dabs(A(ip,iq))
            end do
          end do
          !
          ! If sm = 0., then return.  This is the normal return,
          !  which relies on quadratic convergence to machine
          !  underflow.
          !
          if (sm.eq.0.d0) return
          !
          ! In the first three sweeps carry out the PQ rotation only if
          !  |A_PQ| > tresh, where tresh is some threshold value,
          !  see equation (11.1.25).  Thereafter tresh = 0.
          !
          if (i.lt.4) then
            tresh = 0.2d0*sm/n**2
          else
            tresh = 0.d0
          end if
          !
          do ip=1,n-1
            do iq=ip+1,n
              G = 100.d0*dabs(A(ip,iq))
              !
              ! After four sweeps, skip the rotation if the
              !  off-diagonal element is small.
              !
 	      if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     +            .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
                A(ip,iq) = 0.d0
              else if (dabs(A(ip,iq)).gt.tresh) then
                H = D(iq) - D(ip)
                if (dabs(H)+G.eq.dabs(H)) then
                  !
                  ! T = 1./(2.*theta), equation (11.1.10)
                  !
 	          T =A(ip,iq)/H
 	        else
 	          theta = 0.5d0*H/A(ip,iq)
 	          T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
 	          if (theta.lt.0.d0) T = -T
 	        end if
 	        C = 1.d0/dsqrt(1.d0 + T**2.d0)
 	        S = T*C
 	        tau = S/(1.d0 + C)
 	        H = T*A(ip,iq)
 	        Z(ip) = Z(ip) - H
 	        Z(iq) = Z(iq) + H
 	        D(ip) = D(ip) - H
 	        D(iq) = D(iq) + H
 	        A(ip,iq) = 0.d0
                !
                ! Case of rotations 1 <= J < P
 		!
 	        do j=1,ip-1
 	          G = A(j,ip)
 	          H = A(j,iq)
 	          A(j,ip) = G - S*(H + G*tau)
 	          A(j,iq) = H + S*(G - H*tau)
 	        end do
                !
                ! Case of rotations P < J < Q
                !
 	        do j=ip+1,iq-1
 	          G = A(ip,j)
 	          H = A(j,iq)
 	          A(ip,j) = G - S*(H + G*tau)
 	          A(j,iq) = H + S*(G - H*tau)
 	        end do
                !
                ! Case of rotations Q < J <= N
                !
 	        do j=iq+1,n
                  G = A(ip,j)
 	          H = A(iq,j)
 	          A(ip,j) = G - S*(H + G*tau)
 	          A(iq,j) = H + S*(G - H*tau)
 	        end do
 	        do j = 1,n
 	          G = V(j,ip)
 	          H = V(j,iq)
 	          V(j,ip) = G - S*(H + G*tau)
 	          V(j,iq) = H + S*(G - H*tau)
 	        end do
 	        nrot = nrot + 1
              end if
 	    end do
 	  end do
          !
          ! Update D with the sum of T*A_PQ, and reinitialize Z
          !
 	  do ip=1,n
 	    B(ip) = B(ip) + Z(ip)
 	    D(ip) = B(ip)
 	    Z(ip) = 0.d0
 	  end do
 	end do
 
 
      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in jacobi should never happen'
      istat = 0
 
 
      return
      end subroutine jacobi
 
 !****************************************************************************
 !     The following subroutines calculate the spectral
 !      decomposition of a symmetric 3 by 3 matrix
 !****************************************************************************
 
      subroutine spectral(A,D,V,istat)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      implicit none
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*8 D(3),V(3,3),A(3,3),E(3,3)
 
 
      E = A
      !
      call jacobi(E,3,np,D,V,nrot,istat)
      !call eigsrt(D,V,3,np)
 
 
      return
      end subroutine spectral
 
 
 !     ******************************************************************
 !     This subroutine calculates the determinant of a 3 by 3 matrix [A]
      SUBROUTINE mt_det(A,detF)
      !
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 detF
      DIMENSION  A(3,3)
 
 
      detF = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) +
     *            A(1,3)*A(2,1)*A(3,2) -  A(3,1)*A(2,2)*A(1,3) -
     *  	         A(3,2)*A(2,3)*A(1,1) -  A(3,3)*A(2,1)*A(1,2)
 
      RETURN
      END
      !C
 
      ! This subroutine is to create a nxn identity matrix
      subroutine IDEN(A,n)
      !	implicit none
      	integer n, i, j
      	real*8 A(n,n)
          do i=1,n
              do j=1,n
                  A(i,j)=0.d0
              end do
      	end do
      	do i=1,n
      		A(i,i)=1.d0
      	end do
      end subroutine IDEN
 
 
      ! This subroutine is to create a nxn zero matrix
      subroutine ZEROS(A,n)
      !	implicit none
      	integer n, i, j
      	real*8 A(n,n)
          do i=1,n
              do j=1,n
                  A(i,j)=0.d0
              end do
      	end do
      end subroutine
 
 !
 !  **************************************************************************
 
 
        subroutine inverse(a,c,n)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed
      ! during the calculation
      !===========================================================
      implicit none
      integer n
      real*8 a(n,n), c(n,n)
      real*8 L(n,n), U(n,n), b(n), d(n), x(n)
      real*8 coeff
      integer i, j, k
 
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 alows such operations on matrices
      L=0.0
      U=0.0
      b=0.0
 
      ! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do
 
      ! Step 2: prepare L and U matrices
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do
 
      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
      ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
      ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine inverse
 
 
 !
 !  **************************************************************************
 ! 	This is to calculate the principal value and its corresponding vector for a 3X3 Matrix
 !  D returns eigenvalue; V returns eigenvectors
      subroutine MAT_P(A,n,np,D,V)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     *   sm,tresh,G,T,H,theta,S,C,tau
 
 
      ! Initialize V to the identity matrix
      !
      call IDEN(V,3)
 
 
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
 	     B(ip) = A(ip,ip)
 	     D(ip) = B(ip)
 	     Z(ip) = 0.d0
      end do
 
 
      ! Begin iteration
      !
      nrot = 0
      do i=1,50
          !
          ! Sum off-diagonal elements
          !
          sm = 0.d0
          do ip=1,n-1
            do iq=ip+1,n
 	      sm = sm + dabs(A(ip,iq))
            end do
          end do
          !
          ! If sm = 0., then return.  This is the normal return,
          !  which relies on quadratic convergence to machine
          !  underflow.
          !
          if (sm.eq.0.d0) return
          !
          ! In the first three sweeps carry out the PQ rotation only if
          !  |A_PQ| > tresh, where tresh is some threshold value,
          !  see equation (11.1.25).  Thereafter tresh = 0.
          !
          if (i.lt.4) then
            tresh = 0.2d0*sm/n**2
          else
            tresh = 0.d0
          end if
          !
          do ip=1,n-1
            do iq=ip+1,n
              G = 100.d0*dabs(A(ip,iq))
              !
              ! After four sweeps, skip the rotation if the
              !  off-diagonal element is small.
              !
 	      if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     *       .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
                A(ip,iq) = 0.d0
              else if (dabs(A(ip,iq)).gt.tresh) then
                H = D(iq) - D(ip)
                if (dabs(H)+G.eq.dabs(H)) then
                  !
                  ! T = 1./(2.*theta), equation (11.1.10)
                  !
 	          T =A(ip,iq)/H
 	        else
 	          theta = 0.5d0*H/A(ip,iq)
 	          T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
 	          if (theta.lt.0.d0) T = -T
 	        end if
 	        C = 1.d0/dsqrt(1.d0 + T**2.d0)
 	        S = T*C
 	        tau = S/(1.d0 + C)
 	        H = T*A(ip,iq)
 	        Z(ip) = Z(ip) - H
 	        Z(iq) = Z(iq) + H
 	        D(ip) = D(ip) - H
 	        D(iq) = D(iq) + H
 	        A(ip,iq) = 0.d0
                !
                ! Case of rotations 1 <= J < P
 		!
 	        do j=1,ip-1
 	          G = A(j,ip)
 	          H = A(j,iq)
 	          A(j,ip) = G - S*(H + G*tau)
 	          A(j,iq) = H + S*(G - H*tau)
 	        end do
                !
                ! Case of rotations P < J < Q
                !
 	        do j=ip+1,iq-1
 	          G = A(ip,j)
 	          H = A(j,iq)
 	          A(ip,j) = G - S*(H + G*tau)
 	          A(j,iq) = H + S*(G - H*tau)
 	        end do
                !
                ! Case of rotations Q < J <= N
                !
 	        do j=iq+1,n
                  G = A(ip,j)
 	          H = A(iq,j)
 	          A(ip,j) = G - S*(H + G*tau)
 	          A(iq,j) = H + S*(G - H*tau)
 	        end do
 	        do j = 1,n
 	          G = V(j,ip)
 	          H = V(j,iq)
 	          V(j,ip) = G - S*(H + G*tau)
 	          V(j,iq) = H + S*(G - H*tau)
 	        end do
 	        nrot = nrot + 1
              end if
 	     end do
 	     end do
          !
          ! Update D with the sum of T*A_PQ, and reinitialize Z
          !
 	     do ip=1,n
 	     B(ip) = B(ip) + Z(ip)
 	     D(ip) = B(ip)
 	     Z(ip) = 0.d0
 	     end do
 	     end do
 
 
 
      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in MAT_P should never happen'
      istat = 0
 
 
      return
      end subroutine MAT_P
 
 
      SUBROUTINE DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
 !
      INCLUDE 'ABA_PARAM.INC'
 !
      DIMENSION U(3),TIME(2),COORDS(3)
 !
 !  THE PRIMARY FUNCTION IS F(T)=1.5 SIN(OMEGA T) + .5 COS(OMEGA T)
 !  WHERE OMEGA IS .1 PI AND T THE CURRENT TIME.
 !  THE APPROPRIATE INTEGRALS AND DERIVATIVES ARE COMPUTED ACCORDINGLY.
 !
      OMEGA=2.d0*3.1415926535897932d0
      AOMGT=OMEGA*TIME(1)
      CS=COS(AOMGT)
      SN=SIN(AOMGT)
      AN=0.1
      BN=0.1
      IF(JDOF.EQ.1) THEN
 !	  x direction
          U(1)=AN*(1.d0-CS)
 !		  write(*,*),U(1)
      ELSE IF(JDOF.EQ.2) THEN
 !	  y direction
         U(1)=AN*SN
      ELSE IF(JDOF.EQ.3) THEN
 !	  z direction
         U(1)=0.d0
      ENDIF
      RETURN
      END
