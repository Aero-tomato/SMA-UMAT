!DEC$ FREEFORM
! Open with Notepad++ for format compatibility
!**********************************************************************************XXX
!                                                                                                                   XXX
!                              3D FINITE STRIN SMA UMAT                                       XXX
!                                                                                                                   XXX
!    			COPYRIGHT   LEI XU1, Alexandros Solomou1, Theo Baxevanis2, Dimitris Lagoudas1 	XXX
!           1 Department of Aerospace Engineering, Texas A&M University    XXX
!           2 Department of Mechanical Engineering, University of Houston  XXX
!               Version  (3-D)  20/JAN./2018                                                     XXX
!                                                                                                                   XXX                             
! XXX************************************************************************************************ XXX  
! XXX  *MATERIAL, NAME=SMA                                                                                              XXX                                                                   
! XXX  *DEPVAR                                                                                                                       XXX   
! XXX    26                                                                                                                               XXX   
! XXX     statev(1)           =      x                                   !Martensite Vol. Frac.                           XXX   
! XXX     statev(2)           =      real(flag_fwd)              !FWD finish at 1 or not 0                     XXX   
! XXX     statev(3)           =      real(flag_rev)               !REV  finish at 1 or not 0                      XXX   
! XXX     statev(4)           =      real(transformation)     !Trans. DIR  0-elastic; -1 REV; 1 FWD  XXX   
! XXX     statev(5:10)      =      et(1:6)                          !Trans. Strain                                       XXX   
! XXX     statev(11:16)    =      et_tr(1:6)                     !Trans. Strain at reverse                       XXX   
! XXX     statev(17:22)    =      Lamdat_r(1:6)             !Reverse DIR tensor                              XXX   
! XXX     statev(23)         =      RPLC                                                                                        XXX   
! XXX     statev(24)         =      BOUND_REACHED                                                                  XXX   
! XXX     statev(25)         =      NR_CONVERGENCE                                                                 XXX   
! XXX     statev(26)         =      TOO                                                                                         XXX   
! XXX*******************************************************************************************************************************XXX
! XXX  *User Material, constants=25                                                                                                                                  XXX
! XXX  <Ea>,              <Em>,          <v>,           <alpha>,        <CA>,            <CM>,                <Ms>,         <Mf>                    XXX 
! XXX  <As>,              <Af>,          <Hmax>,         <kt>,          <SigCal>,       <n1>,                 <n2>,         <n3>                    XXX                                             
! XXX  <n4>,             <X_initial>,      <Etd>,  <Tube_Flag>,      <Tube_Axis>,  <Tube_Radius>, <Elastic>,    <Coupling>  XXX                                                         
! XXX   <Model>                                                                                                                                                                  XXX
! XXX********************************************************************************************************************************XXX


SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
                            & NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
INCLUDE 'ABA_PARAM.INC'
CHARACTER*80 CMNAME
DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
                 & PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)				 
REAL*8,ALLOCATABLE ::  PARAM(:),VAR(:)
INTEGER NPARAM, NVAR
REAL*8 variable
NPARAM=31 							!  NUMBER OF  MODEL PARAMETERS, Check subroutine INITIALIZATION
NVAR=45  								!  NUMBER OF  DEFINED DEPENDENT VARIABLES, Check subroutine INITIALIZATION 
ALLOCATE(PARAM(NPARAM))
ALLOCATE(VAR(NVAR))	 
CALL INITIALIZATION(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,DFGRD0,DFGRD1,PARAM,NPARAM,VAR,NVAR)
CALL ELASTIC_PREDICTOR(PARAM,NPARAM,VAR,NVAR,TEMP,DTEMP,STATEV,NSTATV,time)
IF (nint(VAR(5))/=0) THEN
CALL TRANS_CORRECTOR(PARAM,NPARAM,VAR,NVAR,TEMP,DTEMP)
ENDIF
CALL STRESS_UPDATE(PARAM,NPARAM,VAR,NVAR,temp+dtemp,stress)
CALL JACOBIAN_MATRIX(PARAM,NPARAM,VAR,NVAR,stress,temp+dtemp,ddsdde)
CALL SDV_UPDATE(PARAM,NPARAM,VAR,NVAR,STATEV,NSTATV)
RETURN
END SUBROUTINE

!**********************************************************************
!***********************END of MAIN UMAT*************************
!**********************************************************************

SUBROUTINE INITIALIZATION(STRAN,DSTRAN,TEMP,DTEMP,COORDS,DROT,NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,TIME,DTIME,JSTEP,PROPS,NPROPS,STATEV,NSTATV,DFGRD0,DFGRD1,PARAM,NPARAM,VAR,NVAR)
IMPLICIT NONE
!DEFINITION OF ABAQUS PROVIDED VARIABLES
REAL*8 STRAN(6),DSTRAN(6),COORDS(3),DROT(3,3),TIME(2),DFGRD0(3,3),DFGRD1(3,3),STATEV(NSTATV),PROPS(NPROPS)
REAL*8 TEMP,DTEMP,DTIME
INTEGER JSTEP(4)
INTEGER NTENS,NDI,NSHR,NOEL,NPT,LAYER,KSPT,KINC,NSTATV,NPROPS


!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8   Sa,Sm,v,alpha,CA,CM,As,Af,Ms,Mf,Hmax,kt,SigCal,n1,n2,n3,n4,l1,X_Initial,Tube_Radius,kp
REAL*8   rdso, rduo, Dc, Yo, a1, a2, a3
INTEGER EtD,Tube_Flag,Tube_Axis,ELASTIC,COUPLING,MODEL,NLGEOM
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 TOO,RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER transformation,flag_fwd,flag_rev,NR_CONVERGENCE,BOUND_REACHED
INTEGER NVAR
REAL*8 VAR(NVAR)



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& MATLAB GENERATED VARIABLES SPECIFIC FOR THIS SUBROUTINE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!real*8 bs1,bs2,bs3,bs4,bs5,bs6
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
x_initial          =props(18)

Etd                 =nint(props(19))
Tube_Flag       =nint(props(20))
Tube_Axis       =nint(props(21))
Tube_Radius   =props(22)
ELASTIC         =nint(props(23))
Coupling         =nint(props(24))
Model             =nint(props(25))


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
	x =x_initial                       ! Check SMA phase
    xo=x_initial
	too=TEMP
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
CALL SDV_ASSIGN(x,xo,flag_fwd,flag_rev,transformation,et,eto,et_tr,Lamdat_r,rplc,BOUND_REACHED,NR_CONVERGENCE,too,STATEV,NSTATV)
ENDIF

! CALCULATE COMMON MATERIAL PROPERTIES such as: a1 a2 a3, rdso, rduo, Dc, Y
!DEC$ NOFREEFORM  
      t2 = SigCal**2
      t3 = sqrt(t2)
      t7 = kt*t3
      t4 = exp(-t7)
      t5 = ca+cm
      t6 = 1.0D0/t5
      t8 = t4-1.0D0
      t9 = Hmax*t8
      t10 = sa-sm
      t11 = SigCal*t10
      t14 = Hmax*kt*t3*t4
      t12 = t9+t11-t14
      t13 = af-as
      t15 = mf-ms
      t16 = n3+1.0D0
      t17 = 1.0D0/t16
      t18 = t17+1.0D0
      t19 = n1+1.0D0
      t20 = 1.0D0/t19
      t21 = t20+1.0D0
      A0(1,1) = ca*cm*t6*t12*t15*2.0D0
      A0(1,2) = ca*cm*t6*t12*t13*(-2.0D0)
      A0(1,3) = ca*cm*t6*t12*t13*t18*(-1.0D0/2.0D0)-ca*cm*t6*t12*t15*t21
     &*(1.0D0/2.0D0)
      A0(1,4) = ca*cm*t6*t12*2.0D0
      A0(1,5) = ca*cm*t6*t12*(af+ms)
      A0(1,6) = -(t6*t12*(ca-cm))/(t9-t14)
      A0(1,7) = -ca*cm*t6*t12*(af-ms)+ca*cm*t6*t12*t13*t18*(1.0D0/2.0D0)
     &+ca*cm*t6*t12*t15*t21*(1.0D0/2.0D0)
	 
	 
       a1     =A0(1,1)
       a2     =A0(1,2)
       a3     =A0(1,3)
       rdso  =A0(1,4)
       rduo  =A0(1,5)
       Dc     =A0(1,6)
       Yo     =A0(1,7) 

!DEC$ FREEFORM

	  
	  

if(NLGEOM==1)THEN ! 1-> ENABLED , 0--> DISABLED
  FN=DFGRD0
  FN1=DFGRD1
  CALL FIND_LOG(RLOG,FN,FN1)   !Find logarithmic incrementation rotation matrix.
  CALL FIND_E(FN1,e)                  ! Find current  total logarithmic strain  e for n+1 step.
  CALL FIND_E(FN,eo)                 ! Find previous total logarithmic strain eo for  n step.
  CALL ROTATION_SDV(eo,et,eto,et_tr,Lamdat_r,RLOG,NDI,NSHR)  ! RLOG
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
PARAM(24)=TOO     
PARAM(25)=X_UP_BOUND        
PARAM(26)=X_LOW_BOUND
PARAM(27)=RATIO
PARAM(28)=NLGEOM
PARAM(29)=ELASTIC
PARAM(30)=COUPLING
PARAM(31)=MODEL     
   

VAR(1)           =      x 
VAR(2)           =      xo
VAR(3)         =      real(flag_fwd)             
VAR(4)         =      real(flag_rev)             
VAR(5)         =      real(transformation)    
VAR(6:11)    =      et(1:6) 
VAR(12:17)    =      eto(1:6)                                            
VAR(18:23)    =      et_tr(1:6)                  
VAR(24:29)    =      Lamdat_r(1:6)            
VAR(30:35)    =      e(1:6)
VAR(36:41)    =      eo(1:6)
VAR(42)         =      RPLC
VAR(43)         =      real(BOUND_REACHED)
VAR(44)         =      real(NR_CONVERGENCE)
VAR(45)         =      TOO

END SUBROUTINE

SUBROUTINE ROTATION_SDV(eo,et,eto,et_tr,Lamdat_r,RLOG,NDI,NSHR)

IMPLICIT NONE
REAL*8 eo(6),et(6),eto(6),et_tr(6),Lamdat_r(6),RLOG(3,3)

integer nshr,ndi

!CALL ROTATION_PROCEDURE(eo,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(et,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(eto,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(et_tr,RLOG,2,NDI,NSHR)
CALL ROTATION_PROCEDURE(Lamdat_r,RLOG,2,NDI,NSHR)


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

SUBROUTINE SDV_ASSIGN(x,xo,flag_fwd,flag_rev,transformation,et,eto,et_tr,Lamdat_r,rplc,BOUND_REACHED,NR_CONVERGENCE,too,STATEV,NSTATV)
implicit none

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 TOO,RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER transformation,flag_fwd,flag_rev,NR_CONVERGENCE,BOUND_REACHED

integer NSTATV
real*8 STATEV(NSTATV)

     x                             =statev(1)                             
     xo                           =statev(1)                             
     flag_fwd                  =nint(statev(2))                     
     flag_rev                   =nint(statev(3))                    
     transformation         =nint(statev(4))                     
     et(1:6)                    =statev(5:10)                        
     eto(1:6)                  =statev(5:10)                                
	 et_tr(1:6)                =statev(11:16)                      
     Lamdat_r(1:6)         =statev(17:22)                                                  
     RPLC                        =        statev(23)                                                                                      
     BOUND_REACHED     = nint(statev(24))   	                                                                              
     NR_CONVERGENCE   = nint(statev(25))  	                                                                            
     TOO                          =        statev(26)																				
                                         
end
                    
SUBROUTINE PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL) 
IMPLICIT NONE
!DEFINITION OF USER DEFINED MODEL PARAMETERS
REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,kp
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
TOO                  =PARAM(24)
X_UP_BOUND    =PARAM(25)  
X_LOW_BOUND  =PARAM(26)
RATIO               =PARAM(27)
NLGEOM            =nint(PARAM(28))
ELASTIC            =nint(PARAM(29))
COUPLING         =nint(PARAM(30))
MODEL              =nint(PARAM(31)) 

END SUBROUTINE

SUBROUTINE VAR_ASSIGNMENTS(INDEX,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too) 
IMPLICIT NONE

!DEFINITIONS OF ARGUMENTS SPECIFIC FOR THIS SUBROUTINE
INTEGER INDEX


!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 too,RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CODE START @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

IF (INDEX==1)THEN

x                                      =           VAR(1)         
xo                                    =           VAR(2)                 
flag_fwd                           =           nint(VAR(3))       
flag_rev                            =           nint(VAR(4))       
transformation                  =           nint(VAR(5))       
et(1:6)                             =           VAR(6:11)   
eto(1:6)                           =           VAR(12:17)  
et_tr(1:6)                        =           VAR(18:23)  
Lamdat_r(1:6)                 =           VAR(24:29)   
e(1:6)                              =           VAR(30:35)   
eo(1:6)                            =           VAR(36:41)  
rplc                                  =           VAR(42)       
BOUND_REACHED            =           nint(VAR(43))       
NR_CONVERGENCE          =           nint(VAR(44))       
TOO                                =           VAR(45)       

ELSEIF (INDEX==2)THEN


VAR(1)           =      x 
VAR(2)           =      xo
VAR(3)           =      real(flag_fwd)             
VAR(4)           =      real(flag_rev)             
VAR(5)           =      real(transformation)    
VAR(6:11)      =      et(1:6) 
VAR(12:17)    =      eto(1:6)                                            
VAR(18:23)    =      et_tr(1:6)                  
VAR(24:29)    =      Lamdat_r(1:6)            
VAR(30:35)    =      e(1:6)
VAR(36:41)    =      eo(1:6)
VAR(42)         =      RPLC
VAR(43)         =      real(BOUND_REACHED)
VAR(44)         =      real(NR_CONVERGENCE)
VAR(45)         =      TOO

ENDIF

END SUBROUTINE   

subroutine ELASTIC_PREDICTOR(PARAM,NPARAM,VAR,NVAR,TEMP,DTEMP,STATEV,NSTATV,time)
IMPLICIT NONE
!DEFINITION OF ABAQUS PROVIDED VARIABLES
REAL*8 TEMP,DTEMP,STATEV(NSTATV)
INTEGER NSTATV

!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio
integer    NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


REAL*8 Tn, Tn1, time(2)

! variables using in ELASTIC_PREDICTOR_subroutine
REAL*8 PHI_Fwd_Cur,PHI_Fwd_Pre,PHI_Rev_Cur,PHI_Rev_Pre,DelPHI_Fwd,DelPHI_Rev,stress(6)
INTEGER load_DIR, checknan(10), checkinf(10), index1, index2

checknan=0
checkinf=0

CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)
!CALL STRESS_UPDATE(PARAM,NPARAM,VAR,NVAR,temp+dtemp,stress)

if (temp<=0) then
    temp=1.0e-8
endif
Tn=temp                                !previous temp

Tn1=temp+dtemp                  !current temp

CALL PHI_Cal(PARAM,NPARAM,x,e,et,et_tr,Lamdat_r,Tn1,too,PHI_Fwd_Cur,1)
CALL check_nan1(PHI_Fwd_Cur,1,checknan(1))
CALL check_inf1(PHI_Fwd_Cur,1,checkinf(1))
      if(checknan(1)==1)then
      write(*,*)'NAN  in PHI_Fwd_Cur'
      endif
      if(checkinf(1)==1)then
      write(*,*)'INF  in PHI_Fwd_Cur'
      endif

	

CALL PHI_Cal(PARAM,NPARAM,xo,eo,eto,et_tr,Lamdat_r,Tn,too,PHI_Fwd_Pre,1)
CALL check_nan1(PHI_Fwd_Pre,1,checknan(2))
CALL check_inf1(PHI_Fwd_Pre,1,checkinf(2))
      if(checknan(2)==1)then
      write(*,*)'NAN  in PHI_Fwd_Pre'
      endif
      if(checkinf(2)==1)then
      write(*,*)'INF  in PHI_Fwd_Pre'
      endif


CALL PHI_Cal(PARAM,NPARAM,x,e,et,et_tr,Lamdat_r,Tn1,too,PHI_Rev_Cur,-1)
CALL check_nan1(PHI_Rev_Cur,1,checknan(3))
CALL check_inf1(PHI_Rev_Cur,1,checkinf(3))      
      if(checknan(3)==1)then
      write(*,*)'NAN  in PHI_Rev_Cur'
      endif
      if(checkinf(3)==1)then
      write(*,*)'INF  in PHI_Rev_Cur'
      endif      


CALL PHI_Cal(PARAM,NPARAM,xo,eo,eto,et_tr,Lamdat_r,Tn,too,PHI_Rev_Pre,-1)
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
    ! write(*,*)"ep",ep
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
CALL VAR_ASSIGNMENTS(2,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

end subroutine

Subroutine TRANS_CORRECTOR(PARAM,NPARAM,VAR,NVAR,TEMP,DTEMP)
IMPLICIT NONE
!DEFINITION OF ABAQUS PROVIDED VARIABLES
REAL*8 TEMP,DTEMP

!DEFINITION OF  INPUT MODEL PARAMETERS
REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio
integer    NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR



REAL*8 Tn, Tn1

real*8 NR_RE(7), NR_JAC(7,7), NR_JAC_INV(7,7),du(7),NR_JAC1(7,7),Residual(2)
integer Iter, ii, jj,output(4)


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

Tn=temp
Tn1=temp+dtemp

NR_CONVERGENCE=0
Iter=0
NR_RE=0
NR_JAC=0
NR_JAC_INV=0
du=0
residual(1)=10
residual(2)=10


DO WHILE (NR_CONVERGENCE==0 .and. Iter<=100)


 CALL N_R_Residual(PARAM,NPARAM,VAR,NVAR,Tn1,transformation,NR_RE)
 residual(1)=maxval(abs(NR_RE(1:6)))
 residual(2)=abs(NR_RE(7))


 if(residual(1)<1E-8_8 .and. residual(2)<1.0e-6_8) then
    !write(*,*)"iteration",Iter
 	NR_CONVERGENCE=1
	
 else

	CALL N_R_JAC(PARAM,NPARAM,VAR,NVAR,Tn1,transformation,NR_JAC)
    NR_JAC1=NR_JAC
	CALL inverse(NR_JAC1,NR_JAC_INV,7)

	du=-matmul(NR_JAC_INV,NR_RE)
	et=et+du(1:6:1)
	x=x+du(7)


!***************Debug Part***************
output=0

call check_inf1(NR_RE,7,output(1))
call check_inf2(NR_JAC,7,output(2))
call check_nan1(NR_RE,7,output(3))
call check_nan2(NR_JAC,7,output(4))
if(output(1)==1.or.output(2)==1.or.output(3)==1.or.output(4)==1)then
                   !write(*,*)' Total time at beginning of current increment',propsUR(29)  
                   write(*,*)"output",output 
                   write(*,*)"transformation",transformation				   
                   !write(*,*)"kinc",kinc1
                   write(*,*)"iteration",Iter
                   !write(*,*)"Element number",noel1
                   !write(*,*)"Integration number",npt1
                   write(*,*)"Volume Frac.",x
                   write(*,*)"Old Vol. Frac.",xo
                   write(*,*)"temperature",Tn1			   
                   write(*,*)"E********************************************"
				   write(*,*),e
                   write(*,*)"Eo********************************************"
				   write(*,*),eo
                   write(*,*)"ET********************************************"
				   write(*,*),ET
                   write(*,*)"ETo********************************************"
				   write(*,*),ETo                   
                   write(*,*)"du********************************************"
				   write(*,*),du
                   write(*,*)"ET_tr********************************************"
				   write(*,*),ET_tr
                   write(*,*)"lamdat_r********************************************"
				   write(*,*),lamdat_r				   
                       
                   do ii=1,7
                       do jj=1,7
                            write(*,*)"NR_JAC",NR_JAC(ii,jj),ii,jj
                       enddo
                   enddo
                        
                   
                        stop
                   endif
				   
!***************Debug Part***************
if(Iter>=20 .and. residual(1)<1E-4_8 .and. residual(2)<1.0e-2_8) then
 ! if(Iter>=50) then
 	NR_CONVERGENCE=1
	!write(*,*)' Total time at beginning of current increment',propsUR(29)  
	WRITE(*,*)"************delta x or du(7)=",du(7)
	WRITE(*,*)"Hard to Converge in Newton_Raphson,  Iter=",Iter
endif
	iter=iter+1  !
	
!******TRANSFORMATION_ADJUSTMENTS*********
if(X_LOW_BOUND<=x<=X_UP_BOUND)then     ! mixed phase, both fwd and rev can happens
	flag_fwd=0
	flag_rev=0
endif

if(transformation==1 .and. x>=X_UP_BOUND)then
	x=X_UP_BOUND
	NR_CONVERGENCE=1
	flag_fwd=1                                 ! Forward transformation finished
    flag_rev= 0 
endif

if(transformation==1 .and.x<=X_LOW_BOUND)then
	x=X_LOW_BOUND
	flag_fwd=0
	flag_rev=1 
endif

if(transformation==-1.and. x<=X_LOW_BOUND)then
	x=X_LOW_BOUND
	et=0
	NR_CONVERGENCE=1
	flag_fwd=0
	flag_rev=1                                 ! Reverse transformation finished     
  
endif

if(transformation==-1.and. x>=X_UP_BOUND)then
	x=X_UP_BOUND
	flag_fwd=1
	flag_rev=0
endif


endif

if(transformation==1)then
Lamdat_r=et/x ! THOSE VALUES ARE UPDATED ONLY DURING FORWARD TRANSFORMATION
et_tr=et
endif

CALL VAR_ASSIGNMENTS(2,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)


end do
end

Subroutine PHI_Cal(PARAM,NPARAM,x,e,et,et_tr,lamdat_r,Tn,too,phi,index1)

implicit real*8 (t)

!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio
integer    NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

real*8   Tn

! Changing state variables
real*8 e(6),et(6),ep(6),et_tr(6),Lamdat_r(6),backstress(6),stress(6)
real*8 x,xo,xd,xdo,zt,zto,ztd,ztdo

real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6
real*8 lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 phi

CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)

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


   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)
   
 
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
      t12 = Sa-Sm
      t13 = -e1+et1+t10
      t14 = t4*t8*t13*v
      t15 = -e3+et3+t10
      t16 = v-1.0D0
      t17 = -e2+et2+t10
      t18 = t4*t8*t17*v
      t19 = t4*t8*t15*v
      t28 = t4*t8*t16*t17
      t20 = t14+t19-t28
      t26 = t4*t8*t15*t16
      t21 = t14+t18-t26
      t25 = t4*t8*t13*t16
      t22 = t18+t19-t25
      t23 = t12*t22*v
      t24 = t4*t8*t13*v*(1.0D0/2.0D0)
      t27 = t12*t21*v
      t29 = t12*t20*v
      t30 = t4*t8*t17*v*(1.0D0/2.0D0)
      t31 = t4*t8*t15*v*(1.0D0/2.0D0)
      t32 = e4-et4
      t33 = v*2.0D0
      t34 = t33+2.0D0
      t35 = e5-et5
      t36 = v+1.0D0
      t37 = 1.0D0/t36**2
      t38 = 1.0D0/t3**2
      t39 = e6-et6
      t40 = abs(t36)
      t41 = e1*et2*4.0D0
      t42 = e2*et1*4.0D0
      t43 = e1*et3*4.0D0
      t44 = e3*et1*4.0D0
      t45 = e2*et3*4.0D0
      t46 = e3*et2*4.0D0
      t47 = e1**2
      t48 = t47*4.0D0
      t49 = e2**2
      t50 = t49*4.0D0
      t51 = e3**2
      t52 = t51*4.0D0
      t53 = e4**2
      t54 = t53*3.0D0
      t55 = e5**2
      t56 = t55*3.0D0
      t57 = e6**2
      t58 = t57*3.0D0
      t59 = et1**2
      t60 = t59*4.0D0
      t61 = et2**2
      t62 = t61*4.0D0
      t63 = et3**2
      t64 = t63*4.0D0
      t65 = et4**2
      t66 = t65*3.0D0
      t67 = et5**2
      t68 = t67*3.0D0
      t69 = et6**2
      t70 = t69*3.0D0
      t75 = e1*e2*4.0D0
      t76 = e1*e3*4.0D0
      t77 = e2*e3*4.0D0
      t78 = e1*et1*8.0D0
      t79 = e2*et2*8.0D0
      t80 = e3*et3*8.0D0
      t81 = e4*et4*6.0D0
      t82 = e5*et5*6.0D0
      t83 = e6*et6*6.0D0
      t84 = et1*et2*4.0D0
      t85 = et1*et3*4.0D0
      t86 = et2*et3*4.0D0
      t71 = t41+t42+t43+t44+t45+t46+t48+t50+t52+t54+t56+t58+t60+t62+t64+
     &t66+t68+t70-t75-t76-t77-t78-t79-t80-t81-t82-t83-t84-t85-t86
      t72 = abs(t71)
      t73 = abs(t3)
      t74 = t32**2
      t87 = 1.0D0/sqrt(t72)
      t88 = Hmax*(3.0D0/2.0D0)
      t89 = 1.0D0/t40
      t90 = sqrt(t72)
      t91 = 1.0D0/t73
      t95 = kt*t89*t90*t91*(1.0D0/2.0D0)
      t92 = exp(-t95)
      t96 = Hmax*t92*(3.0D0/2.0D0)
      t93 = t88-t96
      t94 = t35**2
      t97 = t39**2
      t98 = 1.0D0/t36
      t99 = e1*2.0D0
      t100 = et1*2.0D0
      t101 = et2*2.0D0
      t102 = e1-e2*2.0D0+e3-et1-et3+t101
      t103 = alpha*too
      t104 = e2*v
      t105 = e3*v
      t106 = et1*v
      t107 = alpha*too*v
      t108 = et3*2.0D0
      t109 = e1+e2-e3*2.0D0-et1-et2+t108
      t110 = e1*v
      t111 = et2*v
      t0 = -Yo-a3-rduo+Tn*rdso+(t24+t30-t4*t8*t15*t16*(1.0D0/2.0D0))*(t2
     &3+t29-t12*t21)+(t24+t31-t4*t8*t16*t17*(1.0D0/2.0D0))*(t23+t27-t12*
     &t20)+(t30+t31-t4*t8*t13*t16*(1.0D0/2.0D0))*(t27+t29-t12*t22)-a1*(-
     &(-x+1.0D0)**n2+x**n1+1.0D0)*(1.0D0/2.0D0)-t12*t34*t37*t38*t74*(1.0
     &D0/8.0D0)-t12*t34*t37*t38*t94*(1.0D0/8.0D0)-t12*t34*t37*t38*t97*(1
     &.0D0/8.0D0)+t37*t38*t40*t73*t74*t87*t93+t37*t38*t40*t73*t87*t93*t9
     &4+t37*t38*t40*t73*t87*t93*t97-Dc*t37*t38*t40*t73*t74*t87*t93-Dc*t3
     &7*t38*t40*t73*t87*t93*t94-Dc*t37*t38*t40*t73*t87*t93*t97+Dc*t4*t20
     &*t40*t73*t87*t93*t98*t102*(2.0D0/3.0D0)+Dc*t4*t21*t40*t73*t87*t93*
     &t98*t109*(2.0D0/3.0D0)+t8*t38*t40*t73*t87*t93*t98*t102*(e2-et2+t10
     &3-t104+t105-t106+t107+t110+t111-Tn*alpha-et3*v-Tn*alpha*v)*(2.0D0/
     &3.0D0)+t8*t38*t40*t73*t87*t93*t98*(e2+e3-et2-et3-t99+t100)*(e1-et1
     &+t103+t104+t105+t106+t107-Tn*alpha-e1*v-et2*v-et3*v-Tn*alpha*v)*(2
     &.0D0/3.0D0)+t8*t38*t40*t73*t87*t93*t98*t109*(e3-et3+t103+t104-t105
     &-t106+t107+t110-t111-Tn*alpha+et3*v-Tn*alpha*v)*(2.0D0/3.0D0)+Dc*t
     &4*t22*t40*t73*t87*t93*t98*(e2+e3-et2-et3-t99+t100)*(2.0D0/3.0D0)



	 
	    phi=t0
 
       elseif (index1==-1)then ! Reverse
      t2 = Sm*x
      t14 = Sa*x
      t3 = Sa+t2-t14
      t4 = 1.0D0/t3
      t5 = v**2
      t6 = t5*2.0D0
      t7 = t6+v-1.0D0
      t8 = 1.0D0/t7
      t9 = alpha*too
      t10 = e2*v
      t11 = e3*v
      t12 = et1*v
      t13 = alpha*too*v
      t15 = e1*v
      t16 = et2*v
      t17 = v+1.0D0
      t18 = 1.0D0/t17
      t19 = e4-et4
      t20 = v*2.0D0
      t21 = t20+2.0D0
      t22 = Sa-Sm
      t23 = e5-et5
      t24 = 1.0D0/t17**2
      t25 = 1.0D0/t3**2
      t26 = e6-et6
      t27 = et3*v
      t29 = Tn*alpha
      t28 = e3-et3+t9+t10-t11-t12+t13+t15-t16+t27-t29-Tn*alpha*v
      t0 = -Yo-a3+rduo-Tn*rdso+a2*(-(-x+1.0D0)**n4+x**n3+1.0D0)*(1.0D0/2
     &.0D0)+lamdat_r2*t4*t8*(e2-et2+t9-t10+t11-t12+t13+t15+t16-Tn*alpha-
     &et3*v-Tn*alpha*v)+lamdat_r3*t4*t8*t28-lamdat_r4*t4*t18*t19*(1.0D0/
     &2.0D0)-lamdat_r5*t4*t18*t23*(1.0D0/2.0D0)-lamdat_r6*t4*t18*t26*(1.
     &0D0/2.0D0)+lamdat_r1*t4*t8*(e1-et1+t9+t10+t11+t12+t13-Tn*alpha-e1*
     &v-et2*v-et3*v-Tn*alpha*v)-t8*t22*t25*t28*(e3-et3+t9-t29)*(1.0D0/2.
     &0D0)+t19**2*t21*t22*t24*t25*(1.0D0/8.0D0)+t21*t22*t23**2*t24*t25*(
     &1.0D0/8.0D0)+t21*t22*t24*t25*t26**2*(1.0D0/8.0D0)-t8*t22*t25*(e1-e
     &t1+t9-Tn*alpha)*(e1-et1+t9+t10+t11+t12+t13-t15-t16-t27-Tn*alpha-Tn
     &*alpha*v)*(1.0D0/2.0D0)-t8*t22*t25*(e2-et2+t9-Tn*alpha)*(e2-et2+t9
     &-t10+t11-t12+t13+t15+t16-t27-Tn*alpha-Tn*alpha*v)*(1.0D0/2.0D0)


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

REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio
integer    NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR

! FEA integer information
integer index1
real*8   Tn


! Changing state variables
real*8 NR_JAC(7,7),A0(7,7)


real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6
real*8 eto1,eto2,eto3,eto4,eto5,eto6,epo1,epo2,epo3,epo4,epo5,epo6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

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
	


   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)  
   

   

   if(abs(x-xo)<=1e-15)then
       x=xo+1.0e-15_8
   end if

   

!DEC$ NOFREEFORM
 
       if(index1==1)then    !Forward
      t2 = v+1.0D0
      t3 = abs(t2)
      t4 = e1*et2*4.0D0
      t5 = e2*et1*4.0D0
      t6 = e1*et3*4.0D0
      t7 = e3*et1*4.0D0
      t8 = e2*et3*4.0D0
      t9 = e3*et2*4.0D0
      t10 = e1**2
      t11 = t10*4.0D0
      t12 = e2**2
      t13 = t12*4.0D0
      t14 = e3**2
      t15 = t14*4.0D0
      t16 = e4**2
      t17 = t16*3.0D0
      t18 = e5**2
      t19 = t18*3.0D0
      t20 = e6**2
      t21 = t20*3.0D0
      t22 = et1**2
      t23 = t22*4.0D0
      t24 = et2**2
      t25 = t24*4.0D0
      t26 = et3**2
      t27 = t26*4.0D0
      t28 = et4**2
      t29 = t28*3.0D0
      t30 = et5**2
      t31 = t30*3.0D0
      t32 = et6**2
      t33 = t32*3.0D0
      t40 = e1*e2*4.0D0
      t41 = e1*e3*4.0D0
      t42 = e2*e3*4.0D0
      t43 = e1*et1*8.0D0
      t44 = e2*et2*8.0D0
      t45 = e3*et3*8.0D0
      t46 = e4*et4*6.0D0
      t47 = e5*et5*6.0D0
      t48 = e6*et6*6.0D0
      t49 = et1*et2*4.0D0
      t50 = et1*et3*4.0D0
      t51 = et2*et3*4.0D0
      t34 = t4+t5+t6+t7+t8+t9+t11+t13+t15+t17+t19+t21+t23+t25+t27+t29+t3
     &1+t33-t40-t41-t42-t43-t44-t45-t46-t47-t48-t49-t50-t51
      t35 = abs(t34)
      t36 = Sm*x
      t39 = Sa*x
      t37 = Sa+t36-t39
      t38 = abs(t37)
      t52 = 1.0D0/t3
      t53 = sqrt(t35)
      t54 = 1.0D0/t38
      t61 = kt*t52*t53*t54*(1.0D0/2.0D0)
      t55 = exp(-t61)
      t56 = x-xo
      t57 = 1.0D0/t2
      t58 = 1.0D0/t37
      t59 = (t34/abs(t34))
      t60 = Hmax*(3.0D0/2.0D0)
      t74 = Hmax*t55*(3.0D0/2.0D0)
      t62 = t60-t74
      t63 = e1*2.0D0
      t64 = et1*2.0D0
      t65 = e1*8.0D0
      t66 = e2*4.0D0
      t67 = e3*4.0D0
      t68 = et1*8.0D0
      t69 = et2*4.0D0
      t70 = et3*4.0D0
      t71 = t65-t66-t67-t68+t69+t70
      t72 = e2+e3-et2-et3-t63+t64
      t73 = 1.0D0/sqrt(t35)
      t75 = 1.0D0/t35
      t76 = 1.0D0/t35**(3.0D0/2.0D0)
      t77 = e1*4.0D0
      t78 = e2*8.0D0
      t79 = et1*4.0D0
      t80 = et2*8.0D0
      t81 = t67-t70+t77-t78-t79+t80
      t82 = t3*t38*t56*t57*t58*t62*t73*(2.0D0/3.0D0)
      t83 = e3*8.0D0
      t84 = et3*8.0D0
      t85 = t66-t69+t77-t79-t83+t84
      t86 = e4*6.0D0
      t87 = et4*6.0D0
      t88 = t86-t87
      t89 = e5*6.0D0
      t90 = et5*6.0D0
      t91 = t89-t90
      t92 = e6*6.0D0
      t93 = et6*6.0D0
      t94 = t92-t93
      t95 = Sa-Sm
      t96 = (t37/abs(t37))
      t97 = et2*2.0D0
      t99 = e2*2.0D0
      t98 = e1+e3-et1-et3+t97-t99
      t100 = 1.0D0/t37**2
      t101 = et3*2.0D0
      t103 = e3*2.0D0
      t102 = e1+e2-et1-et2+t101-t103
      t104 = e4-et4
      t105 = e5-et5
      t106 = e6-et6
      t107 = v**2
      t108 = t107*2.0D0
      t109 = t108+v-1.0D0
      t110 = 1.0D0/t109
      t111 = Tn-too
      t112 = alpha*t111
      t113 = v-1.0D0
      t114 = -e1+et1+t112
      t115 = t58*t110*t114*v*(1.0D0/2.0D0)
      t116 = -e3+et3+t112
      t117 = -e2+et2+t112
      t118 = t58*t95*t110*v
      t119 = t58*t95*t110*t113*v
      t151 = t58*t95*t107*t110
      t120 = t118+t119-t151
      t121 = t58*t110*t117*v*(1.0D0/2.0D0)
      t122 = t58*t110*t116*v*(1.0D0/2.0D0)
      t123 = t58*t110*t117*v
      t124 = t58*t110*t114*v
      t125 = t58*t110*t116*v
      t131 = t58*t110*t113*t116
      t126 = t123+t124-t131
      t130 = t58*t110*t113*t117
      t127 = t124+t125-t130
      t128 = t95*t127*v
      t133 = t58*t110*t113*t114
      t129 = t123+t125-t133
      t132 = t95*t126*v
      t134 = t95*t129*v
      t135 = alpha*too
      t136 = e2*v
      t137 = e3*v
      t138 = et1*v
      t139 = alpha*too*v
      t140 = e1*v
      t141 = et2*v
      t142 = 1.0D0/t2**2
      t143 = t104**2
      t144 = t105**2
      t145 = t106**2
      t146 = et3*v
      t148 = Tn*alpha
      t149 = Tn*alpha*v
      t147 = e3-et3+t135+t136-t137-t138+t139+t140-t141+t146-t148-t149
      t166 = t58*t110*t113*t116*(1.0D0/2.0D0)
      t150 = t115+t121-t166
      t165 = t58*t110*t113*t114*(1.0D0/2.0D0)
      t152 = t121+t122-t165
      t153 = t58*t95*t110*t113
      t154 = t58*t95*t107*t110*2.0D0
      t155 = t153+t154
      t164 = t58*t110*t113*t117*(1.0D0/2.0D0)
      t156 = t115+t122-t164
      t168 = t95*t127
      t157 = t132+t134-t168
      t167 = t95*t126
      t158 = t128+t134-t167
      t159 = t58*t110*t158*v*(1.0D0/2.0D0)
      t170 = t95*t129
      t160 = t128+t132-t170
      t161 = e1-et1+t135+t136+t137+t138+t139-t140-t141-t146-t148-t149
      t162 = e2-et2+t135-t136+t137-t138+t139+t140+t141-t146-t148-t149
      t163 = Dc*t3*t38*t57*t62*t73*t100*t102*t110*v*(2.0D0/3.0D0)
      t169 = t58*t110*t157*v*(1.0D0/2.0D0)
      t171 = t58*t110*t160*v*(1.0D0/2.0D0)
      t172 = Dc*t3*t38*t57*t62*t73*t98*t100*t110*v*(2.0D0/3.0D0)
      t173 = Dc*t3*t38*t57*t62*t73*t100*t110*v*(e2+e3-et2-et3-t63+t64)*(
     &2.0D0/3.0D0)
      t174 = e4*2.0D0
      t175 = et4*2.0D0
      t176 = t174-t175
      t177 = v*2.0D0
      t178 = t177+2.0D0
      t179 = e5*2.0D0
      t180 = et5*2.0D0
      t181 = t179-t180
      t182 = e6*2.0D0
      t183 = et6*2.0D0
      t184 = t182-t183
      t185 = t95*t100*t110*t117*v
      t186 = t95*t100*t110*t114*v
      t187 = t95*t100*t110*t116*v
      t193 = t95*t100*t110*t113*t117
      t188 = t186+t187-t193
      t192 = t95*t100*t110*t113*t114
      t189 = t185+t187-t192
      t190 = t95*t189*v
      t195 = t95*t100*t110*t113*t116
      t191 = t185+t186-t195
      t194 = t95*t188*v
      t196 = t95*t191*v
      t197 = t95*t100*t110*t114*v*(1.0D0/2.0D0)
      t198 = t95*t100*t110*t117*v*(1.0D0/2.0D0)
      t199 = t95*t100*t110*t116*v*(1.0D0/2.0D0)
      t200 = t95**2
      t201 = 1.0D0/t37**3
      A0(1,1) = t3*t38*t56*t57*t58*t62*t73*(-4.0D0/3.0D0)+Hmax*kt*t55*t5
     &6*t57*t58*t59*t71*t72*t75*(1.0D0/4.0D0)-t3*t38*t56*t57*t58*t59*t62
     &*t71*t72*t76*(1.0D0/3.0D0)-1.0D0
      A0(1,2) = t82+t3*t38*t56*t57*t58*t59*t62*t76*t81*(e2+e3-et2-et3-t6
     &3+t64)*(1.0D0/3.0D0)-Hmax*kt*t55*t56*t57*t58*t59*t72*t75*t81*(1.0D
     &0/4.0D0)
      A0(1,3) = t82-Hmax*kt*t55*t56*t57*t58*t59*t75*t85*(e2+e3-et2-et3-t
     &63+t64)*(1.0D0/4.0D0)+t3*t38*t56*t57*t58*t59*t62*t76*t85*(e2+e3-et
     &2-et3-t63+t64)*(1.0D0/3.0D0)
      A0(1,4) = Hmax*kt*t55*t56*t57*t58*t59*t75*t88*(e2+e3-et2-et3-t63+t
     &64)*(1.0D0/4.0D0)-t3*t38*t56*t57*t58*t59*t62*t72*t76*t88*(1.0D0/3.
     &0D0)
      A0(1,5) = Hmax*kt*t55*t56*t57*t58*t59*t75*t91*(e2+e3-et2-et3-t63+t
     &64)*(1.0D0/4.0D0)-t3*t38*t56*t57*t58*t59*t62*t72*t76*t91*(1.0D0/3.
     &0D0)
      A0(1,6) = Hmax*kt*t55*t56*t57*t58*t59*t75*t94*(e2+e3-et2-et3-t63+t
     &64)*(1.0D0/4.0D0)-t3*t38*t56*t57*t58*t59*t62*t72*t76*t94*(1.0D0/3.
     &0D0)
      A0(1,7) = t3*t38*t57*t58*t62*t72*t73*(-2.0D0/3.0D0)-t3*t38*t56*t57
     &*t62*t72*t73*t95*t100*(2.0D0/3.0D0)+t3*t56*t57*t58*t62*t73*t95*t96
     &*(e2+e3-et2-et3-t63+t64)*(2.0D0/3.0D0)-Hmax*kt*t54*t55*t56*t57*t58
     &*t72*t95*t96*(1.0D0/2.0D0)
      A0(2,1) = t82+Hmax*kt*t55*t56*t57*t58*t59*t71*t75*t98*(1.0D0/4.0D0
     &)-t3*t38*t56*t57*t58*t59*t62*t71*t76*t98*(1.0D0/3.0D0)
      A0(2,2) = t3*t38*t56*t57*t58*t62*t73*(-4.0D0/3.0D0)-Hmax*kt*t55*t5
     &6*t57*t58*t59*t75*t81*t98*(1.0D0/4.0D0)+t3*t38*t56*t57*t58*t59*t62
     &*t76*t81*t98*(1.0D0/3.0D0)-1.0D0
      A0(2,3) = t82-Hmax*kt*t55*t56*t57*t58*t59*t75*t85*t98*(1.0D0/4.0D0
     &)+t3*t38*t56*t57*t58*t59*t62*t76*t85*t98*(1.0D0/3.0D0)
      A0(2,4) = Hmax*kt*t55*t56*t57*t58*t59*t75*t88*t98*(1.0D0/4.0D0)-t3
     &*t38*t56*t57*t58*t59*t62*t76*t88*t98*(1.0D0/3.0D0)
      A0(2,5) = Hmax*kt*t55*t56*t57*t58*t59*t75*t91*t98*(1.0D0/4.0D0)-t3
     &*t38*t56*t57*t58*t59*t62*t76*t91*t98*(1.0D0/3.0D0)
      A0(2,6) = Hmax*kt*t55*t56*t57*t58*t59*t75*t94*t98*(1.0D0/4.0D0)-t3
     &*t38*t56*t57*t58*t59*t62*t76*t94*t98*(1.0D0/3.0D0)
      A0(2,7) = t3*t38*t57*t58*t62*t73*t98*(-2.0D0/3.0D0)-t3*t38*t56*t57
     &*t62*t73*t95*t98*t100*(2.0D0/3.0D0)+t3*t56*t57*t58*t62*t73*t95*t96
     &*t98*(2.0D0/3.0D0)-Hmax*kt*t54*t55*t56*t57*t58*t95*t96*t98*(1.0D0/
     &2.0D0)
      A0(3,1) = t82+Hmax*kt*t55*t56*t57*t58*t59*t71*t75*t102*(1.0D0/4.0D
     &0)-t3*t38*t56*t57*t58*t59*t62*t71*t76*t102*(1.0D0/3.0D0)
      A0(3,2) = t82-Hmax*kt*t55*t56*t57*t58*t59*t75*t81*t102*(1.0D0/4.0D
     &0)+t3*t38*t56*t57*t58*t59*t62*t76*t81*t102*(1.0D0/3.0D0)
      A0(3,3) = t3*t38*t56*t57*t58*t62*t73*(-4.0D0/3.0D0)-Hmax*kt*t55*t5
     &6*t57*t58*t59*t75*t85*t102*(1.0D0/4.0D0)+t3*t38*t56*t57*t58*t59*t6
     &2*t76*t85*t102*(1.0D0/3.0D0)-1.0D0
      A0(3,4) = Hmax*kt*t55*t56*t57*t58*t59*t75*t88*t102*(1.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t88*t102*(1.0D0/3.0D0)
      A0(3,5) = Hmax*kt*t55*t56*t57*t58*t59*t75*t91*t102*(1.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t91*t102*(1.0D0/3.0D0)
      A0(3,6) = Hmax*kt*t55*t56*t57*t58*t59*t75*t94*t102*(1.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t94*t102*(1.0D0/3.0D0)
      A0(3,7) = t3*t38*t57*t58*t62*t73*t102*(-2.0D0/3.0D0)-t3*t38*t56*t5
     &7*t62*t73*t95*t100*t102*(2.0D0/3.0D0)+t3*t56*t57*t58*t62*t73*t95*t
     &96*t102*(2.0D0/3.0D0)-Hmax*kt*t54*t55*t56*t57*t58*t95*t96*t102*(1.
     &0D0/2.0D0)
      A0(4,1) = Hmax*kt*t55*t56*t57*t58*t59*t71*t75*t104*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t71*t76*t104
      A0(4,2) = Hmax*kt*t55*t56*t57*t58*t59*t75*t81*t104*(3.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t81*t104
      A0(4,3) = Hmax*kt*t55*t56*t57*t58*t59*t75*t85*t104*(3.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t85*t104
      A0(4,4) = t3*t38*t56*t57*t58*t62*t73*(-2.0D0)-Hmax*kt*t55*t56*t57*
     &t58*t59*t75*t88*t104*(3.0D0/4.0D0)+t3*t38*t56*t57*t58*t59*t62*t76*
     &t88*t104-1.0D0
      A0(4,5) = Hmax*kt*t55*t56*t57*t58*t59*t75*t91*t104*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t76*t91*t104
      A0(4,6) = Hmax*kt*t55*t56*t57*t58*t59*t75*t94*t104*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t76*t94*t104
      A0(4,7) = t3*t38*t57*t58*t62*t73*t104*2.0D0+t3*t38*t56*t57*t62*t73
     &*t95*t100*t104*2.0D0-t3*t56*t57*t58*t62*t73*t95*t96*t104*2.0D0+Hma
     &x*kt*t54*t55*t56*t57*t58*t95*t96*t104*(3.0D0/2.0D0)
      A0(5,1) = Hmax*kt*t55*t56*t57*t58*t59*t71*t75*t105*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t71*t76*t105
      A0(5,2) = Hmax*kt*t55*t56*t57*t58*t59*t75*t81*t105*(3.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t81*t105
      A0(5,3) = Hmax*kt*t55*t56*t57*t58*t59*t75*t85*t105*(3.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t85*t105
      A0(5,4) = Hmax*kt*t55*t56*t57*t58*t59*t75*t88*t105*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t76*t88*t105
      A0(5,5) = t3*t38*t56*t57*t58*t62*t73*(-2.0D0)-Hmax*kt*t55*t56*t57*
     &t58*t59*t75*t91*t105*(3.0D0/4.0D0)+t3*t38*t56*t57*t58*t59*t62*t76*
     &t91*t105-1.0D0
      A0(5,6) = Hmax*kt*t55*t56*t57*t58*t59*t75*t94*t105*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t76*t94*t105
      A0(5,7) = t3*t38*t57*t58*t62*t73*t105*2.0D0+t3*t38*t56*t57*t62*t73
     &*t95*t100*t105*2.0D0-t3*t56*t57*t58*t62*t73*t95*t96*t105*2.0D0+Hma
     &x*kt*t54*t55*t56*t57*t58*t95*t96*t105*(3.0D0/2.0D0)
      A0(6,1) = Hmax*kt*t55*t56*t57*t58*t59*t71*t75*t106*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t71*t76*t106
      A0(6,2) = Hmax*kt*t55*t56*t57*t58*t59*t75*t81*t106*(3.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t81*t106
      A0(6,3) = Hmax*kt*t55*t56*t57*t58*t59*t75*t85*t106*(3.0D0/4.0D0)-t
     &3*t38*t56*t57*t58*t59*t62*t76*t85*t106
      A0(6,4) = Hmax*kt*t55*t56*t57*t58*t59*t75*t88*t106*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t76*t88*t106
      A0(6,5) = Hmax*kt*t55*t56*t57*t58*t59*t75*t91*t106*(-3.0D0/4.0D0)+
     &t3*t38*t56*t57*t58*t59*t62*t76*t91*t106
      A0(6,6) = t3*t38*t56*t57*t58*t62*t73*(-2.0D0)-Hmax*kt*t55*t56*t57*
     &t58*t59*t75*t94*t106*(3.0D0/4.0D0)+t3*t38*t56*t57*t58*t59*t62*t76*
     &t94*t106-1.0D0
      A0(6,7) = t3*t38*t57*t58*t62*t73*t106*2.0D0+t3*t38*t56*t57*t62*t73
     &*t95*t100*t106*2.0D0-t3*t56*t57*t58*t62*t73*t95*t96*t106*2.0D0+Hma
     &x*kt*t54*t55*t56*t57*t58*t95*t96*t106*(3.0D0/2.0D0)
      A0(7,1) = t159+t163+t169+t172-t120*t150-t120*t156+t152*t155-t58*t1
     &10*t113*t160*(1.0D0/2.0D0)+t3*t38*t57*t62*t73*t100*t110*(e1-et1+t1
     &35+t136+t137+t138+t139-Tn*alpha-e1*v-et2*v-et3*v-Tn*alpha*v)*(4.0D
     &0/3.0D0)-Dc*t3*t38*t57*t58*t62*t73*t126*(2.0D0/3.0D0)-Dc*t3*t38*t5
     &7*t58*t62*t73*t127*(2.0D0/3.0D0)+Dc*t3*t38*t57*t58*t62*t73*t129*(4
     &.0D0/3.0D0)-t3*t38*t57*t62*t73*t100*t110*(e2-et2+t135-t136+t137-t1
     &38+t139+t140+t141-Tn*alpha-et3*v-Tn*alpha*v)*(2.0D0/3.0D0)-t3*t38*
     &t57*t62*t73*t100*t110*t147*(2.0D0/3.0D0)-Hmax*kt*t55*t59*t71*t75*t
     &100*t142*t143*(3.0D0/8.0D0)-Hmax*kt*t55*t59*t71*t75*t100*t142*t144
     &*(3.0D0/8.0D0)-Hmax*kt*t55*t59*t71*t75*t100*t142*t145*(3.0D0/8.0D0
     &)+t3*t38*t59*t62*t71*t76*t100*t142*t143*(1.0D0/2.0D0)+t3*t38*t59*t
     &62*t71*t76*t100*t142*t144*(1.0D0/2.0D0)+t3*t38*t59*t62*t71*t76*t10
     &0*t142*t145*(1.0D0/2.0D0)-t3*t38*t57*t62*t73*t98*t100*t110*v*(2.0D
     &0/3.0D0)-t3*t38*t57*t62*t73*t100*t102*t110*v*(2.0D0/3.0D0)+t3*t38*
     &t57*t62*t73*t100*t110*t113*(e2+e3-et2-et3-t63+t64)*(2.0D0/3.0D0)+D
     &c*Hmax*kt*t55*t59*t71*t75*t100*t142*t143*(3.0D0/8.0D0)+Dc*Hmax*kt*
     &t55*t59*t71*t75*t100*t142*t144*(3.0D0/8.0D0)+Dc*Hmax*kt*t55*t59*t7
     &1*t75*t100*t142*t145*(3.0D0/8.0D0)-Dc*t3*t38*t57*t62*t72*t73*t100*
     &t110*t113*(2.0D0/3.0D0)-Dc*t3*t38*t59*t62*t71*t76*t100*t142*t143*(
     &1.0D0/2.0D0)-Dc*t3*t38*t59*t62*t71*t76*t100*t142*t144*(1.0D0/2.0D0
     &)-Dc*t3*t38*t59*t62*t71*t76*t100*t142*t145*(1.0D0/2.0D0)-Hmax*kt*t
     &55*t57*t59*t71*t75*t100*t102*t110*t147*(1.0D0/4.0D0)+Dc*t3*t38*t57
     &*t58*t59*t62*t71*t76*t98*t127*(1.0D0/3.0D0)+Dc*t3*t38*t57*t58*t59*
     &t62*t71*t76*t102*t126*(1.0D0/3.0D0)+t3*t38*t57*t59*t62*t71*t76*t10
     &0*t102*t110*t147*(1.0D0/3.0D0)+t3*t38*t57*t59*t62*t71*t76*t98*t100
     &*t110*t162*(1.0D0/3.0D0)-Hmax*kt*t55*t57*t59*t71*t72*t75*t100*t110
     &*(e1-et1+t135+t136+t137+t138+t139-t140-t141-t146-Tn*alpha-Tn*alpha
     &*v)*(1.0D0/4.0D0)-Hmax*kt*t55*t57*t59*t71*t75*t98*t100*t110*(e2-et
     &2+t135-t136+t137-t138+t139+t140+t141-t146-Tn*alpha-Tn*alpha*v)*(1.
     &0D0/4.0D0)+Dc*t3*t38*t57*t58*t59*t62*t71*t76*t129*(e2+e3-et2-et3-t
     &63+t64)*(1.0D0/3.0D0)+t3*t38*t57*t59*t62*t71*t76*t100*t110*t161*(e
     &2+e3-et2-et3-t63+t64)*(1.0D0/3.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t71
     &*t72*t75*t129*(1.0D0/4.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t71*t75*t98
     &*t127*(1.0D0/4.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t71*t75*t102*t126*(
     &1.0D0/4.0D0)
      A0(7,2) = t159+t163+t171+t173-t120*t150-t120*t152+t155*t156-t58*t1
     &10*t113*t157*(1.0D0/2.0D0)-Dc*t3*t38*t57*t58*t62*t73*t126*(2.0D0/3
     &.0D0)+Dc*t3*t38*t57*t58*t62*t73*t127*(4.0D0/3.0D0)-Dc*t3*t38*t57*t
     &58*t62*t73*t129*(2.0D0/3.0D0)-t3*t38*t57*t62*t73*t100*t110*t147*(2
     &.0D0/3.0D0)-t3*t38*t57*t62*t73*t100*t110*t161*(2.0D0/3.0D0)+t3*t38
     &*t57*t62*t73*t100*t110*t162*(4.0D0/3.0D0)+Hmax*kt*t55*t59*t75*t81*
     &t100*t142*t143*(3.0D0/8.0D0)+Hmax*kt*t55*t59*t75*t81*t100*t142*t14
     &4*(3.0D0/8.0D0)+Hmax*kt*t55*t59*t75*t81*t100*t142*t145*(3.0D0/8.0D
     &0)+t3*t38*t57*t62*t73*t98*t100*t110*t113*(2.0D0/3.0D0)-t3*t38*t59*
     &t62*t76*t81*t100*t142*t143*(1.0D0/2.0D0)-t3*t38*t59*t62*t76*t81*t1
     &00*t142*t144*(1.0D0/2.0D0)-t3*t38*t59*t62*t76*t81*t100*t142*t145*(
     &1.0D0/2.0D0)-t3*t38*t57*t62*t72*t73*t100*t110*v*(2.0D0/3.0D0)-t3*t
     &38*t57*t62*t73*t100*t102*t110*v*(2.0D0/3.0D0)-Dc*Hmax*kt*t55*t59*t
     &75*t81*t100*t142*t143*(3.0D0/8.0D0)-Dc*Hmax*kt*t55*t59*t75*t81*t10
     &0*t142*t144*(3.0D0/8.0D0)-Dc*Hmax*kt*t55*t59*t75*t81*t100*t142*t14
     &5*(3.0D0/8.0D0)-Dc*t3*t38*t57*t62*t73*t98*t100*t110*t113*(2.0D0/3.
     &0D0)+Dc*t3*t38*t59*t62*t76*t81*t100*t142*t143*(1.0D0/2.0D0)+Dc*t3*
     &t38*t59*t62*t76*t81*t100*t142*t144*(1.0D0/2.0D0)+Dc*t3*t38*t59*t62
     &*t76*t81*t100*t142*t145*(1.0D0/2.0D0)+Hmax*kt*t55*t57*t59*t75*t81*
     &t100*t102*t110*t147*(1.0D0/4.0D0)+Hmax*kt*t55*t57*t59*t75*t81*t98*
     &t100*t110*t162*(1.0D0/4.0D0)-Dc*t3*t38*t57*t58*t59*t62*t72*t76*t81
     &*t129*(1.0D0/3.0D0)-Dc*t3*t38*t57*t58*t59*t62*t76*t81*t98*t127*(1.
     &0D0/3.0D0)-Dc*t3*t38*t57*t58*t59*t62*t76*t81*t102*t126*(1.0D0/3.0D
     &0)-t3*t38*t57*t59*t62*t72*t76*t81*t100*t110*t161*(1.0D0/3.0D0)-t3*
     &t38*t57*t59*t62*t76*t81*t100*t102*t110*t147*(1.0D0/3.0D0)-t3*t38*t
     &57*t59*t62*t76*t81*t98*t100*t110*t162*(1.0D0/3.0D0)+Dc*Hmax*kt*t55
     &*t57*t58*t59*t75*t81*t129*(e2+e3-et2-et3-t63+t64)*(1.0D0/4.0D0)+Hm
     &ax*kt*t55*t57*t59*t75*t81*t100*t110*t161*(e2+e3-et2-et3-t63+t64)*(
     &1.0D0/4.0D0)+Dc*Hmax*kt*t55*t57*t58*t59*t75*t81*t98*t127*(1.0D0/4.
     &0D0)+Dc*Hmax*kt*t55*t57*t58*t59*t75*t81*t102*t126*(1.0D0/4.0D0)
      A0(7,3) = t169+t171+t172+t173-t120*t152-t120*t156+t150*t155-t58*t1
     &10*t113*t158*(1.0D0/2.0D0)+Dc*t3*t38*t57*t58*t62*t73*t126*(4.0D0/3
     &.0D0)-Dc*t3*t38*t57*t58*t62*t73*t127*(2.0D0/3.0D0)-Dc*t3*t38*t57*t
     &58*t62*t73*t129*(2.0D0/3.0D0)+t3*t38*t57*t62*t73*t100*t110*t147*(4
     &.0D0/3.0D0)-t3*t38*t57*t62*t73*t100*t110*t161*(2.0D0/3.0D0)-t3*t38
     &*t57*t62*t73*t100*t110*t162*(2.0D0/3.0D0)+Hmax*kt*t55*t59*t75*t85*
     &t100*t142*t143*(3.0D0/8.0D0)+Hmax*kt*t55*t59*t75*t85*t100*t142*t14
     &4*(3.0D0/8.0D0)+Hmax*kt*t55*t59*t75*t85*t100*t142*t145*(3.0D0/8.0D
     &0)+t3*t38*t57*t62*t73*t100*t102*t110*t113*(2.0D0/3.0D0)-t3*t38*t59
     &*t62*t76*t85*t100*t142*t143*(1.0D0/2.0D0)-t3*t38*t59*t62*t76*t85*t
     &100*t142*t144*(1.0D0/2.0D0)-t3*t38*t59*t62*t76*t85*t100*t142*t145*
     &(1.0D0/2.0D0)-t3*t38*t57*t62*t72*t73*t100*t110*v*(2.0D0/3.0D0)-t3*
     &t38*t57*t62*t73*t98*t100*t110*v*(2.0D0/3.0D0)-Dc*Hmax*kt*t55*t59*t
     &75*t85*t100*t142*t143*(3.0D0/8.0D0)-Dc*Hmax*kt*t55*t59*t75*t85*t10
     &0*t142*t144*(3.0D0/8.0D0)-Dc*Hmax*kt*t55*t59*t75*t85*t100*t142*t14
     &5*(3.0D0/8.0D0)-Dc*t3*t38*t57*t62*t73*t100*t102*t110*t113*(2.0D0/3
     &.0D0)+Dc*t3*t38*t59*t62*t76*t85*t100*t142*t143*(1.0D0/2.0D0)+Dc*t3
     &*t38*t59*t62*t76*t85*t100*t142*t144*(1.0D0/2.0D0)+Dc*t3*t38*t59*t6
     &2*t76*t85*t100*t142*t145*(1.0D0/2.0D0)+Hmax*kt*t55*t57*t59*t75*t85
     &*t100*t102*t110*t147*(1.0D0/4.0D0)+Hmax*kt*t55*t57*t59*t75*t85*t98
     &*t100*t110*t162*(1.0D0/4.0D0)-Dc*t3*t38*t57*t58*t59*t62*t72*t76*t8
     &5*t129*(1.0D0/3.0D0)-Dc*t3*t38*t57*t58*t59*t62*t76*t85*t98*t127*(1
     &.0D0/3.0D0)-Dc*t3*t38*t57*t58*t59*t62*t76*t85*t102*t126*(1.0D0/3.0
     &D0)-t3*t38*t57*t59*t62*t72*t76*t85*t100*t110*t161*(1.0D0/3.0D0)-t3
     &*t38*t57*t59*t62*t76*t85*t100*t102*t110*t147*(1.0D0/3.0D0)-t3*t38*
     &t57*t59*t62*t76*t85*t98*t100*t110*t162*(1.0D0/3.0D0)+Dc*Hmax*kt*t5
     &5*t57*t58*t59*t75*t85*t129*(e2+e3-et2-et3-t63+t64)*(1.0D0/4.0D0)+H
     &max*kt*t55*t57*t59*t75*t85*t100*t110*t161*(e2+e3-et2-et3-t63+t64)*
     &(1.0D0/4.0D0)+Dc*Hmax*kt*t55*t57*t58*t59*t75*t85*t98*t127*(1.0D0/4
     &.0D0)+Dc*Hmax*kt*t55*t57*t58*t59*t75*t85*t102*t126*(1.0D0/4.0D0)
      A0(7,4) = t95*t100*t142*t176*t178*(1.0D0/8.0D0)-t3*t38*t62*t73*t10
     &0*t142*t176+Dc*t3*t38*t62*t73*t100*t142*t176-Hmax*kt*t55*t59*t75*t
     &88*t100*t142*t143*(3.0D0/8.0D0)-Hmax*kt*t55*t59*t75*t88*t100*t142*
     &t144*(3.0D0/8.0D0)-Hmax*kt*t55*t59*t75*t88*t100*t142*t145*(3.0D0/8
     &.0D0)+t3*t38*t59*t62*t76*t88*t100*t142*t143*(1.0D0/2.0D0)+t3*t38*t
     &59*t62*t76*t88*t100*t142*t144*(1.0D0/2.0D0)+t3*t38*t59*t62*t76*t88
     &*t100*t142*t145*(1.0D0/2.0D0)+Dc*Hmax*kt*t55*t59*t75*t88*t100*t142
     &*t143*(3.0D0/8.0D0)+Dc*Hmax*kt*t55*t59*t75*t88*t100*t142*t144*(3.0
     &D0/8.0D0)+Dc*Hmax*kt*t55*t59*t75*t88*t100*t142*t145*(3.0D0/8.0D0)-
     &Dc*t3*t38*t59*t62*t76*t88*t100*t142*t143*(1.0D0/2.0D0)-Dc*t3*t38*t
     &59*t62*t76*t88*t100*t142*t144*(1.0D0/2.0D0)-Dc*t3*t38*t59*t62*t76*
     &t88*t100*t142*t145*(1.0D0/2.0D0)-Hmax*kt*t55*t57*t59*t72*t75*t88*t
     &100*t110*t161*(1.0D0/4.0D0)-Hmax*kt*t55*t57*t59*t75*t88*t100*t102*
     &t110*t147*(1.0D0/4.0D0)-Hmax*kt*t55*t57*t59*t75*t88*t98*t100*t110*
     &t162*(1.0D0/4.0D0)+Dc*t3*t38*t57*t58*t59*t62*t76*t88*t98*t127*(1.0
     &D0/3.0D0)+Dc*t3*t38*t57*t58*t59*t62*t76*t88*t102*t126*(1.0D0/3.0D0
     &)+t3*t38*t57*t59*t62*t76*t88*t100*t102*t110*t147*(1.0D0/3.0D0)+t3*
     &t38*t57*t59*t62*t76*t88*t98*t100*t110*t162*(1.0D0/3.0D0)+Dc*t3*t38
     &*t57*t58*t59*t62*t76*t88*t129*(e2+e3-et2-et3-t63+t64)*(1.0D0/3.0D0
     &)+t3*t38*t57*t59*t62*t76*t88*t100*t110*t161*(e2+e3-et2-et3-t63+t64
     &)*(1.0D0/3.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t72*t75*t88*t129*(1.0D0
     &/4.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t75*t88*t98*t127*(1.0D0/4.0D0)-
     &Dc*Hmax*kt*t55*t57*t58*t59*t75*t88*t102*t126*(1.0D0/4.0D0)
      A0(7,5) = t95*t100*t142*t178*t181*(1.0D0/8.0D0)-t3*t38*t62*t73*t10
     &0*t142*t181+Dc*t3*t38*t62*t73*t100*t142*t181-Hmax*kt*t55*t59*t75*t
     &91*t100*t142*t143*(3.0D0/8.0D0)-Hmax*kt*t55*t59*t75*t91*t100*t142*
     &t144*(3.0D0/8.0D0)-Hmax*kt*t55*t59*t75*t91*t100*t142*t145*(3.0D0/8
     &.0D0)+t3*t38*t59*t62*t76*t91*t100*t142*t143*(1.0D0/2.0D0)+t3*t38*t
     &59*t62*t76*t91*t100*t142*t144*(1.0D0/2.0D0)+t3*t38*t59*t62*t76*t91
     &*t100*t142*t145*(1.0D0/2.0D0)+Dc*Hmax*kt*t55*t59*t75*t91*t100*t142
     &*t143*(3.0D0/8.0D0)+Dc*Hmax*kt*t55*t59*t75*t91*t100*t142*t144*(3.0
     &D0/8.0D0)+Dc*Hmax*kt*t55*t59*t75*t91*t100*t142*t145*(3.0D0/8.0D0)-
     &Dc*t3*t38*t59*t62*t76*t91*t100*t142*t143*(1.0D0/2.0D0)-Dc*t3*t38*t
     &59*t62*t76*t91*t100*t142*t144*(1.0D0/2.0D0)-Dc*t3*t38*t59*t62*t76*
     &t91*t100*t142*t145*(1.0D0/2.0D0)-Hmax*kt*t55*t57*t59*t72*t75*t91*t
     &100*t110*t161*(1.0D0/4.0D0)-Hmax*kt*t55*t57*t59*t75*t91*t100*t102*
     &t110*t147*(1.0D0/4.0D0)-Hmax*kt*t55*t57*t59*t75*t91*t98*t100*t110*
     &t162*(1.0D0/4.0D0)+Dc*t3*t38*t57*t58*t59*t62*t76*t91*t98*t127*(1.0
     &D0/3.0D0)+Dc*t3*t38*t57*t58*t59*t62*t76*t91*t102*t126*(1.0D0/3.0D0
     &)+t3*t38*t57*t59*t62*t76*t91*t100*t102*t110*t147*(1.0D0/3.0D0)+t3*
     &t38*t57*t59*t62*t76*t91*t98*t100*t110*t162*(1.0D0/3.0D0)+Dc*t3*t38
     &*t57*t58*t59*t62*t76*t91*t129*(e2+e3-et2-et3-t63+t64)*(1.0D0/3.0D0
     &)+t3*t38*t57*t59*t62*t76*t91*t100*t110*t161*(e2+e3-et2-et3-t63+t64
     &)*(1.0D0/3.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t72*t75*t91*t129*(1.0D0
     &/4.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t75*t91*t98*t127*(1.0D0/4.0D0)-
     &Dc*Hmax*kt*t55*t57*t58*t59*t75*t91*t102*t126*(1.0D0/4.0D0)
      A0(7,6) = t95*t100*t142*t178*t184*(1.0D0/8.0D0)-t3*t38*t62*t73*t10
     &0*t142*t184+Dc*t3*t38*t62*t73*t100*t142*t184-Hmax*kt*t55*t59*t75*t
     &94*t100*t142*t143*(3.0D0/8.0D0)-Hmax*kt*t55*t59*t75*t94*t100*t142*
     &t144*(3.0D0/8.0D0)-Hmax*kt*t55*t59*t75*t94*t100*t142*t145*(3.0D0/8
     &.0D0)+t3*t38*t59*t62*t76*t94*t100*t142*t143*(1.0D0/2.0D0)+t3*t38*t
     &59*t62*t76*t94*t100*t142*t144*(1.0D0/2.0D0)+t3*t38*t59*t62*t76*t94
     &*t100*t142*t145*(1.0D0/2.0D0)+Dc*Hmax*kt*t55*t59*t75*t94*t100*t142
     &*t143*(3.0D0/8.0D0)+Dc*Hmax*kt*t55*t59*t75*t94*t100*t142*t144*(3.0
     &D0/8.0D0)+Dc*Hmax*kt*t55*t59*t75*t94*t100*t142*t145*(3.0D0/8.0D0)-
     &Dc*t3*t38*t59*t62*t76*t94*t100*t142*t143*(1.0D0/2.0D0)-Dc*t3*t38*t
     &59*t62*t76*t94*t100*t142*t144*(1.0D0/2.0D0)-Dc*t3*t38*t59*t62*t76*
     &t94*t100*t142*t145*(1.0D0/2.0D0)-Hmax*kt*t55*t57*t59*t72*t75*t94*t
     &100*t110*t161*(1.0D0/4.0D0)-Hmax*kt*t55*t57*t59*t75*t94*t100*t102*
     &t110*t147*(1.0D0/4.0D0)-Hmax*kt*t55*t57*t59*t75*t94*t98*t100*t110*
     &t162*(1.0D0/4.0D0)+Dc*t3*t38*t57*t58*t59*t62*t76*t94*t98*t127*(1.0
     &D0/3.0D0)+Dc*t3*t38*t57*t58*t59*t62*t76*t94*t102*t126*(1.0D0/3.0D0
     &)+t3*t38*t57*t59*t62*t76*t94*t100*t102*t110*t147*(1.0D0/3.0D0)+t3*
     &t38*t57*t59*t62*t76*t94*t98*t100*t110*t162*(1.0D0/3.0D0)+Dc*t3*t38
     &*t57*t58*t59*t62*t76*t94*t129*(e2+e3-et2-et3-t63+t64)*(1.0D0/3.0D0
     &)+t3*t38*t57*t59*t62*t76*t94*t100*t110*t161*(e2+e3-et2-et3-t63+t64
     &)*(1.0D0/3.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t72*t75*t94*t129*(1.0D0
     &/4.0D0)-Dc*Hmax*kt*t55*t57*t58*t59*t75*t94*t98*t127*(1.0D0/4.0D0)-
     &Dc*Hmax*kt*t55*t57*t58*t59*t75*t94*t102*t126*(1.0D0/4.0D0)
      A0(7,7) = t150*(t190+t194-t95*t191)+t156*(t190+t196-t95*t188)+t152
     &*(t194+t196-t95*t189)+t158*(t197+t198-t95*t100*t110*t113*t116*(1.0
     &D0/2.0D0))+t157*(t197+t199-t95*t100*t110*t113*t117*(1.0D0/2.0D0))+
     &t160*(t198+t199-t95*t100*t110*t113*t114*(1.0D0/2.0D0))-a1*(n2*(-x+
     &1.0D0)**(n2-1.0D0)+n1*x**(n1-1.0D0))*(1.0D0/2.0D0)-t142*t143*t178*
     &t200*t201*(1.0D0/4.0D0)-t142*t144*t178*t200*t201*(1.0D0/4.0D0)-t14
     &2*t145*t178*t200*t201*(1.0D0/4.0D0)-t3*t62*t73*t95*t96*t100*t142*t
     &143-t3*t62*t73*t95*t96*t100*t142*t144-t3*t62*t73*t95*t96*t100*t142
     &*t145+t3*t38*t62*t73*t95*t142*t143*t201*2.0D0+t3*t38*t62*t73*t95*t
     &142*t144*t201*2.0D0+t3*t38*t62*t73*t95*t142*t145*t201*2.0D0+Hmax*k
     &t*t54*t55*t95*t96*t100*t142*t143*(3.0D0/4.0D0)+Hmax*kt*t54*t55*t95
     &*t96*t100*t142*t144*(3.0D0/4.0D0)+Hmax*kt*t54*t55*t95*t96*t100*t14
     &2*t145*(3.0D0/4.0D0)+Dc*t3*t38*t57*t58*t62*t73*t98*t188*(2.0D0/3.0
     &D0)+Dc*t3*t38*t57*t58*t62*t73*t102*t191*(2.0D0/3.0D0)+Dc*t3*t62*t7
     &3*t95*t96*t100*t142*t143+Dc*t3*t62*t73*t95*t96*t100*t142*t144+Dc*t
     &3*t62*t73*t95*t96*t100*t142*t145-Dc*t3*t38*t62*t73*t95*t142*t143*t
     &201*2.0D0-Dc*t3*t38*t62*t73*t95*t142*t144*t201*2.0D0-Dc*t3*t38*t62
     &*t73*t95*t142*t145*t201*2.0D0+Dc*t3*t38*t57*t58*t62*t73*t189*(e2+e
     &3-et2-et3-t63+t64)*(2.0D0/3.0D0)+Dc*t3*t38*t57*t62*t73*t95*t100*t1
     &29*(e2+e3-et2-et3-t63+t64)*(2.0D0/3.0D0)+t3*t38*t57*t62*t73*t95*t1
     &10*t161*t201*(e2+e3-et2-et3-t63+t64)*(4.0D0/3.0D0)-Dc*Hmax*kt*t54*
     &t55*t95*t96*t100*t142*t143*(3.0D0/4.0D0)-Dc*Hmax*kt*t54*t55*t95*t9
     &6*t100*t142*t144*(3.0D0/4.0D0)-Dc*Hmax*kt*t54*t55*t95*t96*t100*t14
     &2*t145*(3.0D0/4.0D0)-Dc*t3*t57*t58*t62*t72*t73*t95*t96*t129*(2.0D0
     &/3.0D0)+Dc*t3*t38*t57*t62*t73*t95*t98*t100*t127*(2.0D0/3.0D0)+Dc*t
     &3*t38*t57*t62*t73*t95*t100*t102*t126*(2.0D0/3.0D0)-Dc*t3*t57*t58*t
     &62*t73*t95*t96*t98*t127*(2.0D0/3.0D0)-Dc*t3*t57*t58*t62*t73*t95*t9
     &6*t102*t126*(2.0D0/3.0D0)-t3*t57*t62*t72*t73*t95*t96*t100*t110*t16
     &1*(2.0D0/3.0D0)-t3*t57*t62*t73*t95*t96*t100*t102*t110*t147*(2.0D0/
     &3.0D0)-t3*t57*t62*t73*t95*t96*t98*t100*t110*t162*(2.0D0/3.0D0)+t3*
     &t38*t57*t62*t73*t95*t102*t110*t147*t201*(4.0D0/3.0D0)+t3*t38*t57*t
     &62*t73*t95*t98*t110*t162*t201*(4.0D0/3.0D0)+Hmax*kt*t54*t55*t57*t9
     &5*t96*t100*t102*t110*t147*(1.0D0/2.0D0)+Hmax*kt*t54*t55*t57*t95*t9
     &6*t98*t100*t110*t162*(1.0D0/2.0D0)+Dc*Hmax*kt*t54*t55*t57*t58*t95*
     &t96*t129*(e2+e3-et2-et3-t63+t64)*(1.0D0/2.0D0)+Hmax*kt*t54*t55*t57
     &*t95*t96*t100*t110*t161*(e2+e3-et2-et3-t63+t64)*(1.0D0/2.0D0)+Dc*H
     &max*kt*t54*t55*t57*t58*t95*t96*t98*t127*(1.0D0/2.0D0)+Dc*Hmax*kt*t
     &54*t55*t57*t58*t95*t96*t102*t126*(1.0D0/2.0D0)


       NR_JAC=A0
       !NR_JAC(13,)=t0
!   **************** New Elemnts *********!**********************************************************************
!   **********************************************************************
   

       elseif (index1==-1)then ! Reverse

      t2 = Sm*x
      t8 = Sa*x
      t3 = Sa+t2-t8
      t4 = v**2
      t5 = t4*2.0D0
      t6 = t5+v-1.0D0
      t7 = 1.0D0/t6
      t9 = 1.0D0/t3
      t10 = Sa-Sm
      t11 = v-1.0D0
      t12 = 1.0D0/t3**2
      t13 = alpha*too
      t14 = e2*v
      t15 = e3*v
      t16 = et1*v
      t17 = alpha*too*v
      t19 = Tn*alpha
      t18 = e2-et2+t13-t19
      t20 = e1*v
      t21 = et2*v
      t22 = e3-et3+t13-t19
      t23 = e1-et1+t13-t19
      t24 = t7*t10*t12*t23*v*(1.0D0/2.0D0)
      t25 = t7*t10*t12*t18*v*(1.0D0/2.0D0)
      t26 = v+1.0D0
      t27 = 1.0D0/t26
      t28 = v*2.0D0
      t29 = t28+2.0D0
      t30 = 1.0D0/t26**2
      t31 = et3*v
      t32 = t10**2
      t33 = 1.0D0/t3**3
      t35 = Tn*alpha*v
      t34 = e3-et3+t13+t14-t15-t16+t17-t19+t20-t21+t31-t35
      t36 = e4-et4
      t37 = e5-et5
      t38 = e6-et6
      A0(1,1) = -1.0D0
      A0(1,7) = lamdat_r1
      A0(2,2) = -1.0D0
      A0(2,7) = lamdat_r2
      A0(3,3) = -1.0D0
      A0(3,7) = lamdat_r3
      A0(4,4) = -1.0D0
      A0(4,7) = lamdat_r4
      A0(5,5) = -1.0D0
      A0(5,7) = lamdat_r5
      A0(6,6) = -1.0D0
      A0(6,7) = lamdat_r6
      A0(7,1) = t25+lamdat_r1*t7*t9*t11-lamdat_r2*t7*t9*v-lamdat_r3*t7*t
     &9*v+t7*t10*t12*(e1-et1+t13+t14+t15+t16+t17-Tn*alpha-e1*v-et2*v-et3
     &*v-Tn*alpha*v)*(1.0D0/2.0D0)-t7*t10*t11*t12*(e1-et1+t13-Tn*alpha)*
     &(1.0D0/2.0D0)+t7*t10*t12*v*(e3-et3+t13-Tn*alpha)*(1.0D0/2.0D0)
      A0(7,2) = t24+t7*t10*t12*(e2-et2+t13-t14+t15-t16+t17+t20+t21-Tn*al
     &pha-et3*v-Tn*alpha*v)*(1.0D0/2.0D0)+lamdat_r2*t7*t9*t11-lamdat_r1*
     &t7*t9*v-lamdat_r3*t7*t9*v-t7*t10*t11*t12*t18*(1.0D0/2.0D0)+t7*t10*
     &t12*t22*v*(1.0D0/2.0D0)
      A0(7,3) = t24+t25+lamdat_r3*t7*t9*t11-lamdat_r1*t7*t9*v-lamdat_r2*
     &t7*t9*v+t7*t10*t12*t34*(1.0D0/2.0D0)-t7*t10*t11*t12*t22*(1.0D0/2.0
     &D0)
      A0(7,4) = lamdat_r4*t9*t27*(1.0D0/2.0D0)-t10*t12*t29*t30*(e4*2.0D0
     &-et4*2.0D0)*(1.0D0/8.0D0)
      A0(7,5) = lamdat_r5*t9*t27*(1.0D0/2.0D0)-t10*t12*t29*t30*(e5*2.0D0
     &-et5*2.0D0)*(1.0D0/8.0D0)
      A0(7,6) = lamdat_r6*t9*t27*(1.0D0/2.0D0)-t10*t12*t29*t30*(e6*2.0D0
     &-et6*2.0D0)*(1.0D0/8.0D0)
      A0(7,7) = a2*(n4*(-x+1.0D0)**(n4-1.0D0)+n3*x**(n3-1.0D0))*(1.0D0/2
     &.0D0)+t29*t30*t32*t33*t36**2*(1.0D0/4.0D0)+t29*t30*t32*t33*t37**2*
     &(1.0D0/4.0D0)+t29*t30*t32*t33*t38**2*(1.0D0/4.0D0)+lamdat_r1*t7*t1
     &0*t12*(e1-et1+t13+t14+t15+t16+t17-t19-t20-t21-t31-t35)+lamdat_r2*t
     &7*t10*t12*(e2-et2+t13-t14+t15-t16+t17-t19+t20+t21-t31-t35)+lamdat_
     &r3*t7*t10*t12*t34-lamdat_r4*t10*t12*t27*t36*(1.0D0/2.0D0)-lamdat_r
     &5*t10*t12*t27*t37*(1.0D0/2.0D0)-lamdat_r6*t10*t12*t27*t38*(1.0D0/2
     &.0D0)-t7*t22*t32*t33*t34-t7*t18*t32*t33*(e2-et2+t13-t14+t15-t16+t1
     &7-t19+t20+t21-t31-Tn*alpha*v)-t7*t23*t32*t33*(e1-et1+t13+t14+t15+t
     &16+t17-t19-t20-t21-t31-Tn*alpha*v)


      NR_JAC=A0


      ENDIF
!DEC$ FREEFORM

end subroutine

subroutine N_R_Residual(PARAM,NPARAM,VAR,NVAR,Tn,index1,NR_RE)

implicit real*8 (t)



!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8  Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,kp
integer NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


! NR Variables
real*8 NR_RE(7,1),A0(7,1)


real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6
real*8 eto1,eto2,eto3,eto4,eto5,eto6,epo1,epo2,epo3,epo4,epo5,epo6
real*8 lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 phi, Tn

integer index1


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

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
	
   

   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)  
   

   

   if(abs(x-xo)<=0)then
       x=xo+1.0e-16_8
   end if

!DEC$ NOFREEFORM   


      if(index1==1)then    !Forward
      t2 = v+1.0D0
      t3 = abs(t2)
      t4 = e1*et2*4.0D0
      t5 = e2*et1*4.0D0
      t6 = e1*et3*4.0D0
      t7 = e3*et1*4.0D0
      t8 = e2*et3*4.0D0
      t9 = e3*et2*4.0D0
      t10 = e1**2
      t11 = t10*4.0D0
      t12 = e2**2
      t13 = t12*4.0D0
      t14 = e3**2
      t15 = t14*4.0D0
      t16 = e4**2
      t17 = t16*3.0D0
      t18 = e5**2
      t19 = t18*3.0D0
      t20 = e6**2
      t21 = t20*3.0D0
      t22 = et1**2
      t23 = t22*4.0D0
      t24 = et2**2
      t25 = t24*4.0D0
      t26 = et3**2
      t27 = t26*4.0D0
      t28 = et4**2
      t29 = t28*3.0D0
      t30 = et5**2
      t31 = t30*3.0D0
      t32 = et6**2
      t33 = t32*3.0D0
      t40 = e1*e2*4.0D0
      t41 = e1*e3*4.0D0
      t42 = e2*e3*4.0D0
      t43 = e1*et1*8.0D0
      t44 = e2*et2*8.0D0
      t45 = e3*et3*8.0D0
      t46 = e4*et4*6.0D0
      t47 = e5*et5*6.0D0
      t48 = e6*et6*6.0D0
      t49 = et1*et2*4.0D0
      t50 = et1*et3*4.0D0
      t51 = et2*et3*4.0D0
      t34 = t4+t5+t6+t7+t8+t9+t11+t13+t15+t17+t19+t21+t23+t25+t27+t29+t3
     &1+t33-t40-t41-t42-t43-t44-t45-t46-t47-t48-t49-t50-t51
      t35 = abs(t34)
      t36 = Sm*x
      t39 = Sa*x
      t37 = Sa+t36-t39
      t38 = abs(t37)
      t52 = 1.0D0/sqrt(t35)
      t53 = Hmax*(3.0D0/2.0D0)
      t54 = 1.0D0/t3
      t55 = sqrt(t35)
      t56 = 1.0D0/t38
      t62 = kt*t54*t55*t56*(1.0D0/2.0D0)
      t57 = exp(-t62)
      t63 = Hmax*t57*(3.0D0/2.0D0)
      t58 = t53-t63
      t59 = x-xo
      t60 = 1.0D0/t2
      t61 = 1.0D0/t37
      t64 = v**2
      t65 = t64*2.0D0
      t66 = t65+v-1.0D0
      t67 = 1.0D0/t66
      t68 = Tn-too
      t69 = alpha*t68
      t70 = Sa-Sm
      t71 = -e1+et1+t69
      t72 = t61*t67*t71*v
      t73 = -e3+et3+t69
      t74 = v-1.0D0
      t75 = -e2+et2+t69
      t76 = t61*t67*t75*v
      t77 = t61*t67*t73*v
      t86 = t61*t67*t74*t75
      t78 = t72+t77-t86
      t84 = t61*t67*t73*t74
      t79 = t72+t76-t84
      t83 = t61*t67*t71*t74
      t80 = t76+t77-t83
      t81 = t70*t80*v
      t82 = t61*t67*t71*v*(1.0D0/2.0D0)
      t85 = t70*t79*v
      t87 = t70*t78*v
      t88 = t61*t67*t75*v*(1.0D0/2.0D0)
      t89 = t61*t67*t73*v*(1.0D0/2.0D0)
      t90 = e4-et4
      t91 = v*2.0D0
      t92 = t91+2.0D0
      t93 = e5-et5
      t94 = 1.0D0/t2**2
      t95 = 1.0D0/t37**2
      t96 = e6-et6
      t97 = t90**2
      t98 = t93**2
      t99 = t96**2
      t100 = et3*2.0D0
      t112 = e3*2.0D0
      t101 = e1+e2-et1-et2+t100-t112
      t102 = et2*2.0D0
      t106 = e2*2.0D0
      t103 = e1+e3-et1-et3+t102-t106
      t104 = e1*2.0D0
      t105 = et1*2.0D0
      t107 = alpha*too
      t108 = e2*v
      t109 = e3*v
      t110 = et1*v
      t111 = alpha*too*v
      t113 = e1*v
      t114 = et2*v
      A0(1,1) = -et1+eto1-t3*t38*t52*t58*t59*t60*t61*(e2+e3-et2-et3-t104
     &+t105)*(2.0D0/3.0D0)
      A0(2,1) = -et2+eto2-t3*t38*t52*t58*t59*t60*t61*t103*(2.0D0/3.0D0)
      A0(3,1) = -et3+eto3-t3*t38*t52*t58*t59*t60*t61*t101*(2.0D0/3.0D0)
      A0(4,1) = -et4+eto4+t3*t38*t52*t58*t59*t60*t61*t90*2.0D0
      A0(5,1) = -et5+eto5+t3*t38*t52*t58*t59*t60*t61*t93*2.0D0
      A0(6,1) = -et6+eto6+t3*t38*t52*t58*t59*t60*t61*t96*2.0D0
      A0(7,1) = -Yo-a3-rduo+Tn*rdso+(t82+t88-t61*t67*t73*t74*(1.0D0/2.0D
     &0))*(t81+t87-t70*t79)+(t82+t89-t61*t67*t74*t75*(1.0D0/2.0D0))*(t81
     &+t85-t70*t78)+(t88+t89-t61*t67*t71*t74*(1.0D0/2.0D0))*(t85+t87-t70
     &*t80)-a1*(-(-x+1.0D0)**n2+x**n1+1.0D0)*(1.0D0/2.0D0)-t70*t92*t94*t
     &95*t97*(1.0D0/8.0D0)-t70*t92*t94*t95*t98*(1.0D0/8.0D0)-t70*t92*t94
     &*t95*t99*(1.0D0/8.0D0)+t3*t38*t52*t58*t94*t95*t97+t3*t38*t52*t58*t
     &94*t95*t98+t3*t38*t52*t58*t94*t95*t99-Dc*t3*t38*t52*t58*t94*t95*t9
     &7-Dc*t3*t38*t52*t58*t94*t95*t98-Dc*t3*t38*t52*t58*t94*t95*t99+Dc*t
     &3*t38*t52*t58*t60*t61*t79*t101*(2.0D0/3.0D0)+Dc*t3*t38*t52*t58*t60
     &*t61*t78*t103*(2.0D0/3.0D0)+t3*t38*t52*t58*t60*t67*t95*t103*(e2-et
     &2+t107-t108+t109-t110+t111+t113+t114-Tn*alpha-et3*v-Tn*alpha*v)*(2
     &.0D0/3.0D0)+t3*t38*t52*t58*t60*t67*t95*(e2+e3-et2-et3-t104+t105)*(
     &e1-et1+t107+t108+t109+t110+t111-Tn*alpha-e1*v-et2*v-et3*v-Tn*alpha
     &*v)*(2.0D0/3.0D0)+t3*t38*t52*t58*t60*t67*t95*t101*(e3-et3+t107+t10
     &8-t109-t110+t111+t113-t114-Tn*alpha+et3*v-Tn*alpha*v)*(2.0D0/3.0D0
     &)+Dc*t3*t38*t52*t58*t60*t61*t80*(e2+e3-et2-et3-t104+t105)*(2.0D0/3
     &.0D0)



 
        NR_RE=A0

       elseif(index1==-1)then ! Reverse
      t2 = x-xo
      t3 = Sm*x
      t15 = Sa*x
      t4 = Sa+t3-t15
      t5 = 1.0D0/t4
      t6 = v**2
      t7 = t6*2.0D0
      t8 = t7+v-1.0D0
      t9 = 1.0D0/t8
      t10 = alpha*too
      t11 = e2*v
      t12 = e3*v
      t13 = et1*v
      t14 = alpha*too*v
      t16 = e1*v
      t17 = et2*v
      t18 = v+1.0D0
      t19 = 1.0D0/t18
      t20 = e4-et4
      t21 = v*2.0D0
      t22 = t21+2.0D0
      t23 = Sa-Sm
      t24 = e5-et5
      t25 = 1.0D0/t18**2
      t26 = 1.0D0/t4**2
      t27 = e6-et6
      t28 = et3*v
      t30 = Tn*alpha
      t29 = e3-et3+t10+t11-t12-t13+t14+t16-t17+t28-t30-Tn*alpha*v
      A0(1,1) = -et1+eto1+lamdat_r1*t2
      A0(2,1) = -et2+eto2+lamdat_r2*t2
      A0(3,1) = -et3+eto3+lamdat_r3*t2
      A0(4,1) = -et4+eto4+lamdat_r4*t2
      A0(5,1) = -et5+eto5+lamdat_r5*t2
      A0(6,1) = -et6+eto6+lamdat_r6*t2
      A0(7,1) = -Yo-a3+rduo-Tn*rdso+a2*(-(-x+1.0D0)**n4+x**n3+1.0D0)*(1.
     &0D0/2.0D0)+lamdat_r2*t5*t9*(e2-et2+t10-t11+t12-t13+t14+t16+t17-Tn*
     &alpha-et3*v-Tn*alpha*v)+lamdat_r3*t5*t9*t29-lamdat_r4*t5*t19*t20*(
     &1.0D0/2.0D0)-lamdat_r5*t5*t19*t24*(1.0D0/2.0D0)-lamdat_r6*t5*t19*t
     &27*(1.0D0/2.0D0)+lamdat_r1*t5*t9*(e1-et1+t10+t11+t12+t13+t14-Tn*al
     &pha-e1*v-et2*v-et3*v-Tn*alpha*v)-t9*t23*t26*t29*(e3-et3+t10-t30)*(
     &1.0D0/2.0D0)+t20**2*t22*t23*t25*t26*(1.0D0/8.0D0)+t22*t23*t24**2*t
     &25*t26*(1.0D0/8.0D0)+t22*t23*t25*t26*t27**2*(1.0D0/8.0D0)-t9*t23*t
     &26*(e1-et1+t10-Tn*alpha)*(e1-et1+t10+t11+t12+t13+t14-t16-t17-t28-T
     &n*alpha-Tn*alpha*v)*(1.0D0/2.0D0)-t9*t23*t26*(e2-et2+t10-Tn*alpha)
     &*(e2-et2+t10-t11+t12-t13+t14+t16+t17-t28-Tn*alpha-Tn*alpha*v)*(1.0
     &D0/2.0D0)
 
 
 	 NR_RE=A0

!DEC$ FREEFORM
      endif

END subroutine

subroutine STRESS_UPDATE(PARAM,NPARAM,VAR,NVAR,Tn,stress)

implicit real*8(t)

!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio
integer    NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR

real*8   Tn,stress(6)


! Changing state variables
real*8 A0(6,1)


real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

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
      t2 = Sm*x
      t14 = Sa*x
      t3 = Sa+t2-t14
      t4 = 1.0D0/t3
      t5 = v**2
      t6 = t5*2.0D0
      t7 = t6+v-1.0D0
      t8 = 1.0D0/t7
      t9 = alpha*too
      t10 = e2*v
      t11 = e3*v
      t12 = et1*v
      t13 = alpha*too*v
      t15 = e1*v
      t16 = et2*v
      t17 = v+1.0D0
      t18 = 1.0D0/t17
      A0(1,1) = -t4*t8*(e1-et1+t9+t10+t11+t12+t13-Tn*alpha-e1*v-et2*v-et3*v-Tn*alpha*v)
      A0(2,1) = -t4*t8*(e2-et2+t9-t10+t11-t12+t13+t15+t16-Tn*alpha-et3*v-Tn*alpha*v)
      A0(3,1) = -t4*t8*(e3-et3+t9+t10-t11-t12+t13+t15-t16-Tn*alpha+et3*v-Tn*alpha*v)
      A0(4,1) = t4*t18*(e4-et4)*(1.0D0/2.0D0)
      A0(5,1) = t4*t18*(e5-et5)*(1.0D0/2.0D0)
      A0(6,1) = t4*t18*(e6-et6)*(1.0D0/2.0D0)


	  stress(1:6)=A0(1:6,1)

end

subroutine JACOBIAN_MATRIX(PARAM,NPARAM,VAR,NVAR,stress,Tn,ddsdde)

implicit real*8(t)

!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio
integer    NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR

real*8   Tn           ! temp+dtemp

! Changing state variables
real*8 stress(6)
real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,stress1,stress2,stress3,stress4,stress5,stress6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 et_tr1,et_tr2,et_tr3,et_tr4,et_tr5,et_tr6,lamdat_r1,lamdat_r2,lamdat_r3,lamdat_r4,lamdat_r5,lamdat_r6

real*8 AAA(6,1),BBB(1,6),LLL(6,6),CCC,A0(6,6),M(6,6),M_inv(6,6),N(6,1),dpds(1,6),dpdx,ddsdde(6,6),S_inv(6,6),AAA_rev(6,1)

CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

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
	


   lamdat_r1=lamdat_r(1)
   lamdat_r2=lamdat_r(2)
   lamdat_r3=lamdat_r(3)
   lamdat_r4=lamdat_r(4)
   lamdat_r5=lamdat_r(5)
   lamdat_r6=lamdat_r(6)  
 
   stress1=stress(1)
   stress2=stress(2)
   stress3=stress(3)
   stress4=stress(4)
   stress5=stress(5)
   stress6=stress(6)



   if(abs(x-xo)<=0)then
       x=xo+1.0e-12_8
   end if







if (transformation == 0) then
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
      t21 = kt*t18
      t19 = exp(-t21)
      t20 = t19-1.0D0
      t36 = stress1*2.0D0
      t22 = stress2+stress3-t36
      t23 = stress2*(1.0D0/3.0D0)
      t24 = stress3*(1.0D0/3.0D0)
      t28 = stress1*(2.0D0/3.0D0)
      t25 = t23+t24-t28
      t26 = x-xo
      t27 = 1.0D0/sqrt(t17)
      t29 = 1.0D0/t17**(3.0D0/2.0D0)
      t40 = stress2*2.0D0
      t30 = stress1+stress3-t40
      t31 = 1.0D0/t17
      t32 = Sa-Sm
      t33 = Hmax*t20*t27*(1.0D0/2.0D0)
      t42 = stress3*2.0D0
      t34 = stress1+stress2-t42
      t39 = t32*x
      t35 = Sa-t39
      t37 = stress1*(1.0D0/3.0D0)
      t41 = stress2*(2.0D0/3.0D0)
      t38 = t24+t37-t41
      t44 = stress3*(2.0D0/3.0D0)
      t43 = t23+t37-t44
      t45 = Hmax*stress4*stress5*t20*t29*9.0D0
      t46 = Hmax*kt*stress4*stress5*t19*t31*9.0D0
      t47 = t45+t46
      t48 = t26*t47
      t49 = v*2.0D0
      t50 = t49+2.0D0
      t51 = t35*t50
      t52 = Hmax*stress4*stress6*t20*t29*9.0D0
      t53 = Hmax*kt*stress4*stress6*t19*t31*9.0D0
      t54 = t52+t53
      t55 = t26*t54
      t56 = Hmax*stress5*stress6*t20*t29*9.0D0
      t57 = Hmax*kt*stress5*stress6*t19*t31*9.0D0
      t58 = t56+t57
      t59 = t26*t58
      A0(1,1) = Sa+t26*(-Hmax*t20*t27+Hmax*t20*t22*t25*t29*(3.0D0/4.0D0)
     &+Hmax*kt*t19*t22*t25*t31*(3.0D0/4.0D0))-t32*x
      A0(1,2) = t26*(t33+Hmax*t20*t25*t29*t30*(3.0D0/4.0D0)+Hmax*kt*t19*
     &t25*t30*t31*(3.0D0/4.0D0))-t35*v
      A0(1,3) = t26*(t33+Hmax*t20*t25*t29*t34*(3.0D0/4.0D0)+Hmax*kt*t19*
     &t25*t31*t34*(3.0D0/4.0D0))-t35*v
      A0(1,4) = -t26*(Hmax*stress4*t20*t25*t29*(9.0D0/2.0D0)+Hmax*kt*str
     &ess4*t19*t25*t31*(9.0D0/2.0D0))
      A0(1,5) = -t26*(Hmax*stress5*t20*t25*t29*(9.0D0/2.0D0)+Hmax*kt*str
     &ess5*t19*t25*t31*(9.0D0/2.0D0))
      A0(1,6) = -t26*(Hmax*stress6*t20*t25*t29*(9.0D0/2.0D0)+Hmax*kt*str
     &ess6*t19*t25*t31*(9.0D0/2.0D0))
      A0(2,1) = t26*(t33+Hmax*t20*t22*t29*t38*(3.0D0/4.0D0)+Hmax*kt*t19*
     &t22*t31*t38*(3.0D0/4.0D0))-t35*v
      A0(2,2) = Sa-t39+t26*(-Hmax*t20*t27+Hmax*t20*t29*t30*t38*(3.0D0/4.
     &0D0)+Hmax*kt*t19*t30*t31*t38*(3.0D0/4.0D0))
      A0(2,3) = t26*(t33+Hmax*t20*t29*t34*t38*(3.0D0/4.0D0)+Hmax*kt*t19*
     &t31*t34*t38*(3.0D0/4.0D0))-t35*v
      A0(2,4) = -t26*(Hmax*stress4*t20*t29*t38*(9.0D0/2.0D0)+Hmax*kt*str
     &ess4*t19*t31*t38*(9.0D0/2.0D0))
      A0(2,5) = -t26*(Hmax*stress5*t20*t29*t38*(9.0D0/2.0D0)+Hmax*kt*str
     &ess5*t19*t31*t38*(9.0D0/2.0D0))
      A0(2,6) = -t26*(Hmax*stress6*t20*t29*t38*(9.0D0/2.0D0)+Hmax*kt*str
     &ess6*t19*t31*t38*(9.0D0/2.0D0))
      A0(3,1) = t26*(t33+Hmax*t20*t22*t29*t43*(3.0D0/4.0D0)+Hmax*kt*t19*
     &t22*t31*t43*(3.0D0/4.0D0))-t35*v
      A0(3,2) = t26*(t33+Hmax*t20*t29*t30*t43*(3.0D0/4.0D0)+Hmax*kt*t19*
     &t30*t31*t43*(3.0D0/4.0D0))-t35*v
      A0(3,3) = Sa-t39+t26*(-Hmax*t20*t27+Hmax*t20*t29*t34*t43*(3.0D0/4.
     &0D0)+Hmax*kt*t19*t31*t34*t43*(3.0D0/4.0D0))
      A0(3,4) = -t26*(Hmax*stress4*t20*t29*t43*(9.0D0/2.0D0)+Hmax*kt*str
     &ess4*t19*t31*t43*(9.0D0/2.0D0))
      A0(3,5) = -t26*(Hmax*stress5*t20*t29*t43*(9.0D0/2.0D0)+Hmax*kt*str
     &ess5*t19*t31*t43*(9.0D0/2.0D0))
      A0(3,6) = -t26*(Hmax*stress6*t20*t29*t43*(9.0D0/2.0D0)+Hmax*kt*str
     &ess6*t19*t31*t43*(9.0D0/2.0D0))
      A0(4,1) = -t26*(Hmax*stress4*t20*t22*t29*(3.0D0/2.0D0)+Hmax*kt*str
     &ess4*t19*t22*t31*(3.0D0/2.0D0))
      A0(4,2) = -t26*(Hmax*stress4*t20*t29*t30*(3.0D0/2.0D0)+Hmax*kt*str
     &ess4*t19*t30*t31*(3.0D0/2.0D0))
      A0(4,3) = -t26*(Hmax*stress4*t20*t29*t34*(3.0D0/2.0D0)+Hmax*kt*str
     &ess4*t19*t31*t34*(3.0D0/2.0D0))
      A0(4,4) = t51+t26*(Hmax*t20*t27*(-3.0D0)+Hmax*t11*t20*t29*9.0D0+Hm
     &ax*kt*t11*t19*t31*9.0D0)
      A0(4,5) = t48
      A0(4,6) = t55
      A0(5,1) = -t26*(Hmax*stress5*t20*t22*t29*(3.0D0/2.0D0)+Hmax*kt*str
     &ess5*t19*t22*t31*(3.0D0/2.0D0))
      A0(5,2) = -t26*(Hmax*stress5*t20*t29*t30*(3.0D0/2.0D0)+Hmax*kt*str
     &ess5*t19*t30*t31*(3.0D0/2.0D0))
      A0(5,3) = -t26*(Hmax*stress5*t20*t29*t34*(3.0D0/2.0D0)+Hmax*kt*str
     &ess5*t19*t31*t34*(3.0D0/2.0D0))
      A0(5,4) = t48
      A0(5,5) = t51+t26*(Hmax*t20*t27*(-3.0D0)+Hmax*t13*t20*t29*9.0D0+Hm
     &ax*kt*t13*t19*t31*9.0D0)
      A0(5,6) = t59
      A0(6,1) = -t26*(Hmax*stress6*t20*t22*t29*(3.0D0/2.0D0)+Hmax*kt*str
     &ess6*t19*t22*t31*(3.0D0/2.0D0))
      A0(6,2) = -t26*(Hmax*stress6*t20*t29*t30*(3.0D0/2.0D0)+Hmax*kt*str
     &ess6*t19*t30*t31*(3.0D0/2.0D0))
      A0(6,3) = -t26*(Hmax*stress6*t20*t29*t34*(3.0D0/2.0D0)+Hmax*kt*str
     &ess6*t19*t31*t34*(3.0D0/2.0D0))
      A0(6,4) = t55
      A0(6,5) = t59
      A0(6,6) = t51+t26*(Hmax*t20*t27*(-3.0D0)+Hmax*t15*t20*t29*9.0D0+Hm
     &ax*kt*t15*t19*t31*9.0D0)

 
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
      t19 = stress3*t2*v
      t20 = sqrt(t18)
      t27 = kt*t20
      t21 = exp(-t27)
      t22 = t21-1.0D0
      t23 = stress3*(1.0D0/3.0D0)
      t24 = 1.0D0/sqrt(t18)
      t25 = stress1*t2*v
      t26 = stress2*t2*v
      t28 = stress1*(1.0D0/3.0D0)
      t29 = stress2*(1.0D0/3.0D0)
      t30 = v*2.0D0
      t31 = t30+2.0D0
      A0(1,1) = t19+t26-stress1*t2+Hmax*t22*t24*(stress1*(-2.0D0/3.0D0)+
     &t23+t29)*(3.0D0/2.0D0)
      A0(2,1) = t19+t25-stress2*t2+Hmax*t22*t24*(stress2*(-2.0D0/3.0D0)+
     &t23+t28)*(3.0D0/2.0D0)
      A0(3,1) = t25+t26-stress3*t2+Hmax*t22*t24*(stress3*(-2.0D0/3.0D0)+
     &t28+t29)*(3.0D0/2.0D0)
      A0(4,1) = -stress4*t2*t31-Hmax*stress4*t22*t24*3.0D0
      A0(5,1) = -stress5*t2*t31-Hmax*stress5*t22*t24*3.0D0
      A0(6,1) = -stress6*t2*t31-Hmax*stress6*t22*t24*3.0D0


 
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
      t23 = kt*t19
      t20 = exp(-t23)
      t21 = t20-1.0D0
      t22 = 1.0D0/sqrt(t18)
      t26 = stress1*2.0D0
      t24 = stress2+stress3-t26
      t25 = 1.0D0/t18**(3.0D0/2.0D0)
      t27 = stress2*(1.0D0/3.0D0)
      t28 = stress3*(1.0D0/3.0D0)
      t31 = stress1*(2.0D0/3.0D0)
      t29 = t27+t28-t31
      t30 = 1.0D0/t18
      t32 = stress1*(1.0D0/3.0D0)
      t35 = stress2*(2.0D0/3.0D0)
      t33 = t28+t32-t35
      t36 = stress3*(2.0D0/3.0D0)
      t34 = t27+t32-t36
      t37 = stress3*t2*v
      t38 = Hmax*stress3*t21*t22*(1.0D0/2.0D0)
      t40 = stress2*2.0D0
      t39 = stress1+stress3-t40
      t41 = stress1*t2*v
      t42 = stress2*t2*v
      t43 = Hmax*stress1*t21*t22*(1.0D0/2.0D0)
      t44 = Hmax*stress2*t21*t22*(1.0D0/2.0D0)
      t46 = stress3*2.0D0
      t45 = stress1+stress2-t46
      t47 = v*2.0D0
      t48 = t47+2.0D0
      A0(1,1) = t37+t38+t42+t44-stress1*t2-Hmax*stress1*t21*t22+Hmax*t21
     &*t22*t29*(3.0D0/2.0D0)+Dc*Hmax*stress1*t21*t22-Dc*Hmax*stress2*t21
     &*t22*(1.0D0/2.0D0)-Dc*Hmax*stress3*t21*t22*(1.0D0/2.0D0)-Dc*Hmax*t
     &21*t22*t29*(3.0D0/2.0D0)-Hmax*t12*t21*t24*t25*(3.0D0/2.0D0)-Hmax*t
     &14*t21*t24*t25*(3.0D0/2.0D0)-Hmax*t16*t21*t24*t25*(3.0D0/2.0D0)+Dc
     &*Hmax*t12*t21*t24*t25*(3.0D0/2.0D0)+Dc*Hmax*t14*t21*t24*t25*(3.0D0
     &/2.0D0)+Dc*Hmax*t16*t21*t24*t25*(3.0D0/2.0D0)-Hmax*kt*t12*t20*t24*
     &t30*(3.0D0/2.0D0)-Hmax*kt*t14*t20*t24*t30*(3.0D0/2.0D0)-Hmax*kt*t1
     &6*t20*t24*t30*(3.0D0/2.0D0)+Hmax*stress1*t21*t24*t25*t29*(3.0D0/4.
     &0D0)+Hmax*stress2*t21*t24*t25*t33*(3.0D0/4.0D0)+Hmax*stress3*t21*t
     &24*t25*t34*(3.0D0/4.0D0)+Dc*Hmax*kt*t12*t20*t24*t30*(3.0D0/2.0D0)+
     &Dc*Hmax*kt*t14*t20*t24*t30*(3.0D0/2.0D0)+Dc*Hmax*kt*t16*t20*t24*t3
     &0*(3.0D0/2.0D0)-Dc*Hmax*stress1*t21*t24*t25*t29*(3.0D0/4.0D0)-Dc*H
     &max*stress2*t21*t24*t25*t33*(3.0D0/4.0D0)-Dc*Hmax*stress3*t21*t24*
     &t25*t34*(3.0D0/4.0D0)+Hmax*kt*stress1*t20*t24*t29*t30*(3.0D0/4.0D0
     &)+Hmax*kt*stress2*t20*t24*t30*t33*(3.0D0/4.0D0)+Hmax*kt*stress3*t2
     &0*t24*t30*t34*(3.0D0/4.0D0)-Dc*Hmax*kt*stress1*t20*t24*t29*t30*(3.
     &0D0/4.0D0)-Dc*Hmax*kt*stress2*t20*t24*t30*t33*(3.0D0/4.0D0)-Dc*Hma
     &x*kt*stress3*t20*t24*t30*t34*(3.0D0/4.0D0)
      A0(1,2) = t37+t38+t41+t43-stress2*t2-Hmax*stress2*t21*t22+Hmax*t21
     &*t22*t33*(3.0D0/2.0D0)-Dc*Hmax*stress1*t21*t22*(1.0D0/2.0D0)+Dc*Hm
     &ax*stress2*t21*t22-Dc*Hmax*stress3*t21*t22*(1.0D0/2.0D0)-Dc*Hmax*t
     &21*t22*t33*(3.0D0/2.0D0)-Hmax*t12*t21*t25*t39*(3.0D0/2.0D0)-Hmax*t
     &14*t21*t25*t39*(3.0D0/2.0D0)-Hmax*t16*t21*t25*t39*(3.0D0/2.0D0)+Dc
     &*Hmax*t12*t21*t25*t39*(3.0D0/2.0D0)+Dc*Hmax*t14*t21*t25*t39*(3.0D0
     &/2.0D0)+Dc*Hmax*t16*t21*t25*t39*(3.0D0/2.0D0)-Hmax*kt*t12*t20*t30*
     &t39*(3.0D0/2.0D0)-Hmax*kt*t14*t20*t30*t39*(3.0D0/2.0D0)-Hmax*kt*t1
     &6*t20*t30*t39*(3.0D0/2.0D0)+Hmax*stress1*t21*t25*t29*t39*(3.0D0/4.
     &0D0)+Hmax*stress2*t21*t25*t33*t39*(3.0D0/4.0D0)+Hmax*stress3*t21*t
     &25*t34*t39*(3.0D0/4.0D0)+Dc*Hmax*kt*t12*t20*t30*t39*(3.0D0/2.0D0)+
     &Dc*Hmax*kt*t14*t20*t30*t39*(3.0D0/2.0D0)+Dc*Hmax*kt*t16*t20*t30*t3
     &9*(3.0D0/2.0D0)-Dc*Hmax*stress1*t21*t25*t29*t39*(3.0D0/4.0D0)-Dc*H
     &max*stress2*t21*t25*t33*t39*(3.0D0/4.0D0)-Dc*Hmax*stress3*t21*t25*
     &t34*t39*(3.0D0/4.0D0)+Hmax*kt*stress1*t20*t29*t30*t39*(3.0D0/4.0D0
     &)+Hmax*kt*stress2*t20*t30*t33*t39*(3.0D0/4.0D0)+Hmax*kt*stress3*t2
     &0*t30*t34*t39*(3.0D0/4.0D0)-Dc*Hmax*kt*stress1*t20*t29*t30*t39*(3.
     &0D0/4.0D0)-Dc*Hmax*kt*stress2*t20*t30*t33*t39*(3.0D0/4.0D0)-Dc*Hma
     &x*kt*stress3*t20*t30*t34*t39*(3.0D0/4.0D0)
      A0(1,3) = t41+t42+t43+t44-stress3*t2-Hmax*stress3*t21*t22+Hmax*t21
     &*t22*t34*(3.0D0/2.0D0)-Dc*Hmax*stress1*t21*t22*(1.0D0/2.0D0)-Dc*Hm
     &ax*stress2*t21*t22*(1.0D0/2.0D0)+Dc*Hmax*stress3*t21*t22-Dc*Hmax*t
     &21*t22*t34*(3.0D0/2.0D0)-Hmax*t12*t21*t25*t45*(3.0D0/2.0D0)-Hmax*t
     &14*t21*t25*t45*(3.0D0/2.0D0)-Hmax*t16*t21*t25*t45*(3.0D0/2.0D0)+Dc
     &*Hmax*t12*t21*t25*t45*(3.0D0/2.0D0)+Dc*Hmax*t14*t21*t25*t45*(3.0D0
     &/2.0D0)+Dc*Hmax*t16*t21*t25*t45*(3.0D0/2.0D0)-Hmax*kt*t12*t20*t30*
     &t45*(3.0D0/2.0D0)-Hmax*kt*t14*t20*t30*t45*(3.0D0/2.0D0)-Hmax*kt*t1
     &6*t20*t30*t45*(3.0D0/2.0D0)+Hmax*stress1*t21*t25*t29*t45*(3.0D0/4.
     &0D0)+Hmax*stress2*t21*t25*t33*t45*(3.0D0/4.0D0)+Hmax*stress3*t21*t
     &25*t34*t45*(3.0D0/4.0D0)+Dc*Hmax*kt*t12*t20*t30*t45*(3.0D0/2.0D0)+
     &Dc*Hmax*kt*t14*t20*t30*t45*(3.0D0/2.0D0)+Dc*Hmax*kt*t16*t20*t30*t4
     &5*(3.0D0/2.0D0)-Dc*Hmax*stress1*t21*t25*t29*t45*(3.0D0/4.0D0)-Dc*H
     &max*stress2*t21*t25*t33*t45*(3.0D0/4.0D0)-Dc*Hmax*stress3*t21*t25*
     &t34*t45*(3.0D0/4.0D0)+Hmax*kt*stress1*t20*t29*t30*t45*(3.0D0/4.0D0
     &)+Hmax*kt*stress2*t20*t30*t33*t45*(3.0D0/4.0D0)+Hmax*kt*stress3*t2
     &0*t30*t34*t45*(3.0D0/4.0D0)-Dc*Hmax*kt*stress1*t20*t29*t30*t45*(3.
     &0D0/4.0D0)-Dc*Hmax*kt*stress2*t20*t30*t33*t45*(3.0D0/4.0D0)-Dc*Hma
     &x*kt*stress3*t20*t30*t34*t45*(3.0D0/4.0D0)
      A0(1,4) = -stress4*t2*t48-Hmax*stress4*t21*t22*6.0D0+Dc*Hmax*stres
     &s4*t21*t22*6.0D0+Hmax*stress4*t12*t21*t25*9.0D0+Hmax*stress4*t14*t
     &21*t25*9.0D0+Hmax*stress4*t16*t21*t25*9.0D0-Dc*Hmax*stress4*t12*t2
     &1*t25*9.0D0-Dc*Hmax*stress4*t14*t21*t25*9.0D0-Dc*Hmax*stress4*t16*
     &t21*t25*9.0D0+Hmax*kt*stress4*t12*t20*t30*9.0D0+Hmax*kt*stress4*t1
     &4*t20*t30*9.0D0+Hmax*kt*stress4*t16*t20*t30*9.0D0-Hmax*stress1*str
     &ess4*t21*t25*t29*(9.0D0/2.0D0)-Hmax*stress2*stress4*t21*t25*t33*(9
     &.0D0/2.0D0)-Hmax*stress3*stress4*t21*t25*t34*(9.0D0/2.0D0)-Dc*Hmax
     &*kt*stress4*t12*t20*t30*9.0D0-Dc*Hmax*kt*stress4*t14*t20*t30*9.0D0
     &-Dc*Hmax*kt*stress4*t16*t20*t30*9.0D0+Dc*Hmax*stress1*stress4*t21*
     &t25*t29*(9.0D0/2.0D0)+Dc*Hmax*stress2*stress4*t21*t25*t33*(9.0D0/2
     &.0D0)+Dc*Hmax*stress3*stress4*t21*t25*t34*(9.0D0/2.0D0)-Hmax*kt*st
     &ress1*stress4*t20*t29*t30*(9.0D0/2.0D0)-Hmax*kt*stress2*stress4*t2
     &0*t30*t33*(9.0D0/2.0D0)-Hmax*kt*stress3*stress4*t20*t30*t34*(9.0D0
     &/2.0D0)+Dc*Hmax*kt*stress1*stress4*t20*t29*t30*(9.0D0/2.0D0)+Dc*Hm
     &ax*kt*stress2*stress4*t20*t30*t33*(9.0D0/2.0D0)+Dc*Hmax*kt*stress3
     &*stress4*t20*t30*t34*(9.0D0/2.0D0)
      A0(1,5) = -stress5*t2*t48-Hmax*stress5*t21*t22*6.0D0+Dc*Hmax*stres
     &s5*t21*t22*6.0D0+Hmax*stress5*t12*t21*t25*9.0D0+Hmax*stress5*t14*t
     &21*t25*9.0D0+Hmax*stress5*t16*t21*t25*9.0D0-Dc*Hmax*stress5*t12*t2
     &1*t25*9.0D0-Dc*Hmax*stress5*t14*t21*t25*9.0D0-Dc*Hmax*stress5*t16*
     &t21*t25*9.0D0+Hmax*kt*stress5*t12*t20*t30*9.0D0+Hmax*kt*stress5*t1
     &4*t20*t30*9.0D0+Hmax*kt*stress5*t16*t20*t30*9.0D0-Hmax*stress1*str
     &ess5*t21*t25*t29*(9.0D0/2.0D0)-Hmax*stress2*stress5*t21*t25*t33*(9
     &.0D0/2.0D0)-Hmax*stress3*stress5*t21*t25*t34*(9.0D0/2.0D0)-Dc*Hmax
     &*kt*stress5*t12*t20*t30*9.0D0-Dc*Hmax*kt*stress5*t14*t20*t30*9.0D0
     &-Dc*Hmax*kt*stress5*t16*t20*t30*9.0D0+Dc*Hmax*stress1*stress5*t21*
     &t25*t29*(9.0D0/2.0D0)+Dc*Hmax*stress2*stress5*t21*t25*t33*(9.0D0/2
     &.0D0)+Dc*Hmax*stress3*stress5*t21*t25*t34*(9.0D0/2.0D0)-Hmax*kt*st
     &ress1*stress5*t20*t29*t30*(9.0D0/2.0D0)-Hmax*kt*stress2*stress5*t2
     &0*t30*t33*(9.0D0/2.0D0)-Hmax*kt*stress3*stress5*t20*t30*t34*(9.0D0
     &/2.0D0)+Dc*Hmax*kt*stress1*stress5*t20*t29*t30*(9.0D0/2.0D0)+Dc*Hm
     &ax*kt*stress2*stress5*t20*t30*t33*(9.0D0/2.0D0)+Dc*Hmax*kt*stress3
     &*stress5*t20*t30*t34*(9.0D0/2.0D0)
      A0(1,6) = -stress6*t2*t48-Hmax*stress6*t21*t22*6.0D0+Dc*Hmax*stres
     &s6*t21*t22*6.0D0+Hmax*stress6*t12*t21*t25*9.0D0+Hmax*stress6*t14*t
     &21*t25*9.0D0+Hmax*stress6*t16*t21*t25*9.0D0-Dc*Hmax*stress6*t12*t2
     &1*t25*9.0D0-Dc*Hmax*stress6*t14*t21*t25*9.0D0-Dc*Hmax*stress6*t16*
     &t21*t25*9.0D0+Hmax*kt*stress6*t12*t20*t30*9.0D0+Hmax*kt*stress6*t1
     &4*t20*t30*9.0D0+Hmax*kt*stress6*t16*t20*t30*9.0D0-Hmax*stress1*str
     &ess6*t21*t25*t29*(9.0D0/2.0D0)-Hmax*stress2*stress6*t21*t25*t33*(9
     &.0D0/2.0D0)-Hmax*stress3*stress6*t21*t25*t34*(9.0D0/2.0D0)-Dc*Hmax
     &*kt*stress6*t12*t20*t30*9.0D0-Dc*Hmax*kt*stress6*t14*t20*t30*9.0D0
     &-Dc*Hmax*kt*stress6*t16*t20*t30*9.0D0+Dc*Hmax*stress1*stress6*t21*
     &t25*t29*(9.0D0/2.0D0)+Dc*Hmax*stress2*stress6*t21*t25*t33*(9.0D0/2
     &.0D0)+Dc*Hmax*stress3*stress6*t21*t25*t34*(9.0D0/2.0D0)-Hmax*kt*st
     &ress1*stress6*t20*t29*t30*(9.0D0/2.0D0)-Hmax*kt*stress2*stress6*t2
     &0*t30*t33*(9.0D0/2.0D0)-Hmax*kt*stress3*stress6*t20*t30*t34*(9.0D0
     &/2.0D0)+Dc*Hmax*kt*stress1*stress6*t20*t29*t30*(9.0D0/2.0D0)+Dc*Hm
     &ax*kt*stress2*stress6*t20*t30*t33*(9.0D0/2.0D0)+Dc*Hmax*kt*stress3
     &*stress6*t20*t30*t34*(9.0D0/2.0D0)


 
       dpds(1,1:6)=A0(1,1:6)
       A0=0
 !********************************************************************************************************
 !********************************************************************************************************
      t0 = a1*(n2*(-x+1.0D0)**(n2-1.0D0)+n1*x**(n1-1.0D0))*(-1
     &.0D0/2.0D0)

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
      t2 = Sa-Sm
      t4 = t2*x
      t3 = Sa-t4
      t5 = v*2.0D0
      t6 = t5+2.0D0
      t7 = t3*t6
      A0(1,1) = t3
      A0(1,2) = -t3*v
      A0(1,3) = -t3*v
      A0(2,1) = -t3*v
      A0(2,2) = t3
      A0(2,3) = -t3*v
      A0(3,1) = -t3*v
      A0(3,2) = -t3*v
      A0(3,3) = t3
      A0(4,4) = t7
      A0(5,5) = t7
      A0(6,6) = t7



 
 	   M_INV=A0
 	   A0=0
 	   CALL inverse(M_INV,M,6)
 !********************************************************************************************************
 !********************************************************************************************************
      t2 = Sa-Sm
      t3 = stress3*t2*v
      t4 = stress1*t2*v
      t5 = stress2*t2*v
      t6 = v*2.0D0
      t7 = t6+2.0D0
      A0(1,1) = lamdat_r1+t3+t5-stress1*t2
      A0(2,1) = lamdat_r2+t3+t4-stress2*t2
      A0(3,1) = lamdat_r3+t4+t5-stress3*t2
      A0(4,1) = lamdat_r4-stress4*t2*t7
      A0(5,1) = lamdat_r5-stress5*t2*t7
      A0(6,1) = lamdat_r6-stress6*t2*t7

	  
       N(1:6,1)=A0(1:6,1)
 	   A0=0
 !********************************************************************************************************
 !********************************************************************************************************
      t2 = Sa-Sm
      t3 = stress3*v
      t4 = stress1*v
      t5 = stress2*v
      t6 = v*2.0D0
      t7 = t6+2.0D0
      A0(1,1) = -lamdat_r1+stress1*t2*(1.0D0/2.0D0)-t2*(-stress1+t3+t5)*
     &(1.0D0/2.0D0)-stress2*t2*v*(1.0D0/2.0D0)-stress3*t2*v*(1.0D0/2.0D0
     &)
      A0(1,2) = -lamdat_r2+stress2*t2*(1.0D0/2.0D0)-t2*(-stress2+t3+t4)*
     &(1.0D0/2.0D0)-stress1*t2*v*(1.0D0/2.0D0)-stress3*t2*v*(1.0D0/2.0D0
     &)
      A0(1,3) = -lamdat_r3+stress3*t2*(1.0D0/2.0D0)-t2*(-stress3+t4+t5)*
     &(1.0D0/2.0D0)-stress1*t2*v*(1.0D0/2.0D0)-stress2*t2*v*(1.0D0/2.0D0
     &)
      A0(1,4) = -lamdat_r4+stress4*t2*t7
      A0(1,5) = -lamdat_r5+stress5*t2*t7
      A0(1,6) = -lamdat_r6+stress6*t2*t7


 
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

REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio
integer    NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


! Changing state variables
real*8 A0(4,1)
real*8 e1,e2,e3,e4,e5,e6,et1,et2,et3,et4,et5,et6,ep1,ep2,ep3,ep4,ep5,ep6,bs1,bs2,bs3,bs4,bs5,bs6
real*8 et_initial1,et_initial2,et_initial3,et_initial4,et_initial5,et_initial6


CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)


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
	  
	  
   ! ep1=ep(1)
   ! ep2=ep(2)
   ! ep3=ep(3)
   ! ep4=ep(4)
   ! ep5=ep(5)
   ! if(maxval(abs(ep))<=1.0e-10)then
       ! ep1=1.0e-10_8
    ! end if



!DEC$ NOFREEFORM
 !******************************************************************************
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
      t54 = sqrt(t53)
      t57 = kt*t34*t54*(1.0D0/2.0D0)
      t55 = exp(-t57)
      t56 = t55-1.0D0
      A0(1,1) = -t56*x
      A0(2,1) = -xdo-t56*x
      A0(3,1) = zto+abs(x-xo)
      A0(4,1) = ztdo+abs(xdo+t56*x)

 
       xd = A0(1,1)
       dxd=A0(2,1) 
       zt  = A0(3,1) 
       ztd=A0(4,1) 
!DEC$ FREEFORM	   
       	  
CALL VAR_ASSIGNMENTS(2,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)

!
!
 end

 
 
subroutine SDV_UPDATE(PARAM,NPARAM,VAR,NVAR,STATEV,NSTATV)
implicit none

!DEFINITION OF  INPUT MODEL PARAMETERS

REAL*8   Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio
integer    NLGEOM,ELASTIC,COUPLING,MODEL
INTEGER NPARAM
REAL*8   PARAM(NPARAM)

!DEFINITION OF USER DEFINED DEPENDENT VARIABLES
REAL*8 x,xo
REAL*8 RPLC
REAL*8 e(6),eo(6),et(6),eto(6),Lamdat_r(6),et_tr(6)
INTEGER flag_fwd,flag_rev,Transformation,Bound_Reached,NR_Convergence
REAL*8 VAR(NVAR)
INTEGER NVAR


INTEGER NSTATV
REAL*8 STATEV(NSTATV)


!CALL PARAM_ASSIGNMENTS(PARAM,NPARAM,Sa,Sm,v,alpha,Ca,Cm,Ms,Mf,As,Af,Hmax,kt,n1,n2,n3,n4,a1,a2,a3,rdso,rduo,Dc,Yo,too,x_Up_Bound,x_Low_Bound,Ratio,NLGEOM,ELASTIC,COUPLING,MODEL)
CALL VAR_ASSIGNMENTS(1,VAR,NVAR,x,xo,flag_fwd,flag_rev,Transformation,et,eto,et_tr,lamdat_r,e,eo,RPLC,Bound_Reached,NR_Convergence,too)


	 statev(1)           =      x
     statev(2)           =      real(flag_fwd)
     statev(3)           =      real(flag_rev)
     statev(4)           =      real(transformation)
     statev(5:10)      =      et(1:6)
     statev(11:16)    =      et_tr(1:6)
	 statev(17:22)    =      Lamdat_r(1:6)
     statev(23)         =      RPLC
	 statev(24)         =      BOUND_REACHED
	 statev(25)         =      NR_CONVERGENCE
	 statev(26)         =      TOO 
	 
end

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
 
 
 ! ****** Exponentional Map to calculate ΔR********
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
