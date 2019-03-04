# SMA-UMAT
A user-defined material subroutine for polycrystalline shape memory alloys under large deformations
! Open with Notepad++ for format compatibility
!*********************************************************************************************XXX
!                                                                                             XXX
!                              3D FINITE STRIN SMA UMAT                                       XXX
!                                                                                             XXX
!    			COPYRIGHT   LEI XU1, Theo Baxevanis2, Dimitris Lagoudas1                            XXX
!           1 Department of Aerospace Engineering, Texas A&M University                       XXX
!           2 Department of Mechanical Engineering, University of Houston                     XXX
!               Version  (3-D)  20/JAN./2018                                                  XXX
!                                                                                             XXX                             
! XXX************************************************************************************************ XXX  
! XXX  *MATERIAL, NAME=SMA                                                                            XXX                                                                   
! XXX  *DEPVAR                                                                                        XXX   
! XXX    26                                                                                           XXX   
! XXX     statev(1)           =      x                                   !Martensite Vol. Frac.       XXX   
! XXX     statev(2)           =      real(flag_fwd)              !FWD finish at 1 or not 0            XXX   
! XXX     statev(3)           =      real(flag_rev)               !REV  finish at 1 or not 0          XXX   
! XXX     statev(4)           =      real(transformation)     !Trans. DIR  0-elastic; -1 REV; 1 FWD   XXX   
! XXX     statev(5:10)      =      et(1:6)                          !Trans. Strain                    XXX   
! XXX     statev(11:16)    =      et_tr(1:6)                     !Trans. Strain at reverse            XXX   
! XXX     statev(17:22)    =      Lamdat_r(1:6)             !Reverse DIR tensor                       XXX   
! XXX     statev(23)         =      RPLC                                                              XXX   
! XXX     statev(24)         =      BOUND_REACHED                                                     XXX   
! XXX     statev(25)         =      NR_CONVERGENCE                                                    XXX   
! XXX     statev(26)         =      TOO                                                               XXX   
! XXX*******************************************************************************************************************************XXX
! XXX  *User Material, constants=25                                                                                                                                  XXX
! XXX  <Ea>,              <Em>,          <v>,           <alpha>,        <CA>,            <CM>,                <Ms>,         <Mf>                    XXX 
! XXX  <As>,              <Af>,          <Hmax>,         <kt>,          <SigCal>,       <n1>,                 <n2>,         <n3>                    XXX                                             
! XXX  <n4>,             <X_initial>,      <Etd>,  <Tube_Flag>,      <Tube_Axis>,  <Tube_Radius>, <Elastic>,    <Coupling>  XXX                                                         
! XXX   <Model>                                                                                                                                                                  XXX
! XXX********************************************************************************************************************************XXX
