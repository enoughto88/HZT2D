!***********************************************************************
!*****              Parameter Setting Module for HZT2D             *****
!*****          Developed by Rei Kawashima (Univ. of Tokyo)        *****
!***********************************************************************

module parameters_mod

   implicit none
   !Program parameter
   character(len=80)           :: TOPDIR    = '/home/kawashima/HZT2D'
   integer         ,parameter  :: ITM       = 100000                    ![-] Maximum number of time steps
   integer         ,parameter  :: ISMP      = 10                        ![-] Sampling interval of distribution
   integer         ,parameter  :: PSMP      = ITM                       ![-] Sampling interval of particles
   double precision,parameter  :: UPS       = 1.0d-6                    ![-] Truncation error
   !Physical parameter
   double precision,parameter  :: PI        = 3.1415926535893d0
   double precision,parameter  :: ECH       = 1.60217733d-19            ![C] Elementary charge
   double precision,parameter  :: KBOLTZ    = 1.3806504d-23             ![m2kgs-1K] Boltzmann constant
   double precision,parameter  :: ME        = 9.1093826d-31             ![kg] Electron mass
   double precision,parameter  :: PHIA      = 300.0d0                   ![V] Anode voltage
   double precision,parameter  :: PHIC      = 0.0d0                     ![V] Cathode voltage
   double precision,parameter  :: TEC       = 12.1d0                    ![eV] Cathode electron temperature
   double precision,parameter  :: SETCFL    = 0.05d0                    ![-] Courant number
   !Grid parameter
   integer         ,parameter  :: NX        = 24                        ![-] x number of cells
   integer         ,parameter  :: NY        = 24                        ![-] y number of cells
   double precision,parameter  :: XL        = 0.050d0                   ![m] Calc. field length in x-direction
   double precision,parameter  :: YL        = 2.0d0*PI*0.040d0          ![m] Calc. field length in y-direction
   double precision,parameter  :: ZL        = 0.015d0                   ![m] Channel length in z-direction
   double precision,parameter  :: DXL       = XL/dble(NX)               ![m] Cell length in x-direction
   double precision,parameter  :: DYL       = YL/dble(NY)               ![m] Cell length in y-direction
   !Particle parameter
   integer         ,parameter  :: NMAX      = 8000000                   ![-] Maximum number of particles
   integer         ,parameter  :: NMINI     = 200000                    ![-] Initial number of particles
   double precision,parameter  :: DTPIC     = 1.0d-8                    ![s] Time step for PIC calculation
   double precision,parameter  :: MACP      = 5.0d10                    ![-] Normal macroparticle size
   double precision,parameter  :: MACI      = 5.0d8                     ![-] Normal macroparticle size
   double precision,parameter  :: NEMIN     = 1.0d16                    ![m-3] Minimum electron number density
   double precision,parameter  :: NUMIN     = 1.0d4                     ![m-1] Minimum total collision frequency
   double precision,parameter  :: MDOT      = 5.0d-6                    ![kg/s] Mass flow rate
   double precision,parameter  :: MI        = 2.196324331d-25           ![kg] Ion mass of xenon
   double precision,parameter  :: TI        = 300.0d0                   ![K] Ion temperature
   double precision,parameter  :: TWALL     = 850.0d0                   ![K] Wall temperature
   integer         ,parameter  :: TRAC      = 0
   integer         ,parameter  :: INPUTP    = 1
   integer         ,parameter  :: INPUTF    = 1
   integer         ,parameter  :: INPUTM    = 1
   integer         ,parameter  :: CALPIC    = 1
   integer         ,parameter  :: CALELE    = 1
   double precision,parameter  :: NESTAR    = 1.0d21
   integer         ,parameter  :: AORDER    = 1

end module parameters_mod


