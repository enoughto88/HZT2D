module global_mod

   use parameters_mod
   implicit none
   !Time Step
   integer              :: it
   integer              :: nstp  = 0
   double precision     :: DTELE = 1.0d-12
   !Average quantities
   double precision     :: inflow
   double precision     :: outflow
   double precision     :: outneu,oution
   double precision     :: acurrent
   double precision     :: momx
   double precision,dimension(1:NXMAX,1:NYMAX)    :: nele_ave           ![m-3] Electron number density
   double precision,dimension(1:NXMAX,1:NYMAX)    :: nneu_ave           ![m-3] Neutral number density
   double precision,dimension(1:NXMAX,1:NYMAX)    :: qion_ave           ![m-3s-1] Ion production rate
   double precision,dimension(1:NXMAX,1:NYMAX)    :: phii_ave           ![V] Space potential
   double precision,dimension(1:NXMAX,1:NYMAX)    :: uele_ave           ![ms-1] x-Velocity
   double precision,dimension(1:NXMAX,1:NYMAX)    :: vele_ave           ![ms-1] y-Velocity
   double precision,dimension(3)                  :: resn_ave           ![-] Residual
   double precision,dimension(1:NXMAX,1:NYMAX)    :: cons = 0.0d0       !Neutral particle consumption
   !Definition of structure

   type particle
      double precision :: xx,yy,vx,vy,vz                                !2D3V position and speed
      double precision :: np                                            !Number of particles in the macroparticle
      double precision :: ip,jp                                         !Cell index of particle position
   end type

end module global_mod


