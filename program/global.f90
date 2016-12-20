module global_mod

   use parameters_mod
   implicit none
   !Time Step
   integer              :: it
   integer              :: nstp  = 0
   integer              :: bound = 0
   double precision     :: DTELE = 1.0d-12
   !Average quantities
   double precision     :: inflow
   double precision     :: outflow
   double precision     :: outneu,oution
   double precision     :: acurrent
   double precision     :: momx
   double precision,dimension(1:NX,1:NY)    :: nele_ave                 ![m-3] Electron number density
   double precision,dimension(1:NX,1:NY)    :: nneu_ave                 ![m-3] Neutral number density
   double precision,dimension(1:NX,1:NY)    :: qion_ave                 ![m-3s-1] Ion production rate
   double precision,dimension(1:NX,1:NY)    :: phii_ave                 ![V] Space potential
   double precision,dimension(1:NX,1:NY)    :: uele_ave                 ![ms-1] x-Velocity
   double precision,dimension(1:NX,1:NY)    :: vele_ave                 ![ms-1] y-Velocity
   double precision,dimension(3)            :: resn_ave                 ![-] Residual
   double precision,dimension(1:NX,1:NY)    :: cons = 0.0d0             !Neutral particle consumption
end module global_mod


