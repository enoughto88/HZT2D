!***********************************************************************
!*****    Electron Fluid Solver Using Hyperbolic System Approach   *****
!*****          Developed by Rei Kawashima (Univ. of Tokyo)        *****
!***********************************************************************

subroutine MassMomentumHES(nele,Tele,phii,uele,vele,qion,Omee,cond,resn,cfl,lpele,time_ele)

   use parameters_mod
   use global_mod
   implicit none
   integer :: nx,ny,k,lpele
   integer :: time_pre,time_post,time_ele
   double precision,dimension(1:NXMAX,1:NYMAX)    :: nele,Tele,qion
   double precision,dimension(1:NXMAX,1:NYMAX),intent(inout) :: phii
   double precision,dimension(1:NXMAX,1:NYMAX),intent(inout) :: uele
   double precision,dimension(1:NXMAX,1:NYMAX),intent(inout) :: vele
   double precision,dimension(1:NXMAX,1:NYMAX)    :: Omee               ![-] Electron Hall parameter
   double precision,dimension(1:NXMAX,1:NYMAX)    :: cond               ![-] Electron conductivity
   double precision,dimension(-1:NXMAX+2,-1:NYMAX+2,3)       :: var     !Ghost cell is used for y-boundary
   double precision,dimension(3)             :: cfl                     !1:x-direction,2:y-direction,3:source
   double precision,dimension(3)             :: resn
   double precision,dimension(1:NXMAX,1:NYMAX,3)   :: rhs
   double precision,dimension(1:NXMAX,1:NYMAX,3)   :: delta
   double precision,dimension(1:NXMAX,1:NYMAX,3)   :: sce
   double precision,dimension(1:NXMAX,1:NYMAX,4)   :: ss
   double precision,dimension(1:NXMAX,1:NYMAX,3)   :: aa1,bb1,cc1,dd1,ee1,ff1,gg1,hh1
   double precision :: maxpot,minpot

   call system_clock(count = time_pre)

   lpele = 0
   if(CALELE.eq.1) then
      do ny = 1,NYMAX
         do nx = 1,NXMAX
            var(nx,ny,1) =-phii(nx,ny)
            var(nx,ny,2) = nele(nx,ny)*uele(nx,ny)/NESTAR
            var(nx,ny,3) = nele(nx,ny)*vele(nx,ny)/NESTAR
         enddo
      enddo
      call HESMagneticTensor(Omee,cond/NESTAR,ss)
      call SourceTerm(nele/NESTAR,Tele,qion/NESTAR,sce)

      do lpele = 1,50000
         !Dirichlet Boundary condition
         do ny = 1,NYMAX
            var(0      ,ny,1) = 2.0d0*(-PHIA)           -var(1 ,ny,1)
            var(-1     ,ny,1) = 2.0d0*var(0,ny,1)       -var(1 ,ny,1)
            var(NXMAX+1,ny,1) = 2.0d0*(-PHIC)           -var(NXMAX,ny,1)
            var(NXMAX+2,ny,1) = 2.0d0*var(NXMAX+1,ny,1) -var(NXMAX,ny,1)
         enddo
         !Neumann Boundary condition
         do ny = 1,NYMAX
            var(0      ,ny,2) = 2.0d0*var(1   ,ny,2)-var(2   ,ny,2)
            var(-1     ,ny,2) = 2.0d0*var(0   ,ny,2)-var(1   ,ny,2)
            var(NXMAX+1,ny,2) = 2.0d0*var(NXMAX  ,ny,2)-var(NXMAX-1,ny,2)
            var(NXMAX+2,ny,2) = 2.0d0*var(NXMAX+1,ny,2)-var(NXMAX  ,ny,2)
         enddo
         !Periodic Boundary condition
         do nx = 1,NXMAX
            var(nx,0      ,1) = var(nx,NYMAX  ,1)
            var(nx,0      ,2) = var(nx,NYMAX  ,2)
            var(nx,0      ,3) = var(nx,NYMAX  ,3)
            var(nx,-1     ,1) = var(nx,NYMAX-1,1)
            var(nx,-1     ,2) = var(nx,NYMAX-1,2)
            var(nx,-1     ,3) = var(nx,NYMAX-1,3)
            var(nx,NYMAX+1,1) = var(nx,1      ,1)
            var(nx,NYMAX+1,2) = var(nx,1      ,2)
            var(nx,NYMAX+1,3) = var(nx,1      ,3)
            var(nx,NYMAX+2,1) = var(nx,2      ,1)
            var(nx,NYMAX+2,2) = var(nx,2      ,2)
            var(nx,NYMAX+2,3) = var(nx,2      ,3)
         enddo

         if(lpele.eq.1) then
            call HESCoefficient1O2D_X(aa1,bb1,cc1,dd1)
            call HESCoefficient1O2D_Y(ee1,ff1,gg1,hh1)
         endif

         !$omp parallel default(none),private(nx,ny),shared(DTELE,rhs,var,sce,aa1,bb1,cc1,dd1,ee1,ff1,gg1,hh1,ss)
         !$omp do
         do ny = 1,NYMAX
            do nx = 1,NXMAX
               rhs(nx,ny,1) =-DTELE*(&
                            aa1(nx,ny,1)*var(nx-1,ny  ,1)+aa1(nx,ny,2)*var(nx  ,ny  ,1)+aa1(nx,ny,3)*var(nx+1,ny  ,1)&
                           +bb1(nx,ny,1)*var(nx-1,ny  ,2)+bb1(nx,ny,2)*var(nx  ,ny  ,2)+bb1(nx,ny,3)*var(nx+1,ny  ,2)&
                           +ee1(nx,ny,1)*var(nx  ,ny-1,1)+ee1(nx,ny,2)*var(nx  ,ny  ,1)+ee1(nx,ny,3)*var(nx  ,ny+1,1)&
                           +ff1(nx,ny,1)*var(nx  ,ny-1,3)+ff1(nx,ny,2)*var(nx  ,ny  ,3)+ff1(nx,ny,3)*var(nx  ,ny+1,3)&
                           +sce(nx,ny,1))

               rhs(nx,ny,2) =-DTELE*(&
                            cc1(nx,ny,1)*var(nx-1,ny  ,1)+cc1(nx,ny,2)*var(nx  ,ny  ,1)+cc1(nx,ny,3)*var(nx+1,ny  ,1)&
                           +dd1(nx,ny,1)*var(nx-1,ny  ,2)+dd1(nx,ny,2)*var(nx  ,ny  ,2)+dd1(nx,ny,3)*var(nx+1,ny  ,2)&
                           +ss(nx,ny,1)*var(nx,ny,2)&
                           +ss(nx,ny,2)*var(nx,ny,3)&
                           +sce(nx,ny,2))

               rhs(nx,ny,3) =-DTELE*(&
                            gg1(nx,ny,1)*var(nx  ,ny-1,1)+gg1(nx,ny,2)*var(nx  ,ny  ,1)+gg1(nx,ny,3)*var(nx  ,ny+1,1)&
                           +hh1(nx,ny,1)*var(nx  ,ny-1,3)+hh1(nx,ny,2)*var(nx  ,ny  ,3)+hh1(nx,ny,3)*var(nx  ,ny+1,3)&
                           +ss(nx,ny,3)*var(nx,ny,2)&
                           +ss(nx,ny,4)*var(nx,ny,3)&
                           +sce(nx,ny,3))
            enddo
         enddo
         !$omp end do
         !$omp end parallel

         call LHSCALC(rhs,delta)
         call CFL_SCE(ss,cfl)

         DTELE = SETCFL/dmax1(cfl(1),cfl(2),cfl(3))*DTELE

         !Update elemental variable
         do ny = 1,NYMAX
            do nx = 1,NXMAX
               do k = 1,3
                  var(nx,ny,k) = var(nx,ny,k)+delta(nx,ny,k)
               enddo
            enddo
         enddo

         !Calculation of normalized difference
         do k = 1,3
            resn(k) = 0.0d0
         enddo
         do ny = 1,NYMAX
            do nx = 1,NXMAX
               resn(1) = resn(1)+(dabs(delta(nx,ny,1))/(dabs(var(nx,ny,1))+1.0d-2))**2.0d0
               resn(2) = resn(2)+(dabs(delta(nx,ny,2))/(dabs(var(nx,ny,2))+1.0d-2))**2.0d0
               resn(3) = resn(3)+(dabs(delta(nx,ny,3))/(dabs(var(nx,ny,3))+1.0d-2))**2.0d0
            enddo
         enddo
         do k = 1,3
            resn(k) = dsqrt(resn(k)/dble(NXMAX)/dble(NYMAX))
         enddo
         if(mod(lpele,5000).eq.0) write(*,'(1I5,3E12.3,2I5)') lpele/5000,resn(1),resn(2),resn(3)
         if(resn(1).lt.UPS .and. resn(2).lt.UPS .and. resn(3).lt.UPS) exit
      enddo


      do ny = 1,NYMAX
         do nx = 1,NXMAX
            phii(nx,ny) =-var(nx,ny,1)
            uele(nx,ny) = var(nx,ny,2)/nele(nx,ny)*NESTAR
            vele(nx,ny) = var(nx,ny,3)/nele(nx,ny)*NESTAR
         enddo
      enddo
   endif
   maxpot = 0.0d0
   minpot = 1.0d6
   do ny = 1,NYMAX
      do nx = 1,NXMAX
         if(phii(nx,ny).gt.maxpot) maxpot = phii(nx,ny)
         if(phii(nx,ny).lt.minpot) minpot = phii(nx,ny)
      enddo
   enddo
   if(maxpot.gt.1.5d0*PHIA .or. minpot.lt.-0.5d0*PHIA) then
      write(*,*) maxpot,minpot
   endif

   if(it.ge.1) then
      do ny = 1,NYMAX
         do nx = 1,NXMAX
            phii_ave(nx,ny) = phii_ave(nx,ny)+phii(nx,ny)/dble(ISMP)
            uele_ave(nx,ny) = uele_ave(nx,ny)+uele(nx,ny)/dble(ISMP)
            vele_ave(nx,ny) = vele_ave(nx,ny)+vele(nx,ny)/dble(ISMP)
         enddo
      enddo
      do ny = 1,NYMAX
         acurrent = acurrent + nele(1,ny)*uele(1,ny)*DYL*ZL
      enddo
      do k = 1,3
         resn_ave(k) = resn_ave(k)+resn(k)
      enddo
   endif


   call system_clock(count = time_post)
   time_ele = time_post - time_pre

   return
endsubroutine


!-----------------------------------------------------------------------

subroutine HESMagneticTensor(Omee,cond,ss)

   use parameters_mod
   implicit none
   integer :: i,j
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)    :: Omee          ![-] Electron Hall parameter
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)    :: cond          ![-] Electron conductivity
   double precision,dimension(1:NXMAX,1:NYMAX,4),intent(out) :: ss

   !Calculation of magnetic tensor
   !$omp parallel default(none),private(i,j),shared(cond,Omee,ss)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
         ss(i,j,1) = 1.0d0/cond(i,j)*1.0d0
         ss(i,j,2) =-1.0d0/cond(i,j)*Omee(i,j)
         ss(i,j,3) = 1.0d0/cond(i,j)*Omee(i,j)
         ss(i,j,4) = 1.0d0/cond(i,j)*1.0d0
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   return
end subroutine




!-----------------------------------------------------------------------

subroutine SourceTerm(ne,te,qion,sce)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)  :: ne,te,qion
   double precision,dimension(1:NXMAX,1:NYMAX,3),intent(out) :: sce

   !Calculation of pressure terms by 2nd-order central differencing
   !$omp parallel default(none),private(i,j),shared(ne,qion,te,sce)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
         sce(i,j,1) =-qion(i,j)
         if     (i.eq.1) then
            sce(i,j,2) = 0.5d0/DXL*(-3.0d0*te(i,j)      +4.0d0*te(i+1,j  )      -te(i+2,j  ))&
                        +0.5d0/DXL*(-3.0d0*dlog(ne(i,j))+4.0d0*dlog(ne(i+1,j  ))-dlog(ne(i+2,j  )))*te(i,j)
         else if(i.ge.2 .and. i.le.NXMAX-1) then
            sce(i,j,2) = 0.5d0/DXL*(te(i+1,j)      -te(i-1,j))&
                        +0.5d0/DXL*(dlog(ne(i+1,j))-dlog(ne(i-1,j)))*te(i,j)
         else if(i.eq.NXMAX) then
            sce(i,j,2) = 0.5d0/DXL*( 3.0d0*te(i,j)      -4.0d0*te(i-1,j  )      +te(i-2,j  ))&
                        +0.5d0/DXL*( 3.0d0*dlog(ne(i,j))-4.0d0*dlog(ne(i-1,j  ))+dlog(ne(i-2,j  )))*te(i,j)
         endif
         if     (j.eq.1) then
            sce(i,j,3) = 0.5d0/DYL*(te(i,j+1)      -te(i,NYMAX))&
                        +0.5d0/DYL*(dlog(ne(i,j+1))-dlog(ne(i,NYMAX)))*te(i,j)
         else if(j.ge.2 .and. j.le.NYMAX-1) then
            sce(i,j,3) = 0.5d0/DYL*(te(i,j+1)      -te(i,j-1))&
                        +0.5d0/DYL*(dlog(ne(i,j+1))-dlog(ne(i,j-1)))*te(i,j)
         else if(j.eq.NYMAX) then
            sce(i,j,3) = 0.5d0/DYL*(te(i,1)        -te(i,j-1))&
                        +0.5d0/DYL*(dlog(ne(i,1))  -dlog(ne(i,j-1)))*te(i,j)
         endif
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   return
end subroutine



!***********************************************************************

subroutine CFL_SCE(ss,cfl)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j
   double precision :: x1
   double precision,dimension(1:NXMAX,1:NYMAX,4),intent(in)  :: ss
   double precision,dimension(3),intent(out)   :: cfl

   !Calculation of CFL and SCE number
   cfl(1) = 1.0d0*DTELE/DXL
   cfl(2) = 1.0d0*DTELE/DYL
   cfl(3) = 0.0d0
   do j = 1,NYMAX
      do i = 1,NXMAX
         x1 = dmax1(dabs(ss(i,j,1)),dabs(ss(i,j,2)),dabs(ss(i,j,3)),dabs(ss(i,j,4)))*DTELE
         if(x1.gt.cfl(3)) then
            cfl(3) = x1
         endif
      enddo
   enddo

   return
endsubroutine

!---------- Linear Expression of HES Coefficients ----------------------

subroutine HESCoefficient2D_X(var,aa,bb,cc,dd)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j
   double precision,dimension(-1:NXMAX+2,-1:NYMAX+2,3)   :: var
   double precision,dimension(0:NXMAX+1,1:NYMAX) :: ss11,ss12,ss21,ss22
   double precision :: db,df,xx
   double precision,dimension(1:NXMAX,1:NYMAX,5)   :: aa,bb,cc,dd


   do j = 1,NYMAX
      do i = 0,NXMAX+1
         db = var(i,j,1)-var(i-1,j,1); df = var(i+1,j,1)-var(i,j,1)
         if(db*df.lt.0.0d0) then
            ss11(i,j) = 0.0d0; ss12(i,j) = 0.0d0
         else if(dabs(db).le.dabs(df)) then
            ss11(i,j) = 1.0d0; ss12(i,j) = 0.0d0
         else if(dabs(db).gt.dabs(df)) then
            ss11(i,j) = 0.0d0; ss12(i,j) = 1.0d0
         endif
         db = var(i,j,2)-var(i-1,j,2); df = var(i+1,j,2)-var(i,j,2)
         if(db*df.lt.0.0d0) then
            ss21(i,j) = 0.0d0; ss22(i,j) = 0.0d0
         else if(dabs(db).le.dabs(df)) then
            ss21(i,j) = 1.0d0; ss22(i,j) = 0.0d0
         else if(dabs(db).gt.dabs(df)) then
            ss21(i,j) = 0.0d0; ss22(i,j) = 1.0d0
         endif
      enddo
   enddo
   xx = 1.0d0/DXL/4.0d0
   !$omp parallel default(none),private(i,j),shared(xx,ss11,ss12,ss21,ss22,aa,bb,cc,dd)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
            aa(i,j,1) = xx*(ss11(i-1,j))
            aa(i,j,2) = xx*(-2.0d0-ss11(i-1,j)+ss12(i-1,j))
            aa(i,j,3) = xx*( 4.0d0-ss11(i+1,j)-ss12(i-1,j))
            aa(i,j,4) = xx*(-2.0d0+ss11(i+1,j)-ss12(i+1,j))
            aa(i,j,5) = xx*(ss12(i+1,j))
            bb(i,j,1) = xx*(ss21(i-1,j))
            bb(i,j,2) = xx*(-2.0d0-ss21(i-1,j)+ss22(i-1,j)-2.0d0*ss21(i  ,j))
            bb(i,j,3) = xx*( 2.0d0*ss21(i  ,j)-ss22(i-1,j)+ss21(i+1,j)-2.0d0*ss22(i  ,j))
            bb(i,j,4) = xx*( 2.0d0-ss21(i+1,j)+ss22(i+1,j)+2.0d0*ss22(i  ,j))
            bb(i,j,5) = xx*(-ss22(i+1,j))
            cc(i,j,1) = xx*(ss11(i-1,j))
            cc(i,j,2) = xx*(-2.0d0-ss11(i-1,j)+ss12(i-1,j)-2.0d0*ss11(i  ,j))
            cc(i,j,3) = xx*( 2.0d0*ss11(i  ,j)-ss12(i-1,j)+ss11(i+1,j)-2.0d0*ss12(i  ,j))
            cc(i,j,4) = xx*( 2.0d0-ss11(i+1,j)+ss12(i+1,j)+2.0d0*ss12(i  ,j))
            cc(i,j,5) = xx*(-ss12(i+1,j))
            dd(i,j,1) = xx*(ss21(i-1,j))
            dd(i,j,2) = xx*(-2.0d0-ss21(i-1,j)+ss22(i-1,j))
            dd(i,j,3) = xx*( 4.0d0-ss21(i+1,j)-ss22(i-1,j))
            dd(i,j,4) = xx*(-2.0d0+ss21(i+1,j)-ss22(i+1,j))
            dd(i,j,5) = xx*(ss22(i+1,j))
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   return
endsubroutine



!---------- Linear Expression of HES Coefficients ----------------------

subroutine HESCoefficient2D_Y(var,ee,ff,gg,hh)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j
   double precision,dimension(-1:NXMAX+2,-1:NYMAX+2,3)   :: var
   double precision,dimension(1:NXMAX,0:NYMAX+1) :: ss11,ss12,ss21,ss22
   double precision :: db,df,xx
   double precision,dimension(1:NXMAX,1:NYMAX,5)   :: ee,ff,gg,hh

   do j = 0,NYMAX+1
      do i = 1,NXMAX
         db = var(i,j,1)-var(i,j-1,1); df = var(i,j+1,1)-var(i,j,1)
         if(db*df.lt.0.0d0) then
            ss11(i,j) = 0.0d0; ss12(i,j) = 0.0d0
         else if(dabs(db).le.dabs(df)) then
            ss11(i,j) = 1.0d0; ss12(i,j) = 0.0d0
         else if(dabs(db).gt.dabs(df)) then
            ss11(i,j) = 0.0d0; ss12(i,j) = 1.0d0
         endif
         db = var(i,j,3)-var(i,j-1,3); df = var(i,j+1,3)-var(i,j,3)
         if(db*df.lt.0.0d0) then
            ss21(i,j) = 0.0d0; ss22(i,j) = 0.0d0
         else if(dabs(db).le.dabs(df)) then
            ss21(i,j) = 1.0d0; ss22(i,j) = 0.0d0
         else if(dabs(db).gt.dabs(df)) then
            ss21(i,j) = 0.0d0; ss22(i,j) = 1.0d0
         endif
      enddo
   enddo
   xx = 1.0d0/DYL/4.0d0
   !$omp parallel default(none),private(i,j),shared(xx,ss11,ss12,ss21,ss22,ee,ff,gg,hh)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
            ee(i,j,1) = xx*(ss11(i,j-1))
            ee(i,j,2) = xx*(-2.0d0-ss11(i,j-1)+ss12(i,j-1))
            ee(i,j,3) = xx*( 4.0d0-ss11(i,j+1)-ss12(i,j-1))
            ee(i,j,4) = xx*(-2.0d0+ss11(i,j+1)-ss12(i,j+1))
            ee(i,j,5) = xx*(ss12(i,j+1))
            ff(i,j,1) = xx*(ss21(i,j-1))
            ff(i,j,2) = xx*(-2.0d0-ss21(i,j-1)+ss22(i,j-1)-2.0d0*ss21(i,j  ))
            ff(i,j,3) = xx*( 2.0d0*ss21(i,j  )-ss22(i,j-1)+ss21(i,j+1)-2.0d0*ss22(i,j  ))
            ff(i,j,4) = xx*( 2.0d0-ss21(i,j+1)+ss22(i,j+1)+2.0d0*ss22(i,j  ))
            ff(i,j,5) = xx*(-ss22(i,j+1))
            gg(i,j,1) = xx*(ss11(i,j-1))
            gg(i,j,2) = xx*(-2.0d0-ss11(i,j-1)+ss12(i,j-1)-2.0d0*ss11(i,j  ))
            gg(i,j,3) = xx*( 2.0d0*ss11(i,j  )-ss12(i,j-1)+ss11(i,j+1)-2.0d0*ss12(i,j  ))
            gg(i,j,4) = xx*( 2.0d0-ss11(i,j+1)+ss12(i,j+1)+2.0d0*ss12(i,j  ))
            gg(i,j,5) = xx*(-ss12(i,j+1))
            hh(i,j,1) = xx*(ss21(i,j-1))
            hh(i,j,2) = xx*(-2.0d0-ss21(i,j-1)+ss22(i,j-1))
            hh(i,j,3) = xx*( 4.0d0-ss21(i,j+1)-ss22(i,j-1))
            hh(i,j,4) = xx*(-2.0d0+ss21(i,j+1)-ss22(i,j+1))
            hh(i,j,5) = xx*(ss22(i,j+1))
      enddo
   enddo
   !$omp end do
   !$omp end parallel


   return
endsubroutine



!---------- Linear Expression of HES Coefficients ----------------------

subroutine HESCoefficient1O2D_X(aa1,bb1,cc1,dd1)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j
   double precision,dimension(1:NXMAX,1:NYMAX,3)   :: aa1,bb1,cc1,dd1

   !$omp parallel default(none),private(i,j),shared(aa1,bb1,cc1,dd1)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
            aa1(i,j,1) =-1.0d0/DXL/2.0d0
            aa1(i,j,2) = 1.0d0/DXL
            aa1(i,j,3) =-1.0d0/DXL/2.0d0
            bb1(i,j,1) =-1.0d0/DXL/2.0d0
            bb1(i,j,2) = 0.0d0
            bb1(i,j,3) = 1.0d0/DXL/2.0d0
            cc1(i,j,1) =-1.0d0/DXL/2.0d0
            cc1(i,j,2) = 0.0d0
            cc1(i,j,3) = 1.0d0/DXL/2.0d0
            dd1(i,j,1) =-1.0d0/DXL/2.0d0
            dd1(i,j,2) = 1.0d0/DXL
            dd1(i,j,3) =-1.0d0/DXL/2.0d0
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   return
endsubroutine



!---------- Linear Expression of HES Coefficients ----------------------

subroutine HESCoefficient1O2D_Y(ee1,ff1,gg1,hh1)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j
   double precision,dimension(1:NXMAX,1:NYMAX,3)   :: ee1,ff1,gg1,hh1


   !$omp parallel default(none),private(i,j),shared(ee1,ff1,gg1,hh1)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
            ee1(i,j,1) =-1.0d0/DYL/2.0d0
            ee1(i,j,2) = 1.0d0/DYL
            ee1(i,j,3) =-1.0d0/DYL/2.0d0
            ff1(i,j,1) =-1.0d0/DYL/2.0d0
            ff1(i,j,2) = 0.0d0
            ff1(i,j,3) = 1.0d0/DYL/2.0d0
            gg1(i,j,1) =-1.0d0/DYL/2.0d0
            gg1(i,j,2) = 0.0d0
            gg1(i,j,3) = 1.0d0/DYL/2.0d0
            hh1(i,j,1) =-1.0d0/DYL/2.0d0
            hh1(i,j,2) = 1.0d0/DYL
            hh1(i,j,3) =-1.0d0/DYL/2.0d0
      enddo
   enddo
   !$omp end do
   !$omp end parallel


   return
endsubroutine




