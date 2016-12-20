!------------------- 1D Vector Solver ----------------------------------

subroutine MassMomentumHES(nele,Tele,phii,uele,vele,qion,Omee,cond,resn,cfl,lpele,time_ele)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j,k,lpele
   integer :: time_pre,time_post,time_ele
   double precision,dimension(1:NX,1:NY)    :: nele,Tele,qion
   double precision,dimension(1:NX,1:NY),intent(inout) :: phii
   double precision,dimension(1:NX,1:NY),intent(inout) :: uele
   double precision,dimension(1:NX,1:NY),intent(inout) :: vele
   double precision,dimension(1:NX,1:NY,3) :: init
   double precision,dimension(1:NX,1:NY)    :: Omee                     ![-] Electron Hall parameter
   double precision,dimension(1:NX,1:NY)    :: cond                     ![-] Electron conductivity
   double precision,dimension(-1:NX+2,-1:NY+2,3)       :: var           !Ghost cell is used for y-boundary
   double precision,dimension(3)             :: cfl                     !1:x-direction,2:y-direction,3:source
   double precision,dimension(3)             :: resn
   double precision,dimension(1:NX,1:NY,3)   :: rhs
   double precision,dimension(1:NX,1:NY,3)   :: delta
   double precision,dimension(1:NX,1:NY,3)   :: sce
   double precision,dimension(1:NX,1:NY,4)   :: ss
   !double precision,dimension(1:NX,1:NY,9)   :: dinv
   double precision,dimension(1:NX,1:NY,5)   :: aa2,bb2,cc2,dd2,ee2,ff2,gg2,hh2
   double precision,dimension(1:NX,1:NY,3)   :: aa1,bb1,cc1,dd1,ee1,ff1,gg1,hh1
   integer :: bb
   double precision :: maxpot,minpot

   call system_clock(count = time_pre)

   bb = bound
   lpele = 0

      do j = 1,NY
         do i = 1,NX
      !      Tele(i,j) = 100.0d0*dsin(4.0d0*PI/dble(NX)*(dble(i)-0.5d0))*dsin(4.0d0*PI/dble(NY)*(dble(j)-0.5d0))
      !      cond(i,j) = NESTAR
      !      qion(i,j) = 0.0d0
      !      Omee(i,j) = 0.0d0
            init(i,j,1) = phii(i,j)
            init(i,j,2) = uele(i,j)
            init(i,j,3) = vele(i,j)
         enddo
      enddo
   if(CALELE.eq.1) then
      do j = 1,NY
         do i = 1+bb,NX
            var(i,j,1) =-phii(i,j)
            var(i,j,2) = nele(i,j)*uele(i,j)/NESTAR
            var(i,j,3) = nele(i,j)*vele(i,j)/NESTAR
         enddo
      enddo
      call HESMagneticTensor(Omee,cond/NESTAR,ss)
      call SourceTerm(nele/NESTAR,Tele,qion/NESTAR,sce)

      do lpele = 1,50000
         !Dirichlet Boundary condition
         do j = 1,NY
            var(0   +bb,j,1) = 2.0d0*(-PHIA)          -var(1 +bb,j,1)
            var(-1  +bb,j,1) = 2.0d0*var(0+bb,j,1)    -var(1 +bb,j,1)
            var(NX+1,j,1) = 2.0d0*(-PHIC)       -var(NX,j,1)
            var(NX+2,j,1) = 2.0d0*var(NX+1,j,1) -var(NX,j,1)
         enddo
         !Neumann Boundary condition
         do j = 1,NY
            var(0   +bb,j,2) = 2.0d0*var(1   +bb,j,2)-var(2   +bb,j,2)
            var(-1  +bb,j,2) = 2.0d0*var(0   +bb,j,2)-var(1   +bb,j,2)
            var(NX+1,j,2) = 2.0d0*var(NX  ,j,2)-var(NX-1,j,2)
            var(NX+2,j,2) = 2.0d0*var(NX+1,j,2)-var(NX  ,j,2)
         enddo
         !Periodic Boundary condition
         do i = 1+bb,NX
            var(i,0   ,1) = var(i,NY  ,1)
            var(i,0   ,2) = var(i,NY  ,2)
            var(i,0   ,3) = var(i,NY  ,3)
            var(i,-1  ,1) = var(i,NY-1,1)
            var(i,-1  ,2) = var(i,NY-1,2)
            var(i,-1  ,3) = var(i,NY-1,3)
            var(i,NY+1,1) = var(i,1   ,1)
            var(i,NY+1,2) = var(i,1   ,2)
            var(i,NY+1,3) = var(i,1   ,3)
            var(i,NY+2,1) = var(i,2   ,1)
            var(i,NY+2,2) = var(i,2   ,2)
            var(i,NY+2,3) = var(i,2   ,3)
         enddo

         if(AORDER.eq.1) then
            if(lpele.eq.1) then
               call HESCoefficient1O2D_X(aa1,bb1,cc1,dd1)
               call HESCoefficient1O2D_Y(ee1,ff1,gg1,hh1)
            endif

            !$omp parallel default(none),private(i,j),shared(DTELE,rhs,var,sce,aa1,bb1,cc1,dd1,ee1,ff1,gg1,hh1,ss,bb)
            !$omp do
            do j = 1,NY
               do i = 1+bb,NX
                  rhs(i,j,1) =-DTELE*(&
                               aa1(i,j,1)*var(i-1,j  ,1)+aa1(i,j,2)*var(i  ,j  ,1)+aa1(i,j,3)*var(i+1,j  ,1)&
                              +bb1(i,j,1)*var(i-1,j  ,2)+bb1(i,j,2)*var(i  ,j  ,2)+bb1(i,j,3)*var(i+1,j  ,2)&
                              +ee1(i,j,1)*var(i  ,j-1,1)+ee1(i,j,2)*var(i  ,j  ,1)+ee1(i,j,3)*var(i  ,j+1,1)&
                              +ff1(i,j,1)*var(i  ,j-1,3)+ff1(i,j,2)*var(i  ,j  ,3)+ff1(i,j,3)*var(i  ,j+1,3)&
                              +sce(i,j,1))

                  rhs(i,j,2) =-DTELE*(&
                               cc1(i,j,1)*var(i-1,j  ,1)+cc1(i,j,2)*var(i  ,j  ,1)+cc1(i,j,3)*var(i+1,j  ,1)&
                              +dd1(i,j,1)*var(i-1,j  ,2)+dd1(i,j,2)*var(i  ,j  ,2)+dd1(i,j,3)*var(i+1,j  ,2)&
                              +ss(i,j,1)*var(i,j,2)&
                              +ss(i,j,2)*var(i,j,3)&
                              +sce(i,j,2))

                  rhs(i,j,3) =-DTELE*(&
                               gg1(i,j,1)*var(i  ,j-1,1)+gg1(i,j,2)*var(i  ,j  ,1)+gg1(i,j,3)*var(i  ,j+1,1)&
                              +hh1(i,j,1)*var(i  ,j-1,3)+hh1(i,j,2)*var(i  ,j  ,3)+hh1(i,j,3)*var(i  ,j+1,3)&
                              +ss(i,j,3)*var(i,j,2)&
                              +ss(i,j,4)*var(i,j,3)&
                              +sce(i,j,3))
               enddo
            enddo
            !$omp end do
            !$omp end parallel
         elseif(AORDER.eq.2) then
            if(lpele.eq.1) then
               call HESCoefficient2D_X(var,aa2,bb2,cc2,dd2)
               call HESCoefficient2D_Y(var,ee2,ff2,gg2,hh2)
            endif

            !$omp parallel default(none),private(i,j),shared(DTELE,rhs,var,sce,aa2,bb2,cc2,dd2,ee2,ff2,gg2,hh2,ss,bb)
            !$omp do
            do j = 1,NY
               do i = 1+bb,NX
                  rhs(i,j,1) =-DTELE*(&
                               aa2(i,j,1)*var(i-2,j  ,1)+aa2(i,j,2)*var(i-1,j  ,1)+aa2(i,j,3)*var(i  ,j  ,1)+aa2(i,j,4)*var(i+1,j  ,1)+aa2(i,j,5)*var(i+2,j  ,1)&
                              +bb2(i,j,1)*var(i-2,j  ,2)+bb2(i,j,2)*var(i-1,j  ,2)+bb2(i,j,3)*var(i  ,j  ,2)+bb2(i,j,4)*var(i+1,j  ,2)+bb2(i,j,5)*var(i+2,j  ,2)&
                              +ee2(i,j,1)*var(i  ,j-2,1)+ee2(i,j,2)*var(i  ,j-1,1)+ee2(i,j,3)*var(i  ,j  ,1)+ee2(i,j,4)*var(i  ,j+1,1)+ee2(i,j,5)*var(i  ,j+2,1)&
                              +ff2(i,j,1)*var(i  ,j-2,3)+ff2(i,j,2)*var(i  ,j-1,3)+ff2(i,j,3)*var(i  ,j  ,3)+ff2(i,j,4)*var(i  ,j+1,3)+ff2(i,j,5)*var(i  ,j+2,3)&
                              +sce(i,j,1))

                  rhs(i,j,2) =-DTELE*(&
                               cc2(i,j,1)*var(i-2,j  ,1)+cc2(i,j,2)*var(i-1,j  ,1)+cc2(i,j,3)*var(i  ,j  ,1)+cc2(i,j,4)*var(i+1,j  ,1)+cc2(i,j,5)*var(i+2,j  ,1)&
                              +dd2(i,j,1)*var(i-2,j  ,2)+dd2(i,j,2)*var(i-1,j  ,2)+dd2(i,j,3)*var(i  ,j  ,2)+dd2(i,j,4)*var(i+1,j  ,2)+dd2(i,j,5)*var(i+2,j  ,2)&
                              +ss(i,j,1)*var(i,j,2)&
                              +ss(i,j,2)*var(i,j,3)&
                              +sce(i,j,2))

                  rhs(i,j,3) =-DTELE*(&
                               gg2(i,j,1)*var(i  ,j-2,1)+gg2(i,j,2)*var(i  ,j-1,1)+gg2(i,j,3)*var(i  ,j  ,1)+gg2(i,j,4)*var(i  ,j+1,1)+gg2(i,j,5)*var(i  ,j+2,1)&
                              +hh2(i,j,1)*var(i  ,j-2,3)+hh2(i,j,2)*var(i  ,j-1,3)+hh2(i,j,3)*var(i  ,j  ,3)+hh2(i,j,4)*var(i  ,j+1,3)+hh2(i,j,5)*var(i  ,j+2,3)&
                              +ss(i,j,3)*var(i,j,2)&
                              +ss(i,j,4)*var(i,j,3)&
                              +sce(i,j,3))
               enddo
            enddo
            !$omp end do
            !$omp end parallel
         else
            write(*,*) 'Error in AORDER...',AORDER
            stop
         endif

         call LHSCALC(ss,rhs,delta)
         call CFL_SCE(ss,cfl)

         DTELE = SETCFL/dmax1(cfl(1),cfl(2),cfl(3))*DTELE

         !Update elemental variable
         do j = 1,NY
            do i = 1+bb,NX
               do k = 1,3
                  var(i,j,k) = var(i,j,k)+delta(i,j,k)
               enddo
            enddo
         enddo

         !Calculation of normalized difference
         do k = 1,3
            resn(k) = 0.0d0
         enddo
         do j = 1,NY
            do i = 1+bb,NX
               resn(1) = resn(1)+(dabs(delta(i,j,1))/(dabs(var(i,j,1))+1.0d-2))**2.0d0
               resn(2) = resn(2)+(dabs(delta(i,j,2))/(dabs(var(i,j,2))+1.0d16*1.0d3/NESTAR))**2.0d0
               resn(3) = resn(3)+(dabs(delta(i,j,3))/(dabs(var(i,j,3))+1.0d16*1.0d3/NESTAR))**2.0d0
            enddo
         enddo
         do k = 1,3
            resn(k) = dsqrt(resn(k)/dble(NX)/dble(NY))
         enddo
         if(mod(lpele,5000).eq.0) write(*,'(1I5,3E12.3,2I5)') lpele/5000,resn(1),resn(2),resn(3)
         if(resn(1).lt.UPS .and. resn(2).lt.UPS .and. resn(3).lt.UPS) exit
      enddo


      do j = 1,NY
         do i = 1+bb,NX
            phii(i,j) =-var(i,j,1)
            uele(i,j) = var(i,j,2)/nele(i,j)*NESTAR
            vele(i,j) = var(i,j,3)/nele(i,j)*NESTAR
         enddo
      enddo
   endif
   maxpot = 0.0d0
   minpot = 1.0d6
   do j = 1,NY
      do i = 1,NX
         if(phii(i,j).gt.maxpot) maxpot = phii(i,j)
         if(phii(i,j).lt.minpot) minpot = phii(i,j)
      enddo
   enddo
   if(maxpot.gt.1.5d0*PHIA .or. minpot.lt.-0.5d0*PHIA) then
      !do j = 1,NY
      !   do i = 1,NX
      !      phii(i,j) = init(i,j,1)
      !      uele(i,j) = init(i,j,2)
      !      vele(i,j) = init(i,j,3)
      !   enddo
      !enddo
      write(*,*) maxpot,minpot
   endif

   if(it.ge.1) then
      do j = 1,NY
         do i = 1,NX
            phii_ave(i,j) = phii_ave(i,j)+phii(i,j)/dble(ISMP)
            uele_ave(i,j) = uele_ave(i,j)+uele(i,j)/dble(ISMP)
            vele_ave(i,j) = vele_ave(i,j)+vele(i,j)/dble(ISMP)
         enddo
      enddo
      do j = 1,NY
         acurrent = acurrent + nele(1+bb,j)*uele(1+bb,j)*DYL*ZL
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
   double precision,dimension(1:NX,1:NY),intent(in)    :: Omee          ![-] Electron Hall parameter
   double precision,dimension(1:NX,1:NY),intent(in)    :: cond          ![-] Electron conductivity
   double precision,dimension(1:NX,1:NY,4),intent(out) :: ss

   !Calculation of magnetic tensor
   !$omp parallel default(none),private(i,j),shared(cond,Omee,ss)
   !$omp do
   do j = 1,NY
      do i = 1,NX
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
   integer :: i,j,bb
   double precision,dimension(1:NX,1:NY),intent(in)  :: ne,te,qion
   double precision,dimension(1:NX,1:NY,3),intent(out) :: sce
   bb = bound

   !Calculation of pressure terms by 2nd-order central differencing
   !$omp parallel default(none),private(i,j),shared(ne,qion,te,sce,bb)
   !$omp do
   do j = 1,NY
      do i = 1+bb,NX
         sce(i,j,1) =-qion(i,j)
         if     (i.eq.1) then
            sce(i,j,2) = 0.5d0/DXL*(-3.0d0*te(i,j)      +4.0d0*te(i+1,j  )      -te(i+2,j  ))&
                        +0.5d0/DXL*(-3.0d0*dlog(ne(i,j))+4.0d0*dlog(ne(i+1,j  ))-dlog(ne(i+2,j  )))*te(i,j)
         else if(i.ge.2 .and. i.le.NX-1) then
            sce(i,j,2) = 0.5d0/DXL*(te(i+1,j)      -te(i-1,j))&
                        +0.5d0/DXL*(dlog(ne(i+1,j))-dlog(ne(i-1,j)))*te(i,j)
         else if(i.eq.NX) then
            sce(i,j,2) = 0.5d0/DXL*( 3.0d0*te(i,j)      -4.0d0*te(i-1,j  )      +te(i-2,j  ))&
                        +0.5d0/DXL*( 3.0d0*dlog(ne(i,j))-4.0d0*dlog(ne(i-1,j  ))+dlog(ne(i-2,j  )))*te(i,j)
         endif
         if     (j.eq.1) then
            !sce(i,j,3) = 0.5d0/DYL*(-3.0d0*te(i,j)      +4.0d0*te(i  ,j+1)      -te(i  ,j+2))&
            !            +0.5d0/DYL*(-3.0d0*dlog(ne(i,j))+4.0d0*dlog(ne(i  ,j+1))-dlog(ne(i  ,j+2)))*te(i,j)
            sce(i,j,3) = 0.5d0/DYL*(te(i,j+1)      -te(i,NY))&
                        +0.5d0/DYL*(dlog(ne(i,j+1))-dlog(ne(i,NY)))*te(i,j)
         else if(j.ge.2 .and. j.le.NY-1) then
            sce(i,j,3) = 0.5d0/DYL*(te(i,j+1)      -te(i,j-1))&
                        +0.5d0/DYL*(dlog(ne(i,j+1))-dlog(ne(i,j-1)))*te(i,j)
         else if(j.eq.NY) then
            !sce(i,j,3) = 0.5d0/DYL*( 3.0d0*te(i,j)      -4.0d0*te(i  ,j-1)      +te(i  ,j-2))&
            !            +0.5d0/DYL*( 3.0d0*dlog(ne(i,j))-4.0d0*dlog(ne(i  ,j-1))+dlog(ne(i  ,j-2)))*te(i,j)
            sce(i,j,3) = 0.5d0/DYL*(te(i,1)        -te(i,j-1))&
                        +0.5d0/DYL*(dlog(ne(i,1))  -dlog(ne(i,j-1)))*te(i,j)
         endif
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   return
end subroutine


!-----------------------------------------------------------------------

subroutine HESDiagonalInverse(NX,NY,dt,aa,bb,cc,dd,ee,ff,gg,hh,ss,dinv)

   implicit none
   integer :: i,j
   integer,intent(in)                                  :: NX,NY
   double precision,dimension(1:NX,1:NY,5),intent(in)  :: aa,bb,cc,dd,ee,ff,gg,hh
   double precision,dimension(1:NX,1:NY,4),intent(in)  :: ss
   double precision,dimension(1:NX,1:NY,9),intent(out) :: dinv
   double precision,dimension(3,3) :: a
   double precision,intent(in)                         :: dt
   double precision              :: det

   !Calculation of switch
   !$omp parallel default(none),private(i,j,a,det),shared(NX,NY,dt,aa,bb,cc,dd,ee,ff,gg,hh,ss,dinv)
   !$omp do
   do j = 1,NY
      do i = 1,NX
         a(1,1) = 1.0d0+dt*(aa(i,j,3)+ee(i,j,3))
         a(1,2) = dt*bb(i,j,3)
         a(1,3) = dt*ff(i,j,3)
         a(2,1) = dt*cc(i,j,3)
         a(2,2) = 1.0d0+dt*(dd(i,j,3))+ss(i,j,1)
         a(2,3) = 0.0d0!dt*ss(i,j,2)
         a(3,1) = dt*gg(i,j,3)
         a(3,2) = 0.0d0!dt*ss(i,j,3)
         a(3,3) = 1.0d0+dt*(hh(i,j,3))+ss(i,j,4)
         det = a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)+a(3,1)*a(1,2)*a(2,3)&
              -a(1,1)*a(3,2)*a(2,3)-a(3,1)*a(2,2)*a(1,3)-a(2,1)*a(1,2)*a(3,3)
         if(dabs(det).lt.1.0d-3) write(*,'(A18,2I5,1E12.3)') 'Warning...zero det',i,j,det;det = 1.0d0
         dinv(i,j,1) = 1.0d0/det*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
         dinv(i,j,2) = 1.0d0/det*(a(1,3)*a(3,2)-a(1,2)*a(3,3))
         dinv(i,j,3) = 1.0d0/det*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
         dinv(i,j,4) = 1.0d0/det*(a(2,3)*a(3,1)-a(2,1)*a(3,3))
         dinv(i,j,5) = 1.0d0/det*(a(1,1)*a(3,3)-a(1,3)*a(3,1))
         dinv(i,j,6) = 1.0d0/det*(a(1,3)*a(2,1)-a(1,1)*a(2,3))
         dinv(i,j,7) = 1.0d0/det*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
         dinv(i,j,8) = 1.0d0/det*(a(1,2)*a(3,1)-a(1,1)*a(3,2))
         dinv(i,j,9) = 1.0d0/det*(a(1,1)*a(2,2)-a(1,2)*a(2,1))
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
   integer :: i,j,bb
   double precision :: x1
   double precision,dimension(1:NX,1:NY,4),intent(in)  :: ss
   double precision,dimension(3),intent(out)   :: cfl

   bb = bound

   !Calculation of CFL and SCE number
   cfl(1) = 1.0d0*DTELE/DXL
   cfl(2) = 1.0d0*DTELE/DYL
   cfl(3) = 0.0d0
   do j = 1,NY
      do i = 1+bb,NX
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
   double precision,dimension(-1:NX+2,-1:NY+2,3)   :: var
   double precision,dimension(0:NX+1,1:NY) :: ss11,ss12,ss21,ss22
   double precision :: db,df,xx
   double precision,dimension(1:NX,1:NY,5)   :: aa,bb,cc,dd


   do j = 1,NY
      do i = 0,NX+1
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
   do j = 1,NY
      do i = 1,NX
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
   double precision,dimension(-1:NX+2,-1:NY+2,3)   :: var
   double precision,dimension(1:NX,0:NY+1) :: ss11,ss12,ss21,ss22
   double precision :: db,df,xx
   double precision,dimension(1:NX,1:NY,5)   :: ee,ff,gg,hh

   do j = 0,NY+1
      do i = 1,NX
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
   do j = 1,NY
      do i = 1,NX
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
   double precision,dimension(1:NX,1:NY,3)   :: aa1,bb1,cc1,dd1

   !$omp parallel default(none),private(i,j),shared(aa1,bb1,cc1,dd1)
   !$omp do
   do j = 1,NY
      do i = 1,NX
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
   double precision,dimension(1:NX,1:NY,3)   :: ee1,ff1,gg1,hh1


   !$omp parallel default(none),private(i,j),shared(ee1,ff1,gg1,hh1)
   !$omp do
   do j = 1,NY
      do i = 1,NX
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




