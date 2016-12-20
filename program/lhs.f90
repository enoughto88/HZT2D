!-------------------Calculation of Right Hand Side----------------------

subroutine LHSCALC(mag,rhs,delta)

   use parameters_mod
   use global_mod
   implicit none
   double precision,dimension(1:NX,1:NY,4)              :: mag
   double precision,dimension(1:NX,1:NY,3)              :: rhs
   double precision,dimension(1:NX,1:NY,3)              :: delta
   !double precision,dimension(1:NX,1:NY,3)              :: delta1

   if(AORDER.eq.0) mag=0.0d0
   !Implicit treatment of source term
   !call EXPLICIT(rhs,delta)
   !call SOURCEIMP(mag,rhs,delta1)
   call LUSGS(rhs,delta)

   return
endsubroutine


!-------------------Explicit method-------------------------------------

subroutine EXPLICIT(rhs,delta)

   use parameters_mod
   implicit none
   integer :: i,j,k
   double precision,dimension(1:NX,1:NY,3),intent(in)   :: rhs
   double precision,dimension(1:NX,1:NY,3),intent(out)  :: delta

   do j = 1,NY
      do i = 1,NX
         do k = 1,3
            delta(i,j,k) = rhs(i,j,k)
         enddo
      enddo
   enddo

   return
endsubroutine


!-------------------Implicit treatment of source term-------------------

subroutine SOURCEIMP(mag,rhs,delta)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j
   double precision,dimension(1:NX,1:NY,4),intent(in)   :: mag
   double precision,dimension(1:NX,1:NY,3),intent(in)   :: rhs
   double precision,dimension(1:NX,1:NY,3),intent(out)  :: delta
   double precision :: det,m1,m2,m3,m4
   integer :: bb

   bb = bound
   do j = 1,NY
      do i = 1+bb,NX
         m1 = mag(i,j,1);m2 = mag(i,j,2)
         m3 = mag(i,j,3);m4 = mag(i,j,4)
         det = DTELE**2.0d0*(m1*m4-m2*m3)+DTELE*(m1+m4)+1.0d0
         if(det.lt.1.0d-3) det = 1.0d-3
         delta(i,j,1) = rhs(i,j,1)
         delta(i,j,2) = 1.0d0/det*(1.0d0+DTELE*m4)*rhs(i,j,2)&
                       -1.0d0/det*DTELE*m2        *rhs(i,j,3)
         delta(i,j,3) =-1.0d0/det*DTELE*m3        *rhs(i,j,2)&
                       +1.0d0/det*(1.0d0+DTELE*m1)*rhs(i,j,3)
      enddo
   enddo

   return
endsubroutine


!-------------------Implicit treatment of x-flux-------------------

subroutine XFLUXIMP(rhs,delta)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j,k
   double precision,dimension(1:NX,1:NY,3),intent(in)   :: rhs
   double precision,dimension(1:NX,1:NY,3),intent(out)  :: delta
   double precision,dimension(1:NX,1:NY,3)  :: deltah
   double precision :: hh
   double precision,dimension(3) :: qq1,qq2

   !LDU-factorization + SGS
   hh = DTELE/DXL
   do j = 1,NY
      !Sweep 1
      do i = 1,NX
         !Calculation of LdQ
         if(i.eq.1) then
            qq1(1) = 0.0d0
            qq1(2) = 0.0d0
            qq1(3) = 0.0d0
         else if(i.ge.2 .and. i.le.NX) then
            qq1(1) =-hh*0.5d0*(rhs(i-1,j,1)+rhs(i-1,j,2))
            qq1(2) =-hh*0.5d0*(rhs(i-1,j,1)+rhs(i-1,j,2))
            qq1(3) = 0.0d0
         endif
         !Calculation of RHS-LdQ
         qq2(1) = rhs(i,j,1)-qq1(1)
         qq2(2) = rhs(i,j,2)-qq1(2)
         qq2(3) = rhs(i,j,3)-qq1(3)
         !Calculation of D^-1
         if(i.eq.1) then
            deltah(i,j,1) = 1.0d0/(1.0d0+hh)*((1.0d0+0.5d0*hh)*qq2(1)+( 0.5d0*hh)     *qq2(2))
            deltah(i,j,2) = 1.0d0/(1.0d0+hh)*(( 0.5d0*hh)     *qq2(1)+(1.0d0+0.5d0*hh)*qq2(2))
            deltah(i,j,3) = qq2(3)
         else if(i.ge.2 .and. i.le.NX-1) then
            deltah(i,j,1) = 1.0d0/(1.0d0+hh)*qq2(1)
            deltah(i,j,2) = 1.0d0/(1.0d0+hh)*qq2(2)
            deltah(i,j,3) = qq2(3)
         else if(i.eq.NX) then
            deltah(i,j,1) = 1.0d0/(1.0d0+hh)*((1.0d0+0.5d0*hh)*qq2(1)+(-0.5d0*hh)     *qq2(2))
            deltah(i,j,2) = 1.0d0/(1.0d0+hh)*((-0.5d0*hh)     *qq2(1)+(1.0d0+0.5d0*hh)*qq2(2))
            deltah(i,j,3) = qq2(3)
         endif
      enddo
      !Sweep 2
      do k = 1,NX
         i = NX-k+1
         !Calculation of UdQ
         if(i.eq.NX) then
            qq1(1) = 0.0d0
            qq1(2) = 0.0d0
            qq1(3) = 0.0d0
         else if(i.ge.1 .and. i.le.NX-1) then
            qq1(1) = hh*0.5d0*(-deltah(i+1,j,1)+deltah(i+1,j,2))
            qq1(2) = hh*0.5d0*( deltah(i+1,j,1)-deltah(i+1,j,2))
            qq1(3) = 0.0d0
         endif
         !Calculation of D^-1
         if(i.eq.1) then
            qq2(1) = 1.0d0/(1.0d0+hh)*((1.0d0+0.5d0*hh)*qq1(1)+( 0.5d0*hh)     *qq1(2))
            qq2(2) = 1.0d0/(1.0d0+hh)*(( 0.5d0*hh)     *qq1(1)+(1.0d0+0.5d0*hh)*qq1(2))
            qq2(3) = qq1(3)
         else if(i.ge.2 .and. i.le.NX-1) then
            qq2(1) = 1.0d0/(1.0d0+hh)*qq1(1)
            qq2(2) = 1.0d0/(1.0d0+hh)*qq1(2)
            qq2(3) = qq1(3)
         else if(i.eq.NX) then
            qq2(1) = 1.0d0/(1.0d0+hh)*((1.0d0+0.5d0*hh)*qq1(1)+(-0.5d0*hh)     *qq1(2))
            qq2(2) = 1.0d0/(1.0d0+hh)*((-0.5d0*hh)     *qq1(1)+(1.0d0+0.5d0*hh)*qq1(2))
            qq2(3) = qq1(3)
         endif
         !Calculation of delta
         delta(i,j,1) = deltah(i,j,1)-qq2(1)
         delta(i,j,2) = deltah(i,j,2)-qq2(2)
         delta(i,j,3) = deltah(i,j,3)-qq2(3)
      enddo
   enddo

   return
endsubroutine


!-------------------Implicit treatment of y-flux-------------------

subroutine YFLUXIMP(rhs,delta)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j,k
   double precision,dimension(1:NX,1:NY,3),intent(in)   :: rhs
   double precision,dimension(1:NX,1:NY,3),intent(out)  :: delta
   double precision,dimension(1:NX,1:NY,3)  :: deltah
   double precision :: hh
   double precision,dimension(3) :: qq1,qq2

   !LDU-factorization + SGS
   hh = DTELE/DYL
   do i = 1,NX
      !Sweep 1
      do j = 1,NY
         !Calculation of LdQ
         if(j.eq.1) then
            qq1(1) = 0.0d0
            qq1(2) = 0.0d0
            qq1(3) = 0.0d0
         else if(j.ge.2 .and. j.le.NY) then
            qq1(1) =-hh*0.5d0*(rhs(i,j-1,1)+rhs(i,j-1,3))
            qq1(2) = 0.0d0
            qq1(3) =-hh*0.5d0*(rhs(i,j-1,1)+rhs(i,j-1,3))
         endif
         !Calculation of RHS-LdQ
         qq2(1) = rhs(i,j,1)-qq1(1)
         qq2(2) = rhs(i,j,2)-qq1(2)
         qq2(3) = rhs(i,j,3)-qq1(3)
         !Calculation of D^-1
         if(j.eq.1) then
            deltah(i,j,1) = 1.0d0/(1.0d0+hh)*((1.0d0+0.5d0*hh)*qq2(1)+( 0.5d0*hh)     *qq2(3))
            deltah(i,j,2) = qq2(2)
            deltah(i,j,3) = 1.0d0/(1.0d0+hh)*(( 0.5d0*hh)     *qq2(1)+(1.0d0+0.5d0*hh)*qq2(3))
         else if(j.ge.2 .and. j.le.NY-1) then
            deltah(i,j,1) = 1.0d0/(1.0d0+hh)*qq2(1)
            deltah(i,j,2) =                  qq2(2)
            deltah(i,j,3) = 1.0d0/(1.0d0+hh)*qq2(3)
         else if(j.eq.NY) then
            deltah(i,j,1) = 1.0d0/(1.0d0+hh)*((1.0d0+0.5d0*hh)*qq2(1)+(-0.5d0*hh)     *qq2(3))
            deltah(i,j,2) = qq2(2)
            deltah(i,j,3) = 1.0d0/(1.0d0+hh)*((-0.5d0*hh)     *qq2(1)+(1.0d0+0.5d0*hh)*qq2(3))
         endif
      enddo
      !Sweep 2
      do k = 1,NY
         j = NY-k+1
         !Calculation of UdQ
         if(j.eq.NY) then
            qq1(1) = 0.0d0
            qq1(2) = 0.0d0
            qq1(3) = 0.0d0
         else if(j.ge.1 .and. j.le.NY-1) then
            qq1(1) = hh*0.5d0*(-deltah(i,j+1,1)+deltah(i,j+1,3))
            qq1(2) = 0.0d0
            qq1(3) = hh*0.5d0*( deltah(i,j+1,1)-deltah(i,j+1,3))
         endif
         !Calculation of D^-1
         if(j.eq.1) then
            qq2(1) = 1.0d0/(1.0d0+hh)*((1.0d0+0.5d0*hh)*qq1(1)+( 0.5d0*hh)     *qq1(3))
            qq2(2) = qq1(2)
            qq2(3) = 1.0d0/(1.0d0+hh)*(( 0.5d0*hh)     *qq1(1)+(1.0d0+0.5d0*hh)*qq1(3))
         else if(j.ge.2 .and. j.le.NY-1) then
            qq2(1) = 1.0d0/(1.0d0+hh)*qq1(1)
            qq2(2) =                  qq1(2)
            qq2(3) = 1.0d0/(1.0d0+hh)*qq1(3)
         else if(j.eq.NY) then
            qq2(1) = 1.0d0/(1.0d0+hh)*((1.0d0+0.5d0*hh)*qq1(1)+(-0.5d0*hh)     *qq1(3))
            qq2(2) = qq1(2)
            qq2(3) = 1.0d0/(1.0d0+hh)*((-0.5d0*hh)     *qq1(1)+(1.0d0+0.5d0*hh)*qq1(3))
         endif
         !Calculation of delta
         delta(i,j,1) = deltah(i,j,1)-qq2(1)
         delta(i,j,2) = deltah(i,j,2)-qq2(2)
         delta(i,j,3) = deltah(i,j,3)-qq2(3)
      enddo
   enddo

   return
endsubroutine


!-------------------LUSGS method for HES--------------------------------

subroutine LUSGSHES2D(NX,NY,dt,aa,bb,cc,dd,ee,ff,gg,hh,dinv,rhs,delta)

   implicit none
   integer :: i,j,k,l
   integer,intent(in)                                  :: NX,NY
   double precision,intent(in)                         :: dt
   double precision,dimension(1:NX,1:NY,5),intent(in)  :: aa,bb,cc,dd,ee,ff,gg,hh
   double precision,dimension(1:NX,1:NY,9),intent(in)  :: dinv
   double precision,dimension(1:NX,1:NY,3),intent(in)  :: rhs
   double precision,dimension(1:NX,1:NY,3),intent(out) :: delta
   double precision,dimension(-1:NX+2,-1:NY+2,3)       :: del1,del2
   double precision :: xx1,xx2,xx3

   do j =-1,NY+2
      do i =-1,NX+2
         do k = 1,3
            del1(i,j,k) = 0.0d0
            del2(i,j,k) = 0.0d0
         enddo
      enddo
   enddo
   !Sweep 1
   do j = 1,NY
      do i = 1,NX
         xx1 = rhs(i,j,1)-dt*(&
               aa(i,j,1)*del1(i-2,j  ,1)+aa(i,j,2)*del1(i-1,j  ,1)&
              +bb(i,j,1)*del1(i-2,j  ,2)+bb(i,j,2)*del1(i-1,j  ,2)&
              +ee(i,j,1)*del1(i  ,j-2,1)+ee(i,j,2)*del1(i  ,j-1,1)&
              +ff(i,j,1)*del1(i  ,j-2,3)+ff(i,j,2)*del1(i  ,j-1,3))
         xx2 = rhs(i,j,2)-dt*(&
               cc(i,j,1)*del1(i-2,j  ,1)+cc(i,j,2)*del1(i-1,j  ,1)&
              +dd(i,j,1)*del1(i-2,j  ,2)+dd(i,j,2)*del1(i-1,j  ,2))
         xx3 = rhs(i,j,3)-dt*(&
               gg(i,j,1)*del1(i  ,j-2,1)+gg(i,j,2)*del1(i  ,j-1,1)&
              +hh(i,j,1)*del1(i  ,j-2,3)+hh(i,j,2)*del1(i  ,j-1,3))
         del1(i,j,1) = dinv(i,j,1)*xx1+dinv(i,j,2)*xx2+dinv(i,j,3)*xx3
         del1(i,j,2) = dinv(i,j,4)*xx1+dinv(i,j,5)*xx2+dinv(i,j,6)*xx3
         del1(i,j,3) = dinv(i,j,7)*xx1+dinv(i,j,8)*xx2+dinv(i,j,9)*xx3
      enddo
   enddo
   !Sweep 2
   do l = 1,NY
      do k = 1,NX
         i = NX-k+1
         j = NY-l+1
         xx1 = dt*(&
               aa(i,j,4)*del2(i+1,j  ,1)+aa(i,j,5)*del2(i+2,j  ,1)&
              +bb(i,j,4)*del2(i+1,j  ,2)+bb(i,j,5)*del2(i+2,j  ,2)&
              +ee(i,j,4)*del2(i  ,j+1,1)+ee(i,j,5)*del2(i  ,j+2,1)&
              +ff(i,j,4)*del2(i  ,j+1,3)+ff(i,j,5)*del2(i  ,j+2,3))
         xx2 = dt*(&
               cc(i,j,4)*del2(i+1,j  ,1)+cc(i,j,5)*del2(i+2,j  ,1)&
              +dd(i,j,4)*del2(i+1,j  ,2)+dd(i,j,5)*del2(i+2,j  ,2))
         xx3 = dt*(&
               gg(i,j,4)*del2(i  ,j+1,1)+gg(i,j,5)*del2(i  ,j+2,1)&
              +hh(i,j,4)*del2(i  ,j+1,3)+hh(i,j,5)*del2(i  ,j+2,3))
         del2(i,j,1) = del1(i,j,1)-(dinv(i,j,1)*xx1+dinv(i,j,2)*xx2+dinv(i,j,3)*xx3)
         del2(i,j,2) = del1(i,j,2)-(dinv(i,j,4)*xx1+dinv(i,j,5)*xx2+dinv(i,j,6)*xx3)
         del2(i,j,3) = del1(i,j,3)-(dinv(i,j,7)*xx1+dinv(i,j,8)*xx2+dinv(i,j,9)*xx3)
      enddo
   enddo
   do j = 1,NY
      do i = 1,NX
         do k = 1,3
            delta(i,j,k) = del2(i,j,k)
         enddo
      enddo
   enddo

   return
endsubroutine




!------------------- LUSGS -------------------

subroutine LUSGS(rhs,delta)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,j,k
   double precision,dimension(1:NX,1:NY,3),intent(in)   :: rhs
   double precision,dimension(1:NX,1:NY,3),intent(out)  :: delta
   double precision,dimension(0:NX+1,0:NY+1,3)          :: del1,del2
   double precision :: hhx,hhy
   integer :: bb

   bb = bound
   hhx = DTELE/DXL/2.0d0
   hhy = DTELE/DYL/2.0d0

   do j = 1,NY
      do i = 1,NX
         do k = 1,3
            delta(i,j,k) = 0.0d0
         enddo
      enddo
   enddo
   do j = 0,NY+1
      do i = 0,NX+1
         do k = 1,3
            del1(i,j,k) = 0.0d0
            del2(i,j,k) = 0.0d0
         enddo
      enddo
   enddo

   !LDU-factorization + SGS
   !Sweep 1
   do j = 1,NY
      do i = 1+bb,NX
         del1(i,j,1) = 1.0d0/(1.0d0+hhx+hhy)*(rhs(i,j,1)-(-hhx*del1(i-1,j,1)-hhx*del1(i-1,j,2)-hhy*del1(i,j-1,1)-hhy*del1(i,j-1,3)))
         del1(i,j,2) = 1.0d0/(1.0d0+hhx)    *(rhs(i,j,2)-(-hhx*del1(i-1,j,1)-hhx*del1(i-1,j,2)                                    ))
         del1(i,j,3) = 1.0d0/(1.0d0+hhy)    *(rhs(i,j,3)-(                                    -hhy*del1(i,j-1,1)-hhy*del1(i,j-1,3)))
      enddo
   enddo
   !Sweep 2
   do j = 1,NY
      do k = 1,NX-bb
         i = NX-k+1
         del2(i,j,1) = del1(i,j,1)-(1.0d0/(1.0d0+hhx+hhy)*(-hhx*del2(i+1,j,1)+hhx*del2(i+1,j,2)-hhy*del2(i,j+1,1)+hhy*del2(i,j+1,3)))
         del2(i,j,2) = del1(i,j,2)-(1.0d0/(1.0d0+hhx)    *( hhx*del2(i+1,j,1)-hhx*del2(i+1,j,2)                                    ))
         del2(i,j,3) = del1(i,j,3)-(1.0d0/(1.0d0+hhy)    *(                                     hhy*del2(i,j+1,1)-hhy*del2(i,j+1,3)))
      enddo
   enddo

   do j = 1,NY
      do i = 1+bb,NX
         do k = 1,3
            delta(i,j,k) = del2(i,j,k)
         enddo
      enddo
   enddo

   return
endsubroutine





