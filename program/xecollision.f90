!***********************************************************************
!***** xecollision.f90: Xenon collision freq. calculator           *****
!***** Main input     : nn,te                                      *****
!***** Main output    : nuc,nui                                    *****
!***********************************************************************

subroutine XenonCollision(Tele,nneu,nele,babs,nm,pic,fcol,fion,qion,Omee,cond)

   use parameters_mod
   use global_mod
   !$ use omp_lib
   implicit none
   integer :: i,j,m,ip,jp
   integer                              ,intent(in)    :: nm
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)    :: nele
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)    :: babs
   double precision,dimension(1:NXMAX,1:NYMAX),intent(inout) :: nneu
   double precision,dimension(NPMAX,11)       ,intent(in)    :: pic
   double precision,dimension(1:NXMAX,1:NYMAX),intent(inout) :: Tele
   double precision,dimension(1:NXMAX,1:NYMAX),intent(inout) :: fcol,fion,qion,Omee,cond
   double precision,dimension(1:NXMAX,1:NYMAX) :: ncon,nacl

   call ElasticColFreq(nneu,Tele,fcol)
   call IonizationColFreq(nneu,Tele,fion)

   !Correction of ionization collision frequency
   do j = 1,NYMAX
      do i = 1,NXMAX
         nacl(i,j) = 0.0d0
         ncon(i,j) = cons(i,j)+nele(i,j)*fion(i,j)*DXL*DYL*ZL*DTPIC
      enddo
   enddo

   !$omp parallel default(none),shared(nm,nacl,pic),private(m,ip,jp)
   !$omp do
   do m = 1,nm
      if(pic(m,7).lt.0.1d0) then
         ip = int(pic(m,8)+0.1d0)
         jp = int(pic(m,9)+0.1d0)
         nacl(ip,jp) = nacl(ip,jp)+pic(m,6)
      endif
   enddo
   !$omp end do
   !$omp end parallel
   do j = 1,NYMAX
      do i = 1,NXMAX
         if(cons(i,j).gt.nacl(i,j)) then
            cons(i,j) = nacl(i,j)
            qion(i,j) = 0.0d0
         elseif(ncon(i,j).gt.nacl(i,j)) then
            qion(i,j)= (nacl(i,j)-cons(i,j))/(DTPIC*DXL*DYL*ZL)
         else
            qion(i,j)= dmax1(nele(i,j)*fion(i,j),dmin1(1.0d20,(nacl(i,j)-cons(i,j))/(DTPIC*DXL*DYL*ZL)))
         endif
      enddo
   enddo

   do j = 1,NYMAX
      do i = 1,NXMAX
         if(fcol(i,j).lt.NUMIN) then
            fcol(i,j) = NUMIN
            cond(i,j) = nele(i,j)*ECH/ME/fcol(i,j)
            Omee(i,j)= dmax1(ECH/ME/fcol(i,j)*babs(i,j),0.0d0)
            !Omee(i,j)= dmin1(-ECH/ME/fcol(i,j)*babs(i,j),0.0d0)
         else
            cond(i,j) = nele(i,j)*ECH/ME/fcol(i,j)
            Omee(i,j)= dmax1(ECH/ME/fcol(i,j)*babs(i,j),0.0d0)
            !Omee(i,j)= dmin1(-ECH/ME/fcol(i,j)*babs(i,j),0.0d0)
         endif
      enddo
   enddo

   if(it.ge.1) then
      do j = 1,NYMAX
         do i = 1,NXMAX
            qion_ave(i,j) = qion_ave(i,j)+qion(i,j)/dble(ISMP)
         enddo
      enddo
   endif

   return
end subroutine


! Calculation of collision frequency using six-order polinomial approximation

subroutine ElasticColFreq(nneu,Tele,fcol)

   use parameters_mod
   !$ use omp_lib
   implicit none
   integer :: i,j
   double precision :: rate_cof,te
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)    :: Tele
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)    :: nneu
   double precision,dimension(1:NXMAX,1:NYMAX),intent(out)   :: fcol
   double precision,dimension(0:6)  :: c

   c(6) =-2.3131961087d-22
   c(5) = 6.3282110910d-20
   c(4) =-6.8214347801d-18
   c(3) = 3.6642194855d-16
   c(2) =-1.0165083960d-14
   c(1) = 1.3495725561d-13
   c(0) =-1.3643694581d-14
   !$omp parallel default(none),shared(c,Tele,nneu,fcol),private(i,j,te,rate_cof)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
         te = Tele(i,j)
         rate_cof = dmax1(c(6)*te**6.0d0+c(5)*te**5.0d0+c(4)*te**4.0d0 &
                         +c(3)*te**3.0d0+c(2)*te**2.0d0+c(1)*te+c(0),0.0d0)
         fcol(i,j) = nneu(i,j)*rate_cof
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   do j = 1,NYMAX
      do i = 1,NXMAX
         if(fcol(i,j).lt.0.0d0) then
            write(*,*) 'Warning...Error in fcol',i,j,fcol(i,j),nneu(i,j)
         endif
      enddo
   enddo

   return
endsubroutine

! ----------------------------------------------------------------------
! Calculation of ionization speed using six-order polinomial approximation

subroutine IonizationColFreq(nneu,Tele,fion)

   use parameters_mod
   !$ use omp_lib
   implicit none
   integer :: i,j
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)    :: Tele
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)    :: nneu
   double precision,dimension(1:NXMAX,1:NYMAX),intent(out)   :: fion
   double precision  :: te,rate_cof
   double precision,dimension(0:6)  :: c

   c(6) = 1.9171338295d-23
   c(5) =-5.1064636650d-21
   c(4) = 5.2708323330d-19
   c(3) =-2.5997663607d-17
   c(2) = 5.6318241551d-16
   c(1) = 7.0937452968d-16
   c(0) =-3.9979581660d-15
   !$omp parallel default(none),shared(c,Tele,nneu,fion),private(i,j,te,rate_cof)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
         te = Tele(i,j)
         rate_cof = dmax1(c(6)*te**6.0d0+c(5)*te**5.0d0+c(4)*te**4.0d0 &
                         +c(3)*te**3.0d0+c(2)*te**2.0d0+c(1)*te+c(0),0.0d0)
         fion(i,j) = nneu(i,j)*rate_cof
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   do j = 1,NYMAX
      do i = 1,NXMAX
         if(fion(i,j).lt.0.0d0) then
            write(*,*) 'Warning...Error in fion',i,j,fion(i,j),nneu(i,j)
         endif
      enddo
   enddo

   return
endsubroutine
! ----------------------------------------------------------------------






