!************ initialcondition.f90 *************************************

subroutine InitialField(Tele,phii,uele,vele)

   use parameters_mod
   use global_mod
   implicit none
   character(len=30) :: dname
   integer :: i,j,ios
   double precision,dimension(1:NXMAX,1:NYMAX),intent(out)   :: Tele
   double precision,dimension(1:NXMAX,1:NYMAX),intent(out)   :: phii,uele,vele
   double precision,dimension(2) :: check

   do j = 1,NYMAX
      do i = 1,NXMAX
         phii(i,j) = PHIA-PHIA*(((dble(i)-0.5d0)*DXL)/XL)
         uele(i,j) = 0.0d0
         vele(i,j) = 0.0d0
      enddo
   enddo

   if(INPUTF.eq.1) then
      dname = '/data/'
      open(16,file=trim(TOPDIR)//trim(dname)//'initialphite.dat',form='formatted',&
         status='old',action='read',position='rewind')
         do j = 1,NYMAX
            do i = 1,NXMAX
               read(16,*,iostat=ios) phii(i,j),Tele(i,j)
                  if(ios.eq.-1) then
                     write(*,*) 'Error...initialphite.dat is wrong...'
                     stop
                  endif
            enddo
         enddo
         read(16,*,iostat=ios) check(1),check(2)
         if(ios.eq.0) then
            write(*,*) 'Error...initialphite.dat is wrong...'
            stop
         endif
      close(16)
      write(*,*) 'Loaded initialphite.dat...'

      open(17,file=trim(TOPDIR)//trim(dname)//'initialfield.dat',form='formatted',&
         status='old',action='read',position='rewind')
         do j = 1,NYMAX
            do i = 1,NXMAX
               read(17,*,iostat=ios) cons(i,j),phii(i,j),uele(i,j),vele(i,j)
                  if(ios.eq.-1) then
                     write(*,*) 'Error...inputfield.dat is wrong...'
                     stop
                  endif
               !Special treatment for flat distribution
               !cons(i,j) = 0.0d0
            enddo
         enddo
         read(17,*,iostat=ios) check(1),check(2)
         if(ios.eq.0) then
            write(*,*) 'Error...inputfield.dat is wrong...'
            stop
         endif
      close(17)
      write(*,*) 'Loaded initialfield.dat...'
   endif

   dname = '/output/datafile/'
   open(unit=12,file=trim(TOPDIR)//trim(dname)//'input.dat',&
      form='formatted',status='unknown')
      do j = 1,NYMAX
         do i = 1,NXMAX
            write(12,'(4E15.5)') &
               Tele(i,j),phii(i,j),uele(i,j),vele(i,j)
         enddo
      enddo
   close(12)
   write(*,*) 'Wrote input.dat...'


   return
endsubroutine



!--------------- Particle information input ----------------------------

subroutine InputParticle(nm,pic,nneu,nele)

   use parameters_mod
   use global_mod
   use mtmod
   character(len=30) :: dname
   integer                              ,intent(inout) :: nm
   double precision,dimension(NPMAX,11)  ,intent(inout) :: pic
   double precision,dimension(1:NXMAX,1:NYMAX),intent(inout) :: nneu
   double precision,dimension(1:NXMAX,1:NYMAX),intent(inout) :: nele
   integer          :: m,ios
   double precision :: aa,bb

   if(INPUTP.eq.0) then
      nm = NMINI
      do m = 1,nm
         pic(m,1)   = grnd()*XL
         pic(m,2)   = grnd()*YL
         aa = dsqrt(-2.0d0*KBOLTZ/MI*TI*dlog(grnd()))
         bb = 2.0d0*PI*grnd()
         pic(m,3)   = aa*dcos(bb)
         pic(m,4)   = aa*dsin(bb)
         aa = dsqrt(-2.0d0*KBOLTZ/MI*TI*dlog(grnd()))
         bb = 2.0d0*PI*grnd()
         pic(m,5)   = aa*dcos(bb)
         pic(m,6)   = MACP
         if(m.lt.int(0.5d0*dble(NMINI))) then
            pic(m,7)   = 0.0d0
         else
            pic(m,7)   = 1.0d0
         endif
         pic(m,8)   = dble(1+int(pic(m,1)/DXL))
         pic(m,9)   = dble(1+int(pic(m,2)/DYL))
         pic(m,10)  = (pic(m,1)-DXL*(pic(m,8)-1.0d0))/DXL
         pic(m,11)  = (pic(m,2)-DYL*(pic(m,9)-1.0d0))/DYL
      enddo
   else if(INPUTP.eq.1) then
      dname = '/data/'
      open(16,file=trim(TOPDIR)//trim(dname)//'initialpic.dat',form='formatted',&
         status='old',action='read',position='rewind')
         do m = 1,NPMAX
            read(16,*,iostat=ios) pic(m,1),pic(m,2),pic(m,3),pic(m,4),pic(m,6),&
                                  pic(m,7),pic(m,8),pic(m,9),pic(m,10),pic(m,11)
            if(ios.lt.0) exit
            !pic(m,3) = pic(m,3)+10.0d0
            !Special treatment of flat distribution for neutrals
            if(pic(m,7).lt.0.1d0) then
               pic(m,2)   = grnd()*YL
               aa = dsqrt(-2.0d0*KBOLTZ/MI*TI*dlog(grnd()))
               bb = 2.0d0*PI*grnd()
               pic(m,4)   = aa*dsin(bb)
               pic(m,5)   = aa*dcos(bb)
            endif
         enddo
      close(16)
      nm = m-1
      write(*,*) 'Loaded initialpic.dat... N. of particle = ',nm
   endif

   !Calculation of initial density
   call Density(nm,pic,nneu,nele)

   return
endsubroutine


!--------------- Magnetic field information input ----------------------

subroutine InputMagneticField(babs,bfnd)

   use parameters_mod
   use global_mod
   implicit none
   character(len=30) :: dname
   integer :: i,j,s,t,ios
   double precision,dimension(1:NXMAX,1:NYMAX)     :: babs                    ![T] Radial magnetic flux density
   double precision,dimension(1:NXMAX+1,1:NYMAX+1) :: bfnd                    ![T] Radial magnetic flux density defined as nodes
   double precision :: check


   if(INPUTM.eq.1) then
      dname = '/data/'
      open(18,file=trim(TOPDIR)//trim(dname)//'SPT100MF.dat',form='formatted',&
         status='old',action='read',position='rewind')
         do j = 1,NYMAX
            do i = 1,NXMAX
               read(18,*,iostat=ios) babs(i,j)
                  if(ios.eq.-1) then
                     write(*,*) 'Error...SPT100MF.dat is wrong...'
                     stop
                  endif
            enddo
         enddo
         read(18,*,iostat=ios) check
         if(ios.eq.0) then
            write(*,*) 'Error...SPT100MF.dat is wrong...'
            stop
         endif
      close(18)
      write(*,*) 'Loaded SPT100MF.dat...'
   endif

   do j = 1,NYMAX+1
      do i = 1,NXMAX+1
         if    (i.eq.1   ) then; s = i+1
         elseif(i.eq.NXMAX+1) then; s = i-1
         else                  ; s = i
         endif
         if    (j.eq.1   ) then; t = j+1
         elseif(j.eq.NYMAX+1) then; t = j-1
         else                  ; t = j
         endif
         bfnd(i,j) = 0.25d0*(babs(s-1,t-1)+babs(s,t-1)+babs(s-1,t)+babs(s,t))
      enddo
   enddo

   return
endsubroutine




