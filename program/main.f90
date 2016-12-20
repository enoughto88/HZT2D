!***********************************************************************
!*****       HZT2D: Axial-Azimuthal 2D Hall Thruster Simulator     *****
!*****          Developed by Rei Kawashima (Univ. of Tokyo)        *****
!*****             Latest Version Distributed at GitHub            *****
!***********************************************************************

program main
   use parameters_mod
   use global_mod
   implicit none
   character(len=6)  :: cname
   character(len=30) :: dname
   integer :: i,j,m,lpele
   integer :: time_pic,time_ele
   double precision,dimension(1:NX,1:NY)     :: nele                    ![m-3]  Electron number density
   double precision,dimension(1:NX,1:NY)     :: nneu                    ![m-3]  Neutral number density
   double precision,dimension(1:NX,1:NY)     :: Tele                    ![eV]   Electron temperature
   double precision,dimension(1:NX,1:NY)     :: fcol                    ![s-1]  Total collision frequency
   double precision,dimension(1:NX,1:NY)     :: fion                    ![s-1]  Ionization collision frequency
   double precision,dimension(1:NX,1:NY)     :: qion                    ![m-3s-1] Ion production rate
   double precision,dimension(1:NX,1:NY)     :: phii                    ![V]    Space potential
   double precision,dimension(1:NX,1:NY)     :: uele                    ![ms-1] x-Velocity
   double precision,dimension(1:NX,1:NY)     :: vele                    ![ms-1] y-Velocity
   double precision,dimension(1:NX,1:NY)     :: babs                    ![T] Radial magnetic flux density
   double precision,dimension(1:NX,1:NY)     :: Omee                    ![-] Electron Hall parameter
   double precision,dimension(1:NX,1:NY)     :: cond                    ![-] Electron conductivity
   double precision,dimension(3)             :: resn = 0.0d0
   double precision,dimension(3)             :: cfl  = 0.0d0
   double precision,dimension(1:NX+1,1:NY+1,2) :: efnd                  ![Vm-1] Electric field defined at nodes
   double precision,dimension(1:NX+1,1:NY+1) :: bfnd                    ![T] Radial magnetic flux density defined as nodes
   integer                                   :: nm                      !N. of particles
   double precision,dimension(NMAX,11)       :: pic                     !Particle information
   double precision :: vxmax,vymax
   integer :: nmi,nmn

   write(*,*) 'Program start...'

   !Initial setups
   call InitialField(Tele,phii,uele,vele)
   call InputParticle(nm,pic,nneu,nele)
   call InputMagneticField(babs,bfnd)

   MARCH: do it = 1,ITM
      if(mod(it,ISMP).eq.0) nstp = nstp+1
      if(mod(it,ISMP).eq.1) then
         momx    = 0.0d0
         inflow  = 0.0d0
         outflow = 0.0d0
         outneu  = 0.0d0
         oution = 0.0d0
         acurrent = 0.0d0
         do i = 1,3
            resn_ave(i) = 0.0d0
         enddo
         do j = 1,NY
            do i = 1,NX
               nele_ave(i,j) = 0.0d0
               nneu_ave(i,j) = 0.0d0
               qion_ave(i,j) = 0.0d0
               phii_ave(i,j) = 0.0d0
            enddo
         enddo
      endif

      call XenonCollision(Tele,nneu,nele,babs,nm,pic,fcol,fion,qion,Omee,cond)
      call MassMomentumHES(nele,Tele,phii,uele,vele,qion,Omee,cond,resn,cfl,lpele,time_ele)
      call ElectricField(phii,efnd)
      call Particle2D3V(nm,pic,efnd,bfnd,qion,nneu,nele,time_pic)

      if(mod(it,ISMP).eq.0) then
         nmn   = 0
         nmi   = 0
         vxmax = 0.0d0
         vymax = 0.0d0
         do m = 1,nm
            if(pic(m,7).lt.0.1d0) nmn = nmn+1
            if(pic(m,7).gt.0.9d0) nmi = nmi+1
            if(dabs(pic(m,3)).gt.vxmax) then
               vxmax = dabs(pic(m,3))
            endif
            if(dabs(pic(m,4)).gt.vymax) then
               vymax = dabs(pic(m,4))
            endif
         enddo
      endif

      call OutputDisplay
      call OutputFile
   enddo MARCH

   write(*,*) 'Program end...'

   stop
contains


!--------------- Linux console output for monitoring -------------------

   subroutine OutputDisplay


      if(mod(it,ISMP).eq.0) then
         write(*,'(A42)')         '*****************************************'
         write(*,'(A8,I10,A8)')  '        ',nstp,'-th step'
         write(*,*)               '---------------- Fluid ------------------'
         write(*,'(A20,E12.3,A)') '         Time step =',DTELE,' [s]'
         write(*,'(A20,I10,A)')   '       N. of Loops =',lpele,''
         write(*,'(A20,E12.3,A)') '     Anode current =',acurrent*ECH/dble(ISMP),' [A]'
         write(*,'(A20,E12.3,A)') '        Residual 1 =',resn_ave(1)/dble(ISMP),' '
         write(*,'(A20,E12.3,A)') '        Residual 2 =',resn_ave(2)/dble(ISMP),' '
         write(*,'(A20,E12.3,A)') '        Residual 3 =',resn_ave(3)/dble(ISMP),' '
         write(*,'(A20,E12.3,A)') '      CFL Number 1 =',cfl(1),' '
         write(*,'(A20,E12.3,A)') '      CFL Number 2 =',cfl(2),' '
         write(*,'(A20,E12.3,A)') '        SCE Number =',cfl(3),' '
         write(*,'(A20,I10,A)')   '          CPU Time =',time_ele,''
         write(*,*)               '---------------- PIC --------------------'
         write(*,'(A20,E12.3,A)') '         Time Step =',DTPIC,' [s]'
         write(*,'(A20,I10,A)')   '   N. of Particles =',nm,''
         write(*,'(A20,I10,A)')   '      \--- Neutral =',nmn,''
         write(*,'(A20,I10,A)')   '      \------- Ion =',nmi,''
         write(*,'(A20,E12.3,A)') '            Inflow =',inflow*MI/DTPIC/dble(ISMP),' [kg s-1]'
         write(*,'(A20,E12.3,A)') '           Outflow =',outflow*MI/DTPIC/dble(ISMP),' [kg s-1]'
         write(*,'(A20,E12.3,A)') '      \--- Neutral =',outneu*MI/DTPIC/dble(ISMP),' [kg s-1]'
         write(*,'(A20,E12.3,A)') '      \------- Ion =',oution*MI/DTPIC/dble(ISMP),' [kg s-1]'
         write(*,'(A20,E12.3,A)') '           x-Force =',momx/DTPIC/dble(ISMP),' [N]'
         write(*,'(A20,E12.3,A)') '           Fastest =',dmax1(vxmax,vymax),''
         write(*,'(A20,E12.3,A)') '        CFL Number =',dmax1(vxmax*DTPIC/DXL,vymax*DTPIC/DYL),''
         write(*,'(A20,I10,A)')   '          CPU Time =',time_pic,''
      endif

      return
   endsubroutine OutputDisplay



!----------------- File output for numerical results -------------------

   subroutine OutputFile

      if(mod(it,ISMP).eq.0) then
         write(cname,'(I6.6)') nstp
         dname = '/output/distribution/'
         open(unit=22,file=trim(TOPDIR)//trim(dname)//'distribution.'//&
            cname//'.dat',form='formatted',status='replace')
            do j = 1,NY
               do i = 1,NX
                  write(22,'(7E15.5)') &
                     nneu_ave(i,j),nele_ave(i,j),Tele(i,j),phii_ave(i,j),qion_ave(i,j),uele_ave(i,j),vele_ave(i,j)
               enddo
            enddo
         close(22)

         open(unit=24,file=trim(TOPDIR)//trim(dname)//'electron.'//&
            cname//'.dat',form='formatted',status='replace')
            do j = 1,NY
               do i = 1,NX
                  write(24,'(5E15.5)') &
                     uele(i,j),cond(i,j),Omee(i,j),fion(i,j),fcol(i,j)
               enddo
            enddo
         close(24)
         !open(unit=24,file=trim(TOPDIR)//trim(dname)//'nodedist.'//&
         !   cname//'.dat',form='formatted',status='replace')
         !   do j = 1,NY+1
         !      do i = 1,NX+1
         !         write(24,'(2E15.5)') &
         !            efnd(i,j,1),efnd(i,j,2)
         !      enddo
         !   enddo
         !close(24)
         dname = '/output/datafile/'
         open(unit=32,file=trim(TOPDIR)//trim(dname)//'performance.dat', &
            form='formatted',status='unknown',position='append')
            write(32,'(1I7,4E15.5)') nstp,outneu*MI/DTPIC/dble(ISMP),oution*MI/DTPIC/dble(ISMP),&
               momx/DTPIC/dble(ISMP),acurrent*ECH/dble(ISMP)
         close(32)

         dname = '/output/datafile/'
         open(unit=23,file=trim(TOPDIR)//trim(dname)//'resn.dat', &
            form='formatted',status='unknown',position='append')
               write(23,'(1I5,1I10,6E15.5)') &
                  nstp,lpele,resn_ave(1)/dble(ISMP),resn_ave(2)/dble(ISMP),resn_ave(3)/dble(ISMP),cfl(1),cfl(2),cfl(3)
         close(23)
      endif

      !Particle and field information for restart
      if(mod(it,PSMP).eq.0) then
         write(cname,'(I6.6)') nstp
         dname = '/output/datafile/'
         open(unit=23,file=trim(TOPDIR)//trim(dname)//'particle.'//&
            cname//'.dat',form='formatted',status='replace')
               do m = 1,nm
                  write(23,'(2E20.10,9E12.3)') &
                     pic(m,1),pic(m,2),pic(m,3),pic(m,4),pic(m,5),&
                     pic(m,6),pic(m,7),pic(m,8),pic(m,9),pic(m,10),&
                     pic(m,11)
               enddo
         close(23)
         write(cname,'(I6.6)') nstp
         dname = '/output/datafile/'
         open(unit=24,file=trim(TOPDIR)//trim(dname)//'field.'//&
            cname//'.dat',form='formatted',status='replace')
            do j = 1,NY
               do i = 1,NX
                  write(24,'(4E15.5)') &
                     cons(i,j),phii(i,j),uele(i,j),vele(i,j)
               enddo
            enddo
         close(24)
      endif

      return
   endsubroutine OutputFile

end program main




!------------------ Electric field calculation -------------------------

subroutine ElectricField(phii,efnd)

   use parameters_mod
   use global_mod
   !$ use omp_lib
   implicit none
   double precision,dimension(1:NX,1:NY),intent(in)          :: phii    ![V]    Space potential
   double precision,dimension(1:NX+1,1:NY+1,2),intent(out)   :: efnd    ![Vm-1] Electric field defined at nodes
   integer :: i,j,s,t


   !$omp parallel default(none),shared(phii,efnd),private(i,j,s,t)
   !$omp do
   do j = 1,NY+1
      do i = 1,NX+1
         if    (i.eq.1   ) then; s = i+1
         elseif(i.eq.NX+1) then; s = i-1
         else                  ; s = i
         endif
         if    (j.eq.1   ) then; t = j+1
         elseif(j.eq.NY+1) then; t = j-1
         else                  ; t = j
         endif
         efnd(i,j,1) =-0.5d0/DXL*(phii(s  ,t-1)-phii(s-1,t-1))&
                      -0.5d0/DXL*(phii(s  ,t  )-phii(s-1,t  ))
         efnd(i,j,2) =-0.5d0/DYL*(phii(s-1,t  )-phii(s-1,t-1))&
                      -0.5d0/DYL*(phii(s  ,t  )-phii(s  ,t-1))
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   return
endsubroutine




