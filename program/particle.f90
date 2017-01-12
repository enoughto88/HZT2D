!***********************************************************************
!***** particle.f90: 2D3V PIC solver                               *****
!***** Main input : ph,qion                                        *****
!***** Calculate motion of heavy particles: nm,pic                 *****
!***** Main output: nn,ne,fwx,fwy                                  *****
!***********************************************************************

subroutine Particle2D3V(nm,pic,efnd,bfnd,qion,nneu,nele,time_pic)

   use parameters_mod
   use global_mod
   use mtmod
   !$ use omp_lib
   implicit none
   integer :: i,j
   integer :: time_pre,time_post,time_pic
   integer                                    ,intent(inout)    :: nm         !Number of particles
   double precision,dimension(NPMAX,11)       ,intent(inout)    :: pic        !Particle information
   double precision,dimension(1:NXMAX,1:NYMAX),intent(in)       :: qion       ![m-3s-1] Electron production rate
   double precision,dimension(1:NXMAX,1:NYMAX),intent(out)      :: nneu       ![m-3] Neutral number density
   double precision,dimension(1:NXMAX,1:NYMAX),intent(out)      :: nele       ![m-3] Electron number density
   double precision,dimension(1:NXMAX+1,1:NYMAX+1,2),intent(in) :: efnd       ![Vm-1] Electric field defined at nodes
   double precision,dimension(1:NXMAX+1,1:NYMAX+1),intent(in)   :: bfnd       ![T] Radial magnetic flux density defined as nodes
   double precision,dimension(NPMAX,2)                     :: efp
   double precision,dimension(NPMAX)                       :: bfp
   double precision,dimension(NPMAX,5)                     :: post


   call system_clock(count = time_pre)

   if(CALPIC.eq.1) then
      call InflowParticle_Differential(nm,pic)
      call Ionization_particle(nm,pic,qion)

      !Determine which cell particles are included
      call ApplyField(nm,pic,efnd,bfnd,efp,bfp)

      !Movement of particles based on Runge-Kutta method
      call Leapfrog_Particle(nm,pic,efp,bfp,post)
      !Treatment of particles on boundary
      call Judge_Particle(nm,pic,post)

      !Density calculation on cells
      call Density(nm,pic,nneu,nele)
   endif

   !Data accumulation for output of average
   if(it.ge.1) then
      do j = 1,NYMAX
         do i = 1,NXMAX
            nneu_ave(i,j) = nneu_ave(i,j)+nneu(i,j)/dble(ISMP)
            nele_ave(i,j) = nele_ave(i,j)+nele(i,j)/dble(ISMP)
         enddo
      enddo
   endif

   !Particle tracking for trajectory check
   if(TRAC.eq.1) then
      call Trajectory(pic)
   endif

   call system_clock(count = time_post)
   time_pic = time_post - time_pre

   return
endsubroutine

!----------------- Particle inflow from anode side ---------------------

subroutine Inflow_particle(nm,pic)

   use mtmod
   use parameters_mod
   use global_mod
   !$ use omp_lib
   implicit none
   integer :: m
   integer                            ,intent(inout) :: nm
   double precision,dimension(NPMAX,11),intent(inout) :: pic
   integer               :: npin
   double precision      :: cc,dd,psiz,vth2,flow

   vth2   = 2.0d0*KBOLTZ/MI*TWALL
   flow   = MDOT/MI*DTPIC
   npin   = int(flow/MACP)+1
   psiz   = flow/dble(npin)
   !Neutral from left side
   do m = 1,npin
      nm = nm+1
      cc = dsqrt(-vth2*dlog(grnd()))
      dd = 2.0d0*PI*grnd()
      pic(nm,3)  = dabs(cc*dsin(dd))+10.0d0
      pic(nm,4)  = cc*dcos(dd)
      cc = dsqrt(-vth2*dlog(grnd()))
      dd = 2.0d0*PI*grnd()
      pic(nm,5)  = cc*dcos(dd)
      pic(nm,1)  = grnd()*DTPIC*pic(nm,3)
      pic(nm,2)  = YL*grnd()
      pic(nm,6)  = psiz
      pic(nm,7)  = 0.0d0
      pic(nm,8)  = dble(1+int(pic(nm,1)/DXL))
      pic(nm,9)  = dble(1+int(pic(nm,2)/DYL))
      pic(nm,10) = (pic(nm,1)-DXL*(pic(nm,8)-1.0d0))/DXL
      pic(nm,11) = (pic(nm,2)-DYL*(pic(nm,9)-1.0d0))/DYL
   enddo
   inflow = inflow+flow

   return
endsubroutine


!***********************************************************************
!*****           Neutral particle inflow from anode side           *****
!*****      Sinusoidal mass flow rate in azimuthal direction       *****
!***********************************************************************

subroutine InflowParticle_Differential(nm,pic)

   use mtmod
   use parameters_mod
   use global_mod
   implicit none
   integer :: m
   integer                            ,intent(inout) :: nm
   double precision,dimension(NPMAX,11),intent(inout) :: pic
   integer               :: npin
   double precision      :: cc,dd,psiz,vth2,flow

   vth2   = 2.0d0*KBOLTZ/MI*TWALL
   flow   = MDOT/MI*DTPIC
   npin   = int(flow/MACP)+1
   psiz   = flow/dble(npin)
   !Neutral from left side
   do m = 1,npin
      nm = nm+1
      cc = dsqrt(-vth2*dlog(grnd()))
      dd = 2.0d0*PI*grnd()
      pic(nm,3)  = dabs(cc*dsin(dd))+10.0d0
      pic(nm,4)  = cc*dcos(dd)
      cc = dsqrt(-vth2*dlog(grnd()))
      dd = 2.0d0*PI*grnd()
      pic(nm,5)  = cc*dcos(dd)
      pic(nm,1)  = grnd()*DTPIC*pic(nm,3)
      pic(nm,2)  = YL*grnd()
      pic(nm,6)  = psiz*(1.0d0-1.0d0*dcos(2.0d0*2.0d0*PI*pic(nm,2)/YL))
      pic(nm,7)  = 0.0d0
      pic(nm,8)  = dble(1+int(pic(nm,1)/DXL))
      pic(nm,9)  = dble(1+int(pic(nm,2)/DYL))
      pic(nm,10) = (pic(nm,1)-DXL*(pic(nm,8)-1.0d0))/DXL
      pic(nm,11) = (pic(nm,2)-DYL*(pic(nm,9)-1.0d0))/DYL
   enddo
   inflow = inflow+flow

   return
endsubroutine




!------------ Particle generation by ionization ------------------------

subroutine Ionization_particle(nm,pic,qion)

   use mtmod
   use parameters_mod
   use global_mod
   !$ use omp_lib
   implicit none
   integer                                ,intent(inout) :: nm
   double precision,dimension(NPMAX,11)    ,intent(inout) :: pic
   double precision,dimension(1:NXMAX,1:NYMAX)  ,intent(in)    :: qion         ![m-3s-1] Electron production rate
   double precision,dimension(1:NXMAX,1:NYMAX)                 :: nacl
   integer         ,dimension(NPMAX)                      :: lcrd
   integer         ,dimension(NPMAX)                      :: npin
   double precision,dimension(NPMAX)                      :: psiz
   integer :: m,k,i,j,ip,jp,ngen

   do j = 1,NYMAX
      do i = 1,NXMAX
         nacl(i,j)   = 0.0d0
         cons(i,j)   = cons(i,j)+qion(i,j)*DXL*DYL*ZL*DTPIC
      enddo
   enddo

   do m = 1,nm
      if(pic(m,7).lt.0.1d0) then
         ip = int(pic(m,8)+0.1d0)
         jp = int(pic(m,9)+0.1d0)
         nacl(ip,jp) = nacl(ip,jp)+pic(m,6)
      endif
   enddo

   do j = 1,NYMAX
      do i = 1,NXMAX
         if(cons(i,j).gt.nacl(i,j)) then
            write(*,*) 'Warning...No neutral particle, program error1',i,j,cons(i,j),qion(i,j),nacl(i,j)
            cons(i,j) = 0.0d0
         else if(cons(i,j).lt.0.0d0) then
            write(*,*) 'Warning...Negative cons',i,j,cons(i,j),qion(i,j)
         endif
      enddo
   enddo

   ngen = 0
   !$omp parallel default(none),shared(lcrd,npin,psiz),private(m)
   !$omp do
   do m = 1,NPMAX
      lcrd(m) = 0
      npin(m) = 0
      psiz(m) = 0.0d0
   enddo
   !$omp end do
   !$omp end parallel
   do m = 1,nm
      if(pic(m,7).lt.0.1d0) then
         ip = int(pic(m,8)+0.1d0)
         jp = int(pic(m,9)+0.1d0)
         if(pic(m,6).le.cons(ip,jp)) then
            pic(m,7)    = 1.0d0
            cons(ip,jp) = cons(ip,jp)-pic(m,6)
         else
            if(cons(ip,jp).gt.MACI) then
               ngen        = ngen+1
               lcrd(ngen)  = m
               psiz(ngen)  = MACI
               pic(m,6)    = pic(m,6)-MACI
               cons(ip,jp) = cons(ip,jp)-MACI
            endif
         endif
      endif
   enddo

   do m = 1,ngen
      nm = nm+1
      k = lcrd(m)
      pic(nm,1)  = pic(k,1)
      pic(nm,2)  = pic(k,2)
      pic(nm,3)  = pic(k,3)
      pic(nm,4)  = pic(k,4)
      pic(nm,5)  = pic(k,5)
      pic(nm,6)  = psiz(m)
      pic(nm,7)  = 1.0d0
      pic(nm,8)  = pic(k,8)
      pic(nm,9)  = pic(k,9)
      pic(nm,10) = pic(k,10)
      pic(nm,11) = pic(k,11)
   enddo


   return
endsubroutine



!------------ Determine which cell particle is included ----------------

subroutine ApplyField(nm,pic,efnd,bfnd,efp,bfp)

   use parameters_mod
   use global_mod
   !$ use omp_lib
   implicit none
   integer                              ,intent(in)       :: nm
   double precision,dimension(NPMAX,11)  ,intent(in)       :: pic        !Particle information
   double precision,dimension(1:NXMAX+1,1:NYMAX+1,2),intent(in) :: efnd       ![Vm-1] Electric field defined at nodes
   double precision,dimension(1:NXMAX+1,1:NYMAX+1),intent(in)   :: bfnd       ![T] Radial magnetic flux density defined as nodes
   double precision,dimension(NPMAX,2),intent(out)         :: efp
   double precision,dimension(NPMAX)  ,intent(out)         :: bfp
   integer :: m,ip,jp
   double precision :: xc,yc
   double precision,dimension(4) :: p


   !$omp parallel default(none),shared(nm,pic,efnd,bfnd,efp,bfp),private(m,ip,jp,xc,yc,p)
   !$omp do
   do m = 1,nm
      if(pic(m,7).lt.0.1d0) then
         efp(m,1) = 0.0d0
         efp(m,2) = 0.0d0
         bfp(m)   = 0.0d0
      else if(pic(m,7).gt.0.9d0) then
         ip = int(pic(m,8)+0.1d0)
         jp = int(pic(m,9)+0.1d0)
         if(ip.le.0 .or. ip.ge.NXMAX+1 .or. &
            jp.le.0 .or. jp.ge.NYMAX+1) then
            write(*,*) 'Warning...Irregular particle position...Belong',m,ip,jp,pic(m,1),pic(m,2)
         else
            xc  = pic(m,10)
            yc  = pic(m,11)
            if(xc.lt.0.0d0 .or. xc.gt.1.0d0) write(*,'(A20,2E12.3,1I5)') 'Error in Belong xc',xc,pic(m,1),ip
            if(yc.lt.0.0d0 .or. yc.gt.1.0d0) write(*,'(A20,2E12.3,1I5)') 'Error in Belong yc',yc,pic(m,2),jp
            p(1) = (1.0d0-xc)*(1.0d0-yc)
            p(2) = xc        *(1.0d0-yc)
            p(3) = (1.0d0-xc)*yc
            p(4) = xc        *yc
            efp(m,1) = efnd(ip  ,jp  ,1)*p(1)+efnd(ip+1,jp  ,1)*p(2) &
                     + efnd(ip  ,jp+1,1)*p(3)+efnd(ip+1,jp+1,1)*p(4)
            efp(m,2) = efnd(ip  ,jp  ,2)*p(1)+efnd(ip+1,jp  ,2)*p(2) &
                     + efnd(ip  ,jp+1,2)*p(3)+efnd(ip+1,jp+1,2)*p(4)
            bfp(m)   = bfnd(ip  ,jp  )  *p(1)+bfnd(ip+1,jp  )  *p(2) &
                     + bfnd(ip  ,jp+1)  *p(3)+bfnd(ip+1,jp+1)  *p(4)
         endif
      endif
   enddo
   !$omp end do
   !$omp end parallel

   return
endsubroutine


!********** Treatment of particles on boundary *************************

subroutine Judge_Particle(nm,pic,post)

   use parameters_mod
   use global_mod
   !$ use omp_lib
   implicit none
   integer                            ,intent(inout)  :: nm
   double precision,dimension(NPMAX,11),intent(inout)  :: pic
   double precision,dimension(NPMAX,5) ,intent(inout)  :: post
   integer                 :: ndel
   integer,dimension(NPMAX) :: lcrd,search
   integer                 :: m,k,l
   double precision        :: dt
   double precision        :: xx1,yy1,xx2,yy2
   double precision        :: vx1,vy1,vz1,vx2,vy2,vz2


   ! Initialization
   ndel = 0
   !$omp parallel default(none),shared(search,lcrd),private(m)
   !$omp do
   do m = 1,NPMAX
      search(m) = 0
      lcrd(m) = 0
   enddo
   !$omp end do
   !$omp end parallel

   JUD: do m = 1,nm
      if(pic(m,6).lt.0.0d0) then
         ndel = ndel+1
         lcrd(ndel) = m
         search(m) = ndel
         write(*, *) 'Irregular deleting...negative density particle',m,pic(m,6)
         cycle JUD
      endif
      if(post(m,1).gt.0.0d0 .and. post(m,1).lt.XL .and. &
         post(m,2).gt.0.0d0 .and. post(m,2).lt.YL) cycle JUD
      xx1 = pic(m,1)
      yy1 = pic(m,2)
      vx1 = pic(m,3)
      vy1 = pic(m,4)
      vz1 = pic(m,5)
      xx2 = post(m,1)
      yy2 = post(m,2)
      vx2 = post(m,3)
      vy2 = post(m,4)
      vz2 = post(m,5)
      dt  = DTPIC
      BCN: do k = 1,10
         if(      xx2.gt.0.0d0 .and. xx2.lt.XL &
            .and. yy2.gt.0.0d0 .and. yy2.lt.YL) then
            exit BCN
         !In the case of x <= 0
         else if(xx2.le.0.0d0) then
            xx2 = 0.0d0-(xx2-0.0d0)
            vx2 = dabs(vx2)
         !In the case of x >= XL
         else if(xx2.ge.XL) then
            ndel = ndel+1
            lcrd(ndel) = m
            search(m) = ndel
            cycle JUD
         !In the case of y <= 0
         else if(yy2.le.0.0d0) then
            yy2 = YL+(yy2-0.0d0)
         !In the case of y >= YL
         else if(yy2.ge.YL) then
            yy2 = 0.0d0+(yy2-YL)
         endif
      enddo BCN
      if(k.eq.10) write(*,*) 'Warning: k=10 at BCN...'
      if(xx2.lt.0.0d0 .or. xx2.gt.XL .or. yy2.lt.0.0d0 .or. yy2.gt.YL) then
         ndel = ndel+1
         lcrd(ndel) = m
         search(m) = ndel
         write(*, *) 'Irregular deleting in BCN...'
         cycle JUD
      endif
      post(m,1) = xx2
      post(m,2) = yy2
      post(m,3) = vx2
      post(m,4) = vy2
      post(m,5) = vz2
   enddo JUD

   !Updating
   !$omp parallel default(none),shared(nm,pic,post),private(m,k)
   !$omp do
   do m = 1,nm
      do k = 1,5
         pic(m,k) = post(m,k)
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   !Delete flow-out particles
   do m = 1,ndel
      k = lcrd(m)
      if(k.lt.1 .or. k.gt.nm) then
         write(*,*) 'Error in particle deleting...',k,nm,ndel,m
         cycle
      endif
      outflow = outflow+pic(k,6)
      if(pic(k,7).lt.0.1d0) then
         outneu = outneu+pic(k,6)
      elseif(pic(k,7).gt.0.9d0) then
         oution = oution+pic(k,6)
      endif
      momx = momx+MI*pic(k,6)*pic(k,3)
      if(search(nm).gt.0) then
         lcrd(search(nm)) = k
         search(k) = search(nm)
      endif
      do l = 1,11
         pic(k,l) = pic(nm,l)
      enddo
      nm = nm-1
   enddo


   return
endsubroutine



!---------- Leap-frog method -------------------------------------------

subroutine Leapfrog_Particle(nm,pic,efp,bfp,post)

   use parameters_mod
   use global_mod
   !$ use omp_lib
   implicit none
   integer                            ,intent(in)  :: nm
   double precision,dimension(NPMAX,11),intent(in)  :: pic
   double precision,dimension(NPMAX,2) ,intent(in)  :: efp
   double precision,dimension(NPMAX)   ,intent(in)  :: bfp
   double precision,dimension(NPMAX,5) ,intent(out) :: post
   integer                                         :: m

   !$omp parallel default(none),shared(nm,pic,efp,bfp,post),private(m)
   !$omp do
   do m = 1,nm
      post(m,1) = pic(m,1)+DTPIC*(pic(m,3)+DTPIC*ECH/MI*efp(m,1)-DTPIC*ECH/MI*pic(m,4)*bfp(m))
      post(m,2) = pic(m,2)+DTPIC*(pic(m,4)+DTPIC*ECH/MI*efp(m,2)+DTPIC*ECH/MI*pic(m,3)*bfp(m))
      post(m,3) = pic(m,3)+DTPIC*ECH/MI*efp(m,1)-DTPIC*ECH/MI*pic(m,4)*bfp(m)
      post(m,4) = pic(m,4)+DTPIC*ECH/MI*efp(m,2)+DTPIC*ECH/MI*pic(m,3)*bfp(m)
      post(m,5) = pic(m,5)
   enddo
   !$omp end do
   !$omp end parallel

   return
endsubroutine



!********** Density calculator of particles ****************************

subroutine Density(nm,pic,nneu,nele)
   use parameters_mod
   use global_mod
   !$ use omp_lib
   implicit none
   integer                            ,intent(in)    :: nm
   double precision,dimension(NPMAX,11),intent(inout) :: pic
   double precision,dimension(1:NXMAX,1:NYMAX),intent(out) :: nneu
   double precision,dimension(1:NXMAX,1:NYMAX),intent(out) :: nele
   double precision,dimension(1:NXMAX+1,1:NYMAX+1)         :: nabdn,nabdi
   double precision,dimension(1:NXMAX+1,1:NYMAX+1)         :: nandn,nandi
   double precision,dimension(1:NXMAX,1:NYMAX)             :: nacln,nacli
   integer          :: i,j,m,ip,jp
   double precision :: xc,yc
   double precision,dimension(4) :: p


   do j = 1,NYMAX+1
      do i = 1,NXMAX+1
         nandn(i,j) = 0.0d0
         nandi(i,j) = 0.0d0
      enddo
   enddo
   !$omp parallel default(none),shared(nm,pic,nandn,nandi),private(m,ip,jp,xc,yc,p)
   !$omp do
   do m = 1,nm
      ip = 1+int(pic(m,1)/DXL)
      jp = 1+int(pic(m,2)/DYL)
      if(pic(m,6).lt.0.0d0) then
         write(*,*) 'Warning...Negative particle ',m,ip,jp,pic(m,6)
      endif
      if(ip.le.0 .or. ip.ge.NXMAX+1 .or. &
         jp.le.0 .or. jp.ge.NYMAX+1) then
         write(*,*) 'Warning...Irregular particle position...density',m,ip,jp,pic(m,1),pic(m,2)
      else
         xc  = (pic(m,1)-DXL*(dble(ip)-1.0d0))/DXL
         yc  = (pic(m,2)-DYL*(dble(jp)-1.0d0))/DYL
         if(xc.lt.0.0d0 .or. xc.gt.1.0d0) write(*,*) 'Error in Density...xc',m,xc,pic(m,1),ip
         if(yc.lt.0.0d0 .or. yc.gt.1.0d0) write(*,*) 'Error in Density...yc',m,yc,pic(m,2),jp
         p(1) = (1.0d0-xc)*(1.0d0-yc)
         p(2) = xc        *(1.0d0-yc)
         p(3) = (1.0d0-xc)*yc
         p(4) = xc        *yc
         !--------------------------------------------------------
         if    (pic(m,7).lt.0.1d0) then
            nandn(ip  ,jp  ) = nandn(ip  ,jp  )+p(1)*pic(m,6)
            nandn(ip+1,jp  ) = nandn(ip+1,jp  )+p(2)*pic(m,6)
            nandn(ip  ,jp+1) = nandn(ip  ,jp+1)+p(3)*pic(m,6)
            nandn(ip+1,jp+1) = nandn(ip+1,jp+1)+p(4)*pic(m,6)
         elseif(pic(m,7).gt.0.9d0) then
            nandi(ip  ,jp  ) = nandi(ip  ,jp  )+p(1)*pic(m,6)
            nandi(ip+1,jp  ) = nandi(ip+1,jp  )+p(2)*pic(m,6)
            nandi(ip  ,jp+1) = nandi(ip  ,jp+1)+p(3)*pic(m,6)
            nandi(ip+1,jp+1) = nandi(ip+1,jp+1)+p(4)*pic(m,6)
         endif
      endif
      pic(m,8)  = dble(ip)
      pic(m,9)  = dble(jp)
      pic(m,10) = xc
      pic(m,11) = yc
   enddo
   !$omp end do
   !$omp end parallel

   do j = 1,NYMAX+1
      do i = 1,NXMAX+1
         if(j.eq.1 ) then
            nabdn(i,j) = nandn(i,1)+nandn(i,NYMAX+1)
            nabdi(i,j) = nandi(i,1)+nandi(i,NYMAX+1)
         else if(j.eq.NYMAX+1) then
            nabdn(i,j) = nandn(i,NYMAX+1)+nandn(i,1)
            nabdi(i,j) = nandi(i,NYMAX+1)+nandi(i,1)
         else
            nabdn(i,j) = nandn(i,j)
            nabdi(i,j) = nandi(i,j)
         endif
      enddo
   enddo
   do j = 1,NYMAX+1
      nabdn(1,j)    = 2.0d0*nabdn(1,j)
      nabdi(1,j)    = 2.0d0*nabdi(1,j)
      nabdn(NXMAX+1,j) = 2.0d0*nabdn(NXMAX+1,j)
      nabdi(NXMAX+1,j) = 2.0d0*nabdi(NXMAX+1,j)
      !nabdn(1,j)    = nabdn(2,j)
      !nabdi(1,j)    = nabdi(2,j)
      !nabdn(NXMAX+1,j) = nabdn(NXMAX,j)
      !nabdi(NXMAX+1,j) = nabdi(NXMAX,j)
   enddo

   !$omp parallel default(none),shared(nandn,nacln,nneu,nandi,nacli,nele),private(i,j)
   !$omp do
   do j = 1,NYMAX
      do i = 1,NXMAX
         nacln(i,j) = 0.25d0*(nabdn(i,j)+nabdn(i+1,j)+nabdn(i,j+1)+nabdn(i+1,j+1))
         nacli(i,j) = 0.25d0*(nabdi(i,j)+nabdi(i+1,j)+nabdi(i,j+1)+nabdi(i+1,j+1))
         nneu(i,j) = nacln(i,j)/DXL/DYL/ZL
         nele(i,j) = dmax1(nacli(i,j)/DXL/DYL/ZL,NEMIN)
      enddo
   enddo
   !$omp end do
   !$omp end parallel



   return
endsubroutine



!--------------- Tracking trajectories of particles --------------------

subroutine Trajectory(pic)

   use parameters_mod
   implicit none
   double precision,dimension(NPMAX,11)  ,intent(inout)  :: pic          !Particle information
   character(len=30) :: dname

   !Generate trajectory information of 50 macroparticles
   dname = '/output/datafile/'
   open(unit=43,file=trim(TOPDIR)//trim(dname)//'trajectory.dat', &
      form='formatted',status='unknown',position='append')
         write(43,'(30E15.5)') &
            pic(int(dble(NMINI)*0.05d0),1),pic(int(dble(NMINI)*0.05d0),2),pic(int(dble(NMINI)*0.05d0),7), &
            pic(int(dble(NMINI)*0.10d0),1),pic(int(dble(NMINI)*0.10d0),2),pic(int(dble(NMINI)*0.10d0),7), &
            pic(int(dble(NMINI)*0.15d0),1),pic(int(dble(NMINI)*0.15d0),2),pic(int(dble(NMINI)*0.15d0),7), &
            pic(int(dble(NMINI)*0.20d0),1),pic(int(dble(NMINI)*0.20d0),2),pic(int(dble(NMINI)*0.20d0),7), &
            pic(int(dble(NMINI)*0.25d0),1),pic(int(dble(NMINI)*0.25d0),2),pic(int(dble(NMINI)*0.25d0),7), &
            pic(int(dble(NMINI)*0.30d0),1),pic(int(dble(NMINI)*0.30d0),2),pic(int(dble(NMINI)*0.30d0),7), &
            pic(int(dble(NMINI)*0.35d0),1),pic(int(dble(NMINI)*0.35d0),2),pic(int(dble(NMINI)*0.35d0),7), &
            pic(int(dble(NMINI)*0.40d0),1),pic(int(dble(NMINI)*0.40d0),2),pic(int(dble(NMINI)*0.40d0),7), &
            pic(int(dble(NMINI)*0.45d0),1),pic(int(dble(NMINI)*0.45d0),2),pic(int(dble(NMINI)*0.45d0),7), &
            pic(int(dble(NMINI)*0.50d0),1),pic(int(dble(NMINI)*0.50d0),2),pic(int(dble(NMINI)*0.50d0),7)
   close(43)

   return
endsubroutine



