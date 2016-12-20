!---------- 2D Time-Dependent Electron Energy Equation Solver ----------
!-------------------- Last Update 2015/11/11 ---------------------------

subroutine EnergySAD2D(nele,uele,vele,Tele,phii,resn)

   use parameters_mod
   use global_mod
   implicit none
   character(len=6) :: cname
   character(len=30) :: dname
   integer :: i,j
   double precision,dimension(1:III,1:JJJ),intent(in)    :: nele        !Electron density
   double precision,dimension(1:III,1:JJJ),intent(in)    :: uele        !x-Electron velocity
   double precision,dimension(1:III,1:JJJ),intent(in)    :: vele        !y-Electron velocity
   !double precision,dimension(1:III,1:JJJ),intent(in)    :: Sene        !Energy source term
   double precision,dimension(1:III,1:JJJ)               :: muex,muey   !Electron mobility
   double precision,dimension(1:III,1:JJJ),intent(inout) :: Tele        !Electron temperature
   double precision,dimension(1:III,1:JJJ),intent(in)    :: phii        !Space potential
   double precision,dimension(-1:III+2,-1:JJJ+2)         :: var         !Variable including BC
   double precision,dimension(-1:III+2,-1:JJJ+2)         :: advx,advy   !Coefficient for advection term including BC
   double precision,dimension(0:III+1,0:JJJ+1)           :: difx,dify   !Coefficient for diffusion term including BC
   double precision,dimension(1:III,1:JJJ)               :: sce         !Source term
   double precision,dimension(1:III,1:JJJ)               :: delta       !Delta of electron internal energy
   double precision,dimension(0:III,1:JJJ)               :: cfl         !CFL number
   double precision,dimension(4),intent(inout)           :: resn
   double precision :: scex,scey

   call Mobility(muex,muey)
   !Internal energy is chosen for variable
   do j = 1,JJJ
      do i = 1,III
         var(i,j)  = 3.0d0/2.0d0*nele(i,j)*Tele(i,j)
         difx(i,j) = 5.0d0/3.0d0*muex(i,j)*Tele(i,j)
         dify(i,j) = 5.0d0/3.0d0*muey(i,j)*Tele(i,j)
         if(i.eq.1) then
            advx(i,j) = 5.0d0/3.0d0*uele(i,j)+muex(i,j)*Tele(i,j)&
               *1.0d0/2.0d0/DXL*(-3.0d0*dlog(nele(i,j))+4.0d0*dlog(nele(i+1,j))-dlog(nele(i+2,j)))
         else if(i.eq.III) then
            advx(i,j) = 5.0d0/3.0d0*uele(i,j)+muex(i,j)*Tele(i,j)&
               *1.0d0/2.0d0/DXL*( 3.0d0*dlog(nele(i,j))-4.0d0*dlog(nele(i-1,j))+dlog(nele(i-2,j)))
         else if(i.ge.2 .and. i.le.III-1) then
            advx(i,j) = 5.0d0/3.0d0*uele(i,j)+muex(i,j)*Tele(i,j)&
               *1.0d0/2.0d0/DXL*(dlog(nele(i+1,j))-dlog(nele(i-1,j)))
         endif
         if(j.eq.1) then
            advy(i,j) = 5.0d0/3.0d0*vele(i,j)+muey(i,j)*Tele(i,j)&
               *1.0d0/2.0d0/DYL*(-3.0d0*dlog(nele(i,j))+4.0d0*dlog(nele(i,j+1))-dlog(nele(i,j+2)))
         else if(j.eq.JJJ) then
            advy(i,j) = 5.0d0/3.0d0*vele(i,j)+muey(i,j)*Tele(i,j)&
               *1.0d0/2.0d0/DYL*( 3.0d0*dlog(nele(i,j))-4.0d0*dlog(nele(i,j-1))+dlog(nele(i,j-2)))
         else if(j.ge.2 .and. j.le.JJJ-1) then
            advy(i,j) = 5.0d0/3.0d0*vele(i,j)+muey(i,j)*Tele(i,j)&
               *1.0d0/2.0d0/DYL*(dlog(nele(i,j+1))-dlog(nele(i,j-1)))
         endif
      enddo
   enddo
   do j = 1,JJJ
      do i = 1,III
         if(i.eq.1) then
            scex  = nele(i,j)*uele(i,j)*1.0d0/2.0d0/DXL*(phii(i+1,j  )-phii(i-1,j  ))
         else if(i.ge.2 .and. i.le.III-1) then
            scex  = nele(i,j)*uele(i,j)*1.0d0/2.0d0/DXL*(phii(i+1,j  )-phii(i-1,j  ))
         else if(i.eq.III) then
            scex  = nele(i,j)*uele(i,j)*1.0d0/2.0d0/DXL*(phii(i+1,j  )-phii(i-1,j  ))
         endif
         if(j.eq.1) then
            scey  = nele(i,j)*vele(i,j)*1.0d0/2.0d0/DYL*(phii(i  ,j+1)-phii(i  ,j-1))
         else if(j.ge.2 .and. j.le.JJJ-1) then
            scey  = nele(i,j)*vele(i,j)*1.0d0/2.0d0/DYL*(phii(i  ,j+1)-phii(i  ,j-1))
         else if(j.eq.JJJ) then
            scey  = nele(i,j)*vele(i,j)*1.0d0/2.0d0/DYL*(phii(i  ,j+1)-phii(i  ,j-1))
         endif
         sce(i,j)  = scex+scey
      enddo
   enddo

   !Boundary condition
   do j = 1,JJJ
      var(0    ,j)  = var(1,j)
      var(-1   ,j)  = 2.0d0*var(0,j)-var(1,j)
      var(III+1,j)  = 2.0d0*(1.5d0*nele(i,j)*TEC)-var(III,j)
      var(III+2,j)  = 2.0d0*var(III+1,j)-var(III,j)
      difx(0    ,j) = difx(1,j)
      difx(III+1,j) = difx(III,j)
      advx(0    ,j) = advx(1  ,j)
      advx(-1   ,j) = 2.0d0*advx(0    ,j)-advx(1  ,j)
      advx(III+1,j) = advx(III,j)
      advx(III+2,j) = 2.0d0*advx(III+1,j)-advx(III,j)
   enddo
   do i = 1,III
      var(i,0    )  = var(i,JJJ)
      var(i,-1   )  = var(i,JJJ-1)
      var(i,JJJ+1)  = var(i,1)
      var(i,JJJ+2)  = var(i,2)
      dify(i,0    ) = dify(i,JJJ)
      dify(i,JJJ+1) = dify(i,1)
      advy(i,0    ) = advy(i,JJJ)
      advy(i,-1   ) = advy(i,JJJ-1)
      advy(i,JJJ+1) = advy(i,1)
      advy(i,JJJ+2) = advy(i,2)
   enddo

   cfl(0,1) = 0.0d0
   do j = 1,JJJ
      do i = 1,III
         cfl(i,j) = dmax1(cfl(i-1,j),dabs(advx(i,j))*DTELE/DXL,dabs(advy(i,j))*DTELE/DYL)
      enddo
      if(j.ne.JJJ) cfl(0,j+1) = cfl(i,j)
   enddo
   !if(cfl(III,JJJ).gt.SETCFL) write(*,'(A16,1E12.3)') 'Too large CFL...',cfl(III,JJJ)

   call AdvectionDiffusion2D(III,JJJ,DTELE,DXL,DYL,var,advx,advy,difx,dify,sce,delta)

   !Update
   do j = 1,JJJ
      do i = 1,III
         var(i,j) = var(i,j)+delta(i,j)
      enddo
   enddo

   do j = 1,JJJ
      do i = 1,III
         Tele(i,j) = var(i,j)/(3.0d0/2.0d0*nele(i,j))
         if(Tele(i,j).lt.0.0d0) then
            write(*,'(A26,2I5,1E12.3)') 'Warning...negative Tele...',i,j,Tele(i,j)
            Tele(i,j) = 0.0d0
         endif
      enddo
   enddo

   resn(4) = 0.0d0
   do j = 1,JJJ
      do i = 1,III
         resn(4) = resn(4)+(dabs(delta(i,j))/(dabs(var(i,j))+1.0d-3))**2.0d0
      enddo
   enddo
   resn(4) = dsqrt(resn(4)/dble(III)/dble(JJJ))

   if(mod(it,ISMP) .eq. 0) then
      !Information for monitoring
      write(*,'(A38)')         '**************************************'
      write(*,'(A10,I10,A8)') '          ',nstp,'-th step'
      write(*,'(A20,E12.3,A)') '         time step =',DTELE,' [s]'
      write(*,'(A20,E12.3,A)') '        residual 4 =',resn(4),' '
      write(*,'(A20,E12.3,A)') '      cfl number 4 =',cfl(III,JJJ),' '

      write(cname,'(I6.6)') nstp
      dname = '/output/distribution/'
      open(unit=22,file=trim(TOPDIR)//trim(dname)//'distribution2.'//&
         cname//'.dat',form='formatted',status='replace')
         do j = 1,JJJ
            do i = 1,III
               write(22,'(2I5,1E15.5)') &
                  i,j,Tele(i,j)
            enddo
         enddo
      close(22)
      dname = '/output/datafile/'
      open(unit=23,file=trim(TOPDIR)//trim(dname)//'resn2.dat', &
         form='formatted',status='unknown',position='append')
            write(23,'(1I5,1E15.5)') &
               nstp,resn(4)
      close(23)

   endif

   return
endsubroutine



!-------- 2D Advection-Diffusion Equation Steady-State Solver ----------
!      Input : NX,NY,dt,dx,dy,u(-1:NX+2,-1:NY+2)
!              C(-1:NX+2,-1:NY+2),D(0:NX+1,0:NY+1),S(1,NX,1,NY)
!      Output: delta(1:NX,1:NY)
!-----------------------------------------------------------------------

subroutine AdvectionDiffusion2D(NX,NY,dt,dx,dy,u,cAx,cAy,cDx,cDy,S,delta)

   implicit none
   integer :: i,j
   integer,intent(in)                                     :: NX,NY      !Number of cells
   double precision,intent(in)                            :: dt         !Time interval
   double precision,intent(in)                            :: dx,dy      !Space interval
   double precision,dimension(-1:NX+2,-1:NY+2),intent(in) :: u          !Variable including BC
   double precision,dimension(-1:NX+2,-1:NY+2),intent(in) :: cAx,cAy    !Advection coefficient including BC
   double precision,dimension(0:NX+1,0:NY+1),intent(in)   :: cDx,cDy    !Diffusion coefficient including BC
   double precision,dimension(1:NX,1:NY),intent(in)       :: S          !Source term
   double precision,dimension(1:NX,1:NY),intent(out)      :: delta      !Delta of variable in implicit form
   double precision,dimension(1:NX,1:NY)                  :: rhs        !Right-hand side in explicit form
   double precision,dimension(1:NX,1:NY,5)                :: aAx,aAy
   double precision,dimension(1:NX,1:NY,3)                :: aDx,aDy
   double precision,dimension(1:NX,1:NY,5)                :: aFx,aFy

   !Coefficient calculation of advection and diffusion terms
   call AdvectionCoefficient2D_X(NX,NY,dt,dx,cAx,u,aAx)
   call AdvectionCoefficient2D_Y(NX,NY,dt,dy,cAy,u,aAy)
   call DiffusionCoefficient2D_X(NX,NY,dt,dx,cDx,aDx)
   call DiffusionCoefficient2D_Y(NX,NY,dt,dy,cDy,aDy)

   !Element calculation in 5-band matrix
   do j = 1,NY
      do i = 1,NX
         aFx(i,j,1) = aAx(i,j,1)
         aFx(i,j,2) = aAx(i,j,2)+aDx(i,j,1)
         aFx(i,j,3) = aAx(i,j,3)+aDx(i,j,2)
         aFx(i,j,4) = aAx(i,j,4)+aDx(i,j,3)
         aFx(i,j,5) = aAx(i,j,5)
         aFy(i,j,1) = aAy(i,j,1)
         aFy(i,j,2) = aAy(i,j,2)+aDy(i,j,1)
         aFy(i,j,3) = aAy(i,j,3)+aDy(i,j,2)
         aFy(i,j,4) = aAy(i,j,4)+aDy(i,j,3)
         aFy(i,j,5) = aAy(i,j,5)
      enddo
   enddo

   do j = 1,NY
      do i = 1,NX
         rhs(i,j) =-(aFx(i,j,1)*u(i-2,j  )+aFx(i,j,2)*u(i-1,j  )&
                    +aFy(i,j,1)*u(i  ,j-2)+aFy(i,j,2)*u(i  ,j-1)&
                    +(aFx(i,j,3)+aFy(i,j,3))*u(i,j)&
                    +aFx(i,j,4)*u(i+1,j  )+aFx(i,j,5)*u(i+2,j  )&
                    +aFy(i,j,4)*u(i  ,j+1)+aFy(i,j,5)*u(i  ,j+2)&
                    )+dt*S(i,j)
      enddo
   enddo

   !Matrix inversion calculation using LU-SGS method
   call LUSGS2D(NX,NY,aFx,aFy,rhs,delta)

   return
endsubroutine



!--------------- Coefficient for 2D x-advection term -------------------
!          Input:  NX,NY,dt,dx,c(-1:NX+2,-1:NY+2),u(-1:NX+2,-1:NY+2)
!          Output: a(1:NX,1:NY,5)
!    dt(d(Cu)/dx) = a(i,j,1)u(i-2,j)+a(i,j,2)u(i-1,j)+a(i,j,3)u(i,j)
!                  +a(i,j,4)u(i+1,j)+a(i,j,5)u(i+2,j)
!-----------------------------------------------------------------------

subroutine AdvectionCoefficient2D_X(NX,NY,dt,dx,c,u,a)

   implicit none
   integer :: i,j
   integer,intent(in)                                     :: NX,NY
   double precision,intent(in)                            :: dt,dx
   double precision,dimension(-1:NX+2,-1:NY+2),intent(in) :: c
   double precision,dimension(-1:NX+2,-1:NY+2),intent(in) :: u
   double precision,dimension(1:NX,1:NY,5),intent(out)    :: a
   double precision,dimension(0:NX+1,1:NY)                :: sc1,sc2,su1,su2
   double precision,dimension(1:NX+1,1:NY)                :: cL,cR,machL,machR
   double precision,dimension(1:NX+1,1:NY)                :: cLp,cRm
   double precision :: db,df,CNORM

   !Determine switch for c
   do j = 1,NY
      do i = 0,NX+1
         db = c(i,j)-c(i-1,j); df = c(i+1,j)-c(i,j)
         if(db*df.lt.0.0d0) then
            sc1(i,j) = 0.0d0; sc2(i,j) = 0.0d0
         else if(dabs(db).le.dabs(df)) then
            sc1(i,j) = 1.0d0; sc2(i,j) = 0.0d0
         else if(dabs(db).gt.dabs(df)) then
            sc1(i,j) = 0.0d0; sc2(i,j) = 1.0d0
         endif
      enddo
   enddo
   !Determine switch for u
   do j = 1,NY
      do i = 0,NX+1
         db = u(i,j)-u(i-1,j); df = u(i+1,j)-u(i,j)
         if(db*df.lt.0.0d0) then
            su1(i,j) = 0.0d0; su2(i,j) = 0.0d0
         else if(dabs(db).le.dabs(df)) then
            su1(i,j) = 1.0d0; su2(i,j) = 0.0d0
         else if(dabs(db).gt.dabs(df)) then
            su1(i,j) = 0.0d0; su2(i,j) = 1.0d0
         endif
      enddo
   enddo
   !Calculate interpolated c
   do j = 1,NY
      do i = 1,NX+1
         cL(i,j) = c(i-1,j)+0.5d0*(sc1(i-1,j)*(c(i-1,j)-c(i-2,j))+sc2(i-1,j)*(c(i  ,j)-c(i-1,j)))
         cR(i,j) = c(i  ,j)-0.5d0*(sc1(i  ,j)*(c(i  ,j)-c(i-1,j))+sc2(i  ,j)*(c(i+1,j)-c(i  ,j)))
      enddo
   enddo
   CNORM = 1.0d0
   !Calculated upwinded c
   do j = 1,NY
      do i = 1,NX+1
         machL(i,j) = cL(i,j)/CNORM
         machR(i,j) = cR(i,j)/CNORM
         !if(dabs(machL(i,j)).ge.1.0d0) then
         !   cLp(i,j) = CNORM*0.50d0*(machL(i,j)+dabs(machL(i,j)))
         !else
         !   cLp(i,j) = CNORM*0.25d0*(machL(i,j)+1.0d0)**2.0d0
         !endif
         !if(dabs(machR(i,j)).ge.1.0d0) then
         !   cRm(i,j) = CNORM*0.50d0*(machR(i,j)-dabs(machR(i,j)))
         !else
         !   cRm(i,j) =-CNORM*0.25d0*(machR(i,j)-1.0d0)**2.0d0
         !endif
         cLp(i,j) = 0.5d0*(cL(i,j)+dabs(cL(i,j)))
         cRm(i,j) = 0.5d0*(cR(i,j)-dabs(cR(i,j)))
      enddo
   enddo
   !Calculate coefficients a
   do j = 1,NY
      do i = 1,NX
         a(i,j,1) = dt/2.0d0/dx*( cLp(i  ,j)*su1(i-1,j))
         a(i,j,2) = dt/2.0d0/dx*(-cLp(i+1,j)*su1(i  ,j)&
                                 -cRm(i  ,j)*su1(i  ,j)&
                                 -cLp(i  ,j)*(2.0d0+su1(i-1,j)-su2(i-1,j)))
         a(i,j,3) = dt/2.0d0/dx*( cLp(i+1,j)*(2.0d0+su1(i  ,j)-su2(i  ,j))&
                                 +cRm(i+1,j)*su1(i+1,j)&
                                 -cLp(i  ,j)*su2(i-1,j)&
                                 -cRm(i  ,j)*(2.0d0-su1(i  ,j)+su2(i  ,j)))
         a(i,j,4) = dt/2.0d0/dx*( cLp(i+1,j)*su2(i  ,j)&
                                 +cRm(i+1,j)*(2.0d0-su1(i+1,j)+su2(i+1,j))&
                                 +cRm(i  ,j)*su2(i  ,j))
         a(i,j,5) = dt/2.0d0/dx*(-cRm(i+1,j)*su2(i+1,j))
      enddo
   enddo

   return
endsubroutine


!--------------- Coefficient for 2D y-advection term -------------------
!       Input:  NX,NY,dt,dy,c(-1:NX+2,-1:NY+2),u(-1:NX+2,-1:NY+2)
!       Output: a(1:NX,1:NY,5)
!    dt(d(Cu)/dy) = a(i,j,1)u(i,j-2)+a(i,j,2)u(i,j-1)+a(i,j,3)u(i,j)
!                  +a(i,j,4)u(i,j+1)+a(i,j,5)u(i,j+2)
!-----------------------------------------------------------------------

subroutine AdvectionCoefficient2D_Y(NX,NY,dt,dy,c,u,a)

   implicit none
   integer :: i,j
   integer,intent(in)                                     :: NX,NY
   double precision,intent(in)                            :: dt,dy
   double precision,dimension(-1:NX+2,-1:NY+2),intent(in) :: c
   double precision,dimension(-1:NX+2,-1:NY+2),intent(in) :: u
   double precision,dimension(1:NX,1:NY,5),intent(out)    :: a
   double precision,dimension(1:NX,0:NY+1)                :: sc1,sc2,su1,su2
   double precision,dimension(1:NX,1:NY+1)                :: cL,cR,machL,machR
   double precision,dimension(1:NX,1:NY+1)                :: cLp,cRm
   double precision :: db,df,CNORM

   !Determine switch for c
   do j = 0,NY+1
      do i = 1,NX
         db = c(i,j)-c(i,j-1); df = c(i,j+1)-c(i,j)
         if(db*df.lt.0.0d0) then
            sc1(i,j) = 0.0d0; sc2(i,j) = 0.0d0
         else if(dabs(db).le.dabs(df)) then
            sc1(i,j) = 1.0d0; sc2(i,j) = 0.0d0
         else if(dabs(db).gt.dabs(df)) then
            sc1(i,j) = 0.0d0; sc2(i,j) = 1.0d0
         endif
      enddo
   enddo
   !Determine switch for u
   do j = 0,NY+1
      do i = 1,NX
         db = u(i,j)-u(i,j-1); df = u(i,j+1)-u(i,j)
         if(db*df.lt.0.0d0) then
            su1(i,j) = 0.0d0; su2(i,j) = 0.0d0
         else if(dabs(db).le.dabs(df)) then
            su1(i,j) = 1.0d0; su2(i,j) = 0.0d0
         else if(dabs(db).gt.dabs(df)) then
            su1(i,j) = 0.0d0; su2(i,j) = 1.0d0
         endif
      enddo
   enddo
   !Calculate interpolated c
   do j = 1,NY+1
      do i = 1,NX
         cL(i,j) = c(i,j-1)+0.5d0*(sc1(i,j-1)*(c(i,j-1)-c(i,j-2))+sc2(i,j-1)*(c(i,j  )-c(i,j-1)))
         cR(i,j) = c(i,j  )-0.5d0*(sc1(i,j  )*(c(i,j  )-c(i,j-1))+sc2(i,j  )*(c(i,j+1)-c(i,j  )))
      enddo
   enddo
   CNORM = 1.0d0
   !Calculated upwinded c
   do j = 1,NY+1
      do i = 1,NX
         machL(i,j) = cL(i,j)/CNORM
         machR(i,j) = cR(i,j)/CNORM
         !if(dabs(machL(i,j)).ge.1.0d0) then
         !   cLp(i,j) = CNORM*0.50d0*(machL(i,j)+dabs(machL(i,j)))
         !else
         !   cLp(i,j) = CNORM*0.25d0*(machL(i,j)+1.0d0)**2.0d0
         !endif
         !if(dabs(machR(i,j)).ge.1.0d0) then
         !   cRm(i,j) = CNORM*0.50d0*(machR(i,j)-dabs(machR(i,j)))
         !else
         !   cRm(i,j) =-CNORM*0.25d0*(machR(i,j)-1.0d0)**2.0d0
         !endif
         cLp(i,j) = 0.5d0*(cL(i,j)+dabs(cL(i,j)))
         cRm(i,j) = 0.5d0*(cR(i,j)-dabs(cR(i,j)))
      enddo
   enddo
   !Calculate coefficients a
   do j = 1,NY
      do i = 1,NX
         a(i,j,1) = dt/2.0d0/dy*( cLp(i,j  )*su1(i,j-1))
         a(i,j,2) = dt/2.0d0/dy*(-cLp(i,j+1)*su1(i,j  )&
                                 -cRm(i,j  )*su1(i,j  )&
                                 -cLp(i,j  )*(2.0d0+su1(i,j-1)-su2(i,j-1)))
         a(i,j,3) = dt/2.0d0/dy*( cLp(i,j+1)*(2.0d0+su1(i,j  )-su2(i,j  ))&
                                 +cRm(i,j+1)*su1(i,j+1)&
                                 -cLp(i,j  )*su2(i,j-1)&
                                 -cRm(i,j  )*(2.0d0-su1(i,j  )+su2(i,j  )))
         a(i,j,4) = dt/2.0d0/dy*( cLp(i,j+1)*su2(i,j  )&
                                 +cRm(i,j+1)*(2.0d0-su1(i,j+1)+su2(i,j+1))&
                                 +cRm(i,j  )*su2(i,j  ))
         a(i,j,5) = dt/2.0d0/dy*(-cRm(i,j+1)*su2(i,j+1))
      enddo
   enddo

   return
endsubroutine


!------------- Coefficient for 2D x-diffusion term ---------------------
!               Input : NX,NY,dt,dx,D(0:NX+1,0:NY+1)
!               Output: a(1:NX,1:NY,3)
!    dt(d/dx(-Dd/dx(u)))
!         = a(i,j,1)u(i-1,j)+a(i,j,2)u(i,j)+a(i,j,3)u(i+1,j)
!
!-----------------------------------------------------------------------

subroutine DiffusionCoefficient2D_X(NX,NY,dt,dx,D,a)

   implicit none
   integer :: i,j
   integer,intent(in)                                   :: NX,NY
   double precision,intent(in)                          :: dt,dx
   double precision,dimension(0:NX+1,0:NY+1),intent(in) :: D
   double precision,dimension(1:NX,1:NY,3),intent(out)  :: a

   !Calculate coefficients a
   do j = 1,NY
      do i = 1,NX
         a(i,j,1) =-dt/2.0d0/dx/dx*(D(i  ,j)+D(i-1,j))
         a(i,j,2) = dt/2.0d0/dx/dx*(D(i+1,j)+2.0d0*D(i  ,j)+D(i-1,j))
         a(i,j,3) =-dt/2.0d0/dx/dx*(D(i+1,j)+D(i  ,j))
      enddo
   enddo

   return
endsubroutine


!------------- Coefficient for 2D y-diffusion term ---------------------
!               Input : NX,NY,dt,dy,D(0:NX+1,0:NY+1)
!               Output: a(1:NX,1:NY,3)
!    dt(d/dy(-Dd/dy(u)))
!         = a(i,j,1)u(i,j-1)+a(i,j,2)u(i,j)+a(i,j,3)u(i,j+1)
!
!-----------------------------------------------------------------------

subroutine DiffusionCoefficient2D_Y(NX,NY,dt,dy,D,a)

   implicit none
   integer :: i,j
   integer,intent(in)                                   :: NX,NY
   double precision,intent(in)                          :: dt,dy
   double precision,dimension(0:NX+1,0:NY+1),intent(in) :: D
   double precision,dimension(1:NX,1:NY,3),intent(out)  :: a

   !Calculate coefficients a
   do j = 1,NY
      do i = 1,NX
         a(i,j,1) =-dt/2.0d0/dy/dy*(D(i,j  )+D(i,j-1))
         a(i,j,2) = dt/2.0d0/dy/dy*(D(i,j+1)+2.0d0*D(i,j  )+D(i,j-1))
         a(i,j,3) =-dt/2.0d0/dy/dy*(D(i,j+1)+D(i,j  ))
      enddo
   enddo

   return
endsubroutine


!--------------- LU-SGS scheme for scalar problem ----------------------
!      Input:  NX,NY,cx(1:NX,1:NY,5),cy(1:NX,1:NY,5),rhs(1:NX,1:NY)
!      Output: delta(1:NX,1:NY)
!           cx(i,j,1)*delta(i-2,j  )+cx(i,j,2)*delta(i-1,j  )
!          +cy(i,j,1)*delta(i  ,j-2)+cy(i,j,2)*delta(i  ,j-1)
!               +(1+cx(i,j,3)+cy(i,j,3))*delta(i  ,j  )
!          +cx(i,j,4)*delta(i+1,j  )+cx(i,j,5)*delta(i+2,j  )
!          +cy(i,j,4)*delta(i  ,j+1)+cy(i,j,5)*delta(i  ,j+2)
!                            = rhs(i,j)
!-----------------------------------------------------------------------

subroutine LUSGS2D(NX,NY,cx,cy,rhs,delta)

   implicit none
   integer :: i,j,k,l
   integer,intent(in)                                 :: NX,NY          !Number of cells
   double precision,dimension(1:NX,1:NY,5),intent(in) :: cx             !x-coefficients in LHS
   double precision,dimension(1:NX,1:NY,5),intent(in) :: cy             !y-coefficients in LHS
   double precision,dimension(1:NX,1:NY),intent(in)   :: rhs            !RHS
   double precision,dimension(1:NX,1:NY),intent(out)  :: delta          !Delta of variable
   double precision,dimension(-1:NX+2,-1:NY+2)        :: dd1,dd2        !Delta in each sweep
   double precision xx,yy

   !Initialization of dd1,dd2
   do j =-1,NY+2
      do i =-1,NX+2
         dd1(i,j) = 0.0d0
         dd2(i,j) = 0.0d0
      enddo
   enddo
   !Sweep1 (1,1) to (NX,NY)
   do j = 1,NY
      do i = 1,NX
         xx = 1.0d0+cx(i,j,3)+cy(i,j,3)
         if(xx.eq.0.0d0) write(*,*) 'Warning in LUSGS...',i,j
         yy = cx(i,j,1)*dd1(i-2,j  )+cx(i,j,2)*dd1(i-1,j  )&
             +cy(i,j,1)*dd1(i  ,j-2)+cy(i,j,2)*dd1(i  ,j-1)
         dd1(i,j) = 1.0d0/xx*(rhs(i,j)-yy)
      enddo
   enddo
   !Sweep2 (NX,NY) to (1,1)
   do l = 1,NY
      do k = 1,NX
         i = NX+1-k
         j = NY+1-l
         xx = 1.0d0+cx(i,j,3)+cy(i,j,3)
         yy = cx(i,j,4)*dd2(i+1,j  )+cx(i,j,5)*dd2(i+2,j  )&
             +cy(i,j,4)*dd2(i  ,j+1)+cy(i,j,5)*dd2(i  ,j+2)
         dd2(i,j) = dd1(i,j)-yy/xx
      enddo
   enddo
   !Copy
   do j = 1,NY
      do i = 1,NX
         delta(i,j) = dd2(i,j)
      enddo
   enddo

   return
endsubroutine






!-----------------------------------------------------------------------

subroutine Mobility(muex,muey)

   use parameters_mod
   implicit none
   integer :: i,j
   double precision,dimension(1:III,1:JJJ),intent(out) :: muex,muey
   double precision :: xx,yy

   !Calculation of magnetic tensor
   !$omp parallel default(none),private(i,j,xx,yy),shared(muex,muey)
   !$omp do
   do j = 1,JJJ
      do i = 1,III
         xx = DXLND*(dble(i)-0.5d0)
         yy = OMEGAE*dexp(-((xx-XLND)/XLND)**2.0d0)
         muex(i,j) = 1.0d0/yy**2.0d0
         muey(i,j) = yy   /yy**2.0d0
      enddo
   enddo
   !$omp end do
   !$omp end parallel

   return
end subroutine



