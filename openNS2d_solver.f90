	subroutine Euler2d_solver(rho,u,v,p,xx,yy) !!

	include 'openNS3d.h'

	integer :: i,j,iter,ns

!	real(kind=OCFD_REAL_KIND) :: T

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,rho
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
!	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) ::

!    T=2.d0
    ns=ceiling(T/dt)
	write(*,*) ns
	do iter=1,ns
		!call TVD_RK3(rho,u,v,p)
		call RK4(rho,u,v,p)

		if (my_id==0) then
			write(*,*) iter
			write(*,*) iter*dt
		endif

	!	if (iter==20000) then
	!		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!		write(filename,"('Taylor',F5.2,'.plt')") iter*dt
	!		call write_data(rho,u,v,p,xx,yy,filename)
	!	endif

		if (mod(iter,Istep_save)==0) then
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			write(filename,"('Taylor',F5.2,'.plt')") iter*dt
			call write_data(rho,u,v,p,xx,yy,filename)
		endif

	end do




	end subroutine








!!!==========================================================================
!!!
!!!==========================================================================
subroutine computeR(rho,u,v,p,R_st,R_nd,R_rd,R_th)
	include 'openNS3d.h'


	integer :: i,j,k

!	real(kind=OCFD_REAL_KIND) ::

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,rho
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st,R_nd,R_rd,R_th
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Epos1,Epos2,Epos3,Epos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Eneg1,Eneg2,Eneg3,Eneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Fpos1,Fpos2,Fpos3,Fpos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Fneg1,Fneg2,Fneg3,Fneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: d2f_Epos1,d2f_Epos2,d2f_Epos3,d2f_Epos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: d2f_Eneg1,d2f_Eneg2,d2f_Eneg3,d2f_Eneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: d2f_Fpos1,d2f_Fpos2,d2f_Fpos3,d2f_Fpos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: d2f_Fneg1,d2f_Fneg2,d2f_Fneg3,d2f_Fneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: df_Epos1,df_Epos2,df_Epos3,df_Epos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: df_Eneg1,df_Eneg2,df_Eneg3,df_Eneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: df_Fpos1,df_Fpos2,df_Fpos3,df_Fpos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: df_Fneg1,df_Fneg2,df_Fneg3,df_Fneg4





!-----------------------------------------------------------------------------------
!CCU45
!-----------------------------------------------------------------------------------
	if (NUM_METHOD_CONV==45) then


	if (npx0 .ne. 1) then
		call check_x2d(rho)
		call check_x2d(u)
		call check_x2d(v)
		call check_x2d(p)

		call VanLeerE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,&
		Eneg1,Eneg2,Eneg3,Eneg4)

		do j=1,ny
			call Du2Dx_PADE4(Epos1(:,j),d2f_Epos1(:,j),Iperiodic_X)
			call Du2Dx_PADE4(Eneg1(:,j),d2f_Eneg1(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind(Epos1(:,j),d2f_Epos1(:,j),df_Epos1(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind(Eneg1(:,j),d2f_Eneg1(:,j),df_Eneg1(:,j),Iperiodic_X)

			call Du2Dx_PADE4(Epos2(:,j),d2f_Epos2(:,j),Iperiodic_X)
			call Du2Dx_PADE4(Eneg2(:,j),d2f_Eneg2(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind(Epos2(:,j),d2f_Epos2(:,j),df_Epos2(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind(Eneg2(:,j),d2f_Eneg2(:,j),df_Eneg2(:,j),Iperiodic_X)

			call Du2Dx_PADE4(Epos3(:,j),d2f_Epos3(:,j),Iperiodic_X)
			call Du2Dx_PADE4(Eneg3(:,j),d2f_Eneg3(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind(Epos3(:,j),d2f_Epos3(:,j),df_Epos3(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind(Eneg3(:,j),d2f_Eneg3(:,j),df_Eneg3(:,j),Iperiodic_X)

			call Du2Dx_PADE4(Epos4(:,j),d2f_Epos4(:,j),Iperiodic_X)
			call Du2Dx_PADE4(Eneg4(:,j),d2f_Eneg4(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind(Epos4(:,j),d2f_Epos4(:,j),df_Epos4(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind(Eneg4(:,j),d2f_Eneg4(:,j),df_Eneg4(:,j),Iperiodic_X)
		end do

	elseif (npx0==1) then

		call VanLeerE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,&
		Eneg1,Eneg2,Eneg3,Eneg4)

		do j=1,ny
			call Du2Dx_PADE4_serial(Epos1(:,j),d2f_Epos1(:,j),Iperiodic_X)
			call Du2Dx_PADE4_serial(Eneg1(:,j),d2f_Eneg1(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind_serial(Epos1(:,j),d2f_Epos1(:,j),df_Epos1(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind_serial(Eneg1(:,j),d2f_Eneg1(:,j),df_Eneg1(:,j),Iperiodic_X)

			call Du2Dx_PADE4_serial(Epos2(:,j),d2f_Epos2(:,j),Iperiodic_X)
			call Du2Dx_PADE4_serial(Eneg2(:,j),d2f_Eneg2(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind_serial(Epos2(:,j),d2f_Epos2(:,j),df_Epos2(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind_serial(Eneg2(:,j),d2f_Eneg2(:,j),df_Eneg2(:,j),Iperiodic_X)

			call Du2Dx_PADE4_serial(Epos3(:,j),d2f_Epos3(:,j),Iperiodic_X)
			call Du2Dx_PADE4_serial(Eneg3(:,j),d2f_Eneg3(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind_serial(Epos3(:,j),d2f_Epos3(:,j),df_Epos3(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind_serial(Eneg3(:,j),d2f_Eneg3(:,j),df_Eneg3(:,j),Iperiodic_X)

			call Du2Dx_PADE4_serial(Epos4(:,j),d2f_Epos4(:,j),Iperiodic_X)
			call Du2Dx_PADE4_serial(Eneg4(:,j),d2f_Eneg4(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind_serial(Epos4(:,j),d2f_Epos4(:,j),df_Epos4(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind_serial(Eneg4(:,j),d2f_Eneg4(:,j),df_Eneg4(:,j),Iperiodic_X)
		end do
	endif

	if (npy0 .ne. 1) then
		call check_y2d(rho)
		call check_y2d(u)
		call check_y2d(v)
		call check_y2d(p)

		call VanLeerF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,&
		Fneg1,Fneg2,Fneg3,Fneg4)

		do i=1,nx
			call Du2Dy_PADE4(Fpos1(i,:),d2f_Fpos1(i,:),Iperiodic_Y)
			call Du2Dy_PADE4(Fneg1(i,:),d2f_Fneg1(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind(Fpos1(i,:),d2f_Fpos1(i,:),df_Fpos1(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind(Fneg1(i,:),d2f_Fneg1(i,:),df_Fneg1(i,:),Iperiodic_Y)

			call Du2Dy_PADE4(Fpos2(i,:),d2f_Fpos2(i,:),Iperiodic_Y)
			call Du2Dy_PADE4(Fneg2(i,:),d2f_Fneg2(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind(Fpos2(i,:),d2f_Fpos2(i,:),df_Fpos2(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind(Fneg2(i,:),d2f_Fneg2(i,:),df_Fneg2(i,:),Iperiodic_Y)

			call Du2Dy_PADE4(Fpos3(i,:),d2f_Fpos3(i,:),Iperiodic_Y)
			call Du2Dy_PADE4(Fneg3(i,:),d2f_Fneg3(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind(Fpos3(i,:),d2f_Fpos3(i,:),df_Fpos3(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind(Fneg3(i,:),d2f_Fneg3(i,:),df_Fneg3(i,:),Iperiodic_Y)

			call Du2Dy_PADE4(Fpos4(i,:),d2f_Fpos4(i,:),Iperiodic_Y)
			call Du2Dy_PADE4(Fneg4(i,:),d2f_Fneg4(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind(Fpos4(i,:),d2f_Fpos4(i,:),df_Fpos4(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind(Fneg4(i,:),d2f_Fneg4(i,:),df_Fneg4(i,:),Iperiodic_Y)
		end do

	elseif (npy0==1) then

		call VanLeerF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,&
		Fneg1,Fneg2,Fneg3,Fneg4)

		do i=1,nx
			call Du2Dy_PADE4_serial(Fpos1(i,:),d2f_Fpos1(i,:),Iperiodic_Y)
			call Du2Dy_PADE4_serial(Fneg1(i,:),d2f_Fneg1(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind_serial(Fpos1(i,:),d2f_Fpos1(i,:),df_Fpos1(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind_serial(Fneg1(i,:),d2f_Fneg1(i,:),df_Fneg1(i,:),Iperiodic_Y)

			call Du2Dy_PADE4_serial(Fpos2(i,:),d2f_Fpos2(i,:),Iperiodic_Y)
			call Du2Dy_PADE4_serial(Fneg2(i,:),d2f_Fneg2(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind_serial(Fpos2(i,:),d2f_Fpos2(i,:),df_Fpos2(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind_serial(Fneg2(i,:),d2f_Fneg2(i,:),df_Fneg2(i,:),Iperiodic_Y)

			call Du2Dy_PADE4_serial(Fpos3(i,:),d2f_Fpos3(i,:),Iperiodic_Y)
			call Du2Dy_PADE4_serial(Fneg3(i,:),d2f_Fneg3(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind_serial(Fpos3(i,:),d2f_Fpos3(i,:),df_Fpos3(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind_serial(Fneg3(i,:),d2f_Fneg3(i,:),df_Fneg3(i,:),Iperiodic_Y)

			call Du2Dy_PADE4_serial(Fpos4(i,:),d2f_Fpos4(i,:),Iperiodic_Y)
			call Du2Dy_PADE4_serial(Fneg4(i,:),d2f_Fneg4(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind_serial(Fpos4(i,:),d2f_Fpos4(i,:),df_Fpos4(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind_serial(Fneg4(i,:),d2f_Fneg4(i,:),df_Fneg4(i,:),Iperiodic_Y)
		end do
	endif
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
	elseif (NUM_METHOD_CONV==5) then
!==================WENO=========================================

	if (npx0 .ne. 1) then
			call check_x2d(rho)
			call check_x2d(u)
			call check_x2d(v)
			call check_x2d(p)



		!call VanLeerE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,&
		!Eneg1,Eneg2,Eneg3,Eneg4)
		!call StegerE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,&
		!Eneg1,Eneg2,Eneg3,Eneg4)
		call LaxE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,&
		Eneg1,Eneg2,Eneg3,Eneg4)


		do j=1,ny
			call du1_weno5(Epos1(:,j),df_Epos1(:,j),nx,hx)
			call du2_weno5(Eneg1(:,j),df_Eneg1(:,j),nx,hx)

			call du1_weno5(Epos2(:,j),df_Epos2(:,j),nx,hx)
			call du2_weno5(Eneg2(:,j),df_Eneg2(:,j),nx,hx)

			call du1_weno5(Epos3(:,j),df_Epos3(:,j),nx,hx)
			call du2_weno5(Eneg3(:,j),df_Eneg3(:,j),nx,hx)

			call du1_weno5(Epos4(:,j),df_Epos4(:,j),nx,hx)
			call du2_weno5(Eneg4(:,j),df_Eneg4(:,j),nx,hx)
		end do

			!call OCFD_DFX_REFBOUND_CHECK_2d(Epos1,df_Epos1,NUM_METHOD_OTH,1)
			!call OCFD_DFX_REFBOUND_CHECK_2d(Eneg1,df_Eneg1,NUM_METHOD_OTH,2)
			!call OCFD_DFX_REFBOUND_CHECK_2d(Epos2,df_Epos2,NUM_METHOD_OTH,1)
			!call OCFD_DFX_REFBOUND_CHECK_2d(Eneg2,df_Eneg2,NUM_METHOD_OTH,2)
			!call OCFD_DFX_REFBOUND_CHECK_2d(Epos3,df_Epos3,NUM_METHOD_OTH,1)
			!call OCFD_DFX_REFBOUND_CHECK_2d(Eneg3,df_Eneg3,NUM_METHOD_OTH,2)
			!call OCFD_DFX_REFBOUND_CHECK_2d(Epos4,df_Epos4,NUM_METHOD_OTH,1)
			!call OCFD_DFX_REFBOUND_CHECK_2d(Eneg4,df_Eneg4,NUM_METHOD_OTH,2)

	elseif (npx0==1) then

		if (Iperiodic_X .eq. 1) then
		rho(-3:0,:)=rho(nx-4:nx-1,:)
		u(-3:0,:)=u(nx-4:nx-1,:)
		v(-3:0,:)=v(nx-4:nx-1,:)
		p(-3:0,:)=p(nx-4:nx-1,:)

		rho(nx+1:nx+4,:)=rho(2:5,:)
		u(nx+1:nx+4,:)=u(2:5,:)
		v(nx+1:nx+4,:)=v(2:5,:)
		p(nx+1:nx+4,:)=p(2:5,:)
		endif

		!call VanLeerE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,&
		!Eneg1,Eneg2,Eneg3,Eneg4)
		!call StegerE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,&
		!Eneg1,Eneg2,Eneg3,Eneg4)
		call LaxE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,&
		Eneg1,Eneg2,Eneg3,Eneg4)


		do j=1,ny
			call du1_weno5(Epos1(:,j),df_Epos1(:,j),nx,hx)
			call du2_weno5(Eneg1(:,j),df_Eneg1(:,j),nx,hx)

			call du1_weno5(Epos2(:,j),df_Epos2(:,j),nx,hx)
			call du2_weno5(Eneg2(:,j),df_Eneg2(:,j),nx,hx)

			call du1_weno5(Epos3(:,j),df_Epos3(:,j),nx,hx)
			call du2_weno5(Eneg3(:,j),df_Eneg3(:,j),nx,hx)

			call du1_weno5(Epos4(:,j),df_Epos4(:,j),nx,hx)
			call du2_weno5(Eneg4(:,j),df_Eneg4(:,j),nx,hx)
		end do

		!	call OCFD_DFX_REFBOUND_CHECK_2d(Epos1,df_Epos1,NUM_METHOD_OTH,1)
		!	call OCFD_DFX_REFBOUND_CHECK_2d(Eneg1,df_Eneg1,NUM_METHOD_OTH,2)
		!	call OCFD_DFX_REFBOUND_CHECK_2d(Epos2,df_Epos2,NUM_METHOD_OTH,1)
		!	call OCFD_DFX_REFBOUND_CHECK_2d(Eneg2,df_Eneg2,NUM_METHOD_OTH,2)
		!	call OCFD_DFX_REFBOUND_CHECK_2d(Epos3,df_Epos3,NUM_METHOD_OTH,1)
		!	call OCFD_DFX_REFBOUND_CHECK_2d(Eneg3,df_Eneg3,NUM_METHOD_OTH,2)
		!	call OCFD_DFX_REFBOUND_CHECK_2d(Epos4,df_Epos4,NUM_METHOD_OTH,1)
		!	call OCFD_DFX_REFBOUND_CHECK_2d(Eneg4,df_Eneg4,NUM_METHOD_OTH,2)
	endif


	if (npy0 .ne. 1) then
		call check_y2d(rho)
		call check_y2d(u)
		call check_y2d(v)
		call check_y2d(p)

		!call VanLeerF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,&
		!Fneg1,Fneg2,Fneg3,Fneg4)
		!call StegerF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,&
		!Fneg1,Fneg2,Fneg3,Fneg4)
		call LaxF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,&
		Fneg1,Fneg2,Fneg3,Fneg4)


		do i=1,nx
			call du1_weno5(Fpos1(i,:),df_Fpos1(i,:),ny,hy)
			call du2_weno5(Fneg1(i,:),df_Fneg1(i,:),ny,hy)

			call du1_weno5(Fpos2(i,:),df_Fpos2(i,:),ny,hy)
			call du2_weno5(Fneg2(i,:),df_Fneg2(i,:),ny,hy)

			call du1_weno5(Fpos3(i,:),df_Fpos3(i,:),ny,hy)
			call du2_weno5(Fneg3(i,:),df_Fneg3(i,:),ny,hy)

			call du1_weno5(Fpos4(i,:),df_Fpos4(i,:),ny,hy)
			call du2_weno5(Fneg4(i,:),df_Fneg4(i,:),ny,hy)
		end do

!			call OCFD_DFY_REFBOUND_CHECK_2d(Fpos1,df_Fpos1,NUM_METHOD_OTH,1)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fneg1,df_Fneg1,NUM_METHOD_OTH,2)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fpos2,df_Fpos2,NUM_METHOD_OTH,1)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fneg2,df_Fneg2,NUM_METHOD_OTH,2)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fpos3,df_Fpos3,NUM_METHOD_OTH,1)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fneg3,df_Fneg3,NUM_METHOD_OTH,2)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fpos4,df_Fpos4,NUM_METHOD_OTH,1)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fneg4,df_Fneg4,NUM_METHOD_OTH,2)

!			call OCFD_DFY_BOUND_CHECK_2d(Fpos1,df_Fpos1,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fneg1,df_Fneg1,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fpos2,df_Fpos2,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fneg2,df_Fneg2,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fpos3,df_Fpos3,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fneg3,df_Fneg3,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fpos4,df_Fpos4,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fneg4,df_Fneg4,NUM_METHOD_OTH)

	elseif (npy0==1) then

!		if (Iperiodic_Y .eq. 1) then
!		Fpos1(:,-3:0)=Fpos1(:,ny-4:ny-1)
!		Fpos2(:,-3:0)=Fpos2(:,ny-4:ny-1)
!		Fpos3(:,-3:0)=Fpos3(:,ny-4:ny-1)
!		Fpos4(:,-3:0)=Fpos4(:,ny-4:ny-1)
!
!		Fpos1(:,ny+1:ny+4)=Fpos1(:,2:5)
!		Fpos2(:,ny+1:ny+4)=Fpos2(:,2:5)
!		Fpos3(:,ny+1:ny+4)=Fpos3(:,2:5)
!		Fpos4(:,ny+1:ny+4)=Fpos4(:,2:5)
!
!
!		Fneg1(:,-3:0)=Fneg1(:,ny-4:ny-1)
!		Fneg2(:,-3:0)=Fneg2(:,ny-4:ny-1)
!		Fneg3(:,-3:0)=Fneg3(:,ny-4:ny-1)
!		Fneg4(:,-3:0)=Fneg4(:,ny-4:ny-1)
!
!		Fneg1(:,ny+1:ny+4)=Fneg1(:,2:5)
!		Fneg2(:,ny+1:ny+4)=Fneg2(:,2:5)
!		Fneg3(:,ny+1:ny+4)=Fneg3(:,2:5)
!		Fneg4(:,ny+1:ny+4)=Fneg4(:,2:5)
!		endif

		!call VanLeerF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,&
		!Fneg1,Fneg2,Fneg3,Fneg4)
		!call StegerF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,&
		!Fneg1,Fneg2,Fneg3,Fneg4)
		call LaxF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,&
		Fneg1,Fneg2,Fneg3,Fneg4)

		do i=1,nx
			call du1_weno5(Fpos1(i,:),df_Fpos1(i,:),ny,hy)
			call du2_weno5(Fneg1(i,:),df_Fneg1(i,:),ny,hy)

			call du1_weno5(Fpos2(i,:),df_Fpos2(i,:),ny,hy)
			call du2_weno5(Fneg2(i,:),df_Fneg2(i,:),ny,hy)

			call du1_weno5(Fpos3(i,:),df_Fpos3(i,:),ny,hy)
			call du2_weno5(Fneg3(i,:),df_Fneg3(i,:),ny,hy)

			call du1_weno5(Fpos4(i,:),df_Fpos4(i,:),ny,hy)
			call du2_weno5(Fneg4(i,:),df_Fneg4(i,:),ny,hy)
		end do

!			call OCFD_DFY_REFBOUND_CHECK_2d(Fpos1,df_Fpos1,NUM_METHOD_OTH,1)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fneg1,df_Fneg1,NUM_METHOD_OTH,2)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fpos2,df_Fpos2,NUM_METHOD_OTH,1)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fneg2,df_Fneg2,NUM_METHOD_OTH,2)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fpos3,df_Fpos3,NUM_METHOD_OTH,1)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fneg3,df_Fneg3,NUM_METHOD_OTH,2)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fpos4,df_Fpos4,NUM_METHOD_OTH,1)
!			call OCFD_DFY_REFBOUND_CHECK_2d(Fneg4,df_Fneg4,NUM_METHOD_OTH,2)


!			call OCFD_DFY_BOUND_CHECK_2d(Fpos1,df_Fpos1,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fneg1,df_Fneg1,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fpos2,df_Fpos2,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fneg2,df_Fneg2,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fpos3,df_Fpos3,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fneg3,df_Fneg3,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fpos4,df_Fpos4,NUM_METHOD_OTH)
!			call OCFD_DFY_BOUND_CHECK_2d(Fneg4,df_Fneg4,NUM_METHOD_OTH)

	endif

!================================================================
	else
	      print*, 'This Numerical Method is not supported'
		  stop
	endif


	R_st=-(df_Epos1+df_Eneg1)-(df_Fpos1+df_Fneg1)
	R_nd=-(df_Epos2+df_Eneg2)-(df_Fpos2+df_Fneg2)
	R_rd=-(df_Epos3+df_Eneg3)-(df_Fpos3+df_Fneg3)
	R_th=-(df_Epos4+df_Eneg4)-(df_Fpos4+df_Fneg4)

	do j=1,ny
		do i=1,nx
			R_rd(i,j)=R_rd(i,j)+rho(i,j)
			R_th(i,j)=R_th(i,j)+rho(i,j)*v(i,j)
		end do
	end do


end subroutine




subroutine RK4(rho,u,v,p)
	include 'openNS3d.h'

	integer :: i,j

	real(kind=OCFD_REAL_KIND),parameter::gama=5.d0/3.d0

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,rho

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st,R_nd,R_rd,R_th
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st1,R_nd1,R_rd1,R_th1
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st2,R_nd2,R_rd2,R_th2
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st3,R_nd3,R_rd3,R_th3

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st,Q_nd,Q_rd,Q_th
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st1,Q_nd1,Q_rd1,Q_th1
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st2,Q_nd2,Q_rd2,Q_th2
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st3,Q_nd3,Q_rd3,Q_th3

!	gama=1.4d0

	call computeR(rho,u,v,p,R_st,R_nd,R_rd,R_th)


	do j=1,ny
		do i=1,nx
			Q_st(i,j)=rho(i,j)
			Q_nd(i,j)=rho(i,j)*u(i,j)
			Q_rd(i,j)=rho(i,j)*v(i,j)
			Q_th(i,j)=p(i,j)/(gama-1.d0)+rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0)
		end do
	end do

	Q_st1=Q_st+(dt/2.d0)*R_st
	Q_nd1=Q_nd+(dt/2.d0)*R_nd
	Q_rd1=Q_rd+(dt/2.d0)*R_rd
	Q_th1=Q_th+(dt/2.d0)*R_th

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st1(i,j)
			u(i,j)=Q_nd1(i,j)/rho(i,j)
			v(i,j)=Q_rd1(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th1(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

	if (npy==0) then
		p(:,-3:1)=1.d0
		rho(:,-3:1)=2.d0
		u(:,-3:1)=0.d0
		v(:,-3:1)=0.d0
	endif

	if (npy==npy0-1) then
		p(:,ny:ny+4)=2.5d0
		rho(:,ny:ny+4)=1.d0
		u(:,ny:ny+4)=0.d0
		v(:,ny:ny+4)=0.d0
	endif


		if (Iperiodic_X .ne. 1) then
			if (npx .eq. 0) then
				do k=1,4
					rho(1-k,:)=rho(k,:)
					u(1-k,:)=-u(k,:)
					v(1-k,:)=v(k,:)
					p(1-k,:)=p(k,:)
				enddo
			endif

			if (npx .eq. npx0-1) then
				do k=nx+1,nx+4
					rho(k,:)=rho(2*nx+1-k,:)
					u(k,:)=-u(2*nx+1-k,:)
					v(k,:)=v(2*nx+1-k,:)
					p(k,:)=p(2*nx+1-k,:)
				enddo
			endif
		endif
!---------------------------------
!	if (Iperiodic_X==0) then
!
!		if (npx==0) then
!
!			u(1,:)=0.d0
!!			v(1,:)=0.d0
!
!		endif
!
!		if (npx==npx0-1) then
!
!			u(nx,:)=0.d0
!!			v(nx,:)=0.d0
!
!		endif
!
!	elseif (Iperiodic_X==1) then
!
!		call check_xbound(u)
!		call check_xbound(v)
!		call check_xbound(p)
!		call check_xbound(rho)
!
!	endif
!-----------------------------------

	call computeR(rho,u,v,p,R_st1,R_nd1,R_rd1,R_th1)


	Q_st2=Q_st+(dt/2.d0)*R_st1
	Q_nd2=Q_nd+(dt/2.d0)*R_nd1
	Q_rd2=Q_rd+(dt/2.d0)*R_rd1
	Q_th2=Q_th+(dt/2.d0)*R_th1

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st2(i,j)
			u(i,j)=Q_nd2(i,j)/rho(i,j)
			v(i,j)=Q_rd2(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th2(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

	if (npy==0) then
		p(:,-3:1)=1.d0
		rho(:,-3:1)=2.d0
		u(:,-3:1)=0.d0
		v(:,-3:1)=0.d0
	endif

	if (npy==npy0-1) then
		p(:,ny:ny+4)=2.5d0
		rho(:,ny:ny+4)=1.d0
		u(:,ny:ny+4)=0.d0
		v(:,ny:ny+4)=0.d0
	endif

	if (Iperiodic_X .ne. 1) then
			if (npx .eq. 0) then
				do k=1,4
					rho(1-k,:)=rho(k,:)
					u(1-k,:)=-u(k,:)
					v(1-k,:)=v(k,:)
					p(1-k,:)=p(k,:)
				enddo
			endif

			if (npx .eq. npx0-1) then
				do k=nx+1,nx+4
					rho(k,:)=rho(2*nx+1-k,:)
					u(k,:)=-u(2*nx+1-k,:)
					v(k,:)=v(2*nx+1-k,:)
					p(k,:)=p(2*nx+1-k,:)
				enddo
			endif
		endif
!---------------------------------
!	if (Iperiodic_X==0) then
!
!		if (npx==0) then
!
!			u(1,:)=0.d0
!!			v(1,:)=0.d0
!
!		endif
!
!		if (npx==npx0-1) then
!
!			u(nx,:)=0.d0
!!			v(nx,:)=0.d0
!
!		endif
!
!	elseif (Iperiodic_X==1) then
!
!		call check_xbound(u)
!		call check_xbound(v)
!		call check_xbound(p)
!		call check_xbound(rho)
!
!	endif
!-----------------------------------

	call computeR(rho,u,v,p,R_st2,R_nd2,R_rd2,R_th2)

	Q_st3=Q_st+dt*R_st2
	Q_nd3=Q_nd+dt*R_nd2
	Q_rd3=Q_rd+dt*R_rd2
	Q_th3=Q_th+dt*R_th2

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st3(i,j)
			u(i,j)=Q_nd3(i,j)/rho(i,j)
			v(i,j)=Q_rd3(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th3(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

	if (npy==0) then
		p(:,-3:1)=1.d0
		rho(:,-3:1)=2.d0
		u(:,-3:1)=0.d0
		v(:,-3:1)=0.d0
	endif

	if (npy==npy0-1) then
		p(:,ny:ny+4)=2.5d0
		rho(:,ny:ny+4)=1.d0
		u(:,ny:ny+4)=0.d0
		v(:,ny:ny+4)=0.d0
	endif

		if (Iperiodic_X .ne. 1) then
			if (npx .eq. 0) then
				do k=1,4
					rho(1-k,:)=rho(k,:)
					u(1-k,:)=-u(k,:)
					v(1-k,:)=v(k,:)
					p(1-k,:)=p(k,:)
				enddo
			endif

			if (npx .eq. npx0-1) then
				do k=nx+1,nx+4
					rho(k,:)=rho(2*nx+1-k,:)
					u(k,:)=-u(2*nx+1-k,:)
					v(k,:)=v(2*nx+1-k,:)
					p(k,:)=p(2*nx+1-k,:)
				enddo
			endif
		endif
!---------------------------------
!	if (Iperiodic_X==0) then
!
!		if (npx==0) then
!
!			u(1,:)=0.d0
!!			v(1,:)=0.d0
!
!		endif
!
!		if (npx==npx0-1) then
!
!			u(nx,:)=0.d0
!!			v(nx,:)=0.d0
!
!		endif
!
!	elseif (Iperiodic_X==1) then
!
!		call check_xbound(u)
!		call check_xbound(v)
!		call check_xbound(p)
!		call check_xbound(rho)
!
!	endif
!-----------------------------------

	call computeR(rho,u,v,p,R_st3,R_nd3,R_rd3,R_th3)

	Q_st=Q_st+(dt/6.d0)*(R_st+2.d0*R_st1+2.d0*R_st2+R_st3)
	Q_nd=Q_nd+(dt/6.d0)*(R_nd+2.d0*R_nd1+2.d0*R_nd2+R_nd3)
	Q_rd=Q_rd+(dt/6.d0)*(R_rd+2.d0*R_rd1+2.d0*R_rd2+R_rd3)
	Q_th=Q_th+(dt/6.d0)*(R_th+2.d0*R_th1+2.d0*R_th2+R_th3)

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st(i,j)
			u(i,j)=Q_nd(i,j)/rho(i,j)
			v(i,j)=Q_rd(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

	if (npy==0) then
		p(:,-3:1)=1.d0
		rho(:,-3:1)=2.d0
		u(:,-3:1)=0.d0
		v(:,-3:1)=0.d0
	endif

	if (npy==npy0-1) then
		p(:,ny:ny+4)=2.5d0
		rho(:,ny:ny+4)=1.d0
		u(:,ny:ny+4)=0.d0
		v(:,ny:ny+4)=0.d0
	endif

		if (Iperiodic_X .ne. 1) then
			if (npx .eq. 0) then
				do k=1,4
					rho(1-k,:)=rho(k,:)
					u(1-k,:)=-u(k,:)
					v(1-k,:)=v(k,:)
					p(1-k,:)=p(k,:)
				enddo
			endif

			if (npx .eq. npx0-1) then
				do k=nx+1,nx+4
					rho(k,:)=rho(2*nx+1-k,:)
					u(k,:)=-u(2*nx+1-k,:)
					v(k,:)=v(2*nx+1-k,:)
					p(k,:)=p(2*nx+1-k,:)
				enddo
			endif
		endif
!---------------------------------
!	if (Iperiodic_X==0) then
 !
!		if (npx==0) then
 !
!			u(1,:)=0.d0
 !
!		endif
 !
!		if (npx==npx0-1) then
 !
!			u(nx,:)=0.d0
 !
!		endif
 !
	!elseif (Iperiodic_X==1) then
 !
	!	call check_xbound(u)
	!	call check_xbound(v)
	!	call check_xbound(p)
	!	call check_xbound(rho)
 !
!	endif
!-----------------------------------

end subroutine RK4



!!!!=============================================================

!!!!==============================================================
subroutine TVD_RK3(rho,u,v,p)
	include 'openNS3d.h'

	integer :: i,j

	real(kind=OCFD_REAL_KIND),parameter::gama=5.d0/3.d0


	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,rho

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st,R_nd,R_rd,R_th
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st1,R_nd1,R_rd1,R_th1
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st2,R_nd2,R_rd2,R_th2
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st3,R_nd3,R_rd3,R_th3

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st,Q_nd,Q_rd,Q_th
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st1,Q_nd1,Q_rd1,Q_th1
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st2,Q_nd2,Q_rd2,Q_th2
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st3,Q_nd3,Q_rd3,Q_th3

!	gama=1.4d0

	call computeR(rho,u,v,p,R_st,R_nd,R_rd,R_th)


	do j=1,ny
		do i=1,nx
			Q_st(i,j)=rho(i,j)
			Q_nd(i,j)=rho(i,j)*u(i,j)
			Q_rd(i,j)=rho(i,j)*v(i,j)
			Q_th(i,j)=p(i,j)/(gama-1.d0)+rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0)
		end do
	end do

	Q_st1=Q_st+dt*R_st
	Q_nd1=Q_nd+dt*R_nd
	Q_rd1=Q_rd+dt*R_rd
	Q_th1=Q_th+dt*R_th

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st1(i,j)
			u(i,j)=Q_nd1(i,j)/rho(i,j)
			v(i,j)=Q_rd1(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th1(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

	if (npy==0) then
		p(:,1:2)=1.d0
		rho(:,1:2)=2.d0
		u(:,1:2)=0.d0
		v(:,1:2)=0.d0
	endif

	if (npy==npy0-1) then
		p(:,ny-1:ny)=2.5d0
		rho(:,ny-1:ny)=1.d0
		u(:,ny-1:ny)=0.d0
		v(:,ny-1:ny)=0.d0
	endif


	call computeR(rho,u,v,p,R_st1,R_nd1,R_rd1,R_th1)


	Q_st2=3.d0/4.d0*Q_st+1.d0/4.d0*Q_st1+(dt/4.d0)*R_st1
	Q_nd2=3.d0/4.d0*Q_nd+1.d0/4.d0*Q_nd1+(dt/4.d0)*R_nd1
	Q_rd2=3.d0/4.d0*Q_rd+1.d0/4.d0*Q_rd1+(dt/4.d0)*R_rd1
	Q_th2=3.d0/4.d0*Q_th+1.d0/4.d0*Q_th1+(dt/4.d0)*R_th1

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st2(i,j)
			u(i,j)=Q_nd2(i,j)/rho(i,j)
			v(i,j)=Q_rd2(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th2(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

	if (npy==0) then
		p(:,1:2)=1.d0
		rho(:,1:2)=2.d0
		u(:,1:2)=0.d0
		v(:,1:2)=0.d0
	endif

	if (npy==npy0-1) then
		p(:,ny-1:ny)=2.5d0
		rho(:,ny-1:ny)=1.d0
		u(:,ny-1:ny)=0.d0
		v(:,ny-1:ny)=0.d0
	endif

	call computeR(rho,u,v,p,R_st2,R_nd2,R_rd2,R_th2)

	Q_st3=Q_st/3.d0+2.d0/3.d0*Q_st2+2.d0/3.d0*dt*R_st2
	Q_nd3=Q_nd/3.d0+2.d0/3.d0*Q_nd2+2.d0/3.d0*dt*R_nd2
	Q_rd3=Q_rd/3.d0+2.d0/3.d0*Q_rd2+2.d0/3.d0*dt*R_rd2
	Q_th3=Q_th/3.d0+2.d0/3.d0*Q_th2+2.d0/3.d0*dt*R_th2

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st3(i,j)
			u(i,j)=Q_nd3(i,j)/rho(i,j)
			v(i,j)=Q_rd3(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th3(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

	if (npy==0) then
		p(:,1:2)=1.d0
		rho(:,1:2)=2.d0
		u(:,1:2)=0.d0
		v(:,1:2)=0.d0
	endif

	if (npy==npy0-1) then
		p(:,ny-1:ny)=2.5d0
		rho(:,ny-1:ny)=1.d0
		u(:,ny-1:ny)=0.d0
		v(:,ny-1:ny)=0.d0
	endif

end subroutine TVD_RK3