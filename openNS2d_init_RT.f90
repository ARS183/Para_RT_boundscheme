
!!!==========================================================================================
!!!Init Taylor Vortex
!!!==========================================================================================

subroutine init_RT(rho,u,v,p,xx,yy) !!
	
	include 'openNS3d.h'
	integer :: i,j,i_global,j_global

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,xx
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: yy,rho,cs
	real(kind=OCFD_REAL_KIND) :: pi,xbeg,ybeg
	integer :: paran  !!!
	

	real(kind=OCFD_REAL_KIND),parameter::gama=5.d0/3.d0
!
	pi=4.d0*datan(1.d0)
	xbeg=0.d0
	ybeg=0.d0

	do i=1,nx
	do j=1,ny
	
		i_global=i_offset(npx)+i-1
		j_global=j_offset(npy)+j-1
		xx(i,j)=xbeg+dble(i_global-1)*hx+0.5d0*hx
		yy(i,j)=ybeg+dble(j_global-1)*hy

		if (yy(i,j) .gt. 0.5d0) then
			rho(i,j)=1.d0
			u(i,j)=0.d0
			p(i,j)=yy(i,j)+1.5d0
			cs(i,j)=dsqrt(gama*p(i,j)/rho(i,j))
			v(i,j)=-0.025d0*cs(i,j)*dcos(8.d0*pi*xx(i,j))
		else
			rho(i,j)=2.d0
			u(i,j)=0.d0
			p(i,j)=2.d0*yy(i,j)+1.d0
			cs(i,j)=dsqrt(gama*p(i,j)/rho(i,j))
			v(i,j)=-0.025d0*cs(i,j)*dcos(8.d0*pi*xx(i,j))
		endif
	enddo
	enddo
	
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
!!
	end subroutine
