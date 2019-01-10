 subroutine LaxE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,Eneg1,Eneg2,Eneg3,Eneg4)
    include 'openNS3d.h'

    integer :: i,j

    real(kind=OCFD_REAL_KIND) :: gama,unit1,unit2,unit3

    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: rho,u,v,p
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Epos1,Epos2,Epos3,Epos4
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Eneg1,Eneg2,Eneg3,Eneg4
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: cs,Me,eng,lamdaE
 

    gama=5.d0/3.d0
    do j=1-LAP,ny+LAP
        do i=1-LAP,nx+LAP
            cs(i,j)=dsqrt(gama*p(i,j)/rho(i,j))
			lamdaE(i,j)=DABS(u(i,j))+cs(i,j)
            eng(i,j)=p(i,j)/(gama-1.d0)+rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0)

            Epos1(i,j)=(rho(i,j)*u(i,j)+lamdaE(i,j)*rho(i,j))/2.d0
            Epos2(i,j)=(rho(i,j)*u(i,j)*u(i,j)+p(i,j)+lamdaE(i,j)*rho(i,j)*u(i,j))/2.d0
            Epos3(i,j)=(rho(i,j)*u(i,j)*v(i,j)+lamdaE(i,j)*rho(i,j)*v(i,j))/2.d0
            Epos4(i,j)=((eng(i,j)+p(i,j))*u(i,j)+lamdaE(i,j)*eng(i,j))/2.d0

            Eneg1(i,j)=(rho(i,j)*u(i,j)-lamdaE(i,j)*rho(i,j))/2.d0
            Eneg2(i,j)=(rho(i,j)*u(i,j)*u(i,j)+p(i,j)-lamdaE(i,j)*rho(i,j)*u(i,j))/2.d0
            Eneg3(i,j)=(rho(i,j)*u(i,j)*v(i,j)-lamdaE(i,j)*rho(i,j)*v(i,j))/2.d0
            Eneg4(i,j)=((eng(i,j)+p(i,j))*u(i,j)-lamdaE(i,j)*eng(i,j))/2.d0

        end do
    end do

end subroutine




subroutine LaxF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,Fneg1,Fneg2,Fneg3,Fneg4)
    include 'openNS3d.h'

    integer :: i,j

    real(kind=OCFD_REAL_KIND) :: gama,unit1,unit2,unit3

    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: rho,u,v,p
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Fpos1,Fpos2,Fpos3,Fpos4
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Fneg1,Fneg2,Fneg3,Fneg4
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: cs,Mf,eng,lamdaF


    gama=5.d0/3.d0
    do j=1-LAP,ny+LAP
        do i=1-LAP,nx+LAP
            cs(i,j)=dsqrt(gama*p(i,j)/rho(i,j))
            lamdaF(i,j)=DABS(v(i,j))+cs(i,j)
            eng(i,j)=p(i,j)/(gama-1.d0)+rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0)
            
			
				Fpos1(i,j)=(rho(i,j)*v(i,j)+lamdaF(i,j)*rho(i,j))/2.d0
                Fpos2(i,j)=(rho(i,j)*v(i,j)*u(i,j)+lamdaF(i,j)*rho(i,j)*u(i,j))/2.d0
				Fpos3(i,j)=(rho(i,j)*v(i,j)*v(i,j)+p(i,j)+lamdaF(i,j)*rho(i,j)*v(i,j))/2.d0
                Fpos4(i,j)=((eng(i,j)+p(i,j))*v(i,j)+lamdaF(i,j)*eng(i,j))/2.d0
				
                Fneg1(i,j)=(rho(i,j)*v(i,j)-lamdaF(i,j)*rho(i,j))/2.d0
                Fneg2(i,j)=(rho(i,j)*v(i,j)*u(i,j)-lamdaF(i,j)*rho(i,j)*u(i,j))/2.d0
				Fneg3(i,j)=(rho(i,j)*v(i,j)*v(i,j)+p(i,j)-lamdaF(i,j)*rho(i,j)*v(i,j))/2.d0
                Fneg4(i,j)=((eng(i,j)+p(i,j))*v(i,j)-lamdaF(i,j)*eng(i,j))/2.d0


        end do
    end do

end subroutine