subroutine StegerE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,Eneg1,Eneg2,Eneg3,Eneg4)
    include 'openNS3d.h'

    integer :: i,j

    real(kind=OCFD_REAL_KIND) :: gama,unit1,unit2,unit3

    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: rho,u,v,p
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Epos1,Epos2,Epos3,Epos4
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Eneg1,Eneg2,Eneg3,Eneg4
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: cs,Me,eng
 

    gama=5.d0/3.d0
    do j=1-LAP,ny+LAP
        do i=1-LAP,nx+LAP
            cs(i,j)=dsqrt(gama*p(i,j)/rho(i,j))
            Me(i,j)=u(i,j)/cs(i,j)
            eng(i,j)=p(i,j)/(gama-1.d0)+rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0)
            
            if (Me(i,j).le.(-1.d0)) then
                Epos1(i,j)=0.d0
                Epos2(i,j)=0.d0
                Epos3(i,j)=0.d0
                Epos4(i,j)=0.d0

                Eneg1(i,j)=rho(i,j)*u(i,j)
                Eneg2(i,j)=rho(i,j)*u(i,j)*u(i,j)+p(i,j)
                Eneg3(i,j)=rho(i,j)*u(i,j)*v(i,j)
                Eneg4(i,j)=(eng(i,j)+p(i,j))*u(i,j)

            else if (Me(i,j).le.0.d0 .and. Me(i,j).gt.(-1.d0)) then
                unit1=0.5d0*rho(i,j)*(u(i,j)+cs(i,j))/gama
                unit2=(gama-1.d0)/gama*rho(i,j)*u(i,j)
                unit3=0.5d0*rho(i,j)*(u(i,j)-cs(i,j))/gama

                Epos1(i,j)=unit1
                Epos2(i,j)=unit1*(u(i,j)+cs(i,j))
                Epos3(i,j)=unit1*v(i,j)
                Epos4(i,j)=unit1*(cs(i,j)*cs(i,j)/(gama-1.d0)+cs(i,j)*u(i,j) &
						  +0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0))

                Eneg1(i,j)=unit2+unit3
                Eneg2(i,j)=unit2*u(i,j)+unit3*(u(i,j)-cs(i,j))
                Eneg3(i,j)=unit2*v(i,j)+unit3*v(i,j)
                Eneg4(i,j)=unit2*0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0)+unit3*(cs(i,j) &
						  *cs(i,j)/(gama-1.d0)+0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0)-u(i,j)*cs(i,j))
				
			else if (Me(i,j).le.1.d0 .and. Me(i,j).gt.0.d0) then
				unit1=0.5d0*rho(i,j)*(u(i,j)+cs(i,j))/gama
                unit2=(gama-1.d0)/gama*rho(i,j)*u(i,j)
                unit3=0.5d0*rho(i,j)*(u(i,j)-cs(i,j))/gama
				
				Eneg1(i,j)=unit3
                Eneg2(i,j)=unit3*(u(i,j)-cs(i,j))
                Eneg3(i,j)=unit3*v(i,j)
                Eneg4(i,j)=unit3*(cs(i,j)*cs(i,j)/(gama-1.d0)-cs(i,j)*u(i,j) &
						  +0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0))

                Epos1(i,j)=unit1+unit2
                Epos2(i,j)=unit2*u(i,j)+unit1*(u(i,j)+cs(i,j))
                Epos3(i,j)=unit2*v(i,j)+unit1*v(i,j)
                Epos4(i,j)=unit2*0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0)+unit1*(cs(i,j) &
						  *cs(i,j)/(gama-1.d0)+0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0)+u(i,j)*cs(i,j))
			
			
			else
				Epos1(i,j)=rho(i,j)*u(i,j)
                Epos2(i,j)=rho(i,j)*u(i,j)*u(i,j)+p(i,j)
                Epos3(i,j)=rho(i,j)*u(i,j)*v(i,j)
                Epos4(i,j)=(eng(i,j)+p(i,j))*u(i,j)

                Eneg1(i,j)=0.d0
                Eneg2(i,j)=0.d0
                Eneg3(i,j)=0.d0
                Eneg4(i,j)=0.d0
            
            end if
        end do
    end do

end subroutine




subroutine StegerF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,Fneg1,Fneg2,Fneg3,Fneg4)
    include 'openNS3d.h'

    integer :: i,j

    real(kind=OCFD_REAL_KIND) :: gama,unit1,unit2,unit3

    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: rho,u,v,p
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Fpos1,Fpos2,Fpos3,Fpos4
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Fneg1,Fneg2,Fneg3,Fneg4
    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: cs,Mf,eng


    gama=5.d0/3.d0
    do j=1-LAP,ny+LAP
        do i=1-LAP,nx+LAP
            cs(i,j)=dsqrt(gama*p(i,j)/rho(i,j))
            Mf(i,j)=v(i,j)/cs(i,j)
            eng(i,j)=p(i,j)/(gama-1.d0)+rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0)
			
            if (Mf(i,j).le.(-1.d0)) then
                Fpos1(i,j)=0.d0
                Fpos3(i,j)=0.d0
                Fpos2(i,j)=0.d0
                Fpos4(i,j)=0.d0

                Fneg1(i,j)=rho(i,j)*v(i,j)
                Fneg3(i,j)=rho(i,j)*v(i,j)*v(i,j)+p(i,j)
                Fneg2(i,j)=rho(i,j)*v(i,j)*u(i,j)
                Fneg4(i,j)=(eng(i,j)+p(i,j))*v(i,j)

            else if (Mf(i,j).le.0.d0 .and. Mf(i,j).gt.(-1.d0)) then
                unit1=0.5d0*rho(i,j)*(v(i,j)+cs(i,j))/gama
                unit2=(gama-1.d0)/gama*rho(i,j)*v(i,j)
                unit3=0.5d0*rho(i,j)*(v(i,j)-cs(i,j))/gama

                Fpos1(i,j)=unit1
                Fpos3(i,j)=unit1*(v(i,j)+cs(i,j))
                Fpos2(i,j)=unit1*u(i,j)
                Fpos4(i,j)=unit1*(cs(i,j)*cs(i,j)/(gama-1.d0)+cs(i,j)*v(i,j) &
						  +0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0))

                Fneg1(i,j)=unit2+unit3
                Fneg3(i,j)=unit2*v(i,j)+unit3*(v(i,j)-cs(i,j))
                Fneg2(i,j)=unit2*u(i,j)+unit3*u(i,j)
                Fneg4(i,j)=unit2*0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0)+unit3*(cs(i,j) &
						  *cs(i,j)/(gama-1.d0)+0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0)-v(i,j)*cs(i,j))
				
			else if (Mf(i,j).le.1.d0 .and. Mf(i,j).gt.0.d0) then
				unit1=0.5d0*rho(i,j)*(v(i,j)+cs(i,j))/gama
                unit2=(gama-1.d0)/gama*rho(i,j)*v(i,j)
                unit3=0.5d0*rho(i,j)*(v(i,j)-cs(i,j))/gama
				
				Fneg1(i,j)=unit3
                Fneg3(i,j)=unit3*(v(i,j)-cs(i,j))
                Fneg2(i,j)=unit3*u(i,j)
                Fneg4(i,j)=unit3*(cs(i,j)*cs(i,j)/(gama-1)-cs(i,j)*v(i,j) &
						  +0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0))

                Fpos1(i,j)=unit1+unit2
                Fpos3(i,j)=unit2*v(i,j)+unit1*(v(i,j)+cs(i,j))
                Fpos2(i,j)=unit2*u(i,j)+unit1*u(i,j)
                Fpos4(i,j)=unit2*0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0)+unit1*(cs(i,j) &
						  *cs(i,j)/(gama-1.d0)+0.5d0*(u(i,j)**2.d0+v(i,j)**2.d0)+v(i,j)*cs(i,j))
			
			
			else
				Fpos1(i,j)=rho(i,j)*v(i,j)
                Fpos3(i,j)=rho(i,j)*v(i,j)*v(i,j)+p(i,j)
                Fpos2(i,j)=rho(i,j)*v(i,j)*u(i,j)
                Fpos4(i,j)=(eng(i,j)+p(i,j))*v(i,j)

                Fneg1(i,j)=0.d0
                Fneg3(i,j)=0.d0
                Fneg2(i,j)=0.d0
                Fneg4(i,j)=0.d0
            
            end if
        end do
    end do

end subroutine