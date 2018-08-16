	subroutine OCFD_DF_AECDS(u,f,n,h,NUM_METHOD)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,NUM_METHOD

	if (NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_2nd) then
	call OCFD_DF_SCD(u,f,n,h,2)
	goto 100
	
	else if(NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_4th) then
	call OCFD_DF_SCD(u,f,n,h,4)
	goto 100

	else if(NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_6th) then
	call OCFD_DF_SCD(u,f,n,h,6)
	goto 100


	else
	      print*, 'This Numerical Method is not supported'
		  stop
	endif


100	continue


	end subroutine

!================================
!=======f1:dx---upwind
!=======f2:dx---donwwind
!=======s :dxx
	subroutine OCFD_DF_CONV(u,f1,f2,s,n,h,NUM_METHOD)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f1(1-LAP:n+LAP),f2(1-LAP:n+LAP),s(1-LAP:n+LAP)
	integer n,NUM_METHOD

	if (NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_o2nd) then
	call OCFD_DF_SCD_o(u,f1,f2,s,n,h,4)
	goto 100
	
	else if(NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_o4th) then
	call OCFD_DF_SCD_o(u,f1,f2,s,n,h,6)
	goto 100

	else if(NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_o6th) then
	call OCFD_DF_SCD_o(u,f1,f2,s,n,h,8)
	goto 100
	else if(NUM_METHOD.eq.OCFD_NUMERICAL_WENO_5th)    then
	call OCFD_DF_WENO_5th(u,f1,f2,s,n,h)
	goto 100
	else if(NUM_METHOD.eq.OCFD_NUMERICAL_WENO_7th)    then
	call OCFD_DF_WENO_7th(u,f1,f2,s,n,h)
	goto 100	
	


	else
	      print*, 'This Numerical Method is not supported'
		  stop     
	endif


100	continue


	end subroutine



!=================================

	subroutine DF_center(u,f,n,k)
	include 'openNS3d.h'
	integer n,i,k
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP) 

	
	do i=1-LAP+k,n+LAP-k
	f(i)=0.5d0*(u(i+1)-u(i-1))
	enddo
	end subroutine

!==============================

		
	subroutine DF_forward(u,f,n,k)
	include 'openNS3d.h'
	integer n,i,k
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP) 

	
	do i=1-LAP+k,n+LAP-k
	f(i)=(u(i+1)-u(i))
	enddo
	end subroutine
!===============================

	subroutine DF_backward(u,f,n,k)
	include 'openNS3d.h'		
	integer n,i,k
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP) 

	
	do i=1-LAP+k,n+LAP-k
	f(i)=(u(i)-u(i-1))
	enddo
	end subroutine

!================================
	subroutine D2F_center(u,f,n,k)
	include 'openNS3d.h'		
	integer n,i,k
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP) 

	
	do i=1-LAP+k,n+LAP-k
	f(i)=(u(i+1)+u(i-1)-2.d0*u(i))
	enddo
	end subroutine

!================================
    	subroutine D2F_center_2d(u,df2,n,order)
    	include 'openNS3d.h'
    	integer n,i,k,order
    	real(kind=OCFD_REAL_KIND),dimension(1-LAP:n+LAP)::u
    	real(kind=OCFD_REAL_KIND),dimension(1-LAP:n+LAP,LAP) :: df2
    
    	call D2F_center(u,df2(:,1),n,1)
    	do k=2,order/2
    	call D2F_center(df2(:,k-1),df2(:,k),n,k)
    	enddo
    
    
    	end subroutine
!======================================

	subroutine OCFD_DF_SCD(u,f,n,h,order)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	integer n,i,order,k
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP) 
    	real(kind=OCFD_REAL_KIND),dimension(1-LAP:n+LAP,LAP) :: df1,df2	
	real(kind=OCFD_REAL_KIND),dimension(3) :: a 

   
    	a(1)=1.d0
    	a(2)=-dble(1)/dble(6)
    	a(3)=dble(1)/dble(30)
	f=0.d0

	
    	call DF_center(u,df1(:,1),n,1)
	call D2F_center_2d(u,df2,n,order)
    	do k=2,order/2
	call DF_center(df2(:,k-1),df1(:,k),n,k)
	enddo

	do i=1,order/2
	f=f+a(i)*df1(:,i)
	enddo
    	f=f/h

	end subroutine

!============================
	subroutine OCFD_D2F_SCD(u,s,n,h,order)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,h2
	integer n,i,order,k
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),s(1-LAP:n+LAP)
    	real(kind=OCFD_REAL_KIND),dimension(1-LAP:n+LAP,LAP) :: df2	
	real(kind=OCFD_REAL_KIND),dimension(4) :: b

   
    	
	b(1)=1.d0
    	b(2)=-dble(1)/dble(12)
    	b(3)=dble(1)/dble(90)
	b(4)=-dble(1)/dble(560)
	
	s=0.d0
	h2=h*h

	
	call D2F_center_2d(u,df2,n,order)

	
    	do i=1,order/2
    	s=s+b(i)*df2(:,i)
    	enddo
    	s=s/h2	

	

	end subroutine
!============================
	subroutine OCFD_DF_SCD_o(u,f1,f2,s,n,h,order)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,h2
	integer n,i,order,k
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f1(1-LAP:n+LAP),f2(1-LAP:n+LAP),s(1-LAP:n+LAP)
    	real(kind=OCFD_REAL_KIND),dimension(1-LAP:n+LAP,LAP) :: df1,df2	
	real(kind=OCFD_REAL_KIND),dimension(4) :: a,b,c

   
    	a(1)=1.d0
    	a(2)=-dble(1)/dble(6)
    	a(3)=dble(1)/dble(30)
	a(4)=-dble(1)/dble(140)
	b(1)=1.d0
    	b(2)=-dble(1)/dble(12)
    	b(3)=dble(1)/dble(90)
	b(4)=-dble(1)/dble(560)
	c(1)=0.d0
	c(2)=0.d0
	c(3)=0.0536d0
	c(4)=-0.0153d0
	f1=0.d0
	f2=0.d0
	s=0.d0
	h2=h*h

	
    	call DF_center(u,df1(:,1),n,1)
	call D2F_center_2d(u,df2,n,order)

	
    	do i=1,order/2
    	s=s+b(i)*df2(:,i)
    	enddo
    	s=s/h2	

	f1=df1(:,1)
	f2=df1(:,1)
    	do k=2,order/2-1
	call DF_center(df2(:,k-1),df1(:,k),n,k)
	f1=f1+a(k)*df1(:,k)
	f2=f2+a(k)*df1(:,k)	
	enddo
	
	k=order/2
	call DF_backward(df2(:,k-1),df1(:,k),n,k)
	f1=f1+c(k)*df1(:,k)
	call DF_forward(df2(:,k-1),df1(:,k),n,k)
	f2=f2+c(k)*df1(:,k)

	f1=f1/h
	f2=f2/h
	
	end subroutine
	

!===============WENO===============================
	subroutine OCFD_DF_WENO_5th(u,fx1,fx2,s,n,h)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,h2
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),fx1(1-LAP:n+LAP),fx2(1-LAP:n+LAP),s(1-LAP:n+LAP)
	integer n,i
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:n+LAP,LAP) ::df2	
	real(kind=OCFD_REAL_KIND),dimension(4) :: b

	b(1)=1.d0
    	b(2)=-dble(1)/dble(12)
    	b(3)=dble(1)/dble(90)
	b(4)=-dble(1)/dble(560)	
	s=0.d0
	h2=h*h
	call D2F_center_2d(u,df2,n,4)
    	do i=1,2
    	s=s+b(i)*df2(:,i)
    	enddo
    	s=s/h2	


	 
	 do i=1-LAP,n+LAP
	 fx1(i)=u(i)
	 fx2(i)=u(i)
	 end do	

         call dx1_weno5(fx1,n,h)
	 call dx2_weno5(fx2,n,h)
	  

	    
	 return
	end subroutine


	subroutine OCFD_DF_WENO_7th(u,fx1,fx2,s,n,h)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,h2
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),fx1(1-LAP:n+LAP),fx2(1-LAP:n+LAP),s(1-LAP:n+LAP)
	integer n,i
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:n+LAP,LAP) ::df2	
	real(kind=OCFD_REAL_KIND),dimension(4) :: b

	b(1)=1.d0
    	b(2)=-dble(1)/dble(12)
    	b(3)=dble(1)/dble(90)
	b(4)=-dble(1)/dble(560)	
	s=0.d0
	h2=h*h
	call D2F_center_2d(u,df2,n,4)
    	do i=1,2
    	s=s+b(i)*df2(:,i)
    	enddo
    	s=s/h2	


	 
	 do i=1-LAP,n+LAP
	 fx1(i)=u(i)
	 fx2(i)=u(i)
	 end do	

         call dx1_weno7(fx1,n,h)
	 call dx2_weno7(fx2,n,h)
	  

	    
	 return
	end subroutine



	subroutine dx1_weno5(v,n,h)
	include 'openNS3d.h'
      	integer n,k
      	real(kind=OCFD_REAL_KIND):: v(1-LAP:n+LAP),hj(1-LAP:n+LAP),h
      	real(kind=OCFD_REAL_KIND):: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23

     ep=1.d-6
	 C03=3.d0/10.d0
   	 C13=3.d0/5.d0
	 C23=1.d0/10.d0

      do k=0,n
         S0=13.d0/12.d0*(v(k)-2.d0*v(k+1)+v(k+2))**2+   &
             1.d0/4.d0*(3.d0*v(k)-4.d0*v(k+1)+v(k+2))**2

         S1=13.d0/12.d0*(v(k-1)-2.d0*v(k)+v(k+1))**2+   &
            1.d0/4.d0*(v(k-1)-v(k+1))**2

         S2=13.d0/12.d0*(v(k-2)-2.d0*v(k-1)+v(k))**2+   &
           1.d0/4.d0*(v(k-2)-4.d0*v(k-1)+3.d0*v(k))**2

     a0=C03/((ep+S0)**2)
	 a1=C13/((ep+S1)**2)
	 a2=C23/((ep+S2)**2)

	 W0=a0/(a0+a1+a2)
     	 W1=a1/(a0+a1+a2)
	 W2=a2/(a0+a1+a2)

	 q03=1.d0/3.d0*v(k)+5.d0/6.d0*v(k+1)-1.d0/6.d0*v(k+2)
	 q13=-1.d0/6.d0*v(k-1)+5.d0/6.d0*v(k)+1.d0/3.d0*v(k+1)
	 q23=1.d0/3.d0*v(k-2)-7.d0/6.d0*v(k-1)+11.d0/6.d0*v(k)

	 hj(k)=W0*q03+W1*q13+W2*q23
	 
       enddo

       do k=1,n
	 v(k)=(hj(k)-hj(k-1))/h
       enddo
	 return
       end  subroutine



      subroutine dx2_weno5(v,n,h)
	  include 'openNS3d.h'
      	 integer n,k
     	 real(kind=OCFD_REAL_KIND):: v(1-LAP:n+LAP),hj(1-LAP:n+LAP),h
      	 real(kind=OCFD_REAL_KIND):: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23

	 ep=1.d-6
	 
	 C03=3.d0/10.d0
   	 C13=3.d0/5.d0
	 C23=1.d0/10.d0

       do k=1,n+1

         S0=13.d0/12.d0*(v(k)-2.d0*v(k-1)+v(k-2))**2+  &
           1.d0/4.d0*(3.d0*v(k)-4.d0*v(k-1)+v(k-2))**2
       	 	 
         S1=13.d0/12.d0*(v(k+1)-2.d0*v(k)+v(k-1))**2+  &
             1.d0/4.d0*(v(k+1)-v(k-1))**2

         S2=13.d0/12.d0*(v(k+2)-2.d0*v(k+1)+v(k))**2+   &
            1.d0/4.d0*(v(k+2)-4.d0*v(k+1)+3.d0*v(k))**2

     a0=C03/((ep+S0)**2)
	 a1=C13/((ep+S1)**2)
	 a2=C23/((ep+S2)**2)

	 W0=a0/(a0+a1+a2)
     W1=a1/(a0+a1+a2)
	 W2=a2/(a0+a1+a2)

	 q03=1.d0/3.d0*v(k)+5.d0/6.d0*v(k-1)-1.d0/6.d0*v(k-2)
	 q13=-1.d0/6.d0*v(k+1)+5.d0/6.d0*v(k)+1.d0/3.d0*v(k-1)
	 q23=1.d0/3.d0*v(k+2)-7.d0/6.d0*v(k+1)+11.d0/6.d0*v(k)

	 hj(k-1)=W0*q03+W1*q13+W2*q23
	 enddo

       do k=1,n
	   v(k)=(hj(k)-hj(k-1))/h
       enddo	
       return

	 end subroutine


	subroutine dx1_weno7(v,n,h)
	  include 'openNS3d.h'
     	 integer n,k
	  real(kind=OCFD_REAL_KIND):: v(1-LAP:n+LAP),hj(1-LAP:n+LAP),h
      real(kind=OCFD_REAL_KIND):: ep,C0,C1,C2,C3,S0,S1,S2,S3,  &
            s10,s11,s12,s13,s20,s21,s22,s23,s30,s31,s32,s33,  &
            a0,a1,a2,a3,W0,W1,W2,W3,q0,q1,q2,q3

     ep=1.d-7
	 C0=1.d0/35.d0
   	 C1=12.d0/35.d0
	 C2=18.d0/35.d0
	 C3=4.d0/35.d0

      do k=0,n
! 1  ?��????        
      S10=-2.d0/6.d0*v(k-3)+9.d0/6.d0*v(k-2)-18.d0/6.d0*v(k-1) +11.d0/6.d0*v(k)
      S11=1.d0/6.d0*v(k-2)-6.d0/6.d0*v(k-1)+3.d0/6.d0*v(k)    +2.d0/6.d0*v(k+1)
      S12=-2.d0/6.d0*v(k-1)-3.d0/6.d0*v(k)+6.d0/6.d0*v(k+1)-1.d0/6.d0*v(k+2)
      S13=-11.d0/6.d0*v(k)+18.d0/6.d0*v(k+1)-9.d0/6.d0*v(k+2) +2.d0/6.d0*v(k+3)

         S20=-v(k-3)+4.d0*v(k-2)-5.d0*v(k-1)+2.d0*v(k)              ! 2 ?��????
         S21=             v(k-1)-2.d0*v(k)  +v(k+1)         
         S22=             v(k)  -2.d0*v(k+1)+v(k+2)         
         S23=2.d0*v(k)-5.d0*v(k+1)+4.d0*v(k+2)-1.d0*v(k+3)         

         S30=-v(k-3)+3.d0*v(k-2)-3.d0*v(k-1)+v(k)                   ! 3 ?��????                
         S31=-v(k-2)+3.d0*v(k-1)-3.d0*v(k  )+v(k+1)                 
         S32=-v(k-1)+3.d0*v(k  )-3.d0*v(k+1)+v(k+2)                 
         S33=-v(k  )+3.d0*v(k+1)-3.d0*v(k+2)+v(k+3)                 

       S0=S10*S10+13.d0/12.d0*S20*S20+1043.d0/960.d0*S30*S30 +1.d0/12.d0*S10*S30
       S1=S11*S11+13.d0/12.d0*S21*S21+1043.d0/960.d0*S31*S31 +1.d0/12.d0*S11*S31
       S2=S12*S12+13.d0/12.d0*S22*S22+1043.d0/960.d0*S32*S32 +1.d0/12.d0*S12*S32
       S3=S13*S13+13.d0/12.d0*S23*S23+1043.d0/960.d0*S33*S33 +1.d0/12.d0*S13*S33


     a0=C0/((ep+S0)**2)
	 a1=C1/((ep+S1)**2)
	 a2=C2/((ep+S2)**2)
	 a3=C3/((ep+S3)**2)


	 W0=a0/(a0+a1+a2+a3)
     W1=a1/(a0+a1+a2+a3)
	 W2=a2/(a0+a1+a2+a3)
	 W3=a3/(a0+a1+a2+a3)

!  4?��??��???????????
	 q0=-3.d0/12.d0*v(k-3)+13.d0/12.d0*v(k-2)-23.d0/12.d0*v(k-1) +25.d0/12.d0*v(k)
	 q1=1.d0/12.d0*v(k-2)-5.d0/12.d0*v(k-1)+13.d0/12.d0*v(k)  +3.d0/12.d0*v(k+1)
	 q2=-1.d0/12.d0*v(k-1)+7.d0/12.d0*v(k)+7.d0/12.d0*v(k+1)  -1.d0/12.d0*v(k+2)
	 q3=3.d0/12.d0*v(k)+13.d0/12.d0*v(k+1)-5.d0/12.d0*v(k+2)  +1.d0/12.d0*v(k+3)

!  ??4??4?��??��?????����????1??7?��??��?????
	 hj(k)=W0*q0+W1*q1+W2*q2+W3*q3
	 
       enddo

       do k=1,n
	 v(k)=(hj(k)-hj(k-1))/h
       enddo
	 return
	 end
!-------------------------------------------------
!----- ????????????

      subroutine dx2_weno7(v,n,h)
	  include 'openNS3d.h'
      integer n,k
	  real(kind=OCFD_REAL_KIND):: v(1-LAP:n+LAP),hj(1-LAP:n+LAP),h
      real(kind=OCFD_REAL_KIND):: ep,C0,C1,C2,C3,S0,S1,S2,S3,   &
            s10,s11,s12,s13,s20,s21,s22,s23,s30,s31,s32,s33,   &
            a0,a1,a2,a3,W0,W1,W2,W3,q0,q1,q2,q3


     ep=1.d-7
	 C0=1.d0/35.d0
   	 C1=12.d0/35.d0
	 C2=18.d0/35.d0
	 C3=4.d0/35.d0

      do k=1,n+1
! 1  ?��????
         S10=-2.d0/6.d0*v(k+3)+9.d0/6.d0*v(k+2)-18.d0/6.d0*v(k+1)  +11.d0/6.d0*v(k)
         S11=1.d0/6.d0*v(k+2)-6.d0/6.d0*v(k+1)+3.d0/6.d0*v(k)      +2.d0/6.d0*v(k-1)
         S12=-2.d0/6.d0*v(k+1)-3.d0/6.d0*v(k)+6.d0/6.d0*v(k-1)     -1.d0/6.d0*v(k-2)
         S13=-11.d0/6.d0*v(k)+18.d0/6.d0*v(k-1)-9.d0/6.d0*v(k-2)   +2.d0/6.d0*v(k-3)

         S20=-v(k+3)+4.d0*v(k+2)-5.d0*v(k+1)+2.d0*v(k)              ! 2 ?��????
         S21=             v(k+1)-2.d0*v(k)  +v(k-1)         
         S22=             v(k)  -2.d0*v(k-1)+v(k-2)         
         S23=2.d0*v(k)-5.d0*v(k-1)+4.d0*v(k-2)-1.d0*v(k-3)         

         S30=-v(k+3)+3.d0*v(k+2)-3.d0*v(k+1)+v(k)                   ! 3 ?��????                
         S31=-v(k+2)+3.d0*v(k+1)-3.d0*v(k  )+v(k-1)                 
         S32=-v(k+1)+3.d0*v(k  )-3.d0*v(k-1)+v(k-2)                 
         S33=-v(k  )+3.d0*v(k-1)-3.d0*v(k-2)+v(k-3)                 

       S0=S10*S10+13.d0/12.d0*S20*S20  +1043.d0/960.d0*S30*S30 +1.d0/12.d0*S10*S30
       S1=S11*S11+13.d0/12.d0*S21*S21  +1043.d0/960.d0*S31*S31 +1.d0/12.d0*S11*S31
       S2=S12*S12+13.d0/12.d0*S22*S22  +1043.d0/960.d0*S32*S32 +1.d0/12.d0*S12*S32
       S3=S13*S13+13.d0/12.d0*S23*S23  +1043.d0/960.d0*S33*S33 +1.d0/12.d0*S13*S33

     a0=C0/((ep+S0)**2)
	 a1=C1/((ep+S1)**2)
	 a2=C2/((ep+S2)**2)
	 a3=C3/((ep+S3)**2)


	 W0=a0/(a0+a1+a2+a3)
     W1=a1/(a0+a1+a2+a3)
	 W2=a2/(a0+a1+a2+a3)
	 W3=a3/(a0+a1+a2+a3)

!  4?��??��???????????
	 q0=-3.d0/12.d0*v(k+3)+13.d0/12.d0*v(k+2)-23.d0/12.d0*v(k+1) +25.d0/12.d0*v(k)
	 q1=1.d0/12.d0*v(k+2)-5.d0/12.d0*v(k+1)+13.d0/12.d0*v(k)     +3.d0/12.d0*v(k-1)
	 q2=-1.d0/12.d0*v(k+1)+7.d0/12.d0*v(k)+7.d0/12.d0*v(k-1)    -1.d0/12.d0*v(k-2)
	 q3=3.d0/12.d0*v(k)+13.d0/12.d0*v(k-1)-5.d0/12.d0*v(k-2)    +1.d0/12.d0*v(k-3)

!  ??4??4?��??��?????����????1??7?��??��?????
	 hj(k-1)=W0*q0+W1*q1+W2*q2+W3*q3
	 
       enddo

       do k=1,n
	   v(k)=(hj(k)-hj(k-1))/h
       enddo
	 return
	 end



!--by�½�ǿ 2018��
!!!==========================UCC45==========================!!!
subroutine DF_UCC45_P(u,s,f,n,h,flag)
	include 'openNS3d.h'		
	integer n,i,flag
	real(kind=OCFD_REAL_KIND)::h,sig,pi
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::d(1-LAP:n+LAP),s(1-LAP:n+LAP)

	sig=0.54395d0
	pi=4.d0*datan(1.d0)

	!��һ��ģ�f(1)�ڵ��ñ�����ǰ�ȵ��ñ߽纯���������

	if (flag==1) then
		do i=3,n
			d(i)=(-3.d0*sig*u(i-2)+(-42.d0-2.d0*sig)*u(i-1)+(48.d0-12.d0*sig)*u(i) &
			+(-6.d0+18.d0*sig)*u(i+1)-sig*u(i+2))/12.d0/h+(1.d0-sig)*h*s(i)
			f(i)=(d(i)-(1.d0+sig)*f(i-1))/(2.d0+sig)
		enddo
	elseif (flag==2) then
		do i=2,n-2
			d(i)=(-3.d0*sig*u(i-2)+(-42.d0-2.d0*sig)*u(i-1)+(48.d0-12.d0*sig)*u(i) &
			+(-6.d0+18.d0*sig)*u(i+1)-sig*u(i+2))/12.d0/h+(1.d0-sig)*h*s(i)
			f(i)=(d(i)-(1.d0+sig)*f(i-1))/(2.d0+sig)
		enddo
	elseif (flag==0) then
		do i=2,n
			d(i)=(-3.d0*sig*u(i-2)+(-42.d0-2.d0*sig)*u(i-1)+(48.d0-12.d0*sig)*u(i) &
			+(-6.d0+18.d0*sig)*u(i+1)-sig*u(i+2))/12.d0/h+(1.d0-sig)*h*s(i)
			f(i)=(d(i)-(1.d0+sig)*f(i-1))/(2.d0+sig)
		enddo
	elseif (flag==12) then
		do i=3,n-2
			d(i)=(-3.d0*sig*u(i-2)+(-42.d0-2.d0*sig)*u(i-1)+(48.d0-12.d0*sig)*u(i) &
			+(-6.d0+18.d0*sig)*u(i+1)-sig*u(i+2))/12.d0/h+(1.d0-sig)*h*s(i)
			f(i)=(d(i)-(1.d0+sig)*f(i-1))/(2.d0+sig)
		enddo
	endif
	
end subroutine



subroutine DF_UCC45_M(u,s,f,n,h,flag)
	include 'openNS3d.h'		
	integer n,i,flag
	real(kind=OCFD_REAL_KIND)::h,sig,pi
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::d(1-LAP:n+LAP),s(1-LAP:n+LAP)

	sig=0.54395d0
	pi=4.d0*datan(1.d0)

	!��һ��ģ�f(1)�ڵ��ñ�����ǰ�ȵ��ñ߽纯���������

	if (flag==1) then
		do i=n-1,3,-1
			d(i)=(3.d0*sig*u(i+2)+(42.d0+2.d0*sig)*u(i+1)+(-48.d0+12.d0*sig)*u(i) &
			+(6.d0-18.d0*sig)*u(i-1)+sig*u(i-2))/12.d0/h-(1.d0-sig)*h*s(i)
			f(i)=(d(i)-(1.d0+sig)*f(i+1))/(2.d0+sig)
		enddo
	elseif (flag==2) then
		do i=n-2,1,-1
			d(i)=(3.d0*sig*u(i+2)+(42.d0+2.d0*sig)*u(i+1)+(-48.d0+12.d0*sig)*u(i) &
			+(6.d0-18.d0*sig)*u(i-1)+sig*u(i-2))/12.d0/h-(1.d0-sig)*h*s(i)
			f(i)=(d(i)-(1.d0+sig)*f(i+1))/(2.d0+sig)
		enddo
	elseif (flag==0) then
		do i=n-1,1,-1
			d(i)=(3.d0*sig*u(i+2)+(42.d0+2.d0*sig)*u(i+1)+(-48.d0+12.d0*sig)*u(i) &
			+(6.d0-18.d0*sig)*u(i-1)+sig*u(i-2))/12.d0/h-(1.d0-sig)*h*s(i)
			f(i)=(d(i)-(1.d0+sig)*f(i+1))/(2.d0+sig)
		enddo
	elseif (flag==12) then   !chuan xing
		do i=n-2,3,-1
			d(i)=(3.d0*sig*u(i+2)+(42.d0+2.d0*sig)*u(i+1)+(-48.d0+12.d0*sig)*u(i) &
			+(6.d0-18.d0*sig)*u(i-1)+sig*u(i-2))/12.d0/h-(1.d0-sig)*h*s(i)
			f(i)=(d(i)-(1.d0+sig)*f(i+1))/(2.d0+sig)
		enddo
	endif
	
end subroutine
!-------------------------------------------------------------
!
!
!-------------------------periodic----------------------------
subroutine DF_UCC45_P_circle(u,s,f,n,h)
	include 'openNS3d.h'		
	integer n,i
	real(kind=OCFD_REAL_KIND)::h,sig,pi
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::d(1-LAP:n+LAP),s(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::center(1-LAP:n+LAP),up(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::low(1-LAP:n+LAP)

	sig=0.54395d0
	pi=4.d0*datan(1.d0)

	!��һ��ģ�f(1)�ڵ��ñ�����ǰ�ȵ��ñ߽纯���������


		do i=3,n-2
			d(i)=(-3.d0*sig*u(i-2)+(-42.d0-2.d0*sig)*u(i-1)+(48.d0-12.d0*sig)*u(i) &
			+(-6.d0+18.d0*sig)*u(i+1)-sig*u(i+2))/12.d0/h+(1.d0-sig)*h*s(i)
		enddo

		d(1)=(-3.d0*sig*u(n-2)+(-42.d0-2.d0*sig)*u(n-1)+(48.d0-12.d0*sig)*u(1) &
		+(-6.d0+18.d0*sig)*u(2)-sig*u(3))/12.d0/h+(1.d0-sig)*h*s(1)

		d(2)=(-3.d0*sig*u(n-1)+(-42.d0-2.d0*sig)*u(1)+(48.d0-12.d0*sig)*u(2) &
		+(-6.d0+18.d0*sig)*u(3)-sig*u(4))/12.d0/h+(1.d0-sig)*h*s(2)

		d(n-1)=(-3.d0*sig*u(n-3)+(-42.d0-2.d0*sig)*u(n-2)+(48.d0-12.d0*sig)*u(n-1) &
		+(-6.d0+18.d0*sig)*u(n)-sig*u(2))/12.d0/h+(1.d0-sig)*h*s(n-1)

		d(n)=(-3.d0*sig*u(n-2)+(-42.d0-2.d0*sig)*u(n-1)+(48.d0-12.d0*sig)*u(n) &
		+(-6.d0+18.d0*sig)*u(2)-sig*u(3))/12.d0/h+(1.d0-sig)*h*s(n)

		do i=1,n
			center(i)=2.d0+sig
			up(i)=0.d0
			low(i)=1.d0+sig
		enddo

		call trid_circle(1,n,n,center,up,low,d)

		f(1:n)=d(1:n)
	
end subroutine



subroutine DF_UCC45_M_circle(u,s,f,n,h)
	include 'openNS3d.h'		
	integer n,i
	real(kind=OCFD_REAL_KIND)::h,sig,pi
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::d(1-LAP:n+LAP),s(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::center(1-LAP:n+LAP),up(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::low(1-LAP:n+LAP)

	sig=0.54395d0
	pi=4.d0*datan(1.d0)

	!��һ��ģ�f(1)�ڵ��ñ�����ǰ�ȵ��ñ߽纯���������


		do i=n-2,3,-1
			d(i)=(3.d0*sig*u(i+2)+(42.d0+2.d0*sig)*u(i+1)+(-48.d0+12.d0*sig)*u(i) &
			+(6.d0-18.d0*sig)*u(i-1)+sig*u(i-2))/12.d0/h-(1.d0-sig)*h*s(i)
		enddo

		d(n)=(3.d0*sig*u(3)+(42.d0+2.d0*sig)*u(2)+(-48.d0+12.d0*sig)*u(n) &
		+(6.d0-18.d0*sig)*u(n-1)+sig*u(n-2))/12.d0/h-(1.d0-sig)*h*s(n)
		
		d(n-1)=(3.d0*sig*u(2)+(42.d0+2.d0*sig)*u(n)+(-48.d0+12.d0*sig)*u(n-1) &
		+(6.d0-18.d0*sig)*u(n-2)+sig*u(n-3))/12.d0/h-(1.d0-sig)*h*s(n-1)
	
		d(2)=(3.d0*sig*u(4)+(42.d0+2.d0*sig)*u(3)+(-48.d0+12.d0*sig)*u(2) &
		+(6.d0-18.d0*sig)*u(1)+sig*u(n-1))/12.d0/h-(1.d0-sig)*h*s(2)

		d(1)=(3.d0*sig*u(3)+(42.d0+2.d0*sig)*u(2)+(-48.d0+12.d0*sig)*u(1) &
		+(6.d0-18.d0*sig)*u(n-1)+sig*u(n-2))/12.d0/h-(1.d0-sig)*h*s(1)

		do i=1,n
			center(i)=2.d0+sig
			low(i)=0.d0
			up(i)=1.d0+sig
		enddo

		call trid_circle(1,n,n,center,up,low,d)

		f(1:n)=d(1:n)


end subroutine

!=============================================================




!---------------------------------------------------------------------
!���׵���
!---------------------------------------------------------------------
subroutine D2F_PADE4(u,s,n,h)
	include 'openNS3d.h'		
	integer n,i
	real(kind=OCFD_REAL_KIND)::h,h2,pi,pi2
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::right(1-LAP:n+LAP),s(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::center(1-LAP:n+LAP),up(1-LAP:n+LAP),low(1-LAP:n+LAP)
	
	pi=4.d0*datan(1.d0)
	pi2=pi*pi
	h2=h*h

!��һ��ģ�f(1)�ڵ��ñ�����ǰ�ȵ��ñ߽纯���������

	do i=2,n-1
		center(i)=5.d0/6.d0
		up(i)=1.d0/12.d0
		low(i)=1.d0/12.d0
		right(i)=(u(i+1)+u(i-1)-2.d0*u(i))/h2
	enddo
	right(2)=right(2)-low(2)*s(1)
	right(n-1)=right(n-1)-up(n-1)*s(n)	

	call trid(2,n-1,n,center,up,low,right)

	s(2:n-1)=right(2:n-1)
end subroutine


!-----------------���ڱ߽�---------------------��
subroutine D2F_PADE4_circle(u,s,n,h)
	include 'openNS3d.h'		
	integer n,i
	real(kind=OCFD_REAL_KIND)::h,h2,pi,pi2
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::right(1-LAP:n+LAP),s(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::center(1-LAP:n+LAP),up(1-LAP:n+LAP),low(1-LAP:n+LAP)
	
	pi=4.d0*datan(1.d0)
	pi2=pi*pi
	h2=h*h

!��һ��ģ�f(1)�ڵ��ñ�����ǰ�ȵ��ñ߽纯���������

	
		center=5.d0/6.d0
		up=1.d0/12.d0
		low=1.d0/12.d0

	right(1)=(u(2)+u(n-1)-2.d0*u(1))/h2
	do i=2,n-1
		right(i)=(u(i+1)+u(i-1)-2.d0*u(i))/h2
	enddo
	right(n)=(u(2)+u(n-1)-2.d0*u(n))/h2	

	call trid_circle(1,n,n,center,up,low,right)

	s(1:n)=right(1:n)
end subroutine
!-----------------------------------------------



!--------------------------------------------------------------------
!���׵���D6F_Center
!--------------------------------------------------------------------	
subroutine D6F_Center(u,d6f,n,h,flag)
	include 'openNS3d.h'		
	integer n,i,flag
	real(kind=OCFD_REAL_KIND)::h,h2,h4,h6
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),d6f(1-LAP:n+LAP)

	h2=h*h
	h4=h2*h2
	h6=h4*h2

	if (flag==1) then
		!������߽߱�
		do i=4,n
			d6f(i)=(u(i-3)-6.d0*u(i-2)+15.d0*u(i-1)-20.d0*u(i) &
			+15.d0*u(i+1)-6.d0*u(i+2)+u(i+3))/h6
		enddo
	elseif (flag==2) then
		!�����ұ߽߱�
		do i=1,n-3
			d6f(i)=(u(i-3)-6.d0*u(i-2)+15.d0*u(i-1)-20.d0*u(i) &
			+15.d0*u(i+1)-6.d0*u(i+2)+u(i+3))/h6
		enddo
	elseif (flag==12) then
		!����
		do i=4,n-3
			d6f(i)=(u(i-3)-6.d0*u(i-2)+15.d0*u(i-1)-20.d0*u(i) &
			+15.d0*u(i+1)-6.d0*u(i+2)+u(i+3))/h6
		enddo
	elseif (flag==0) then
		!�ڲ�����
		do i=1,n
			d6f(i)=(u(i-3)-6.d0*u(i-2)+15.d0*u(i-1)-20.d0*u(i) &
			+15.d0*u(i+1)-6.d0*u(i+2)+u(i+3))/h6
		enddo
	endif

end subroutine

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!


!!===========׷�Ϸ�===========!!
subroutine trid(n0,nn,n,a,b,c,d)
	include 'openNS3d.h'		
	integer k,i,n0,n,nn
	real(kind=OCFD_REAL_KIND)::a(1-LAP:n+LAP),b(1-LAP:n+LAP),c(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::d(1-LAP:n+LAP),p(1-LAP:n+LAP),q(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::y(1-LAP:n+LAP)

	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	!c   a---center  b---up   c---low     d---right
	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	p(n0)=a(n0)
	q(n0)=b(n0)/a(n0)
		  
	do k=n0+1,nn-1
		p(k)=a(k)-c(k)*q(k-1)
		q(k)=b(k)/p(k)
	end do
		  
	p(nn)=a(nn)-c(nn)*q(nn-1)
	y(n0)=d(n0)/p(n0)
		  
	do k=n0+1,nn
		y(k)=(d(k)-c(k)*y(k-1))/p(k)
	end do
	 d(nn)=y(nn)
		  
	do  k=nn-1,n0,-1
		d(k)=y(k)-q(k)*d(k+1)
	end do	  
	
	return
end	subroutine
    

!===ѭ�����Խǵ�׷�Ϸ�===
!subroutine trid_circle_incorrect(n0,n,center,up,low,right)
	!  !
 !     ! a:center b:up c:low d:right e:upright-corner f:downleft-corner
	!include 'openNS3d.h'
	!real(kind=OCFD_REAL_KIND)::urc,dlc,w1,w2,det
	!real(kind=OCFD_REAL_KIND)::center(1-LAP:n+LAP),up(1-LAP:n+LAP)
	!real(kind=OCFD_REAL_KIND)::low(1-LAP:n+LAP),u(1-LAP:n+LAP),v(1-LAP:n+LAP)
	!real(kind=OCFD_REAL_KIND)::right(1-LAP:n+LAP),x(1-LAP:n+LAP)
	!integer :: n0,n
 !
	!  urc=low(1)
	!  dlc=up(n)
	!  
	!  call trid(n0,n,center,up,low,right)
	!  x(n0:n)=right(n0:n)
 !
	!  right(n0)=1.d0
	!  right(n0+1:n)=0.d0
	!  call trid(n0,n,center,up,low,right)
	!  u(n0:n)=right(n0:n)
 !
	!  right(n)=1.d0
	!  right(n0:n-1)=0.d0
	!  call trid(n0,n,center,up,low,right)
	!  v(n0:n)=right(n0:n)
 !
	!  det=u(1)*v(n)-(v(1)+1.d0/dlc)*(u(n)+1.d0/urc)
	!  w1=(v(n)*x(1)-(v(1)+1.d0/dlc)*x(n))/det
	!  w2=(-(u(n)+1.d0/urc)*x(1)+u(1)*x(n))/det
 !
	!  x=x-w1*u-w2*v
 !
	!  right(n0:n)=x(n0:n)
 !
	!  return
 !
	!  end subroutine


	
!===ѭ�����Խǵ�׷�Ϸ�===
subroutine trid_circle(n0,nn,n,b,c,a,d)
	include 'openNS3d.h'
	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	!c   b---center  c---up   a---low     d---right
	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer :: n0,n,nn
	real(kind=OCFD_REAL_KIND)::s
	real(kind=OCFD_REAL_KIND)::b(1-LAP:n+LAP),c(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::a(1-LAP:n+LAP),uu(1-LAP:n+LAP),vv(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::d(1-LAP:n+LAP),p(1-LAP:n+LAP),q(1-LAP:n+LAP)

	real(kind=OCFD_REAL_KIND)::r(1-LAP:n+LAP),w(1-LAP:n+LAP)

	integer::j

	p(n0)=b(n0)
	q(n0)=c(n0)/p(n0)
	r(n0)=a(n0)/b(n0)

    do j=n0+1,nn-2
        p(j)=b(j)-a(j)*q(j-1)
        q(j)=c(j)/p(j)
		r(j)=-a(j)*r(j-1)/p(j)
    end do
    p(nn-1)=b(nn-1)-a(nn-1)*q(nn-2)
    r(nn-1)=(c(nn-1)-a(nn-1)*r(nn-2))/p(nn-1)
	p(nn)=b(nn)-a(nn)*r(nn-1)
	s=c(nn)/p(nn)

	w(nn-1)=r(nn-1)
    do j=nn-2,n0,-1 
        w(j)=r(j)-q(j)*w(j+1)
    end do
	
	uu(n0)=d(n0)/p(n0)
	do j=n0+1,nn
        uu(j)=(d(j)-a(j)*uu(j-1))/p(j)
    end do

	vv(nn)=uu(nn)
	vv(nn-1)=uu(nn-1)
	do j=nn-2,n0,-1
        vv(j)=uu(j)-q(j)*vv(j+1)
    end do

	d(nn)=(vv(nn)-s*vv(n0))/(1.d0-w(n0)*s)
	do j=nn-1,n0,-1
        d(j)=vv(j)-w(j)*d(nn)
	end do

end



!===================WENO5======================
subroutine du1_weno5(v,df,n,h)
	include 'openNS3d.h'
      	integer n,k
      	real(kind=OCFD_REAL_KIND):: v(1-LAP:n+LAP),df(1-LAP:n+LAP),hj(1-LAP:n+LAP),h
      	real(kind=OCFD_REAL_KIND):: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23

     ep=1.d-6
	 C03=3.d0/10.d0
   	 C13=3.d0/5.d0
	 C23=1.d0/10.d0

      do k=0,n
         S0=13.d0/12.d0*(v(k)-2.d0*v(k+1)+v(k+2))**2+   &
             1.d0/4.d0*(3.d0*v(k)-4.d0*v(k+1)+v(k+2))**2

         S1=13.d0/12.d0*(v(k-1)-2.d0*v(k)+v(k+1))**2+   &
            1.d0/4.d0*(v(k-1)-v(k+1))**2

         S2=13.d0/12.d0*(v(k-2)-2.d0*v(k-1)+v(k))**2+   &
           1.d0/4.d0*(v(k-2)-4.d0*v(k-1)+3.d0*v(k))**2

     a0=C03/((ep+S0)**2)
	 a1=C13/((ep+S1)**2)
	 a2=C23/((ep+S2)**2)

	 W0=a0/(a0+a1+a2)
     	 W1=a1/(a0+a1+a2)
	 W2=a2/(a0+a1+a2)

	 q03=1.d0/3.d0*v(k)+5.d0/6.d0*v(k+1)-1.d0/6.d0*v(k+2)
	 q13=-1.d0/6.d0*v(k-1)+5.d0/6.d0*v(k)+1.d0/3.d0*v(k+1)
	 q23=1.d0/3.d0*v(k-2)-7.d0/6.d0*v(k-1)+11.d0/6.d0*v(k)

	 hj(k)=W0*q03+W1*q13+W2*q23
	 
       enddo

       do k=1,n
	 df(k)=(hj(k)-hj(k-1))/h
       enddo
	 return
       end  subroutine



      subroutine du2_weno5(v,df,n,h)
	  include 'openNS3d.h'
      	 integer n,k
     	 real(kind=OCFD_REAL_KIND):: v(1-LAP:n+LAP),hj(1-LAP:n+LAP),df(1-LAP:n+LAP),h
      	 real(kind=OCFD_REAL_KIND):: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23

	 ep=1.d-6
	 
	 C03=3.d0/10.d0
   	 C13=3.d0/5.d0
	 C23=1.d0/10.d0

       do k=1,n+1

         S0=13.d0/12.d0*(v(k)-2.d0*v(k-1)+v(k-2))**2+  &
           1.d0/4.d0*(3.d0*v(k)-4.d0*v(k-1)+v(k-2))**2
       	 	 
         S1=13.d0/12.d0*(v(k+1)-2.d0*v(k)+v(k-1))**2+  &
             1.d0/4.d0*(v(k+1)-v(k-1))**2

         S2=13.d0/12.d0*(v(k+2)-2.d0*v(k+1)+v(k))**2+   &
            1.d0/4.d0*(v(k+2)-4.d0*v(k+1)+3.d0*v(k))**2

     a0=C03/((ep+S0)**2)
	 a1=C13/((ep+S1)**2)
	 a2=C23/((ep+S2)**2)

	 W0=a0/(a0+a1+a2)
     W1=a1/(a0+a1+a2)
	 W2=a2/(a0+a1+a2)

	 q03=1.d0/3.d0*v(k)+5.d0/6.d0*v(k-1)-1.d0/6.d0*v(k-2)
	 q13=-1.d0/6.d0*v(k+1)+5.d0/6.d0*v(k)+1.d0/3.d0*v(k-1)
	 q23=1.d0/3.d0*v(k+2)-7.d0/6.d0*v(k+1)+11.d0/6.d0*v(k)

	 hj(k-1)=W0*q03+W1*q13+W2*q23
	 enddo

       do k=1,n
	   df(k)=(hj(k)-hj(k-1))/h
       enddo	
       return

	 end subroutine