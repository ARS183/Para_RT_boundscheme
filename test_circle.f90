program test_circle

	implicit none
	include 'openNS3d.h'

    integer::i

    integer,parameter::nn=5

    real(kind=OCFD_REAL_KIND),allocatable,dimension(:) :: right,s,center,low,up

	allocate(right(1-LAP:nn+LAP),s(1-LAP:nn+LAP))
	allocate(center(1-LAP:nn+LAP),up(1-LAP:nn+LAP),low(1-LAP:nn+LAP))

!更一般的，f(1)在调用本函数前先调用边界函数计算出来

	do i=1,nn
		center(i)=2.d0
		up(i)=1.d0
		low(i)=0.d0
	enddo
	
    right(1)=4.d0
    right(2)=7.d0
    right(3)=10.d0
    right(4)=13.d0
    right(5)=11.d0


	call trid_circle2(1,nn,center,up,low,right)

	s(1:nn)=right(1:nn)

    do i=1,nn
        write(*,*) s(i)
    enddo

end program