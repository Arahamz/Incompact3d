!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-RT.f90
!!!      AUTHOR: Arash Hamzehloo
!!! DESCRIPTION: This module describes a single-mode 3D Rayleigh–Taylor instability.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module rayleigh_taylor

  USE decomp_2d
  USE variables
  USE param
  use var, only : dx

  IMPLICIT NONE

  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_rt, boundary_conditions_rt, postprocess_rt

contains


  subroutine init_rt (rho1,ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    real(mytype) :: x, y, z, r, x0, y0, z0, yy, eta, y1, x1
    integer :: k,j,i,is, ii, code

    if (iscalar==1) then
       phi1(:,:,:,1) = zero
    endif

    if (iin.eq.0) then !empty domain
       if (nrank==0) write(*,*) "Empty initial domain!"
       ux1=zero; uy1=zero; uz1=zero
    endif

    x0 = zero
    y0 = -yly / four
    z0 = zero

    eta = xlx/(xsize(1) - one)

    if (iin.eq.1) then

       do k=1,xsize(3)
          z = real(k + xstart(3) - 2, mytype) * dz - half * zlz
          do j=1,xsize(2)
             y = real(j + xstart(2) - 2, mytype) * dy - half * yly
             do i=1,xsize(1)
                x = real(i + xstart(1) - 2, mytype) * dx - half * xlx


                ux1(i,j,k) = u1
                uy1(i,j,k) = zero
                uz1(i,j,k) = zero

                rho1(i,j,k,1) = one

                y1 = -0.05_mytype*(cos(2*pi*x)+cos(2*pi*z))

                phi1(i, j, k, 1) = half*(one+tanh((y1 - y)/(sqrt(two)*eta)))

             enddo
          enddo
       enddo

    endif

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_rt

  subroutine boundary_conditions_rt (ux,uy,uz,phi,ep)

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: x, y, z, x0, y0, z0, cx!, inflow!, outflow
    integer :: ierr
    integer :: i, j, k, is

    ! x0 = zero
    ! y0 = yly / four
    ! z0 = zero
    !
    ! do k=1,xsize(3)
    !    do j=1,xsize(2)
    !      y = real(j + xstart(2) - 2, mytype) * dy
    !     if (y.le.y0) then
    !       bxx1(j,k)=(three/two)*(one+sin(two*pi*seven*(itime*dt)))*(two*y - (y*y))
    !       phi(1,j,k,1) = -one
    !     else
    !       bxx1(j,k)=(three/two)*(one+sin(two*pi*seven*(itime*dt)))
    !       phi(1,j,k,1) = one
    !     endif
    !       bxy1(j,k)=zero
    !       bxz1(j,k)=zero
    !
    !       inflow = inflow + bxx1(j,k)
    !    enddo
    ! enddo
    !
    !
    ! call MPI_ALLREDUCE(inflow,outflow,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    ! outflow = outflow / ny / nz
    !
    ! do k=1,xsize(3)
    !    do j=1,xsize(2)
    !
    !          bxxn(j,k)=ux(xsize(1),j,k)-dt*outflow*(ux(xsize(1),j,k)-ux(xsize(1)-1,j,k))
    !          bxyn(j,k)=uy(xsize(1),j,k)-dt*outflow*(uy(xsize(1),j,k)-uy(xsize(1)-1,j,k))
    !          bxzn(j,k)=uz(xsize(1),j,k)-dt*outflow*(uz(xsize(1),j,k)-uz(xsize(1)-1,j,k))
    !
    !         phi(xsize(1),j,k,1) = phi(xsize(1),j,k,1) - dt * outflow * (phi(xsize(1),j,k,1) - phi(xsize(1)-1,j,k,1))
    !
    !    enddo
    ! enddo

    !call inflow (phi)

    !call outflow (ux,uy,uz,phi)


    if (iscalar.ne.0) then
       if (nclxS1.eq.1) then
          i = 1
          phi(i,:,:,:) = phi(i+1,:,:,:)
       endif
       if (nclxSn.eq.1) then
          i = xsize(1)
          phi(i,:,:,:) = phi(i - 1,:,:,:)

       endif
       if (itimescheme.ne.7) then
          if ((nclyS1.eq.1).and.(xstart(2).eq.1)) then
             !! Generate a hot patch on bottom boundary
             phi(:,1,:,:) = phi(:,2,:,:)
          endif
          if ((nclySn.eq.1).and.(xend(2).eq.ny)) THEN
             phi(:,xsize(2),:,:) = phi(:,xsize(2)-1,:,:)
          endif
       endif

       if ((nclzS1.eq.1).and.(xstart(3).eq.1)) then
          i = 1
          phi(:,:,i,:) = phi(:,:,i+1,:)
       endif
       if ((nclzSn.eq.1).and.(xend(3).eq.nz)) THEN
          phi(:,:,xsize(3),:) = phi(:,:,xsize(3)-1,:)
       endif
    endif

    ! IF (nclx1.EQ.2) THEN
    !    bxx1(:,:) = zero
    !    bxy1(:,:) = zero
    !    bxz1(:,:) = zero
    ! ENDIF
    ! IF (nclxn.EQ.2) THEN
    !    bxxn(:,:) = zero
    !    bxyn(:,:) = zero
    !    bxzn(:,:) = zero
    ! ENDIF
    !
    ! IF (ncly1.EQ.2) THEN
    !    byx1(:,:) = zero
    !    byy1(:,:) = u2
    !    byz1(:,:) = zero
    ! ENDIF
    ! IF (nclyn.EQ.2) THEN
    !    byxn(:,:) = zero
    !    byyn(:,:) = u2
    !    byzn(:,:) = zero
    ! ENDIF
    !
    ! IF (nclz1.EQ.2) THEN
    !    bzx1(:,:) = zero
    !    bzy1(:,:) = zero
    !    bzz1(:,:) = zero
    ! ENDIF
    ! IF (nclzn.EQ.2) THEN
    !    bzxn(:,:) = zero
    !    bzyn(:,:) = zero
    !    bzzn(:,:) = zero
    ! ENDIF

  end subroutine boundary_conditions_rt

  subroutine inflow (phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    integer  :: j,k,is
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)
    do k=1,xsize(3)
       do j=1,xsize(2)
          bxx1(j,k)=one+bxo(j,k)*inflow_noise
          bxy1(j,k)=zero+byo(j,k)*inflow_noise
          bxz1(j,k)=zero+bzo(j,k)*inflow_noise
       enddo
    enddo

    if (iscalar.eq.1) then
       do is=1, numscalar
          do k=1,xsize(3)
             do j=1,xsize(2)
                phi(1,j,k,is)=cp(is)
             enddo
          enddo
       enddo
    endif

    return
  end subroutine inflow

  subroutine outflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609.
    uxmin=1609.
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
          if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
       enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    if (u1.eq.zero) then
       cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
    elseif (u1.eq.one) then
       cx=uxmax1*gdt(itr)*udx
    elseif (u1.eq.two) then
       cx=u2*gdt(itr)*udx    !works better
    else
       stop
    endif

    do k=1,xsize(3)
       do j=1,xsize(2)
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
       enddo
    enddo

    if (iscalar==1) then
       if (u2.eq.zero) then
          cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
       elseif (u2.eq.one) then
          cx=uxmax1*gdt(itr)*udx
       elseif (u2.eq.two) then
          cx=u2*gdt(itr)*udx    !works better
       else
          stop
       endif

       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
          enddo
       enddo
    endif

    if (nrank==0) write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

    return
  end subroutine outflow

  subroutine postprocess_rt(ux1,uy1,uz1,phi1,ep1)

    USE decomp_2d
    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

  end subroutine postprocess_rt

end module rayleigh_taylor
