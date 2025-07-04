! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  use mod_param     , only:is_impdiff,is_impdiff_1d
#if defined(_LES)
  use mod_param     , only:eps
#endif
  use mod_types
  implicit none
  private
  public chkdt
#if defined(_LES)
  public chkdt_les
#endif
  contains
#if defined(_LES)
  subroutine chkdt_les(n,dl,dzci,dzfi,visc,visct,u,v,w,dtmax)
    !
    ! compute maximum allowed time step, refer to Pieter Wesseling (P200)
    ! for the stability conditions of the advective and diffusion terms
    ! the eddy viscosity term is taken into account in the calculation of dt.
    ! It is acceptable not to consider it, since it is larger than both
    ! the viscous and advective terms. In WRLES, the viscous term is commonly
    ! implicitly treated, so the eddy viscosity term does not influence dtmax,
    ! even when the grid is very fine
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:,0:,0:) :: visct,u,v,w
    real(rp), intent(out) :: dtmax
    real(rp) :: dxi,dyi,dzi
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dti
    real(rp) :: dtidx,dtidy,dtidz,dtid,dl2i
    real(rp) :: viscx,viscy,viscz
    integer :: i,j,k
    !
    dxi = 1./dl(1)
    dyi = 1./dl(2)
    dzi = 1./dl(3)
    dl2i  = dxi*dxi+dyi*dyi
    !
    dti = 0.
    dti = 0.
    !$acc data copy(dti,dtid) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(max:dti,dtid) async(1) &
    !$acc private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) &
    !$acc private(viscx,viscy,viscz,dtidx,dtidy,dtidz)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ux = abs(u(i,j,k))
          vx = 0.25*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          uy = 0.25*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          uz = 0.25*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          dti = max(dti,dtix,dtiy,dtiz)
          !
          viscx = 0.5_rp*(visct(i,j,k)+visct(i+1,j  ,k  ))
          viscy = 0.5_rp*(visct(i,j,k)+visct(i  ,j+1,k  ))
          viscz = 0.5_rp*(visct(i,j,k)+visct(i  ,j  ,k+1))
          dtidx = viscx*(dl2i+dzfi(k)*dzfi(k))
          dtidy = viscy*(dl2i+dzfi(k)*dzfi(k))
          dtidz = viscz*(dl2i+dzci(k)*dzci(k))
          !
          if(is_impdiff .and. .not.is_impdiff_1d) then
            ! nothing to do
          else
            dtidx = dtidx + visc*dl2i
            dtidy = dtidy + visc*dl2i
            dtidz = dtidz + visc*dl2i
            if(.not.is_impdiff_1d) then
              dtidx = dtidx + visc*(dzfi(k)*dzfi(k))
              dtidy = dtidy + visc*(dzfi(k)*dzfi(k))
              dtidz = dtidz + visc*(dzci(k)*dzci(k))
            endif
          endif
          dtid  = max(dtid,dtidx,dtidy,dtidz)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    if(dti  == 0._rp) dti  = 1._rp
    if(dtid == 0._rp) dtid = eps
    dtmax = min(0.4125_rp/dtid,1.732_rp/dti) ! viscous CFL could be 1.5
    call MPI_ALLREDUCE(MPI_IN_PLACE,dtmax,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    if(dti < epsilon(0._rp)) dti = 1.
  end subroutine chkdt_les
#endif
  subroutine chkdt(n,dl,dzci,dzfi,visc,alpha,u,v,w,dtmax)
    !
    ! computes maximum allowed time step
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in) :: visc,alpha
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out) :: dtmax
    real(rp) :: dxi,dyi,dzi
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dti
    integer :: i,j,k
    real(rp), save :: dlmin
    logical , save :: is_first = .true.
    !
    dxi = 1./dl(1)
    dyi = 1./dl(2)
    dzi = 1./dl(3)
    if(is_first) then ! calculate dlmin only once
      is_first = .false.
      dlmin = minval(dl(1:2))
      if(.not.is_impdiff_1d) then
        dlmin = min(dlmin,minval(1./dzfi))
      end if
      call MPI_ALLREDUCE(MPI_IN_PLACE,dlmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    end if
    !
    dti = 0.
    !$acc data copy(dti) async(1)
    !$acc parallel loop collapse(3) default(present) private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) reduction(max:dti) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) REDUCTION(max:dti)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ux = abs(u(i,j,k))
          vx = 0.25*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          uy = 0.25*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          uz = 0.25*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          dti = max(dti,dtix,dtiy,dtiz)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,dti,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(dti < epsilon(0._rp)) dti = 1.
    if(is_impdiff .and. .not.is_impdiff_1d) then
      dtmax = sqrt(3.)/dti
    else
      dtmax = min(1.65/12./max(visc,alpha)*dlmin**2,sqrt(3.)/dti)
    end if
  end subroutine chkdt
end module mod_chkdt
