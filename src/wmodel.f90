! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
#if defined(_LES)
module mod_wmodel
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan,is_finite => ieee_is_finite
  use mpi
  use mod_types
  use mod_typedef, only: bound
  use mod_param, only: kap_log,b_log,eps
  implicit none
  private
  public updt_wallmodelbc
  contains
    !
  subroutine updt_wallmodelbc(n,is_bound,lwm,l,dl,zc,zf,dzc,dzf,visc,h,index_wm,u,v,w, &
                              bcu,bcv,bcw,bcu_mag,bcv_mag,bcw_mag)
    !
    ! update wall model bc wall by wall
    ! 
    implicit none
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:1,3) :: is_bound
    integer , intent(in), dimension(0:1,3) :: lwm,index_wm
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc,zf,dzc,dzf
    real(rp), intent(in) :: visc,h
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    type(bound), intent(inout) :: bcu,bcv,bcw
    type(bound), intent(in   ) :: bcu_mag,bcv_mag,bcw_mag
    real(rp) :: wei,coef,uh,vh,wh,u1,u2,v1,v2,w1,w2,u_mag,v_mag,w_mag,tauw(2)
    integer  :: nh,i,j,k,i1,i2,j1,j2,k1,k2
    !
    nh = 1
    !
    if(is_bound(0,1).and.lwm(0,1)/=0) then
      call cmpt_wallmodelbc(n,0,1,nh,lwm(0,1),l,dl,zc,zf,dzc,dzf,visc,h,index_wm(0,1),v,w, &
                            bcv%x,bcw%x,bcv_mag%x,bcw_mag%x)
    end if
    if(is_bound(1,1).and.lwm(1,1)/=0) then
      call cmpt_wallmodelbc(n,1,1,nh,lwm(1,1),l,dl,zc,zf,dzc,dzf,visc,h,index_wm(1,1),v,w, &
                            bcv%x,bcw%x,bcv_mag%x,bcw_mag%x)
    end if
    if(is_bound(0,2).and.lwm(0,2)/=0) then
      call cmpt_wallmodelbc(n,0,2,nh,lwm(0,2),l,dl,zc,zf,dzc,dzf,visc,h,index_wm(0,2),u,w, &
                            bcu%y,bcw%y,bcu_mag%y,bcw_mag%y)
    end if
    if(is_bound(1,2).and.lwm(1,2)/=0) then
      call cmpt_wallmodelbc(n,1,2,nh,lwm(1,2),l,dl,zc,zf,dzc,dzf,visc,h,index_wm(1,2),u,w, &
                            bcu%y,bcw%y,bcu_mag%y,bcw_mag%y)
    end if
    if(is_bound(0,3).and.lwm(0,3)/=0) then
      call cmpt_wallmodelbc(n,0,3,nh,lwm(0,3),l,dl,zc,zf,dzc,dzf,visc,h,index_wm(0,3),u,v, &
                            bcu%z,bcv%z,bcu_mag%z,bcv_mag%z)
    end if
    if(is_bound(1,3).and.lwm(1,3)/=0) then
      call cmpt_wallmodelbc(n,1,3,nh,lwm(1,3),l,dl,zc,zf,dzc,dzf,visc,h,index_wm(1,3),u,v, &
                            bcu%z,bcv%z,bcu_mag%z,bcv_mag%z)
    end if
  end subroutine updt_wallmodelbc
  !
  subroutine cmpt_wallmodelbc(n,ibound,idir,nh,mtype,l,dl,zc,zf,dzc,dzf,visc,h,index,vel1,vel2, &
                              bcvel1,bcvel2,bcvel1_mag,bcvel2_mag)
    !
    ! compute wall model bc of a wall
    !
    ! wall-parallel velocity at ghost cells is used only for computing the viscous terms,
    ! including its left- and right-hand sides. It is not used for convective terms, or
    ! the correction procedure. In WMLES, a wall is a no-slip wall with corrected wall stress.
    ! When wall stress is required for computing viscous terms, the wall is regarded as a
    ! Neumann bc. When wall velocity is required for computing work, it is regarded as a
    ! no-slip wall. When filtering/strain rate is required for computing eddy viscosity,
    ! the wall is a slip wall, with the wall velocity extrapolated from the interior. Hence,
    ! a ghost point can have three different values.
    !
    ! index 0 must be calculated for the right/front/top walls (chkdt), but not necessary
    ! for the opposite walls. However, index 0 for the left/back/bottom walls is necessary
    ! when a subgrid model is involved (cmpt_dw_plus). Hence, the best practice is to
    ! include index 0.
    !
    ! The singularity at the corner for cavity flow affects the calculation of time step.
    ! However, the velocity is small at the corner, so it does not impose a limitation on
    ! the time step. Hence, no special treatment is necessary at the corner.
    ! The ghost point outside of the singularity is essentially not used for computing
    ! u and w at the corners, since u and w at the corners are set as zero by the
    ! no-penetration boundary condition. When eddy viscosity (strain rate) is used,
    ! its calculation uses the ghost point outside of the singularity. However,
    ! when wall model is used, one-sided difference is used in the calculation of
    ! of wall gradients (strain rate), so the ghost point outside of the singularity is
    ! not used. If van Driest damping is used, the ghost point outside of the singularity
    ! is used, but the damping is not important at the corner. In summary,
    ! WMLES: the van Driest damping and time step calculation use the ghost point
    ! outside of the singularity.
    ! WRLES: the van Driest damping, time step calculation and strain rate calculation
    ! use the ghost point outside of the singularity.
    ! DNS: time step calculation uses the ghost point outside of the singularity.
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: ibound,idir,nh,mtype
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc,zf,dzc,dzf
    real(rp), intent(in) :: visc,h
    integer , intent(in) :: index
    real(rp), intent(in), dimension(0:,0:,0:), target :: vel1,vel2
    real(rp), intent(inout), dimension(1-nh:,1-nh:,0:), target :: bcvel1,bcvel2
    real(rp), intent(in), dimension(1-nh:,1-nh:,0:), target :: bcvel1_mag,bcvel2_mag
    real(rp), pointer, dimension(:,:,:) :: u,v,w
    real(rp), pointer, dimension(:,:,:) :: bcu,bcv,bcw
    real(rp), pointer, dimension(:,:,:) :: bcu_mag,bcv_mag,bcw_mag
    real(rp) :: sgn,wei,coef,uh,vh,wh,u1,u2,v1,v2,w1,w2,u_mag,v_mag,w_mag,tauw(2),visci
    integer  :: i,j,k,i1,i2,j1,j2,k1,k2
    !
    visci = 1._rp/visc
    !
    select case(idir)
    case(1)
      if(ibound==0) then
        i2 = index
        i1 = index-1
        coef = (h-(i1-0.5_rp)*dl(1))/dl(1)
        sgn = 1._rp
      else
        i2 = index
        i1 = index+1
        coef = (h-(n(1)-i1+0.5_rp)*dl(1))/dl(1)
        sgn = -1._rp
      end if
      v       => vel1
      w       => vel2
      bcv     => bcvel1
      bcw     => bcvel2
      bcv_mag => bcvel1_mag
      bcw_mag => bcvel2_mag
      !$acc parallel loop collapse(2) default(present)  private(v1,v2,w1,w2,v_mag,w_mag,vh,wh,tauw) async(1)
      do k = 1,n(3)
        do j = 0,n(2)
          v1    = v(i1,j,k)
          v2    = v(i2,j,k)
          w1    = 0.25_rp*(w(i1,j,k) + w(i1,j+1,k) + w(i1,j,k-1) + w(i1,j+1,k-1))
          w2    = 0.25_rp*(w(i2,j,k) + w(i2,j+1,k) + w(i2,j,k-1) + w(i2,j+1,k-1))
          v_mag = bcv_mag(j,k,ibound)
          w_mag = 0.25_rp*(bcw_mag(j,k  ,ibound) + bcw_mag(j+1,k  ,ibound) + &
                           bcw_mag(j,k-1,ibound) + bcw_mag(j+1,k-1,ibound))
          vh    = vel_relative(v1,v2,coef,v_mag)
          wh    = vel_relative(w1,w2,coef,w_mag)
          call wallmodel(mtype,vh,wh,h,l(1),visc,tauw)
          bcv(j,k,ibound) = sgn*visci*tauw(1)
        end do
      end do
      !$acc parallel loop collapse(2) default(present)  private(v1,v2,w1,w2,v_mag,w_mag,vh,wh,tauw,wei) async(1)
      do k = 0,n(3)
        do j = 1,n(2)
          wei   = (zf(k)-zc(k))/dzc(k)
          v1    = 0.5_rp*((1._rp-wei)*(v(i1,j-1,k) + v(i1,j,k)) + wei*(v(i1,j-1,k+1) + v(i1,j,k+1)))
          v2    = 0.5_rp*((1._rp-wei)*(v(i2,j-1,k) + v(i2,j,k)) + wei*(v(i2,j-1,k+1) + v(i2,j,k+1)))
          w1    = w(i1,j,k)
          w2    = w(i2,j,k)
          v_mag = 0.5_rp*((1._rp-wei)*(bcv_mag(j-1,k  ,ibound) + bcv_mag(j,k  ,ibound)) + &
                                 wei *(bcv_mag(j-1,k+1,ibound) + bcv_mag(j,k+1,ibound)))
          w_mag = bcw_mag(j,k,ibound)
          vh    = vel_relative(v1,v2,coef,v_mag)
          wh    = vel_relative(w1,w2,coef,w_mag)
          call wallmodel(mtype,vh,wh,h,l(1),visc,tauw)
          bcw(j,k,ibound) = sgn*visci*tauw(2)
        end do
      end do
    case(2)
      if(ibound==0) then
        j2 = index
        j1 = index-1
        coef = (h-(j1-0.5)*dl(2))/dl(2)
        sgn = 1._rp
      else
        j2 = index
        j1 = index+1
        coef = (h-(n(2)-j1+0.5)*dl(2))/dl(2)
        sgn = -1._rp
      end if
      u       => vel1
      w       => vel2
      bcu     => bcvel1
      bcw     => bcvel2
      bcu_mag => bcvel1_mag
      bcw_mag => bcvel2_mag
      !$acc parallel loop collapse(2) default(present)  private(u1,u2,w1,w2,u_mag,w_mag,uh,wh,tauw) async(1)
      do k = 1,n(3)
        do i = 0,n(1)
          u1    = u(i,j1,k)
          u2    = u(i,j2,k)
          w1    = 0.25_rp*(w(i,j1,k) + w(i+1,j1,k) + w(i,j1,k-1) + w(i+1,j1,k-1))
          w2    = 0.25_rp*(w(i,j2,k) + w(i+1,j2,k) + w(i,j2,k-1) + w(i+1,j2,k-1))
          u_mag = bcu_mag(i,k,ibound)
          w_mag = 0.25_rp*(bcw_mag(i,k  ,ibound) + bcw_mag(i+1,k  ,ibound) + &
                           bcw_mag(i,k-1,ibound) + bcw_mag(i+1,k-1,ibound))
          uh    = vel_relative(u1,u2,coef,u_mag)
          wh    = vel_relative(w1,w2,coef,w_mag)
          call wallmodel(mtype,uh,wh,h,l(2),visc,tauw)
          bcu(i,k,ibound) = sgn*visci*tauw(1)
        end do
      end do
      !$acc parallel loop collapse(2) default(present)  private(u1,u2,w1,w2,u_mag,w_mag,uh,wh,tauw,wei) async(1)
      do k = 0,n(3)
        do i = 1,n(1)
          wei   = (zf(k)-zc(k))/dzc(k)
          u1    = 0.5_rp*((1._rp-wei)*(u(i-1,j1,k) + u(i,j1,k)) + wei*(u(i-1,j1,k+1) + u(i,j1,k+1)))
          u2    = 0.5_rp*((1._rp-wei)*(u(i-1,j2,k) + u(i,j2,k)) + wei*(u(i-1,j2,k+1) + u(i,j2,k+1)))
          w1    = w(i,j1,k)
          w2    = w(i,j2,k)
          u_mag = 0.5_rp*((1._rp-wei)*(bcu_mag(i-1,k  ,ibound) + bcu_mag(i,k  ,ibound)) + &
                                 wei *(bcu_mag(i-1,k+1,ibound) + bcu_mag(i,k+1,ibound)))
          w_mag = bcw_mag(i,k,ibound)
          uh    = vel_relative(u1,u2,coef,u_mag)
          wh    = vel_relative(w1,w2,coef,w_mag)
          call wallmodel(mtype,uh,wh,h,l(2),visc,tauw)
          bcw(i,k,ibound) = sgn*visci*tauw(2)
        end do
      end do
    case(3)
      if(ibound==0) then
        k2 = index
        k1 = index-1
        coef = (h-zc(k1))/dzc(k1)
        sgn = 1._rp
      else
        k2 = index
        k1 = index+1
        coef = (h-(l(3)-zc(k1)))/(dzc(k2))
        sgn = -1._rp
      end if
      u       => vel1
      v       => vel2
      bcu     => bcvel1
      bcv     => bcvel2
      bcu_mag => bcvel1_mag
      bcv_mag => bcvel2_mag
      !$acc parallel loop collapse(2) default(present)  private(u1,u2,v1,v2,u_mag,v_mag,uh,vh,tauw) async(1)
      do j = 1,n(2)
        do i = 0,n(1)
          u1    = u(i,j,k1)
          u2    = u(i,j,k2)
          v1    = 0.25_rp*(v(i,j,k1) + v(i+1,j,k1) + v(i,j-1,k1) + v(i+1,j-1,k1))
          v2    = 0.25_rp*(v(i,j,k2) + v(i+1,j,k2) + v(i,j-1,k2) + v(i+1,j-1,k2))
          u_mag = bcu_mag(i,j,ibound)
          v_mag = 0.25_rp*(bcv_mag(i,j  ,ibound) + bcv_mag(i+1,j  ,ibound) + &
                           bcv_mag(i,j-1,ibound) + bcv_mag(i+1,j-1,ibound))
          uh    = vel_relative(u1,u2,coef,u_mag)
          vh    = vel_relative(v1,v2,coef,v_mag)
          call wallmodel(mtype,uh,vh,h,l(3),visc,tauw)
          bcu(i,j,ibound) = sgn*visci*tauw(1)
        end do
      end do
      !$acc parallel loop collapse(2) default(present)  private(u1,u2,v1,v2,u_mag,v_mag,uh,vh,tauw) async(1)
      do j = 0,n(2)
        do i = 1,n(1)
          u1    = 0.25_rp*(u(i-1,j,k1) + u(i,j,k1) + u(i-1,j+1,k1) + u(i,j+1,k1))
          u2    = 0.25_rp*(u(i-1,j,k2) + u(i,j,k2) + u(i-1,j+1,k2) + u(i,j+1,k2))
          v1    = v(i,j,k1)
          v2    = v(i,j,k2)
          u_mag = 0.25_rp*(bcu_mag(i-1,j  ,ibound) + bcu_mag(i,j  ,ibound) + &
                           bcu_mag(i-1,j+1,ibound) + bcu_mag(i,j+1,ibound))
          v_mag = bcv_mag(i,j,ibound)
          uh    = vel_relative(u1,u2,coef,u_mag)
          vh    = vel_relative(v1,v2,coef,v_mag)
          call wallmodel(mtype,uh,vh,h,l(3),visc,tauw)
          bcv(i,j,ibound) = sgn*visci*tauw(2)
        end do
      end do
    end select
  end subroutine cmpt_wallmodelbc
  !
  function vel_relative(v1,v2,coef,bcv_mag)
    !
    ! compute relative velocity to a wall
    !
    implicit none
    real(rp), intent(in) :: v1,v2,coef,bcv_mag
    real(rp) :: vel_relative
    !
    !$acc routine seq
    vel_relative = (1._rp-coef)*v1 + coef*v2
    vel_relative = vel_relative - bcv_mag
  end function vel_relative
  !
  subroutine wallmodel(mtype,uh,vh,h,l1d,visc,tauw)
    !
    ! Newton-Raphson, 3~7 iters. It is bad to initialize with an assumed linear profile
    ! below the first cell center, because its slope is quite uncertain due to the
    ! erroneous velocity at the first cell center. We do not save tauw_tot from previous
    ! step, so those values are not directly available. At zero velocity,
    ! utau=visc/h*exp(-kap_log*b_log). It is necessary to use abs to avoid negative values
    ! of utau. Note that arrays (u/v/w) are initialized as zero (ghost points), which
    ! makes a reasonable guess at end points. The convergence criterion is
    ! |tauw/tauw_old-1.0|<1.0e-4, corresponding to |utau/utau_old-1.0|<0.5e-4.
    ! For a channel, it is efficient to use the computed utau available at the nearest
    ! grid point, which helps reduce ~50% of the iterations. However, "if" statements
    ! must be introduced at special points in square duct/cavity.
    !
    implicit none
    integer, parameter :: WM_LAM = -1, &
                          WM_LOG =  1
    integer, intent(in)  :: mtype
    real(rp), intent(in) :: uh,vh,h,l1d,visc
    real(rp), intent(out), dimension(2) :: tauw
    real(rp) :: upar,utau,f,fp,conv,utau_old,tauw_tot
    real(rp) :: umax,del
    !
    !$acc routine seq
    select case(mtype)
    case(WM_LOG)
      conv = 1._rp
      upar = sqrt(uh*uh+vh*vh)
      utau = max(sqrt(upar/h*visc),visc/h*exp(-kap_log*b_log))
      do while(conv>0.5e-4_rp)
        utau_old = utau
        f  = upar/utau-1._rp/kap_log*log(h*utau/visc)-b_log
        fp = -1._rp/utau*(upar/utau+1._rp/kap_log)
        utau = abs(utau-f/fp)
        conv = abs(utau/utau_old-1._rp)
      end do
      tauw_tot = utau*utau
      tauw(1) = tauw_tot*uh/(upar+eps)
      tauw(2) = tauw_tot*vh/(upar+eps)
    case(WM_LAM)
      upar = sqrt(uh*uh+vh*vh)
      del  = 0.5_rp*l1d
      umax = upar/(h/del*(2._rp-h/del))
      tauw_tot = 2._rp/del*umax*visc
      tauw(1) = tauw_tot*uh/(upar+eps)
      tauw(2) = tauw_tot*vh/(upar+eps)
    end select
  end subroutine wallmodel
end module mod_wmodel
#endif