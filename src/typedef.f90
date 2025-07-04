! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
#if defined(_LES)
module mod_typedef
  use mod_types, only: rp,sp,dp,i8,MPI_REAL_RP
  type bound
  real(rp), allocatable, dimension(:,:,:) :: x
  real(rp), allocatable, dimension(:,:,:) :: y
  real(rp), allocatable, dimension(:,:,:) :: z
  end type bound
end module mod_typedef
#endif