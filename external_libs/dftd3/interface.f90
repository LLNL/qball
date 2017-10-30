!
! Copyright (C) 2017 Xavier Andrade
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

!This is a bit ugly, as we do not export the object but we keep a
!global variable inside a module. This is to avoid the complication of
!passing a derived datatype between Fortran and C.

module dftd3_data
  use dftd3_api
  
  implicit none

  type(dftd3_calc) :: calc
  
end module dftd3_data

! ---------------------------------------

subroutine f90_dftd3_init()
  use dftd3_data

  implicit none

  type(dftd3_input) :: input
  
  call dftd3_init(calc, input)
  call dftd3_set_functional(calc, func = 'pbe', version = 4, tz = .false.)
  
end subroutine f90_dftd3_init

! ---------------------------------------

subroutine f90_dftd3_end()
  implicit none
  
end subroutine f90_dftd3_end

! ---------------------------------------

subroutine f90_dftd3_pbc_dispersion(natoms, coords, izp, latvecs, disp, grads, stress)
  use dftd3_data

  implicit none

  integer(4), intent(in)  :: natoms
  real(8),    intent(in)  :: coords(1:3, 1:natoms)
  integer(4), intent(in)  :: izp(1:natoms)
  real(8),    intent(in)  :: latvecs(1:3, 1:3)
  real(8),    intent(out) :: disp
  real(8),    intent(out) :: grads(1:3, 1:natoms)
  real(8),    intent(out) :: stress(1:3, 1:3)

  call dftd3_pbc_dispersion(calc, coords, izp, latvecs, disp, grads, stress)
  
end subroutine f90_dftd3_pbc_dispersion
