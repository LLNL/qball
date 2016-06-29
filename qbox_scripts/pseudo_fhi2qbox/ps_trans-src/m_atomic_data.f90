module m_atomic_data
  use m_dft_op_types_def
  save
  private
  integer(kind=INTG), parameter :: max_atom = 200

  type atomic_data
     character(len=2) :: name
     real(kind=DP) :: mass
  end type atomic_data

  type(atomic_data) ::  atom_data(max_atom)



  public get_mass_from_table, init_mass_data
contains
  function get_mass_from_table(atm)
    implicit none
    character(len=*), intent(in) :: atm
    integer(kind=INTG) i
    real(kind=DP) get_mass_from_table
    do i = 1, ubound(atom_data,1)
       if(atom_data(i)%name == atm) then
!          print *,'atom ',atm,' ',atom_data(i)%mass
          get_mass_from_table = atom_data(i)%mass
       end if
    end do
  end function get_mass_from_table
  subroutine init_mass_data()
    implicit none
    atom_data(1) = atomic_data('H', 1.0079400d0)
    atom_data(2) = atomic_data('He', 4.0002602d0)
    atom_data(3) = atomic_data('Li', 6.9410000d0)
    atom_data(4) = atomic_data('Be', 9.0121800d0)
    atom_data(5) = atomic_data('B',  10.8220000d0)
    atom_data(6) = atomic_data('C', 12.0110000d0)
    atom_data(7) = atomic_data('N', 14.0076000d0)
    atom_data(8) = atomic_data('O', 15.9994000d0)
    atom_data(9) = atomic_data('F', 18.9984030d0)
    atom_data(10) = atomic_data('Ne', 20.1790000d0)
    atom_data(11) = atomic_data('Na', 22.9897700d0)
    atom_data(12) = atomic_data('Mg', 24.3050000d0)
    atom_data(13) = atomic_data('Al', 26.9815400d0)
    atom_data(14) = atomic_data('Si', 28.0860000d0)
    atom_data(15) = atomic_data('P', 30.9737600d0)
    atom_data(16) = atomic_data('S', 32.0660000d0)
    atom_data(17) = atomic_data('Cl', 35.4530000d0)
    atom_data(18) = atomic_data('Ar', 39.9480000d0)
    atom_data(19) = atomic_data('K', 39.0983000d0)
    atom_data(20) = atomic_data('Ca', 40.0780000d0)
    atom_data(21) = atomic_data('Sc', 44.9559100d0)
    atom_data(22) = atomic_data('Ti', 47.8800000d0)
    atom_data(23) = atomic_data('V', 47.8800000d0)
    atom_data(24) = atomic_data('Cr', 51.9961000d0)
    atom_data(25) = atomic_data('Mn', 54.9380000d0)
    atom_data(26) = atomic_data('Fe', 55.8470000d0)
    atom_data(27) = atomic_data('Co', 58.9332000d0)
    atom_data(28) = atomic_data('Ni', 58.6900000d0)
    atom_data(29) = atomic_data('Cu', 63.5460000d0)
    atom_data(30) = atomic_data('Zn', 65.3900000d0)
    atom_data(31) = atomic_data('Ga', 69.7230000d0)
    atom_data(32) = atomic_data('Ge', 72.5900000d0)
    atom_data(33) = atomic_data('As', 74.9216000d0)
    atom_data(34) = atomic_data('Se', 78.9600000d0)
    atom_data(35) = atomic_data('Br', 79.9040000d0)
    atom_data(36) = atomic_data('Kr', 83.8000000d0)
    atom_data(37) = atomic_data('Rb', 85.4678000d0)
    atom_data(38) = atomic_data('Sr', 87.6200000d0)
    atom_data(39) = atomic_data('Y', 88.9059000d0)
    atom_data(40) = atomic_data('Zr', 91.2240000d0)
    atom_data(41) = atomic_data('Nb', 92.9064000d0)
    atom_data(42) = atomic_data('Mo', 95.9400000d0)
    !  atom_data(43) = atomic_data('Tc', 0.0000000d0)
    atom_data(44) = atomic_data('Ru',101.0700000d0)
    atom_data(45) = atomic_data('Rh', 102.9055000d0)
    atom_data(46) = atomic_data('Pd', 106.4200000d0)
    atom_data(47) = atomic_data('Ag', 107.8682000d0)
    atom_data(48) = atomic_data('Cd', 112.4100000d0)
    atom_data(49) = atomic_data('In', 114.8200000d0)
    atom_data(50) = atomic_data('Sn', 118.7100000d0)
    atom_data(51) = atomic_data('Sb', 121.7500000d0)
    atom_data(52) = atomic_data('Te', 127.6000000d0)
    atom_data(53) = atomic_data('I', 126.9045000d0)
    atom_data(54) = atomic_data('Xe', 131.2900000d0)
    atom_data(55) = atomic_data('Mu', 0.1134940d0)
    atom_data(56) = atomic_data('BM', 28.0860000d0)
    atom_data(101) = atomic_data('D', 2.01400d0) !from BUTSURIGAKUJITEN p892
    atom_data(102) = atomic_data('T', 3.01605d0) !from BUTSURIGAKUJITEN p892
  end subroutine init_mass_data
end module m_atomic_data
