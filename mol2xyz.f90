program mol2xyz

  use fundamental_constants
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use mol2_class
  use file_functions
  use string_functions, only: real
  use crystallography_class

  implicit none

  type (mol2_object)      :: mol2

  type(crystallography_object) :: xtal

  character(len=6), allocatable :: atomlabels(:)
  real, allocatable :: mycoords(:,:)

  integer :: i, j

  character(len=100) ::  fname

  if ( .not. have_args() ) then
     write(stderr,*) "Must provide mol2 file(s) as command line argument(s)"
     stop
  end if

  do while (have_args())

     ! Grab the file name from the command line
     fname = next_arg()
     
     mol2 = fname

     allocate(mycoords(3,bond_num(mol2)), atomlabels(atom_num(mol2)))

     mycoords = coords(mol2) 
     atomlabels = labels(mol2)

     write(*,'(I0)') atom_num(mol2)
     do i = 1,atom_num(mol2)
        write(*,'(A,3F10.4)') atomlabels(i), mycoords(:,i)
     end do

     deallocate(mycoords, atomlabels)

  end do

end program mol2xyz

