program catmol2

  use mol2_class
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use iso_varying_string
  use file_functions, only: stderr
  use string_functions, only: join

  implicit none

  !! $Log: catmol2.f90,v $
  !! Revision 1.1  2006/05/18 01:30:26  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: catmol2.f90,v 1.1 2006/05/18 01:30:26 aidan Exp aidan $"

  type (mol2_object)              :: all, input

  integer :: error
  
  type (varying_string) ::  myoptions(1)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'

  ! This call parses the command line arguments for command line options
  call get_options(myoptions, error)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (error > 0) then
!     write(stderr,*) 'ERROR! Unknown options: ',join(bad_options()," ")
     write(stderr,*) 'ERROR! Unknown options: '
     call usage
     STOP
  end if

  ! Check if we just want to print the usage
  if (option_exists('help')) then
     call usage
     STOP
  end if

  if ( .not. have_args() ) then
     write(stderr,*) "Must provide mol2 filenames as command line arguments"
     call usage
     stop
  end if

  all = char(next_arg())

  do while (have_args())

     input = char(next_arg())

     call add(all, input)

  end do

  call print(all)

contains

  subroutine usage
    
    write(stderr,*)
    write(stderr,*) 'catmol2 takes the atoms from a series of mol2 files and places them '
    write(stderr,*) 'in a single mol2 file, which is printed to standard output '
    write(stderr,*)
    write(stderr,*) 'Usage: catmol2 [OPTIONS] <mol2 file> <mol2 file> ...'
    write(stderr,*)
    write(stderr,*) '  Options:'
    write(stderr,*)
    write(stderr,*) '  --help          - print this message'
    write(stderr,*)

  end subroutine usage

end program catmol2
