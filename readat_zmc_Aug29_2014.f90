!------------------------------------------------------------------------------!
!          Program proper starts here                                          !
!------------------------------------------------------------------------------!

subroutine readat(LLL,isite,cat,inum,xxx,yyy,zzz)

  !------------------------------------------------------------------------------!
  !          Use these modules -- requires compile flag -lmodules                !
  !------------------------------------------------------------------------------!

  use cmdline_arguments
  Use file_functions, only: open, exists, read_buffer, stderr, stdout
  use iso_varying_string
  use string_functions
  use variable_array

  implicit none

  !------------------------------------------------------------------------------!
  !          These variables are for dovetailing with diffuse.f                  !
  !------------------------------------------------------------------------------!

  real, intent(out)                               :: xxx,yyy,zzz
  INTEGER,intent(out)                             :: inum
  INTEGER,dimension(3),intent(in)                 :: LLL
  CHARACTER(len=4),intent(in)                     :: cat
  INTEGER,intent(in)                              :: isite

  !------------------------------------------------------------------------------!
  !          These variables are internal to readat, though some must be saved   !
  !------------------------------------------------------------------------------!

  real, allocatable :: x(:,:,:,:), y(:,:,:,:), z(:,:,:,:)
  CHARACTER(len=2)                                :: ty1, ty2
  character(len=2),allocatable ::occc(:,:,:,:)
  INTEGER                                         :: itest,ia,ib,ic
  integer, allocatable :: look_up(:,:,:,:)
  integer :: sitenum

  !------------------------------------------------------------------------------!
  !                      End of declarations                                     !
  !------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------!
  !                      Need to save a couple of arrays which we set up on      !
  !       first pass through the readat.f
  !------------------------------------------------------------------------------!

  save  x, y, z,itest, occc, look_up

  !------------------------------------------------------------------------------!
  !      If positions already established, skip over this part                   !
  !------------------------------------------------------------------------------!

  IF (itest.eq.1234) goto 657

  call getcoords()

657 itest = 1234

  inum = 0
  ia=LLL(1)
  ib=LLL(2)
  ic=LLL(3)

  sitenum = look_up(ia,ib,ic,isite) 

  if(sitenum.eq.0) goto 156  

  ty1 = occc(ia,ib,ic,sitenum)
  ty2 = cat(1:2)

  if(ty1.ne.ty2) goto 156
  ! Now, only if we get through all these checks do we reset inum to 1 and give the xyz's

  xxx = x(ia,ib,ic,sitenum)
  yyy = y(ia,ib,ic,sitenum)
  zzz = z(ia,ib,ic,sitenum)
  inum = 1

156 return


contains

  !------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------!

  subroutine getcoords()

    type (varying_string) :: myoptions(3)
    integer :: ierr, dunit, simsize(3),possible_sites,occupied_sites,i,j,ia,ib,ic, id
    character(len=50)  :: cver
    logical            :: quiet

    character(len=100)                   :: outname 

    !------------------------------------------------------------------------------!
    ! These are our accepted command line options (see subroutine usage for
    ! an explanation)
    !------------------------------------------------------------------------------!

    myoptions(1) = 'help'
    myoptions(2) = 'quiet' 
    myoptions(3) = 'version'

    ! Version string.  Update when change program, maintain its length
    cver = 'August 29 2014                                    '


    ! This call parses the command line arguments for command line options

    call get_options(myoptions, ierr)

    ! Check we weren't passed duff options -- spit the dummy if we were
    if (ierr > 0) then
       write(stderr,*) 'ERROR! Unknown option(s): ',join(bad_options()," ")
       call usage()
       stop
    end if
    if((option_exists('version')).or.(option_exists('help')))then
       if (option_exists('version')) then
          write(stdout,*)'----------------------------------------------------------------------------'
          write(stdout,*)'ReadatZMC version: ',cver
          write(stdout,*)'----------------------------------------------------------------------------'
       end if
       if (option_exists('help')) call usage()
       stop
    end if

    quiet = .FALSE.
    if (option_exists('quiet')) quiet = .TRUE.

    !------------------------------------------------------------------------------!
    !                      Initialise keyword list.                                !
    !------------------------------------------------------------------------------!


    if ( .not. have_args() ) then
       write(stderr,*)'----------------------------------------------------------------------------'
       write(stderr,*)'Must provide filename for reading in the crystal.'
       write(stderr,*)'It is probably a filename of the form forename.diffuse.'
       write(stderr,*)'Program exiting.'             
       write(stderr,*)'----------------------------------------------------------------------------'
       stop
    end if   ! arg test
    outname  = next_arg()  

    !dunit = freeunit()
    !open(unit=dunit, file=outname)

    dunit = open(outname,status='old')

    read(dunit,*)simsize(1),simsize(2),simsize(3),possible_sites,occupied_sites
    if(.not.quiet) then
       write(stdout,*)'----------------------------------------------------------------------------'
       write(stdout,*)'Reading in MC crystal from file ',trim(outname)
       write(stdout,*)'Crystal dimensions are: ',simsize
       write(stdout,*)'There are ',possible_sites,' sites in the unit cell.'
       write(stdout,*)'Of which at most ',occupied_sites,' can be occupied at once.'
       write(stdout,*)'----------------------------------------------------------------------------'
    end if
    allocate(x(simsize(1),simsize(2),simsize(3),occupied_sites),stat=ierr)
    if (ierr.ne.0) call allerr('x(simsize(1),simsize(2),simsize(3),occupied_sites)    ')
    allocate(y(simsize(1),simsize(2),simsize(3),occupied_sites),stat=ierr)
    if (ierr.ne.0) call allerr('y(simsize(1),simsize(2),simsize(3),occupied_sites)    ')
    allocate(z(simsize(1),simsize(2),simsize(3),occupied_sites),stat=ierr)
    if (ierr.ne.0) call allerr('z(simsize(1),simsize(2),simsize(3),occupied_sites)    ')

    allocate(occc(simsize(1),simsize(2),simsize(3),occupied_sites),stat=ierr)
    if (ierr.ne.0) call allerr('occc(simsize(1),simsize(2),simsize(3),occupied_sites)')

    allocate(look_up(simsize(1),simsize(2),simsize(3),possible_sites),stat=ierr)
    if (ierr.ne.0) call allerr('look_up(....)          ')

    look_up = 0


    do i = 1, simsize(1)*simsize(2)*simsize(3)
       read(dunit,*)ia,ib,ic
       do j = 1, possible_sites
          read(dunit,*)id,look_up(ia,ib,ic,id)
       end do
    end do

    do i = 1, simsize(1)*simsize(2)*simsize(3)
       read(dunit,*)ia,ib,ic
       do j = 1, maxval(look_up(ia,ib,ic,:))
          read(dunit,*)x(ia,ib,ic,j),y(ia,ib,ic,j),z(ia,ib,ic,j),occc(ia,ib,ic,j)
       end do
    end do

    close(dunit)

    return 

  end subroutine getcoords

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine allerr(xyzfile)          

    character(len=*)                             :: xyzfile

    write(stderr,*)'----------------------------------------------------------------------------'
    write(stderr,*) 'Problem allocating memory to ',trim(xyzfile)
    write(stderr,*) 'Do you have enough for the simulation size you want?'
    write(stderr,*) 'Exiting.'
    write(stderr,*)'----------------------------------------------------------------------------'

    stop

  end subroutine allerr

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine usage()  

   use file_functions, only: stderr

    write(stderr,*)   
    write(stderr,*)'|---------------------------------------------------------------------|'
    write(stderr,*)'| Usage:                                                              |'
    write(stderr,*)'|                                                                     |'
    write(stderr,*)'| DZMC [--option_1] [--option_2] ... [--option_n] filename            |'
    write(stderr,*)'|                                                                     |'
    write(stderr,*)'| filename is a file of atom coordinates put out by the --diffuse     |'
    write(stderr,*)'| option of ZMC.                                                      |'
    write(stderr,*)'|                                                                     |'
    write(stderr,*)'| --quiet   Sends nothing to the screen except error messages.        |'
    write(stderr,*)'|                                                                     |'
    write(stderr,*)'| --help    Prints this information and exits.                        |'
    write(stderr,*)'|                                                                     |'
    write(stderr,*)'| --version Prints version information and exits.                     |'
    write(stderr,*)'|                                                                     |'
    write(stderr,*)'|---------------------------------------------------------------------|'

  end subroutine usage

  !------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------!

end subroutine readat

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

