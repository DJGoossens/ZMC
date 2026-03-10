
!  This is ZMC, Monte Carlo program for modelling diffuse scattering from 
!  crystals.  Version is in string variable cver.
!
!  Detailed development history is in file zmc_documentation.txt
!
!  Written by D.J.Goossens, goossens@rsc.anu.edu.au
!                           darren.goossens@gmail.com
!  
!  Use and development by others is encouraaged.
!
!  Uses module library provided by A.P.Heerdegen, aidan@rsc.anu.edu.au
!  (modified in an ad hoc way to compile with newer gfortran in 2026 by DJG)
!
!  All uses of 'write(6,*)' are for debugging and should usually be commented
!  out or not present at all. Program uses print or write(stdout/stderr...); 
!  the write(6..) format search for.
!
!  Learned on FORTRAN 77. Though this program is Fortran 95, roughly, a lot of 
!  my coding habits are older, and those of someone who never learned formally.
!
!  I am sure much could be done better, but the algorithm seems to work...
!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! subroutines are declared in two modules at the front of the program
!------------------------------------------------------------------------------!

module divs
  ! Integer Parameters
  integer, parameter:: div_max = 50    !max number of bins for iacc/irej/iup
end module divs

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module zmc_subroutines

  use crystallography_class
  use file_functions, only: open, exists, stderr, stdout, read_buffer, stdin
  use string_functions
  use keyword_class, only: exists, get_value, has_value, initialise, keyword_object, assignment(=),read, num_value
  use zmatrix_class, only: assignment(=),zmatrix_object,load_zmatrix,copy,as_xyz,as_cartesian,print,parameter,labels,num

  implicit none

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

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

  subroutine file_does_not_exist(filename)

    character(len=100), intent(in) ::   filename

    write(stderr,*)'----------------------------------------------------------------------------'
    write(stderr,*)'File does not exist:    ',trim(filename)
    write(stderr,*)'Program exiting.'
    write(stderr,*)'----------------------------------------------------------------------------'
    stop

  end subroutine file_does_not_exist

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  subroutine file_does_exist(filename)

    character(len=100), intent(in) ::   filename
    character(len=100) ::   fi2
    character(len=1) :: response

    fi2 = ''
    fi2 = trim(filename)
    write(stdout,*)'File already exists:    ',trim(fi2)
482 response = ''
    write(stdout,'(a21)',advance='no')' Overwrite? (y or n) '
    read(stdin,*)response
    if((response.eq.'n').or.(response.eq.'N')) stop
    if((response.ne.'y').and.(response.ne.'Y')) goto 482
    write(stdout,*)'----------------------------------------------------------------------------'
    return

  end subroutine file_does_exist

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine fatal_keyword(keywrd, buffer)

    character(len=*),intent(in) :: keywrd
    character(len=*),intent(in) :: buffer

    write(stderr,*)'----------------------------------------------------------------------------'
    write(stderr,*)'Problem reading value of ',trim(keywrd)
    write(stderr,*)'Offending line is: ',trim(buffer)
    write(stderr,*)'Program exiting.'
    write(stderr,*)'----------------------------------------------------------------------------'
    stop

    return

  end  subroutine fatal_keyword

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_optional_integer(value_given,keyword,key,quiet)

    character(len=500)                :: buffer
    integer, intent(out)              :: value_given
    integer                           :: ierr
    character(len=*),intent(in)       :: keyword
    type (keyword_object)             :: key
    logical, intent(in)               :: quiet 


    value_given = -9999
    if (exists(key, keyword)) then
       if (has_value(key, keyword)) then
          ! Grab 1 integer 
          buffer=' '
          buffer = get_value(key,keyword)
          read(buffer,*,iostat=ierr)value_given
          if(ierr.ne.0)then
             write(stderr,*)'-----------------------------------------------'
             write(stderr,*)'Illegal input in field ',keyword
             write(stderr,*)'Offending input is: ',trim(buffer)
             write(stderr,*)'Program exiting.'
             write(stderr,*)'-----------------------------------------------'
             stop
          else
             !              make sure they are acceptable numbers.
             if(value_given.lt.1)then
                if (.not. quiet)write(stdout,*)'Value for ',keyword,' is zero or negative.'
                if (.not. quiet)write(stdout,*)'Is this OK?'
                if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
             end if
             if (.not. quiet)write(stdout,*)'Value for ',keyword,' specified explicitly is: ',value_given
             if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
          end if
       else
          if (.not. quiet)write(stdout,*)keyword,' keyword given but it has no value!'
          if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
       end if
    else
       if (.not. quiet)write(stdout,*)'Value for ',keyword,' not specified explicitly.'
       if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if

    return

  end subroutine get_optional_integer

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------


  subroutine get_simsize(simsize, key,quiet)

    character(len=500)                 :: buffer
    integer,dimension(3), intent(out)  :: simsize
    integer                            :: ierr
    type (keyword_object)              :: key
    logical, intent(in)                :: quiet 
    !      
    if (exists(key, 'CRYSTAL')) then
       if (has_value(key, 'CRYSTAL')) then
          buffer=' '
          buffer = get_value(key,'CRYSTAL')
          read(buffer,*,iostat=ierr)simsize(1),simsize(2),simsize(3)
          if(ierr.ne.0)then
             write(stderr,*)'-----------------------------------------------'
             write(stderr,*)'Illegal input in field CRYSTAL'
             write(stderr,*)'Offending input is: ',trim(buffer)
             write(stderr,*)'Program exiting.'
             write(stderr,*)'-----------------------------------------------'
             stop
          else
             !              make sure they are acceptable numbers.
             if(minval(simsize).lt.1)then
                write(stderr,*)'-----------------------------------------------'
                write(stderr,*)'One or more model crystal dimensions invalid.'
                write(stderr,*)'Should give three integer, all larger than zero.'
                write(stderr,*)'-----------------------------------------------'
                stop
             end if
          end if
       else
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'One or more model crystal dimensions invalid.'
          write(stderr,*)'Should give three integer, all larger than zero.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if
    else
       write(stderr,*)'-----------------------------------------------'
       write(stderr,*)'No CRYSTAL keyword given.'
       write(stderr,*)'Model CRYSTAL size in unit cells must be given.'
       write(stderr,*)'-----------------------------------------------'
       stop
    end if
    if (.not. quiet)write(stdout,'(" Crystal size is: ",3i6)')simsize
    if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'

    return

  end subroutine get_simsize

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_cell(xtal,key,quiet)

    character(len=500)     ::  buffer
    type (crystallography_object),intent(out)                           :: xtal
    integer :: ierr, i
    real :: cellpar(6)
    type (keyword_object)              :: key
    logical, intent(in)                :: quiet 

    if (exists(key, 'CELL')) then
       if (has_value(key, 'CELL')) then
          buffer=' '
          buffer = get_value(key,'CELL')
          read(buffer,*,iostat=ierr)(cellpar(i),i=1,6)
          if(ierr.ne.0)then
             write(stderr,*)'-----------------------------------------------'
             write(stderr,*)'Illegal input in field CELL'
             write(stderr,*)'Offending input is: ',trim(buffer)
             write(stderr,*)'Program exiting.'
             write(stderr,*)'-----------------------------------------------'
             stop
          end if
       end if
    else
       write(stderr,*)'-----------------------------------------------'
       write(stderr,*)'No unit cell specified!'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'-----------------------------------------------'
    end if

    xtal = cellpar(:)
    if (.not. quiet)write(stdout,'(1x,"Cell parameters are: ",6f8.4)')(cellpar(i),i=1,6)
    if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'

    return

  end  subroutine get_cell

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_total_zm_files(total_num_zm,value_given,key,quiet)

    integer, intent(in)   ::  value_given
    integer, intent(out)  ::  total_num_zm
    type (keyword_object)              :: key
    logical, intent(in)                :: quiet 

    ! See how many zmatrix keywords there are in the keyword file:
    total_num_zm = num_value(key, 'ZMATFILE')
    if (total_num_zm.lt.1) then
       write(stderr,*)'-----------------------------------------------'
       write(stderr,*)'At least one z-matrix must be specified.'
       write(stderr,*)'-----------------------------------------------'
       stop
    end if
    if(value_given.ne.-9999) then
       ! Do comparison
       if(value_given.ne.total_num_zm) then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'NUMZMATS specified ',value_given,' z-matrices to be supplied,'
          write(stderr,*)'but ',total_num_zm,' z-matrix filename(s) given.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if
    end if
    if (.not. quiet)write(stdout,'(" Number of z-matrix files to read: ",i5)')total_num_zm

    return

  end subroutine get_total_zm_files

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_zmats(zmat,total_num_zm,num_at,key,quiet)

    character(len=100)                 ::  zname
    character(len=500)                 ::  buffer
    integer, intent(in)                ::  total_num_zm
    type (zmatrix_object),intent(out)  :: zmat(:)
    integer                            :: i,ierr
    integer, intent(out)               :: num_at(:)
    integer, allocatable               :: zindex(:)
    type (keyword_object)              :: key
    logical, intent(in)                :: quiet 

    allocate(zindex(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('zindex(total_num_zm) ')
    if(num_value(key,'ZMATFILE').lt.1) then
       write(stderr,*)'-----------------------------------------------'
       write(stderr,*)'Must give at least one z-matrix file.'
       write(stderr,*)'Specify using "ZMATFILE" keyword.'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'-----------------------------------------------'
       stop
    end if
    do i = 1,total_num_zm
       buffer = ' '
       buffer = get_value(key,'ZMATFILE',i)
       zname = ' '
       read(buffer,*,iostat=ierr)zindex(i), zname
       if(ierr.ne.0)then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'Illegal input in field ','ZMATFILE'
          write(stderr,*)'Offending input is: ',trim(buffer)
          write(stderr,*)'Program exiting.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if
       zmat(zindex(i)) = load_zmatrix(zname)
       if (.not. quiet)write(stdout,'(1x,"Reading z-matrix ",i3," of ",i3," from : ",a34)')i,total_num_zm,trim(zname)
       if (.not. quiet)write(stdout,*)'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
       if (.not. quiet)call print(stdout,zmat(zindex(i)))
       if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
       num_at(zindex(i)) = num(zmat(zindex(i)))
    end do
    do i = 1,total_num_zm
       if(.not.(any(zindex.eq.i))) then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'z-matrix indexing illegal.  Check input file.'
          write(stderr,*)'Program exiting.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if
    end do

    deallocate(zindex)

    return

  end subroutine get_zmats

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_ninternal(nintern,total_num_zm,key,quiet)

    character(len=100)     ::  keyword
    character(len=500)     ::  buffer
    character(len=4)       ::  zmat, text
    integer                :: i,ierr, n
    integer, intent(out)   :: nintern(:)
    integer, intent(in)    :: total_num_zm
    integer, allocatable   :: zmat_num(:)
    type (keyword_object)              :: key
    logical, intent(in)                :: quiet 

    allocate(zmat_num(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('zmat_num(total_num_zm) ')

    ! intern() and nintern() are already sized.
    ! First, look for any instances of NUMINTERNAL

    zmat = 'ZMAT'
    keyword = 'NUMINTERNAL'
    nintern = -9999
    if (exists(key, keyword)) then
       n = num_value(key, keyword)
       do i = 1,n
          buffer=' '
          buffer = get_value(key,keyword,i)
          read(buffer,*,iostat=ierr)text,zmat_num(i),nintern(zmat_num(i))
          if(ierr.ne.0)then
             call fatal_keyword(keyword,buffer)
          else
             text = .ucase.text
             !              make sure they are acceptable numbers.
             if (zmat.ne.text) call fatal_keyword(keyword,buffer)
             if((zmat_num(i).lt.1).or.(nintern(zmat_num(i)).lt.0))call fatal_keyword(keyword,buffer)
             if (.not. quiet)write(stdout,'(" Number of internal d.f. for z-matrix ",i5)')zmat_num(i)
             if (.not. quiet)write(stdout,'(" specified explicitly is: ",i5)')nintern(zmat_num(i))
             if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
          end if
       end do ! i = 1,n

    end if

    if (minval(nintern).eq.-9999) then
       if (.not. quiet)write(stdout,*)'Value for ',trim(keyword),' not specified explicitly'
       if (.not. quiet)write(stdout,*)'for one or more z-matrices.'
       if (.not. quiet)write(stdout,*)'Will try to work it out from other inputs.'
       if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if
    deallocate(zmat_num)

    return

  end subroutine get_ninternal

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_internals(nintern,total_num_zm,intern,key,quiet, num_at)

    character(len=100)     ::  keyword
    character(len=500)     ::  buffer
    character(len=8)       :: int_type , zmat, inte, text1, text2
    integer                :: i,ierr, n, zmat_num,internal_num,atom, int_type_num
    integer, intent(inout) :: nintern(:)
    integer, intent(inout) :: intern(:,:,:)
    integer, allocatable   :: counter(:)
    integer, intent(in)    :: total_num_zm
    type (keyword_object)              :: key
    logical, intent(in)                :: quiet 
    integer, intent(in)                :: num_at(:)

    ! First, work out how many internals per z-matrix, assuming it has
    ! not yet been specified.

    allocate(counter(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('counter(total_num_zm) ')
    counter = 0
    zmat = ''
    inte = ''
    zmat = 'ZMAT'
    inte = 'INT'
    keyword = 'INTERNAL'
    n = num_value(key, keyword)
    do i = 1,n

       ! read over occurances of INTERNAL and (1) see if any refer to  z-matrix that
       ! does not exist, and (2) see that the correct number of internals have been
       ! specified if nintern() is not filled.  If nintern() is not filled, fill 
       ! it and make sure it is filled.

       buffer=' '
       buffer = get_value(key,keyword,i)
       read(buffer,*,iostat=ierr)text1,zmat_num,text2,internal_num,int_type,atom
       if(ierr.ne.0) then
          call fatal_keyword(keyword,buffer)
          text2 = .ucase. text2
          text1 = .ucase. text1
          if((trim(text2).ne.trim(inte)).or.(trim(text1).ne.trim(zmat)))  then
             write(stderr,*)'-----------------------------------------------'
             write(stderr,*)'Entry in ',trim(keyword),' field in input file is invalid.'
             write(stderr,*)'Problem is in ',trim(text1),' or ',trim(text2)
             write(stderr,*)'z-matrix number given is invalid.'
             write(stderr,*)'-----------------------------------------------'
             stop
          end if
       else
          if((zmat_num.lt.1).or.(zmat_num.gt.total_num_zm))call fatal_keyword(keyword,buffer)
          counter(zmat_num) = counter(zmat_num) + 1
       end if ! if ierr
    end do

    do zmat_num = 1,total_num_zm
       if (nintern(zmat_num).eq.-9999) then
          nintern(zmat_num) = counter(zmat_num)
       else
          if (nintern(zmat_num).ne.counter(zmat_num))then
             write(stderr,*)'-----------------------------------------------'
             write(stderr,*)'Number of internal degrees of freedom defined '
             write(stderr,*)'for z-matrix ',zmat_num,' disagrees with the   '
             write(stderr,*)'number expected based on value of NUMINTERNAL.'
             write(stderr,*)'-----------------------------------------------'
             stop
          end if
       end if
    end do

    ! Now actually fill intern()

    do i = 1,n
       buffer=' '
       buffer = get_value(key,keyword,i)
       read(buffer,*,iostat=ierr)text1,zmat_num,text2,internal_num,int_type,atom
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       if((zmat_num.lt.1).or.(zmat_num.gt.total_num_zm))call fatal_keyword(keyword,buffer)
       if((internal_num.lt.1).or.(internal_num.gt.nintern(zmat_num)))call fatal_keyword(keyword,buffer)
       if((atom.lt.1).or.(atom.gt.num_at(zmat_num)))call fatal_keyword(keyword,buffer)

       ! Process the text.
       ! Make sure the field IS text.

       int_type_num = 0
       int_type = .ucase.int_type
       if(.isletter.trim(int_type)) then
          if (trim(int_type).eq.'LENGTH') int_type_num = 1
          if (trim(int_type).eq.'ANGLE') int_type_num = 2
          if (trim(int_type).eq.'DIHEDRAL') int_type_num = 3
          if (int_type_num.eq.0) call fatal_keyword(keyword,buffer)
          intern(internal_num,2,zmat_num) = int_type_num
          intern(internal_num,1,zmat_num) = atom
       else
          call fatal_keyword(keyword,buffer)
       end if

       if(.not.quiet) then
          write(stdout,'(" For atom ",i5," in zmat ",i5," the ",a10," can vary.")')atom,zmat_num,trim(.lcase.int_type)
          if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'

       end if

    end do

    deallocate(counter)

    return

  end subroutine  get_internals

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  get_num_loc_occ(num_loc_given, num_loc_qxyz, num_loc_occ, total_num_zm, mol_max, key,quiet, simsize)

    character(len=100)     ::  occname
    character(len=500)     ::  buffer
    integer                :: i,ierr, ia,ib,ic,il,iz,im, ddunit, loc_max
    integer, intent(in)    :: total_num_zm,num_loc_given, num_loc_qxyz
    integer, intent(out)   :: num_loc_occ, mol_max
    type (keyword_object)  :: key
    logical, intent(in)    :: quiet 
    logical                :: finished
    integer,dimension(3), intent(in)   :: simsize

    buffer = get_value(key,'OCCFILE')
    read(buffer,*,iostat=ierr)occname
    if(ierr.ne.0)then
       write(stderr,*)'-----------------------------------------------'
       write(stderr,*)'Illegal input in field ','OCCFILE'
       write(stderr,*)'Offending input is: ',trim(buffer)
       write(stderr,*)'Have you specified a valid occupancy file?'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'-----------------------------------------------'
       stop
    end if

    ddunit = open(occname, status='old') 

    if(ddunit.eq.-1) call file_does_not_exist(trim(occname))

    ! Read until end of file, looking at the location and molecule fields and
    ! noting down the maxvals in each.

    loc_max = -9999
    mol_max = -9999
    i = 0
    do 
       call read_buffer(ddunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
       read(buffer,*,end=13,iostat=ierr)ia,ib,ic,il,iz,im
       i = i + 1
       if(ierr.ne.0) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Error in reading occupancy structure file, line ',i
          write(stderr,*)'of file ',trim(occname)
          write(stderr,*)'Offending line is: ',trim(buffer)
          write(stderr,*)'Program exiting.'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if (il.gt.loc_max)loc_max = il
       if(im.gt.mol_max) mol_max = im

       ! zocc contains which zmatrix is on which location.
       ! mocc contains the number of the occurance of that zmatrix ON THAT LOCATION,
       ! which is to say it gets reset to zero when you move to a new location

       if((ia.gt.simsize(1)).or.(ia.lt.1)) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'a dimension in occupancy file is invalid.  Value: ',ia
          write(stderr,*)'Program exiting. Code AOC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if((ib.gt.simsize(2)).or.(ib.lt.1)) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'b dimension in occupancy file is invalid.    Value: ',ib
          write(stderr,*)'Program exiting. Code BOC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if((ic.gt.simsize(3)).or.(ic.lt.1)) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'c dimension in occupancy file is invalid.    Value: ',ic
          write(stderr,*)'Program exiting. Code COC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if((iz.gt.total_num_zm).or.(iz.lt.1)) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'z-matrix number in occupancy file is invalid.  Value: ',iz
          write(stderr,*)'Program exiting. Code ZOC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if(im.lt.1) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Molecule number in occupancy file is invalid.  Value: ',im
          write(stderr,*)'Program exiting. Code MOC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if(il.lt.1) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Location number in occupancy file is invalid.  Value: ',il
          write(stderr,*)'Program exiting. Code LNO1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if

    end do

13  close(ddunit)

    num_loc_occ = loc_max
    if (num_loc_given.ne.-9999) then
       ! It is meaningful to do a test.   
       if(num_loc_occ.ne.num_loc_given) then       
          write(stderr,*)'----------------------------------------------------------------------------'      
          write(stderr,*)'Too many or not enough or wrong locations specified in occupancy file.'
          write(stderr,*)'Or NUMLOCS is wrong. num_loc_given, loc_max',num_loc_given,loc_max       
          write(stderr,*)'Program exiting.'      
          write(stderr,*)'----------------------------------------------------------------------------'
          stop   
       end if
    end if

    ! Same test for num_loc_qxyzz

    if (num_loc_qxyz.ne.-9999) then
       ! It is meaningful to do a test.   
       if(num_loc_occ.ne.num_loc_qxyz) then  
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Too many or not enough or wrong locations specified in occupancy file.' 
          write(stderr,*)'Or NUMLOCS is wrong. num_loc_occ, num_loc_qxyz',num_loc_occ, num_loc_qxyz       
          write(stderr,*)'Program exiting.'             
          write(stderr,*)'----------------------------------------------------------------------------' 
          stop               
       end if
    end if

    if (num_loc_occ.ne.-9999) then
       if(.not. quiet)      &
            write(stdout,'(" Interrogation of occupancy file suggests there are",i5," locations in unit cell.")')num_loc_occ
       if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if

    return

  end subroutine  get_num_loc_occ

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_num_loc_qxyz(num_loc_given,num_loc_qxyz,total_num_zm,key,quiet)

    character(len=100)     ::  qxyzname
    character(len=500)     ::  buffer
    integer                :: i,ierr, total_num_qxyz
    integer, intent(in)    :: total_num_zm,num_loc_given
    integer, intent(out)   :: num_loc_qxyz
    integer, allocatable   :: zindex(:)
    integer                :: l,j,dummy1, dummy2,dunit
    type (keyword_object)  :: key
    logical, intent(in)    :: quiet 
    logical                :: finished

    num_loc_qxyz = -9999
    allocate(zindex(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('zindex(total_num_zm) ')

    total_num_qxyz = 0
    total_num_qxyz = num_value(key,'QXYZFILE')
    if (total_num_qxyz.eq.0) then
       if(.not.quiet)write(stdout,*)'No QXYZFILE given, this may be a problem!'
    end if
    if (total_num_qxyz.ne.total_num_zm) then
       if(.not.quiet)write(stdout,*)'Number of QXYZFILEs not same as number of z-matrices.'
       if(.not.quiet)write(stdout,*)'This may be a problem!'
    end if
    do i = 1,total_num_qxyz
       buffer = ' '
       buffer = get_value(key,'QXYZFILE',i)
       read(buffer,*,iostat=ierr)zindex(i), qxyzname
       if(ierr.ne.0)then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'Illegal input in field ','QXYZFILE'
          write(stderr,*)'Offending input is: ',trim(buffer)
          write(stderr,*)'Program exiting.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if

       ! OK, we have the name, let's see if we can open the file...

       dunit = open(qxyzname, status='old')
       if(dunit.eq.-1) call file_does_not_exist(trim(qxyzname))

       ! OK, the file is open.  Now, we have to read until we get an empty line.

       do
          call read_buffer(dunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
          read(buffer,*,end=12,iostat=ierr)dummy1,l,dummy2,j
          if(ierr.ne.0) then
             write(stderr,*)'----------------------------------------------------------------------------'
             write(stderr,*)'Problem reading lines in .qxyz file: ',trim(qxyzname)
             write(stderr,*)'Offending line is: ',trim(buffer)
             write(stderr,*)'Program exiting.'
             write(stderr,*)'----------------------------------------------------------------------------'
             stop
          end if

          if (l.gt.num_loc_qxyz)num_loc_qxyz = l

       end do ! read through file

12     close(dunit)

    end do ! i=1,zm

    if (num_loc_given.ne.-9999) then
       ! It is meaningful to do a test.
       if(num_loc_qxyz.eq.num_loc_given) then
       else
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Too many or not enough or wrong locations specified in qxyz files.'
          write(stderr,*)'Or NUMLOCS is wrong. num_loc_given, num_loc_qxyz',num_loc_given,num_loc_qxyz
          write(stderr,*)'Program exiting.'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       ! num_loc_given was not given, so num_loc_qxyz is now our best guess
       ! anyway.       
    end if

    do i = 1,total_num_zm
       if(.not.(any(zindex.eq.i))) then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'QXYZ file zmatrix numbers illegal.  Check input file.'
          write(stderr,*)'Program exiting.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if

    end do
    if (num_loc_qxyz.ne.-9999) then
       if(.not. quiet)write(stdout,'(" Interrogation of qxyz files suggests there are",i5," locations in unit cell.")')num_loc_qxyz
    end if

    deallocate(zindex)

    return


  end subroutine get_num_loc_qxyz

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  get_zocc_mocc(zocc,mocc,simsize,num_loc,total_num_zm,key,quiet)

    integer,intent(inout)                     :: zocc(:,:,:,:),mocc(:,:,:,:)
    type (keyword_object)  :: key
    logical, intent(in)    :: quiet 
    integer,dimension(3), intent(in)   :: simsize
    integer, intent(in) :: num_loc,total_num_zm

    if (exists(key, 'OCCFILE')) then
       call  read_zocc_mocc(zocc,mocc,simsize,num_loc,total_num_zm,key,quiet)
    else
       call work_out_zocc_mocc(zocc,mocc,simsize,total_num_zm,num_loc,key,quiet)
    end if

    return

  end subroutine  get_zocc_mocc

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  read_zocc_mocc(zocc,mocc,simsize,num_loc,total_num_zm,key,quiet)

    character(len=100)     ::  occname
    character(len=500)     ::  buffer
    integer                :: i,ierr, ia,ib,ic,il, ddunit
    integer,intent(out)    :: zocc(:,:,:,:),mocc(:,:,:,:)
    logical                :: finished
    type (keyword_object)  :: key
    logical, intent(in)    :: quiet 
    integer , intent(in)   :: num_loc,total_num_zm
    integer,dimension(3), intent(in)   :: simsize

    buffer = get_value(key,'OCCFILE')
    occname = ''
    read(buffer,*,iostat=ierr)occname
    if(ierr.ne.0)then
       write(stderr,*)'-----------------------------------------------'
       write(stderr,*)'Illegal input in field ','OCCFILE'
       write(stderr,*)'Offending input is: ',trim(buffer)
       write(stderr,*)'Have you specified a valid occupancy file?'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'-----------------------------------------------'
       stop
    end if

    ddunit = open(occname, status='old')

    if(ddunit.eq.-1) call file_does_not_exist(trim(occname))

    do i = 1, simsize(1)*simsize(2)*simsize(3)*num_loc
       call read_buffer(ddunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
       read(buffer,*,iostat=ierr)ia,ib,ic,il,zocc(ia,ib,ic,il),mocc(ia,ib,ic,il)
       if(ierr.ne.0) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Error in reading occupancy structure file, line ',i
          write(stderr,*)'of file ',trim(occname)
          write(stderr,*)'Offending line is: ',trim(buffer)
          write(stderr,*)'Program exiting.'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if

       ! zocc contains which zmatrix is on which location.
       ! mocc contains the number of the occurance of that zmatrix ON THAT
       ! LOCATION,
       ! which is to say it gets reset to zero when you move to a new location

       if((ia.gt.simsize(1)).or.(ia.lt.1)) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'a dimension in occupancy file is invalid.  Value: ',ia
          write(stderr,*)'Program exiting. Code AOC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if((ib.gt.simsize(2)).or.(ib.lt.1)) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'b dimension in occupancy file is invalid.    Value: ',ib
          write(stderr,*)'Program exiting. Code BOC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if((ic.gt.simsize(3)).or.(ic.lt.1)) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'c dimension in occupancy file is invalid.    Value: ',ic
          write(stderr,*)'Program exiting. Code COC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if((zocc(ia,ib,ic,il).gt.total_num_zm).or.(zocc(ia,ib,ic,il).lt.1)) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'z-matrix number in occupancy file is invalid.  Value: ',zocc(ia,ib,ic,il)
          write(stderr,*)'Program exiting. Code ZOC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if(mocc(ia,ib,ic,il).lt.1) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Molecule number in occupancy file is invalid.  Value: ',mocc(ia,ib,ic,il)
          write(stderr,*)'Program exiting. Code MOC1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
       if(il.lt.1) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Location number in occupancy file is invalid.  Value: ',il
          write(stderr,*)'Program exiting. Code LNO1'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if

    end do
    close(ddunit)

    if (.not. quiet)write(stdout,'(1x,"Occupancies read in from         : ",a30)')trim(occname)
    if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'

    if(total_num_zm.ne.maxval(zocc)) then
       write(stderr,*)'----------------------------------------------------------------------------'
       write(stderr,*)'Mismatch between z-matrix numbers in occupancy file and in input files'
       write(stderr,*)'Program exiting. Code TMZ3'
       write(stderr,*)'----------------------------------------------------------------------------'
       stop
    end if

    return

  end  subroutine read_zocc_mocc

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine work_out_zocc_mocc(zocc,mocc,simsize,total_num_zm,num_loc,key,quiet)

    character(len=500)     ::  buffer
    character(len=100)     ::  qxyzname
    integer                :: i,ierr,  total_num_qxyz, j,k
    integer, allocatable   :: zindex(:),zm_on_lock(:),mol_on_lock(:),loc_counter(:)
    integer                :: l,dummy1, dummy2, dunit
    integer,intent(out)    :: zocc(:,:,:,:),mocc(:,:,:,:)
    integer, intent(in)    :: total_num_zm, num_loc
    logical                :: finished
    type (keyword_object)  :: key
    logical, intent(in)    :: quiet 
    integer,dimension(3), intent(in)   :: simsize

    allocate(zindex(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('zindex(total_num_zm) ')
    allocate(mol_on_lock(num_loc),stat=ierr)
    if (ierr.ne.0) call allerr('mol_on_lock(num_loc) ')
    allocate(zm_on_lock(num_loc),stat=ierr)
    if (ierr.ne.0) call allerr('zm_on_lock(num_loc) ')
    allocate(loc_counter(num_loc),stat=ierr)
    if (ierr.ne.0) call allerr('loc_counter(num_loc) ')

    total_num_qxyz = 0
    total_num_qxyz = num_value(key,'QXYZFILE')
    if (total_num_qxyz.eq.0) then
       if(.not.quiet)write(stdout,*)'No QXYZFILE given, this may be a problem!'
       if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if
    if (total_num_qxyz.ne.total_num_zm) then
       if(.not.quiet)write(stdout,*)'Number of QXYZFILEs not same as number of z-matrices.'
       if(.not.quiet)write(stdout,*)'This may be a problem!'
       if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if

    loc_counter = 0

    do i = 1,total_num_qxyz
       buffer = get_value(key,'QXYZFILE',i)
       read(buffer,*,iostat=ierr)zindex(i), qxyzname
       if(ierr.ne.0)then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'Illegal input in field ','QXYZFILE'
          write(stderr,*)'Offending input is: ',trim(buffer)
          write(stderr,*)'Program exiting.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if

       ! OK, we have the name, let's see if we can open the file...

       dunit = open(qxyzname, status='old')

       if(dunit.eq.-1) call file_does_not_exist(trim(qxyzname))

       ! OK, the file is open.  Now, we have to read until we get an empty line.

       do
          call read_buffer(dunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
          read(buffer,*,end=16,iostat=ierr)dummy1,l,dummy2,j
          if(ierr.ne.0) then
             write(stderr,*)'----------------------------------------------------------------------------'
             write(stderr,*)'Problem reading lines in .qxyz file: ',trim(qxyzname)
             write(stderr,*)'Offending line is: ',trim(buffer)
             write(stderr,*)'Program exiting.'
             write(stderr,*)'----------------------------------------------------------------------------'
             stop
          end if
          loc_counter(l) = loc_counter(l) + 1
          mol_on_lock(l) = j
          zm_on_lock(l) = zindex(i)
       end do ! read through file

16     close(dunit)

    end do ! i=1,zm

    if((maxval(loc_counter).gt.1).or.(maxval(mol_on_lock).gt.1))then
       ! One loc can have more than one thing on it and we need an occfile.
       write(stderr,*)'----------------------------------------------------------------------------'
       write(stderr,*)'Must supply an occupancy structure; at least'
       write(stderr,*)'one location can have more than one "thing" on it.'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'----------------------------------------------------------------------------'
       stop
    else
       ! Replicate mol_on_lock() and zm_on_lock() across crystal
       do i = 1,simsize(1)
          do j = 1,simsize(2)
             do k = 1,simsize(3)
                do l = 1,num_loc
                   mocc(i,j,k,l) = mol_on_lock(l)
                   zocc(i,j,k,l) = zm_on_lock(l)
                end do
             end do
          end do
       end do
    end if
    if(.not. quiet)write(stdout,*)'Occupancy structure worked out from QXYZ files.'
    if(.not. quiet)write(stdout,*)'Unless there is an error this means there is no occupancy disorder.'
    if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'

    deallocate(zindex)
    deallocate(mol_on_lock)
    deallocate(zm_on_lock)
    deallocate(loc_counter)

    return

  end subroutine work_out_zocc_mocc

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_what_goes_where(loc_counter,zm_on_loc,mol_on_loc,which_mol,total_num_zm,num_loc,quiet,key)

    character(len=500)     ::  buffer
    character(len=100)     ::  qxyzname
    integer                :: i,ierr,  total_num_qxyz, j
    integer, allocatable   :: zindex(:),old_val_a(:)
    integer,intent(out)    :: zm_on_loc(:,:),mol_on_loc(:,:),loc_counter(:),which_mol(:,:,:)
    integer                :: l,dummy1, dummy2,dunit
    logical                :: finished
    integer, intent(in)    :: total_num_zm, num_loc
    type (keyword_object)  :: key
    logical, intent(in)    :: quiet 

    allocate(zindex(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('zindex(total_num_zm) ')
    allocate(old_val_a(num_loc),stat=ierr)
    if (ierr.ne.0) call allerr('old_val_a(num_loc) ')

    total_num_qxyz = 0
    total_num_qxyz = num_value(key,'QXYZFILE')
    if (total_num_qxyz.eq.0) then
       if(.not.quiet)write(stdout,*)'No QXYZFILE given, this may be a problem!'
       if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if
    if (total_num_qxyz.ne.total_num_zm) then
       if(.not.quiet)write(stdout,*)'Number of QXYZFILEs not same as number of z-matrices.'
       if(.not.quiet)write(stdout,*)'This may be a problem!'
       if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if

    loc_counter = 0
    zm_on_loc = 0
    mol_on_loc = 0

    do i = 1,total_num_qxyz
       buffer = get_value(key,'QXYZFILE',i)
       read(buffer,*,iostat=ierr)zindex(i), qxyzname
       if(ierr.ne.0) then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'Illegal input in field ','QXYZFILE'
          write(stderr,*)'Offending input is: ',trim(buffer)
          write(stderr,*)'Program exiting.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if

       ! OK, we have the name, let's see if we can open the file...

       dunit = open(qxyzname, status='old')

       if(dunit.eq.-1) call file_does_not_exist(trim(qxyzname))

       ! OK, the file is open.  Now, we have to read until we get an empty line.

       old_val_a = loc_counter
       do
          call read_buffer(dunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
          read(buffer,*,end=26,iostat=ierr)dummy1,l,dummy2,j
          if(ierr.ne.0) then
             write(stderr,*)'----------------------------------------------------------------------------'
             write(stderr,*)'Problem reading lines in .qxyz file: ',trim(qxyzname)
             write(stderr,*)'Offending line is: ',trim(buffer)
             write(stderr,*)'Program exiting.'
             write(stderr,*)'----------------------------------------------------------------------------'
             stop
          end if

          if (old_val_a(l).eq.loc_counter(l)) loc_counter(l) = loc_counter(l) + 1

          zm_on_loc(l,loc_counter(l)) = zindex(i)

          mol_on_loc(l,zindex(i)) =  mol_on_loc(l,zindex(i)) + 1

          which_mol(l,zindex(i),mol_on_loc(l,zindex(i))) = j

       end do ! read through file

26     close(dunit)

    end do ! i=1,zm

    if(.not.quiet) then
       do l = 1,num_loc
          write(stdout,'(" Location",i5," can have",i5," type(s) of z-matrix(ices).")')l,loc_counter(l)
          write(stdout,*)'These are types',(zm_on_loc(l,j),j=1,loc_counter(l))
          do i = 1,loc_counter(l)
             write(stdout,'(" Type",i5," is in",i5," orientation(s) on that location.")')zm_on_loc(l,i),mol_on_loc(l,zm_on_loc(l,i))
             write(stdout,'(" These are orientation(s): ")',advance='no')
             do j = 1, mol_on_loc(l,zm_on_loc(l,i))
                write(stdout,'(i4,1x)',advance='no')which_mol(l,zm_on_loc(l,i),j) 
             end do
             write(stdout,*)
          end do
       end do

       write(stdout,*)'----------------------------------------------------------------------------'
    end if

    deallocate(zindex)
    deallocate(old_val_a)

    return

  end subroutine get_what_goes_where

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

end module zmc_subroutines

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

module varmod

  use crystallography_class
  use rotmatrix_class
  use quaternion_class
  use iso_varying_string
  use keyword_class, only: exists, get_value, has_value, initialise, keyword_object, assignment(=),read, num_value
  use zmatrix_class, only: assignment(=),zmatrix_object,load_zmatrix,copy,as_xyz,as_cartesian,print,parameter,labels,num

  implicit none

  ! Characters

  character(len=100)       :: outname, oname, header, inname, keyfile
  character(len=107)       :: fname2,newname
  character(len=6)         ::inline

  ! Integers

  ! variable names prefixed by o usually mean origin, as in origin of a contact
  ! vector -- when we are looking at the energy of a single molecule, we
  ! want to sum over all the contact vectors attached to that molecule,
  ! so that molecule is considered as being at the centre of a sort of network
  ! of vectors all of which originate on an atom of that molecule.
  ! Hence it is the 'origin' molecule..
  !
  ! Names prefixed by d mean destination, so they describe where that contact 
  ! vector connects to.  For example, a line in the contact vector file would
  ! look like
  !
  !  ol  oz  om  oat  da  db  dc  dl  dz  dm  dat  length  type
  !
  ! Where the contact vector's origin in on atom oat of the om'th occurance
  ! of z-matrix number oz on location (site if you like) ol in the unit cell.
  !
  ! The destination of the contact vector is atom dat on the dm'th instance of 
  ! z-matrix number dz on location dl in the unit cell da along a, db along b and dc 
  ! along c (where da = 1 means move one a lattice vector further along) relative to
  ! the unit cell containing the origin molecule.

  integer                                  :: dunit, mcell, i1,i2,ierr
  !                                        dunit = unit no. for reads & writes
  !                                        mcell = total num. of mol. in the crystal
  !                                        i1,i2 = random number seeds
  !                                        ierr is for error trapping

  integer                                  :: typ, ntyp,ncy,mtyp, icy, ic2, ntyp_given
  !                                        typ,ntyp,mtyp relate to num. of vec types, 
  !                                        ncy, ic2 and icy -> number of MC cycles

  integer                                  :: total_num_zm ,k,m,iz,il,ol
  integer                                  :: given_total_num_zm 
  !                                        num_zm = number of z-matrices
  !                                        k,m,iz = indexing variables
  !                                        iz used for zmatrix
  !                                        given_ is value specified
  !                                        optionally explicitly.

  integer, allocatable                     :: num_at(:),nintern(:)
  !                                        num_at() = number of atoms in a given zmat
  !                                        nintern() = number of internal degrees of freedom

  real, allocatable                        :: invals_ge(:)
  !                                        working array used in get_energy()

  integer, dimension(3)                    :: simsize
  !                                        dimensions of simulated crystal in unit cells

  integer                                  :: iat, ninspr, maxval_mocc, ninspr_given
  !                                        iat = counter for atom number
  !                                        ninspr = total number of internal spring constants

  integer,allocatable                      :: conmol(:,:,:)
  !                                        number of contact vectors for given molecule
  !                                        of a given type on a given location.

  integer,allocatable                      :: da(:,:,:,:),db(:,:,:,:),dc(:,:,:,:)
  integer,allocatable                      :: dz(:,:,:,:),dm(:,:,:,:),dl(:,:,:,:)
  integer,allocatable                      :: ty(:,:,:,:),oat(:,:,:,:),dat(:,:,:,:)
  !                                        contact vector descriptor arrays 
  !                                        da, db, dc = unit cell translations
  !                                        dz = destination z-matrix
  !                                        dm = destination molecule
  !                                        ty = vector type
  !                                        oat = origin atom for contact vector
  !                                        dat = destination atom.

  integer                                  :: oa, ob, oc, oz, om, nin
  integer,allocatable                      :: irej(:), iacc(:), iup(:)
  !                                        o? identify origin molecule for MC
  !                                        irej etc are diagnostic.  

  integer,allocatable                      :: intern(:,:,:)
  !                                        describes the internal variables --
  !                                        which atom and what variable (bond length,
  !                                        bond angle or dihedral angle) given the
  !                                        ordinal number of the internal variable and
  !                                        which type of molecule (which z-matrix)

  integer,allocatable                      :: wrap(:,:)
  !                                        boundary conditions for MC calculations

  integer, allocatable                     :: insprtype(:,:)
  !                                        type for each internal spring.  
  !                                        i.e., which spring constant to apply (from insprcon).

  integer, allocatable                     :: reject(:,:,:)
  !                                        acceptance/rejection statistics for each variable 
  !                                        as a function of increment size.

  integer                                  :: ibin,iivar
  !                                        variables to use in the examination of 
  !                                        acceptance/rejection ratios

  integer                                  :: adjcy,badjcy
  !                                        variable increments will be adjusted every adjcy 
  !                                        cycles to keep about 50% rejections
  !                                        spring constants will be globally scaled every badjcy
  !                                        cycles to get the desired Biso (see B_adjust)

  integer                                  :: num_loc, num_loc_given, num_loc_qxyz, num_loc_occ, mol_max, con_test
  !                                        number of locations in the unit cell
  !                                        there are a number of ways of thinking about a
  !                                        location, which was a name chosen (rather than site
  !                                        or position) exactly because it is free from
  !                                        crystallographic connotations.  A location is a region
  !                                        of the unit cell which can be occupied by 1 or more
  !                                        MUTUALLY EXCLUSIVE objects.  If there is
  !                                        no occupancy aspect to a problem, then it is pretty
  !                                        much the same as molecular site.
  !                                        _given is the number specified by keyword in the
  !                                        input file; if omitted, it would be
  !                                        worked out from other input files
  !                                        if given, it is worked out and compared to given value for
  !                                        error checking.

  integer,allocatable                      :: zocc(:,:,:,:),mocc(:,:,:,:)
  !                                        arrays of occupancies, z-matrix and molecular 
  !                                        number on each location. 

  integer,allocatable                      :: crosstype(:,:)
  !                                        which insprtype to use for this cross-term in internal d.f.

  integer,allocatable                      :: cross_terms(:,:,:)
  !                                        which internals to use in a cross term calculation

  integer, allocatable                     :: ncross(:), ncross_given(:)
  integer, allocatable                     :: how_many_types(:),which_types(:,:), instances_types(:,:)
  integer, allocatable                     :: which_mol(:,:,:)
  !                                        book-keeping arrays saying what goes where.

  ! Reals

  real, dimension(6)                       :: cellpar
  !                                        cell parameters read into this
  !                                        1:3 are a, b, c in Angstrom
  !                                        4:6 are alpha, beta, gamma in degrees. However
  !                                        all internal manipulation of angles use radians.

  real                                     ::temp, rdiv_max, one_on_temp
  !                                        temp is MC temperature.
  !                                        rdiv_max is real(div_max)

  real*8                                   :: edif,new_e,old_e
  !                                        old_e is energy before modifying configuration
  !                                        new_e is after and edif is difference

  real,allocatable                         :: translate(:,:,:,:)
  !                                        coordinates of origin atom of z-matrices 
  !                                        in average structure

  real,allocatable                         :: xyzcrystal(:,:,:,:,:)
  !                                        stored origin coordinates for all molecules in model crystal

  real, allocatable                        :: tempin(:)
  !                                       temporary storage of internal variables of molecule being MC'd

  real, dimension(3)                       :: rsimsize, tempxyz
  !                                        simsize as reals, temporary storage of molecule origin 
  !                                        coordinates for mol. being MC'd

  real, dimension(3,3)                     :: celltrans
  !                                        unit cell translation vectors for 100, 010, 001

  real,allocatable                         :: carts(:,:,:,:,:,:)
  !                                        Cartesian positions of all atoms in simulation

  real,allocatable                         :: tempcarts(:,:)
  !                                        temporary storage of Cartesians for the selected molecule

  real,allocatable                         ::length(:,:,:,:)
  !                                        average lengths of contact vectors, 
  !                                        the d0 in E = Sum_i{Ei(di-d0i)^2}

  real                                     :: rvar(1), B_aim
  !                                        random variable; Bfactor we are aiming for, if we are aiming!

  real, allocatable                        :: sprcon(:),size_effect(:)
  !                                        sprcon = 'spring constants' for contact vectors.  
  !                                        the Ei in E = Sum_i{Ei(di-d0i)^2}
  !                                        size-effect gives the size-effect parameters for 
  !                                        each sprcon, no size-effect on  internals
  !                                        size-effect works by replacing d0i with
  !                                        (d0i+size_effect(i))

  real,allocatable                         :: widths(:,:)
  !                                        controls sizes of random increments used on
  !                                        xyz and q in MC
  !                                        In other words, when we modify a molecule's
  !                                        coordinates, we generate shifts in the origin
  !                                        position (xyz) and in the 4 components of the
  !                                        quaternion.
  !                                        We generate a random number rvar and then get
  !                                        the shift from (rvar-0.5)*width(i), so that
  !                                        for a given variable, i, the shift is in the range
  !                                         -(width(i)/2) to +(width(i)/2).  These shifts, along
  !                                        with the shifts on the internal variables, if any,
  !                                        are assembled into the vector tweak() (see below).

  real,allocatable                         :: tweak(:)
  !                                        tweak is the vector of increments added to molecular 
  !                                        variables before calculating new_e

  real,allocatable                         :: initw(:,:), mcw(:,:)
  !                                        initw controls the initial widths of variable 
  !                                        distributions for all variables, internals AND
  !                                        externals (external = xyz and q)

  real,allocatable                         :: inwidths(:,:)
  !                                        controls sizes of increments used on 
  !                                        internal d.f. in MC. see widths().

  real,allocatable                         :: incrystal(:,:,:,:,:)
  !                                        values of the variables for the internal 
  !                                        degrees of freedom for every molecule

  real, allocatable                        :: insprcon(:), cross_spr(:)
  !                                        spring constants for internal springs.

  real                                     :: radjcy,bradjcy
  !                                        real(adjcy) and real(badjcy) 

  ! Custom Types

  type (quaternion), allocatable           :: quatern(:,:,:)
  !                                        quaternion for each molecule in the 
  !                                        average unit cell

  type (quaternion), allocatable           :: qcrystal(:,:,:,:)
  !                                        quaternion for each molecule in the 
  !                                        model crystal

  type (crystallography_object)            :: xtal
  !                                        contains cell parameters (see cellpar)

  type (quaternion)                        :: tempq
  !                                        temporary storage for quaternion of 
  !                                        molecule being tweaked in MC

  type (keyword_object),save               :: key2
  type (varying_string),save               :: mykeywords(41)
  type (varying_string),save               :: keyword
  !                                        related to keywords used for input

  integer                                  :: i

  type (rotmatrix)                         :: rmat
  !                                        rotation matrix

  type (zmatrix_object),allocatable        :: zmat(:)
  !                                        z-matrix(ices) of the molecules in the simulation

  type (zmatrix_object),save               :: z1
  !                                        one z-matrix

  real, allocatable                        :: new_coords(:,:)
  !                                        new_coords is a 3xN array, effectively.
  !                                        Cartesian coordinates of atoms in a given z-matrix

  character(len=500)                       :: buffer
  character(len=50)                        :: cver
  logical                                  :: finished,quiet
  type (varying_string),save               :: myoptions(20)

  !------------------------------------------------------------------------------!
  !                      End of declarations                                     !
  !------------------------------------------------------------------------------!

end module varmod

module energymod

  use varmod

contains

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_energy(ol,om,oz,oa,ob,oc,energy)

    !----------------------------------------------------------------------------------------------
    ! Get the energy of a molecule in the crystal; sum over contact vectors, 
    ! internal degrees of freedom and cross terms between internal d.f.
    !----------------------------------------------------------------------------------------------

    integer                                 :: ia2, ib2, ic2, typ,k,nin,ddm,ddz
    integer,intent(in)                      :: oa,ob,oc,om,oz,ol

    real*8,intent(out)                      :: energy
    real                                    :: x2,y2,z2,dist,s_eff
    real,dimension(3)                       :: xyz

    type (zmatrix_object)                   :: tempzmat

    !----------------------------------------------------------------------------------------------
    ! Loop over contacts to get energy                                 
    ! Make sure that the molecule we are going to exists (occupancies)
    !----------------------------------------------------------------------------------------------

    energy = 0.0

    do k = 1,conmol(ol,om,oz)
       ia2 = wrap(oa+da(ol,oz,om,k),1)
       ib2 = wrap(ob+db(ol,oz,om,k),2)
       ic2 = wrap(oc+dc(ol,oz,om,k),3)
       ddz = zocc(ia2,ib2,ic2,dl(ol,oz,om,k))
       ddm = mocc(ia2,ib2,ic2,dl(ol,oz,om,k))
       if( (ddz.eq.dz(ol,oz,om,k)).and. (ddm.eq.dm(ol,oz,om,k)) ) then

          !------------------------------------------------------------------------------!
          ! If the molecule at the end of the contact vector is not the one that is
          ! there, then we don't add on the energy due to that vector, reasonably enough.
          ! The size effect is implemented here as an additive length. 
          !------------------------------------------------------------------------------!

          x2 = carts(1,oa,ob,oc,ol,oat(ol,oz,om,k))
          y2 = carts(2,oa,ob,oc,ol,oat(ol,oz,om,k))
          z2 = carts(3,oa,ob,oc,ol,oat(ol,oz,om,k))
          xyz(1) = carts(1,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
          xyz(2) = carts(2,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
          xyz(3) = carts(3,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
          xyz = xyz+da(ol,oz,om,k)*celltrans(1,:)+db(ol,oz,om,k)*celltrans(2,:)+  &
               dc(ol,oz,om,k)*celltrans(3,:)
          typ = ty(ol,oz,om,k)
          dist = sqrt((xyz(1)-x2)*(xyz(1)-x2)+(xyz(2)-y2)*(xyz(2)-y2)+(xyz(3)-z2)*(xyz(3)-z2))
          s_eff = size_effect(typ)
          energy = energy + (dist-(length(ol,oz,om,k)+s_eff))*        &
               (dist-(length(ol,oz,om,k)+s_eff))* sprcon(typ)
       end if
    end do

    !------------------------------------------------------------------------------!
    ! Internal energy:
    !------------------------------------------------------------------------------!

    call copy(tempzmat,zmat(oz))
    nin = nintern(oz)
    invals_ge(1:nin)= parameter(tempzmat,intern(1:nin,1,oz),intern(1:nin,2,oz))

    do k = 1, nin
       x2 =  invals_ge(k)
       ! x2 is the average value stored in the z matrix
       y2 = incrystal(oa,ob,oc,ol,k)
       ! y2 in the current value of that internal d.f.
       ia2 = insprtype(k,oz)
       ! ia2 is the type of internal spring
       z2 = insprcon(ia2)
       ! z2 is the spring constant of that internal spring
       energy = energy+(x2-y2)*(x2-y2)*z2
    end do

    !------------------------------------------------------------------------------!
    ! Here is where we get the contribution from the cross terms of the internals.
    !------------------------------------------------------------------------------!

    nin = ncross(oz)
    do k = 1, nin 
       ! Which two internals to use for cross term?
       ia2 = cross_terms(oz,k,1)
       ib2 = cross_terms(oz,k,2)
       ! Then get the deviations from ideal values...
       x2 =  invals_ge(ia2)-incrystal(oa,ob,oc,ol,ia2)
       y2 =  invals_ge(ib2)-incrystal(oa,ob,oc,ol,ib2)
       ! Then get the energy contribution
       energy = energy+x2*y2*cross_spr(crosstype(k,oz))
    end do

    return

  end subroutine get_energy

end module energymod

module mainmod

  use divs
  use zmc_subroutines
  use energymod

  implicit none

contains

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine replicator()

    integer                               :: i,j,k,l,m,n,p,z

    real                                  :: rvar(1),rvar2

    real,dimension(4)                     :: ddd

    !----------------------------------------------------------------------------------------------
    !          Replicate Molecule external coords across the model crystal
    !          putting in the randomness described by initw.
    !----------------------------------------------------------------------------------------------

    do i=1,simsize(1)
       do j=1,simsize(2)
          do k=1,simsize(3)
             do l = 1,num_loc
                z = zocc(i,j,k,l)
                m = mocc(i,j,k,l)
                do p = 1,4
                   call rannum(rvar,1)
                   ddd(p) = 2.0*(rvar(1)-.5)*initw(3+p,z)
                end do
                qcrystal(i,j,k,l)= quatern(l,z,m)+ddd(1:4)
                do n = 1,3
                   call rannum(rvar,1)
                   rvar2 = 2.0*(rvar(1)-0.5)*initw(n,z)
                   xyzcrystal(i,j,k,l,n)=translate(l,z,m,n)+rvar2
                end do
             end do
          end do
       end do
    end do

    return

  end subroutine replicator

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine getcarts()

    integer                               :: ia,ib,ic,i,j,iz,inn,il

    real,allocatable                      :: new_coords(:,:)

    type (zmatrix_object)                 :: z1
    type (rotmatrix)                      :: rmat

    !----------------------------------------------------------------------------------------------
    !          Now just fill up the carts array.  This needs updating the z-matrix
    !          for each molecule, then getting its rmat and so on                
    !----------------------------------------------------------------------------------------------

    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             do il = 1,num_loc
                iz = zocc(ia,ib,ic,il)
                allocate(new_coords(3,(num_at(iz))),stat=ierr)
                if (ierr.ne.0) call allerr('new_coords(3,new_coords(3,(num_at(iz)))    1                 ')
                call copy(z1,zmat(iz))
                inn = nintern(iz)
                incrystal(ia,ib,ic,il,1:inn)=                           &
                     parameter(z1,intern(1:inn,1,iz),intern(1:inn,2,iz), & 
                     incrystal(ia,ib,ic,il,1:inn))
                new_coords = as_xyz(z1)
                ! Make a rotation matrix out of this quaternion
                rmat = as_rotmatrix(qcrystal(ia,ib,ic,il))
                ! Pre-rotate the coordinates using this rotation matrix
                call rotate(rmat,new_coords)
                do j = 1,num_at(iz)
                   do i = 1,3
                      carts(i,ia,ib,ic,il,j)=new_coords(i,j)+          &  
                           xyzcrystal(ia,ib,ic,il,i) 
                   end do
                end do
                deallocate(new_coords)
             end do
          end do
       end do
    end do

    return

  end subroutine getcarts

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine getcon(oname)

    !----------------------------------------------------------------------------------------------
    ! Read in the contact vectors and put them into a bunch of arrays.
    !----------------------------------------------------------------------------------------------

    character(len=100),intent(in)            :: oname
    character(len=100)                       :: header2
    character(len=500) :: buffer

    integer                                  :: dunit, typ,ierr
    integer                                  :: k,m,iz, iiz
    integer                                  :: ia, ib, ic, iat, iiat, mm,il,iil
    integer,allocatable                      :: newa(:,:)
    integer                                  :: newb,icount

    real                                     ::dist

    logical                                  :: finished

    !----------------------------------------------------------------------------------------------
    ! Open the file with contact vectors in it 
    !----------------------------------------------------------------------------------------------

    allocate(newa(total_num_zm,num_loc),stat=ierr)
    if (ierr.ne.0) call allerr('newa(total_num_zm,num_loc)                          ')

    dunit = open(oname, status = 'old') 

    if(dunit.eq.-1) call file_does_not_exist(trim(oname))

    read(dunit,'(a100)',iostat=ierr)header2
    if(ierr.ne.0) then
       write(stderr,*)'----------------------------------------------------------------------------'
       write(stderr,*)'Invalid header in contact vector list.'
       write(stderr,*)'Offending line is: ',trim(buffer)
       write(stderr,*)'----------------------------------------------------------------------------'
       write(stderr,*)'Exiting.'
       stop
    end if

    if (.not. quiet)write(stdout,'(1x,"Header of contact vector file ",a23, " is... ")')trim(oname)
    if (.not. quiet)write(stdout,'(a100)')header2

    !-------------------------------------------------------------------------------------------!
    ! Columns in input file are: location number, origin z-matrix type, origin molecule number,
    ! origin atom number, relative position of destination mol. along a axis, b axis, c axis,  
    ! destination location, dest. z-matrix type, dest. molecule number, dest. atom number,    
    ! contact vector length 
    !-------------------------------------------------------------------------------------------!

    conmol=0
    newa = 0
    icount = 0

    ! Pre-read to get some numbers.

    do 
       call read_buffer(dunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
       read(buffer,*,end=1,iostat=ierr)il,iz,m,iat,ia,ib,ic,iil,iiz,mm,iiat,dist,typ
       icount = icount +1
       if(ierr.ne.0) then
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Invalid field in contact vector list.'
          write(stderr,*)'Ignoring header & comment lines, it is line ',icount
          write(stderr,*)'Line is: ',trim(buffer)
          write(stderr,*)'----------------------------------------------------------------------------'
          write(stderr,*)'Exiting.'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if

       if((il.lt.1).or.(iil.lt.1))then
          write(stderr,*)'Invalid location number specified.  Is it negative?',il,iil
          write(stderr,*)'Program exiting. Code LMAX1'
          stop
       end if

       if((iz.lt.1).or.(iiz.lt.1))then
          write(stderr,*)'Invalid z-matrix number specified.  Is it negative?',iz,iiz
          write(stderr,*)'Program exiting. Code ZMAX1'
          stop
       end if

       if((m.lt.1).or.(mm.lt.1))then
          write(stderr,*)'Invalid molecule number (instance of zmat on location) specified:',m,mm
          write(stderr,*)'  Is it negative?          Program exiting. Code MMAX1'
          stop
       end if

       if((iat.lt.1).or.(iiat.lt.1))then
          write(stderr,*)'Invalid atom number specified. Is it negative?',iat,iiat  
          write(stderr,*)'    Program exiting. Code AMAX1'
          stop
       end if

       if(typ.lt.1)then
          write(stderr,*)'Invalid contact vector type. Is it negative?',typ
          write(stderr,*)'    Program exiting. Code TMAX1'
          stop
       end if

       if(m.gt.newa(iz,il))newa(iz,il)=m
       conmol(il,m,iz) = conmol(il,m,iz)+1
    end do

1   rewind(dunit)


    allocate(da(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('da(num_loc,total_num_zm,maxval(mocc),maxval(conmol))')
    allocate(db(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('db(num_loc,total_num_zm,maxval(mocc),maxval(conmol))')
    allocate(dc(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('dc(num_loc,total_num_zm,maxval(mocc),maxval(conmol))')
    allocate(dz(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('dz(num_loc,total_num_zm,maxval(mocc),maxval(conmol))')
    allocate(dm(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('dm(num_loc,total_num_zm,maxval(mocc),maxval(conmol))')
    allocate(ty(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('ty(num_loc,total_num_zm,maxval(mocc),maxval(conmol))')
    allocate(dl(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('dl(num_loc,total_num_zm,maxval(mocc),maxval(conmol))')
    allocate(oat(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('oat()                                               ')
    allocate(dat(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('dat()                                               ')
    allocate(length(num_loc,total_num_zm,maxval(mocc),maxval(conmol)),stat=ierr)
    if (ierr.ne.0) call allerr('length()                                            ')

    read(dunit,'(a100)')header2
    newb = 0
    do il = 1,num_loc
       do iz = 1,total_num_zm
          do m = 1,newa(iz,il)
             newb=newb+conmol(il,m,iz)
          end do
       end do
    end do

    if (.not. quiet)write(stdout,'(1x,"Contact vector file has ",i6," rows.")')newb
    if(icount.ne.newb) then
       write(stderr,*)'Error in counting rows of file'
       write(stderr,*)icount,' vs ',newb,'.  Error code NEWB'
    end if

    conmol=0

    ty = 0

    do  k = 1,newb
       call read_buffer(dunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
       read(buffer,*)il,iz,m,iat,ia,ib,ic,iil,iiz,mm,iiat,dist,typ
       conmol(il,m,iz) = conmol(il,m,iz)+1
       da(il,iz,m,conmol(il,m,iz))=ia
       db(il,iz,m,conmol(il,m,iz))=ib
       dc(il,iz,m,conmol(il,m,iz))=ic
       oat(il,iz,m,conmol(il,m,iz))=iat
       dat(il,iz,m,conmol(il,m,iz))=iiat
       dz(il,iz,m,conmol(il,m,iz))=iiz
       dm(il,iz,m,conmol(il,m,iz))=mm
       ty(il,iz,m,conmol(il,m,iz))=typ
       length(il,iz,m,conmol(il,m,iz))=dist
       dl(il,iz,m,conmol(il,m,iz))=iil
    end do

    close(dunit)
    ntyp = maxval(ty)
    deallocate(newa)
    if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

    return

  end subroutine getcon

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine replinternal()

    !----------------------------------------------------------------------------------------------
    ! Create the initial state of the system regarding internal degrees of freedom, based
    ! on the information given in the MC input file.
    !----------------------------------------------------------------------------------------------

    integer                                :: iz, i , ia,ib,ic,il,inn

    real,allocatable                       :: invals(:,:)
    real                                   :: r(1)

    type (zmatrix_object)                  ::  testzmat

    !----------------------------------------------------------------------------------------------
    ! First, get the values from zmat and put into the array invals
    ! Remember that the angles are stored natively in radians.    
    ! Only convert to degrees for input/output                   
    !----------------------------------------------------------------------------------------------

    allocate(invals(maxval(num_at),total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('invals(maxval(num_at),total_num_zm)                 ')

    do iz = 1,total_num_zm
       call copy(testzmat,zmat(iz))
       i = nintern(iz)
       invals(1:i,iz)= parameter(testzmat,intern(1:i,1,iz),intern(1:i,2,iz))
    end do

    !----------------------------------------------------------------------------------------------
    ! Replicate internal d.f. across crystal with variations governed by initw
    !----------------------------------------------------------------------------------------------

    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             do il = 1,num_loc
                iz = zocc(ia,ib,ic,il)
                do inn = 1,nintern(iz)
                   call rannum(r,1)
                   r(1) = (r(1)-0.5)*initw(7+inn,iz)
                   incrystal(ia,ib,ic,il,inn)=invals(inn,iz)+r(1)
                end do
             end do
          end do
       end do
    end do
    deallocate(invals)

    return

  end subroutine replinternal

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine gettweak(oz)

    !----------------------------------------------------------------------------------------------
    ! Generate the vector of changes to apply to a molecule during an MC cycle.
    !----------------------------------------------------------------------------------------------

    integer,intent(in)                     :: oz
    integer                                :: i

    real                                   :: rvar(1)

    !----------------------------------------------------------------------------------------------
    ! Generate the vector of random increments, governed by widths and inwidths
    !----------------------------------------------------------------------------------------------

    tweak = 0.

    do i = 1,7
       call rannum(rvar,1)
       tweak(i) = (rvar(1)-.5)*widths(i,oz)
    end do

    do i = 1,nintern(oz)
       call rannum(rvar,1)
       tweak(i+7)=(rvar(1)-.5)*inwidths(i,oz)
    end do

    return

  end subroutine gettweak



  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine writecrystal(outname,header)

    !----------------------------------------------------------------------------------------------
    !  Output the model crystal file to outname
    !----------------------------------------------------------------------------------------------

    character(len=100),intent(in)               :: outname, header

    integer                                     :: ia,ib,ic,iz,i,dunit,il

    logical                                     :: improp

    !----------------------------------------------------------------------------------------------
    ! Loop over crystal writing out the xyz, q and internals for each molecule 
    !----------------------------------------------------------------------------------------------

    dunit = open(outname)
    write(dunit,*)trim(header)
    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             do il = 1,num_loc
                iz = zocc(ia,ib,ic,il)
                improp = improper(qcrystal(ia,ib,ic,il))
                write(dunit,'(4i3," ",7(f8.5," "),L3," ")', advance = 'no')ia,ib,ic,il,  &
                     xyzcrystal(ia,ib,ic,il,:),&
                     as_array(qcrystal(ia,ib,ic,il)),improp
                do i = 1,nintern(iz)
                   write(dunit,'(f8.5," ")', advance='no')incrystal(ia,ib,ic,il,i)
                end do
                write(dunit,*)
             end do
          end do
       end do
    end do

    close(dunit)

    return

  end subroutine writecrystal


  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------!
  ! NOTE that B_adjust may want to be customised for the problem being tackled.             
  ! This version ought to work in general but uses ALL atoms which might (a) be slow and   
  ! (b) simply be not necessary and (c) some might be dummy atoms whose Biso is not really
  ! relevant.  At the very least you need to set B_aim (in MC input file)  to be the     
  ! B-factor you want to aim for.  Recall that B = 8*(pi^2)*U.                          
  ! The core of this routine is stolen from diffuse.f [1], so it should give Biso      
  ! which are then repeated when you run diffuse.                                     
  ! ([1] B.D.Butler and T.R.Welberry, J. Appl. Cryst., 25 (1992) 391-399))           
  !-----------------------------------------------------------------------------------------!

  subroutine  B_adjust(B_aim)

    integer                                 :: iz, il
    integer                                 :: ia,ib,ic,im,iat

    real                                    :: B_aim

    !-----------------------------------------------------------------------------------------!
    ! Internal to subroutine, same as in diffuse.f; unchanged
    !-----------------------------------------------------------------------------------------!

    real*8,dimension(3)   :: xyz
    real                  :: rpts,tot0,Biso,n,Biso_x,Biso_y,Biso_z,xpi2,Btot
    real*8                :: sumx,sumy,sumz,sumua,sumub,sumuc,sumaa,sumab,sumac

    xpi2=(4.*atan(1.))**2
    tot0  = 0.0
    Btot = 0.

    do iz = 1, total_num_zm
       do im = 1, maxval(mocc)
          do il = 1, num_loc
             do iat = 1,num_at(iz)

                !Below is from diffuse.f

                n=0.0
                sumx=0.d0
                sumy=0.d0
                sumz=0.d0
                sumua=0.d0
                sumub=0.d0
                sumuc=0.d0

                ! Loop over atoms
                ! Loop over all cells and do the sums ...

                do  ic=1,simsize(3)
                   do  ib=1,simsize(2)
                      do  ia=1,simsize(1)
                         if((mocc(ia,ib,ic,il).eq.im).and.(zocc(ia,ib,ic,il).eq.iz)) then
                            xyz(:)=carts(:,ia,ib,ic,il,iat)
                            n=n+1.0
                            sumx=sumx+xyz(1)
                            sumy=sumy+xyz(2)
                            sumz=sumz+xyz(3)
                            sumua=sumua+xyz(1)**2
                            sumub=sumub+xyz(2)**2
                            sumuc=sumuc+xyz(3)**2
                         end if
                      end do
                   end do
                end do

                ! Normalize, compute Debye factor, and get outa here ...

                if(n.gt.0) then
                   rpts=1./real(n)
                   xyz(1)=sumx*rpts
                   xyz(2)=sumy*rpts
                   xyz(3)=sumz*rpts
                   sumua=sumua*rpts
                   sumub=sumub*rpts
                   sumuc=sumuc*rpts
                   sumaa=xyz(1)**2
                   sumab=xyz(2)**2
                   sumac=xyz(3)**2
                   Biso_x=8.*xpi2*(sumua-sumaa)
                   Biso_y=8.*xpi2*(sumub-sumab)
                   Biso_z=8.*xpi2*(sumuc-sumac)
                   Biso=(Biso_x+Biso_y+Biso_z)/3.
                   tot0 = tot0 +1.0
                   Btot = Btot+Biso

                   ! Below is not from diffuse.f

                end if
             end do
          end do
       end do
    end do

    Btot          = Btot/tot0
    rpts          = Btot/(Btot+B_aim)
    !!    if (rpts.lt.0.1) rpts=0.1
    !!   if (rpts.gt.0.9) rpts=0.9

    if (.not. quiet)write(stdout,*)'B_aim, B_current :',B_aim,Btot

    sprcon = sprcon*(0.5+rpts)
    insprcon = insprcon*(0.5+rpts)
    cross_spr = cross_spr*(0.5+rpts)

    return

  end subroutine B_adjust

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine incupdate() 

    !---------------------------------------------------------------------------------------------!
    ! This changes the widths of the increments put on variables in an MC cycle.   This   
    ! is to optimise the acceptance/rejection ratio.  It is a (fairly) reliable routine  
    ! but can act up under extreme circumstances.                                       
    ! I am not sure how serious is the speed penalty of doing this; it may well be best to use it
    ! to get a 'reasonable' set of increments then fix them for most of the time. 
    !---------------------------------------------------------------------------------------------!

    integer                                 :: iz, iivar

    real,allocatable                        :: rejgrad(:), gradnum(:),gradtot(:)
    real                                    :: tratio, tcyc
    real,allocatable                        :: rejarray(:,:)

    !-------------------------------------------------------------------------!
    ! Step 1                                                                  
    ! We get the gradient of the rejection ratio with the size of the        
    ! increments being put on each variable.                                
    ! This means we loop over all z-matrices and all variables positioning 
    ! each z-matrix.  Put the gradients in rejarray(); note that the gradient
    ! is just a subtraction because in rise/run all have the same 'run' -- all
    ! use the same number of bins, div_max from divs                         
    !-------------------------------------------------------------------------!

    allocate(rejgrad(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('rejgrad(total_num_zm)                               ')
    allocate(gradnum(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('gradnum(total_num_zm)                               ')
    allocate(gradtot(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('gradtot(total_num_zm)                               ')
    allocate(rejarray(7+maxval(num_at),total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('rejarray(7+maxval(num_at),total_num_zm)             ')

    gradtot = 0.
    gradnum  = 0.

    do iz = 1,total_num_zm

       !-------------------------------------------------------------------------!
       ! 20 March 2007:  If not enough cycles done on that z-mat, skip over it...
       !-------------------------------------------------------------------------!

       tcyc = real(iacc(iz)+iup(iz)+irej(iz))
       if(tcyc.gt.400.0) then
          do iivar = 1,7           
             if(widths(iivar,iz).gt.0.00) then
                rejgrad(iz) = (reject(iz,iivar,div_max)-reject(iz,iivar,1))
                rejarray(iivar,iz) = rejgrad(iz)
                gradtot(iz) = gradtot(iz)+rejgrad(iz)
                gradnum(iz) = gradnum(iz)+1.
             end if
          end do

          do iivar = 7+1, 7+nintern(iz)          
             if(inwidths(iivar-7,iz).gt.0.0000) then
                rejgrad(iz) = (reject(iz,iivar,div_max)-reject(iz,iivar,1))
                rejarray(iivar,iz) = rejgrad(iz)
                gradtot(iz) = gradtot(iz)+rejgrad(iz)
                gradnum(iz) = gradnum(iz)+1.
             end if
          end do

          !---------------------------------------------------------------------------------!
          ! Step 2                                                                    
          ! We want all the variables to have the same impact on the acc/rej ratio.  
          ! on other words, we want all the values in rejarray() to be the same.    
          ! We also want the total rejection fraction to be around 50%.            
          ! The first of these considerations scales the widths() and inwidths() relative
          ! to each other, the second scales them all up and down together.             
          !                                                                            
          ! Get average gradient -- reuse rejgrad                                     
          ! Get the total rejection ratio                                            
          !---------------------------------------------------------------------------------!

          tratio =  tcyc/real(irej(iz))

          !-----------------------------------------------------------------------!
          ! We scale the widths according to the ratio of the gradient            
          ! to the average, average/grad, such                                   
          ! that variables with big gradients will have their effects reduced:  
          ! But we need some safeguards to stop things blowing up.             
          ! If a width is less than 0.000005, we assume it is zero and is meant
          ! to be zero and it is not changed.                                 
          !-----------------------------------------------------------------------!

          rejgrad(iz) = gradtot(iz)/gradnum(iz)
          if(gradnum(iz).lt.1)return

          do iivar = 1,7           
             if(widths(iivar,iz).gt.0.000005) then

                ! Reuse gradtot

                gradtot = rejgrad(iz)/rejarray(iivar,iz)*0.5*tratio
                if (gradtot(iz).lt.0.80)gradtot(iz) = 0.80
                if (gradtot(iz).gt.1.25)gradtot(iz) = 1.25
                widths(iivar,iz) = widths(iivar,iz)*gradtot(iz)
             end if
          end do

          do iivar = 1, nintern(iz)  
             if(inwidths(iivar,iz).gt.0.000005) then
                gradtot(iz) = rejgrad(iz)/rejarray(iivar+7,iz)*0.5*tratio
                if (gradtot(iz).lt.0.8)gradtot(iz) = 0.8
                if (gradtot(iz).gt.1.25)gradtot(iz) = 1.25
                inwidths(iivar,iz) = inwidths(iivar,iz)*gradtot(iz)
             end if
          end do
       end if
    end do

    ! Now we do the below in the main body to make sure it always happens
    iacc = 0
    iup = 0
    irej = 0
    reject = 0


    deallocate(rejgrad)
    deallocate(rejarray)
    deallocate(gradnum)
    deallocate(gradtot)

    return

  end subroutine incupdate

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_info(outname,outf)

    !----------------------------------------------------------------------------------------------
    !  This prints out some useful guff about the end simulation.  Histograms, mostly.
    !----------------------------------------------------------------------------------------------

    character(len=1)                              :: outf
    character(len=100),intent(in)                 :: outname
    character(len=107)                            :: xyzfile

    integer                                       :: ia2, ib2, ic2, typ,k,ddm,ddz,nbin,ibin
    integer                                       :: oa,ob,oc,oz,ol,om,stdou2,i,j
    integer,allocatable                           :: histogram(:,:)

    real                                          :: x2,y2,z2,dist
    real,dimension(3)                             :: xyz
    real,allocatable                              :: invals(:) 
    double precision,allocatable                  :: total_len(:),total_sq(:)
    double precision,allocatable                  :: shortest(:),longest(:)
    integer*8 ,allocatable                        :: counter(:)

    double precision,allocatable            :: count_in(:,:),total_val(:,:),total_val_sq(:,:)
    real,allocatable                        :: short_in(:,:), long_in(:,:)
    integer,allocatable                     :: histo_in(:,:,:)

    allocate(invals(maxval(nintern)),stat=ierr)
    if (ierr.ne.0) call allerr('invals(maxval(nintern))                             ')
    allocate(counter(ntyp),stat=ierr)
    if (ierr.ne.0) call allerr('counter(ntyp)                                       ')
    allocate(shortest(ntyp),stat=ierr)
    if (ierr.ne.0) call allerr('invals(shortest(ntyp)                               ')
    allocate(longest(ntyp),stat=ierr)
    if (ierr.ne.0) call allerr('longest(ntyp)                                       ')
    allocate(total_len(ntyp),stat=ierr)
    if (ierr.ne.0) call allerr('total_len(ntyp)                                     ')
    allocate(total_sq(ntyp),stat=ierr)
    if (ierr.ne.0) call allerr('total_sq(ntyp)                                      ')
    allocate(histogram(ntyp,25),stat=ierr)
    if (ierr.ne.0) call allerr('histogram(ntyp,25)                                  ')
    allocate(count_in(total_num_zm,maxval(num_at)),stat=ierr)
    if (ierr.ne.0) call allerr('count_in(total_num_zm,maxval(num_at))               ')
    allocate(total_val(total_num_zm,maxval(num_at)),stat=ierr)
    if (ierr.ne.0) call allerr('total_val(total_num_zm,maxval(num_at))              ')
    allocate(total_val_sq(total_num_zm,maxval(num_at)),stat=ierr)
    if (ierr.ne.0) call allerr('total_val_sq(total_num_zm,maxval(num_at))           ')
    allocate(short_in(total_num_zm,maxval(num_at)),stat=ierr)
    if (ierr.ne.0) call allerr('short_in(total_num_zm,maxval(num_at))               ')
    allocate(long_in(total_num_zm,maxval(num_at)),stat=ierr)
    if (ierr.ne.0) call allerr('long_in(total_num_zm,maxval(num_at))                ')
    allocate(histo_in(total_num_zm,maxval(num_at),25),stat=ierr)
    if (ierr.ne.0) call allerr('histo_in(total_num_zm,maxval(num_at),25)            ')

    if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
    if (.not. quiet)write(stdout,*)'Histogram generation routine'
    if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

    stdou2 = stdout

    if((outf.eq.'y').or.(outf.eq.'Y'))then
       call extension(outname,'summary',xyzfile)
       if (.not. quiet)write(stdout,*)'Writing histograms to ',trim(xyzfile)
       stdou2 = open(xyzfile)
    else
       if (.not. quiet)write(stdout,*)'Writing histograms to stdout'
    end if
    do i = 1,ntyp
       counter(i) = 0
       total_len(i) = 0.0  
       total_sq(i) = 0.0  
       shortest(i) = 1000.0
       longest(i) = 0.000
    end do
    nbin = 25
    do i = 1,ntyp
       do j = 1,nbin
          histogram(i,j) = 0
       end do
    end do
    total_num_zm=0

    do oa = 3,simsize(1)-2
       do ob = 3,simsize(2)-2
          do oc = 3,simsize(3)-2
             do ol = 1,num_loc
                om = mocc(oa,ob,oc,ol)
                oz = zocc(oa,ob,oc,ol)
                if (total_num_zm.lt.zocc(oa,ob,oc,ol))total_num_zm=zocc(oa,ob,oc,ol)
                do k = 1,conmol(ol,om,oz)
                   ia2 = wrap(oa+da(ol,oz,om,k),1)
                   ib2 = wrap(ob+db(ol,oz,om,k),2)
                   ic2 = wrap(oc+dc(ol,oz,om,k),3)
                   ddz = zocc(ia2,ib2,ic2,dl(ol,oz,om,k))
                   ddm = mocc(ia2,ib2,ic2,dl(ol,oz,om,k))
                   if( (ddz.eq.dz(ol,oz,om,k)).and. (ddm.eq.dm(ol,oz,om,k)) ) then
                      x2 = carts(1,oa,ob,oc,ol,oat(ol,oz,om,k))
                      y2 = carts(2,oa,ob,oc,ol,oat(ol,oz,om,k))
                      z2 = carts(3,oa,ob,oc,ol,oat(ol,oz,om,k))
                      xyz(1) = carts(1,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
                      xyz(2) = carts(2,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
                      xyz(3) = carts(3,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
                      xyz = xyz+da(ol,oz,om,k)*celltrans(1,:)+db(ol,oz,om,k)*celltrans(2,:)+  &
                           dc(ol,oz,om,k)*celltrans(3,:)
                      typ = ty(ol,oz,om,k)
                      dist = sqrt((xyz(1)-x2)*(xyz(1)-x2)+(xyz(2)-y2)*(xyz(2)-y2)+(xyz(3)-z2)*(xyz(3)-z2))
                      counter(typ) = counter(typ) + 1
                      total_len(typ) = total_len(typ) + dist 
                      total_sq(typ) = total_sq(typ) + dist*dist 
                      if (longest(typ).lt.dist) longest(typ) = dist
                      if (shortest(typ).gt.dist)shortest(typ) = dist
                   end if
                end do
             end do
          end do
       end do
    end do

    !----------------------------------------------------------------------------------------------
    ! Output info about the contact vector lengths (max, min, std dev...).
    !----------------------------------------------------------------------------------------------

    if (.not. quiet)write(stdout,*)'Writing out global averages for contact vectors...'
    write(stdou2,*)'----------------------------------------------------------------------------'
    write(stdou2,*)'Average contact vector lengths across the whole simulation...'
    write(stdou2,*)'Contact type, number, average length, shortest, longest, stdev'
    do k = 1,ntyp
       write(stdou2,'(i3,i8,4f8.4)')k,counter(k),total_len(k)/real(counter(k)),shortest(k),longest(k),             &
            sqrt(total_sq(k)/real(counter(k))-(total_len(k)/real(counter(k)))**2)
    end do

    !----------------------------------------------------------------------------------------------
    ! And generate some histograms...
    !----------------------------------------------------------------------------------------------

    do oa = 3,simsize(1)-2
       do ob = 3,simsize(2)-2
          do oc = 3,simsize(3)-2
             do ol = 1,num_loc
                om = mocc(oa,ob,oc,ol)
                oz = zocc(oa,ob,oc,ol)
                do k = 1,conmol(ol,om,oz)
                   ia2 = wrap(oa+da(ol,oz,om,k),1)
                   ib2 = wrap(ob+db(ol,oz,om,k),2)
                   ic2 = wrap(oc+dc(ol,oz,om,k),3)
                   ddz = zocc(ia2,ib2,ic2,dl(ol,oz,om,k))
                   ! ddz is the zmatrix which is actually  
                   ! there at the end of the contact vector
                   ddm = mocc(ia2,ib2,ic2,dl(ol,oz,om,k))
                   if( (ddz.eq.dz(ol,oz,om,k)).and. (ddm.eq.dm(ol,oz,om,k)) ) then
                      ! ddm is the occurance of that zmatrix  
                      ! which is actually there at the end of 
                      ! the contact vector.                   
                      x2 = carts(1,oa,ob,oc,ol,oat(ol,oz,om,k))
                      y2 = carts(2,oa,ob,oc,ol,oat(ol,oz,om,k))
                      z2 = carts(3,oa,ob,oc,ol,oat(ol,oz,om,k))
                      xyz(1) = carts(1,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
                      xyz(2) = carts(2,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
                      xyz(3) = carts(3,ia2,ib2,ic2,dl(ol,oz,om,k),dat(ol,oz,om,k))
                      xyz = xyz+da(ol,oz,om,k)*celltrans(1,:)+db(ol,oz,om,k)*celltrans(2,:)+  &
                           dc(ol,oz,om,k)*celltrans(3,:)
                      typ = ty(ol,oz,om,k)
                      dist = sqrt((xyz(1)-x2)*(xyz(1)-x2)+(xyz(2)-y2)*(xyz(2)-y2)+(xyz(3)-z2)*(xyz(3)-z2))
                      do ibin = 1, nbin
                         if( ( dist.gt.(shortest(typ)+(real(ibin-1))*((longest(typ)-shortest(typ))/real(nbin))) )    &
                              .and.(dist.lt.(shortest(typ)+(real( ibin ))*((longest(typ)-shortest(typ))/real(nbin)))))    &
                              histogram(typ,ibin)=histogram(typ,ibin)+1
                      end do
                   end if
                end do
             end do
          end do
       end do
    end do

    if (.not. quiet)write(stdout,*)'Writing out histograms for contact vectors...'
    do k = 1,ntyp
       write(stdou2,*)'----------------------------------------------------------------------------'
       write(stdou2,*)'Histogram for contact vector ',k
       write(stdou2,*)'----------------------------------------------------------------------------'
       do ibin = 1,nbin
          write(stdou2,'(1x,i3,f8.4,i15)')ibin,shortest(k)+(real(ibin)-0.5)*((longest(k)-shortest(k))/real(nbin)),histogram(k,ibin)
       end do
       write(stdou2,*)'----------------------------------------------------------------------------'
    end do

    !----------------------------------------------------------------------------------------------
    !  Then do the same for the internals....
    !----------------------------------------------------------------------------------------------

    count_in = 0.0
    total_val = 0.0  
    total_val_sq = 0.0  
    short_in = 1000.0
    long_in = -10000.000
    nbin = 25
    histo_in = 0

    if (.not. quiet)write(stdout,*)'Number of z-matrices is:',total_num_zm

    !----------------------------------------------------------------------------------------------
    ! Find min and max of each internal d.f.
    !----------------------------------------------------------------------------------------------

    do oa = 1,simsize(1)
       do ob = 1,simsize(1)
          do oc = 1,simsize(1)
             do ol = 1,num_loc
                oz = zocc(oa,ob,oc,ol)

                do typ = 1,nintern(oz)

                   count_in(oz,typ) =    count_in(oz,typ) + 1.0    
                   total_val(oz,typ) = total_val(oz,typ) + incrystal(oa,ob,oc,ol,typ)
                   total_val_sq(oz,typ) = total_val_sq(oz,typ) + (incrystal(oa,ob,oc,ol,typ))**2.0
                   if(short_in(oz,typ).gt.incrystal(oa,ob,oc,ol,typ))short_in(oz,typ)=incrystal(oa,ob,oc,ol,typ)
                   if(long_in(oz,typ).lt.incrystal(oa,ob,oc,ol,typ))long_in(oz,typ)=incrystal(oa,ob,oc,ol,typ)


                end do
             end do
          end do
       end do
    end do

    !----------------------------------------------------------------------------------------------
    ! Info about internals
    !----------------------------------------------------------------------------------------------

    if (.not. quiet)write(stdout,*)'Writing out global averages for internal degrees of freedom...'
    write(stdou2,*)'----------------------------------------------------------------------------'
    write(stdou2,*)'Average values of internals across the whole simulation...'
    write(stdou2,*)'zmat, deg. f., number, average, biggest, smallest, stdev'
    do oz = 1,total_num_zm
       do k = 1,nintern(oz)
          write(stdou2,'(i3,i3,f8.0,5f8.4)')oz,k,count_in(oz,k),total_val(oz,k)/count_in(oz,k),    &
               long_in(oz,k), short_in(oz,k),            &
               sqrt(total_val_sq(oz,k)/count_in(oz,k)-(total_val(oz,k)/count_in(oz,k))**2.0)
       end do
    end do

    !----------------------------------------------------------------------------------------------
    ! And again, generate some histograms...
    !----------------------------------------------------------------------------------------------

    do oa = 1,simsize(1)
       do ob = 1,simsize(1)
          do oc = 1,simsize(1)
             do ol = 1,num_loc
                oz = zocc(oa,ob,oc,ol)
                do typ = 1,nintern(oz)

                   dist = incrystal(oa,ob,oc,ol,typ)
                   do ibin = 1, nbin
                      if( ( dist.gt.(short_in(oz,typ)+(real(ibin-1))*((long_in(oz,typ)-short_in(oz,typ))/real(nbin))) )    &
                           .and.(dist.lt.(short_in(oz,typ)+(real(ibin))*((long_in(oz,typ)-short_in(oz,typ))/real(nbin)))))    &
                           histo_in(oz,typ,ibin)=histo_in(oz,typ,ibin)+1
                   end do
                end do

             end do
          end do
       end do
    end do

    if (.not. quiet)write(stdout,*)'Writing out histograms for internal degrees of freedom...'
    do oz = 1,total_num_zm
       do k = 1,nintern(oz)
          write(stdou2,*)'----------------------------------------------------------------------------'
          write(stdou2,*)'Histogram for internal ',k
          write(stdou2,*)'Of z-matrix            ',oz
          write(stdou2,*)'----------------------------------------------------------------------------'
          do ibin = 1,nbin
             write(stdou2,'(1x,i3,f8.4,i15)')ibin,short_in(oz,k)+      &
                  (real(ibin)-0.5)*((long_in(oz,k)-short_in(oz,k))/real(nbin)),histo_in(oz,k,ibin)
          end do
          write(stdou2,*)'----------------------------------------------------------------------------'
       end do
       write(stdou2,*)'----------------------------------------------------------------------------'
    end do

    if((outf.eq.'y').or.(outf.eq.'Y'))then
       close(stdou2)
    end if

    deallocate(count_in)
    deallocate(total_val)
    deallocate(total_val_sq)
    deallocate(short_in)
    deallocate(long_in)
    deallocate(histo_in)
    deallocate(invals)
    deallocate(counter)
    deallocate(shortest)
    deallocate(longest)
    deallocate(total_len)
    deallocate(total_sq)
    deallocate(histogram)
    return

  end subroutine get_info

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  readcrystal(inname,outname)

    character(len=100),intent(in)                       :: outname, inname
    character(len=100)                                  :: header

    integer                                             :: ia,ib,ic,il,ierr,im
    integer                                             :: iia,iib,iic,iil,iz,i,dunit
    real,allocatable                                    :: interns(:)

    real,dimension(3)                                   :: trans
    real,dimension(4)                                   :: ddd

    logical                                             :: improp

    type (quaternion)                                   :: q

    allocate(interns(maxval(nintern)),stat=ierr)
    if (ierr.ne.0) call allerr('interns(maxval(nintern))                            ')
    if (.not. quiet)write(stdout,*)'Model crystal to be read in from ',trim(inname)
    if(inname.eq.outname) then
       if (.not. quiet)write(stdout,*)'Note: If performing MC, file will be overwritten unless number of MC cycles is zero!!!!!'
    end if
    if (.not. quiet)write(stdout,*)'File header is:                    '
    dunit = open(inname)

    read(dunit,'(a100)')header
    if (.not. quiet)write(stdout,*)header

    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             do il = 1,num_loc
                read(dunit,'(4i3,1x,7(f8.5,1x),L3,1x)', advance = 'no',iostat=ierr)iia,iib,iic,iil,  &
                     trans(:),ddd(:),improp
                if(ierr.ne.0) then
                   write(stderr,*)'----------------------------------------------------------------------------'
                   write(stderr,*)'Problem reading in model crystal file:',trim(inname)
                   write(stderr,*)'Last good line reads:'
                   write(stderr,'(1x,4i3,1x,7(f8.5,1x),L3,1x)')iia,iib,iic,iil,&
                        trans(:),ddd(:),improp
                   write(stderr,*)'----------------------------------------------------------------------------'
                   write(stderr,*)'Exiting.'
                   write(stderr,*)'----------------------------------------------------------------------------'
                   stop
                end if
                q = as_quaternion(ddd)
                iz = zocc(iia,iib,iic,iil)
                im = mocc(iia,iib,iic,iil)
                improp = improper(q,improp)
                do i = 1,nintern(iz)
                   read(dunit,'(f8.5,1x)', advance='no',iostat=ierr)interns(i)
                   if(ierr.ne.0) then
                      write(stderr,*)'----------------------------------------------------------------------------'
                      write(stderr,*)'Problem reading in model crystal file:',trim(inname)
                      write(stderr,*)'Problem line begins:'
                      write(stderr,'(1x,4i3,1x,7(f8.5,1x),L3,1x)')iia,iib,iic,iil,&
                           trans(:),ddd(:),improp
                      write(stderr,*)'----------------------------------------------------------------------------'
                      write(stderr,*)'Exiting.'
                      write(stderr,*)'----------------------------------------------------------------------------'
                      stop
                   end if
                end do

                !----------------------------------------------------------------
                ! Now we have read in a single row of the file. Put it into arrays
                !----------------------------------------------------------------

                qcrystal(iia,iib,iic,iil)=q
                do i = 1,3
                   xyzcrystal(iia,iib,iic,iil,i)=trans(i)
                end do
                do i = 1,nintern(iz)
                   incrystal(iia,iib,iic,iil,i)=interns(i)
                end do

                !----------------------------------------------------------------
                ! Go to next row
                !----------------------------------------------------------------

                read(dunit,*)
             end do
          end do
       end do
    end do

    if(.not.quiet)write(stdout,*)'Crystal read in from ',inname

    deallocate(interns)
    close(dunit)

  end subroutine readcrystal

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine write_xyz(xyzfile,xtal,addcelltrans)

    !----------------------------------------------------------------
    ! Write out xyz coordinates of all atoms to some vast ASCII file.
    !----------------------------------------------------------------

    character(len=107),intent(in)                       :: xyzfile
    character(len=1), intent(in)                        :: addcelltrans

    integer                                             :: ia,ib,ic,il,im,iz,j,dunit,i

    real, dimension(3)                                  :: cvec,cvec2

    type (crystallography_object),intent(in)            :: xtal

    if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
    if (.not. quiet)write(stdout,*)'Writing Cartesian coordinates to ',trim(xyzfile)
    if((addcelltrans.eq.'y').or.(addcelltrans.eq.'Y')) then
       if (.not. quiet)write(stdout,*)'Coordinates will include unit cell translations'
    else
       if (.not. quiet)write(stdout,*)'Coordinates will not include unit cell translations'
    end if

    dunit = open(xyzfile)

    write(dunit,*)'ia  ib  ic  il  iz  im   iat         x         y         z'

    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             cvec(:) = as_cartesian(xtal,real((/(ia-1),(ib-1),(ic-1)/)))
             do il = 1,num_loc
                iz = zocc(ia,ib,ic,il)
                im = mocc(ia,ib,ic,il)
                do j = 1,num_at(iz)
                   do i = 1,3
                      if((addcelltrans.eq.'y').or.(addcelltrans.eq.'Y')) then
                         cvec2(i)=cvec(i)+carts(i,ia,ib,ic,il,j)
                      else
                         cvec2(i)=carts(i,ia,ib,ic,il,j)
                      end if
                   end do
                   write(dunit,'(6i4,i6,3f11.5)')ia,ib,ic,il,iz,im,j,(cvec2(i),i=1,3)
                end do
             end do
          end do
       end do
    end do

    close(dunit)

    return

  end subroutine write_xyz

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine write_frac(xyzfile,xtal)

    !----------------------------------------------------------------
    ! Write out xyz coordinates of all atoms to some vast ASCII file. 
    !----------------------------------------------------------------

    character(len=107),intent(in)                       :: xyzfile

    integer                                             :: ia,ib,ic,il,im,iz,j,dunit,i

    real,dimension(3)                                   :: cvec,fvec

    type (crystallography_object),intent(in)            :: xtal

    if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
    if (.not. quiet)write(stdout,*)'Writing fractional coordinates to ',trim(xyzfile)
    dunit = open(xyzfile)

    write(dunit,*)'ia  ib  ic  il  iz  im   iat         x         y         z'

    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             do il = 1,num_loc
                iz = zocc(ia,ib,ic,il)
                im = mocc(ia,ib,ic,il)
                do j = 1,num_at(iz)
                   do i = 1,3
                      cvec(i)=carts(i,ia,ib,ic,il,j)
                   end do
                   fvec(:)=as_fractional(xtal,cvec(:))
                   write(dunit,'(6i4,i6,3f10.5)')ia,ib,ic,il,iz,im,j,(fvec(i),i=1,3)
                end do
             end do
          end do
       end do
    end do

    close(dunit)

    return

  end subroutine write_frac

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine peanut(outname,oora)

    !----------------------------------------------------------------
    ! Produces correlation diagrams which plot correlation in the displacements of 
    ! molecular origins or average positions across atoms of contacting molecules 
    ! as a function of direction of displacement in each plane ab, ac, bc...     
    ! Hence there will be three diagrams for for each contact vector type -- note
    ! that there will be redundancy as some contact vectors will join the same pairs
    ! of molecules as others. 
    !----------------------------------------------------------------

    character(len=1), intent(in)                  :: oora
    character(len=4)                              :: kip
    character(len=112)               ,dimension(3):: xyzfile
    character(len=100), intent(in)                :: outname

    integer                                       :: ia2, ib2, ic2, typ,k,ddm,ddz,otyp
    integer                                       :: oa,ob,oc,om,oz,ol,icorr_angle,stdou2
    integer                                       :: stdou1,stdou3
    integer                                       :: ia,ib,ic,im,iz,il,i,j,iflag

    double precision, dimension(3)                :: ovec,dvec,avg
    double precision,allocatable                  :: globeave(:,:,:,:)
    double precision,allocatable                  :: globeavecount(:,:,:)
    real,allocatable                              :: working_array(:,:,:,:,:)
    real                                          :: avgcount,corr_angle
    real                                          :: unit_x,unit_y,ovec_proj
    real                                          :: dvec_proj,pi,corr1,corr2
    real                                          :: corr3,corr 
    double precision,dimension(360,3)             :: sumxy360,sumxx360
    double precision                              :: n360
    double precision,dimension(360,3)             :: sumyy360,sumx360,sumy360

    !----------------------------------------------------------------
    ! First, to do the peanut diagrams we are going to use the position of the origin atom
    ! or the average of all atoms as the position of the molecule -- depends on command line.
    !----------------------------------------------------------------

    allocate(globeave(3,num_loc,total_num_zm,maxval(mocc)),stat=ierr)
    if (ierr.ne.0) call allerr('globeave(3,num_loc,total_num_zm,maxval(mocc))       ')
    allocate(globeavecount(num_loc,total_num_zm,maxval(mocc)),stat=ierr)
    if (ierr.ne.0) call allerr('globeavecount(num_loc,total_num_zm,maxval(mocc))    ')
    allocate(working_array(3,simsize(1),simsize(2),simsize(3),num_loc),stat=ierr)
    if (ierr.ne.0) call allerr('working_array()                                     ')

    pi = 3.141592654

    if (.not. quiet) then
       write(stdout,*)'----------------------------------------------------------------------------'
       write(stdout,*)'Displacement correlation extraction subroutine'
       write(stdout,*)'Note that the planes are Cartesian, and do NOT relate directly to the'
       write(stdout,*)'real-space lattice planes in systems of less than orthorhombic symmetry.'
       write(stdout,*)'xy plane should be OK regardless.'
       write(stdout,*)
       write(stdout,*)'Files will have names of the form r_NUM_OA_XY_outname.out'
       write(stdout,*)'where NUM is the contact vector type and OA is o or a for origin or'
       write(stdout,*)'average. XY is the plane. (xy, yz, zx).  outname is the name specified in MC file'
       write(stdout,*)
       write(stdout,*)'To plot these numbers, fire up gnuplot and, if you are in the same directory'
       write(stdout,*)'as the output files, type (at the prompt, "gnuplot>") something like:'
       write(stdout,*)'gnuplot> set angles degrees; set polar; plot "r_01_o_xy_outname.out" '
       write(stdout,*)'(Note the quotes around the filename).'
    end if

    globeave = 0.0
    globeavecount = 0.0
    if ((oora.eq.'o').or.(oora.eq.'O')) then

       !----------------------------------------------------------------
       ! Extract origin and put it in working array.  Need to get the 'global' 
       ! average to use as 'origin' of displacements.
       !----------------------------------------------------------------

       do ia = 1,simsize(1)
          do ib = 1,simsize(2)
             do ic = 1,simsize(3)
                do il = 1,num_loc
                   iz = zocc(ia,ib,ic,il)
                   im = mocc(ia,ib,ic,il)
                   globeavecount(il,iz,im) = globeavecount(il,iz,im)+1.0
                   do j = 1,3
                      working_array(j,ia,ib,ic,il) = carts(j,ia,ib,ic,il,1)
                      globeave(j,il,iz,im) = globeave(j,il,iz,im)+ carts(j,ia,ib,ic,il,1) 
                   end do
                end do
             end do
          end do
       end do
    else

       !----------------------------------------------------------------
       ! Extract average and put that in the working array.  Note: we ARE using dummy atoms!!!
       !----------------------------------------------------------------

       do ia = 1,simsize(1)
          do ib = 1,simsize(2)
             do ic = 1,simsize(3)
                do il = 1,num_loc
                   iz = zocc(ia,ib,ic,il)
                   im = mocc(ia,ib,ic,il)
                   avg = 0.0
                   avgcount = 0.0
                   do i = 1,num_at(iz)
                      globeavecount(il,iz,im) = globeavecount(il,iz,im)+1.0
                      avgcount = avgcount+1.0
                      do j = 1,3
                         avg(j)=avg(j)+carts(j,ia,ib,ic,il,i)  
                         globeave(j,il,iz,im) = globeave(j,il,iz,im)+ carts(j,ia,ib,ic,il,i) 
                      end do
                   end do
                   do j = 1,3
                      working_array(j,ia,ib,ic,il) = avg(j)/avgcount
                   end do
                end do
             end do
          end do
       end do
    end if
    do il = 1,num_loc
       do iz = 1,total_num_zm
          do im= 1,maxval(mocc)
             do j = 1,3
                if(globeavecount(il,iz,im).gt.0.5) then
                   globeave(j,il,iz,im)=globeave(j,il,iz,im)/globeavecount(il,iz,im)
                else
                   globeave(j,il,iz,im)=99999.9
                end if
             end do
          end do
       end do
    end do

    if (.not. quiet)write(stdout,*)'Working arrays filled.'

    !----------------------------------------------------------------
    ! OK, now we should have a single atom's coordinates in the array working-array(), which will 
    ! be the positions of the origin atom OR the average of all in the molecule, including dummies
    ! We also have the global average position of the origin of each molecule in globeave() 
    !
    ! First, we select an origin molecule.  We'll sum over                                
    ! an inner subset of the crystal in order to remove the need for wrapping arrays.    
    !
    ! Cycle over contact vector types --- slow and clumsy but                          
    ! conceptually simple.                                                            
    !----------------------------------------------------------------

    do otyp = 1,ntyp
       n360 = 0.0
       do icorr_angle = 1,360
          do j = 1,3
             sumxy360(icorr_angle,j) = 0.0
             sumxx360(icorr_angle,j) = 0.0
             sumyy360(icorr_angle,j) = 0.0
             sumx360(icorr_angle,j) = 0.0
             sumy360(icorr_angle,j) = 0.0
          end do
       end do
       do oa = 2,(simsize(1)-2),2
          do ob = 2,(simsize(2)-2),2
             do oc = 2,(simsize(3)-2),2
                do ol = 1,num_loc
                   iflag = 0
                   oz = zocc(oa,ob,oc,ol)
                   om = mocc(oa,ob,oc,ol)

                   !-------------------------------------------------
                   ! So the origin molecule is at oa ob oc ol oz om
                   ! Now, get the first contact vector coming off this molecule
                   !-------------------------------------------------

                   do k = 1,conmol(ol,om,oz)
                      if (iflag.eq.0) then
                         typ = ty(ol,oz,om,k)
                         if(typ.eq.otyp) then

                            !----------------------------------
                            ! Go ahead and do the calculation.
                            !----------------------------------

                            ia2 = oa+da(ol,oz,om,k)
                            ib2 = ob+db(ol,oz,om,k)
                            ic2 = oc+dc(ol,oz,om,k)
                            if((oa.eq.ia2).and.(ob.eq.ib2).and.(oc.eq.ic2).and.   &
                                  (ol.eq.dl(ol,oz,om,k)))then

                               !----------------------------------
                               ! do nothing
                               !----------------------------------

                            else
                               ddz = zocc(ia2,ib2,ic2,dl(ol,oz,om,k))

                               !----------------------------------
                               ! ddz is the zmatrix which is actually
                               ! there at the end of the contact vector
                               !----------------------------------

                               ddm = mocc(ia2,ib2,ic2,dl(ol,oz,om,k))
                               if( (ddz.eq.dz(ol,oz,om,k)).and. (ddm.eq.dm(ol,oz,om,k)) ) then

                                  !----------------------------------
                                  ! Then this contact vector is done
                                  ! for this molecule and should move on
                                  !----------------------------------

                                  iflag = 1

                                  !----------------------------------
                                  ! Origin and dest. coords are 
                                  !----------------------------------

                                  do j=1,3
                                     ovec(j)=working_array(j,oa,ob,oc,ol)
                                     dvec(j)=working_array(j,ia2,ib2,ic2,dl(ol,oz,om,k))
                                  end do

                                  !----------------------------------
                                  ! Then subtract from these their 'average' positions:
                                  !----------------------------------

                                  do j=1,3
                                     ovec(j)=ovec(j)-globeave(j,ol,oz,om)
                                     dvec(j)=dvec(j)-globeave(j,dl(ol,oz,om,k),ddz,ddm)
                                     if(globeave(j,ol,oz,om).ge.99999.0) then
                                        write(stderr,*)'Problem with global average.  Code GLOB1'
                                        STOP
                                     end if
                                     if(globeave(j,dl(ol,oz,om,k),ddz,ddm).ge.99999.0) then
                                        write(stderr,*)'Problem with global average.  Code GLOB2'
                                        STOP
                                     end if
                                  end do

                                  !----------------------------------
                                  ! So now we have two vectors displaced relative to their
                                  ! respective 'origins'.                                
                                  ! Now decompose them.                                 
                                  ! Use a Cartesian set of coordinates                 
                                  ! along x, y and z.  NOTE that this is NOT          
                                  ! the crystal axes in a low symmetry system.       
                                  ! Of course in fact they are already decomposed   
                                  ! as they are stored as Cartesian.               
                                  ! Consider each plane one at a time.            
                                  !----------------------------------

                                  !----------------------------------
                                  ! (1) The xy plane.
                                  !----------------------------------

                                  n360=n360+1.0
                                  do icorr_angle = 1,360
                                     corr_angle=real(icorr_angle)*pi/180.0 

                                     !----------------------------------
                                     ! Components of unit vector in direction corr_angle is:
                                     !----------------------------------

                                     unit_x = cos(corr_angle)
                                     unit_y = sin(corr_angle)

                                     !----------------------------------
                                     ! Length of projection of displacement onto this direction is:
                                     ! (Note, must keep the sign....)
                                     !----------------------------------

                                     ovec_proj = ovec(1)*unit_x+ovec(2)*unit_y
                                     dvec_proj = dvec(1)*unit_x+dvec(2)*unit_y

                                     !----------------------------------
                                     ! This is the pair of numbers whose
                                     ! correlation s to be calculated.
                                     !----------------------------------

                                     sumxy360(icorr_angle,1) = sumxy360(icorr_angle,1) + ovec_proj*dvec_proj
                                     sumxx360(icorr_angle,1) = sumxx360(icorr_angle,1) + ovec_proj**2
                                     sumyy360(icorr_angle,1) = sumyy360(icorr_angle,1) + dvec_proj**2
                                     sumx360(icorr_angle,1)  = sumx360(icorr_angle,1)  + ovec_proj
                                     sumy360(icorr_angle,1)  = sumy360(icorr_angle,1)  + dvec_proj

                                  end do

                                  !----------------------------------
                                  ! (2) The yz plane.
                                  !----------------------------------

                                  do icorr_angle = 1,360
                                     corr_angle=real(icorr_angle)*pi/180.0 
                                     unit_x = cos(corr_angle)
                                     unit_y = sin(corr_angle)
                                     ovec_proj = ovec(2)*unit_x+ovec(3)*unit_y
                                     dvec_proj = dvec(2)*unit_x+dvec(3)*unit_y
                                     sumxy360(icorr_angle,2)=sumxy360(icorr_angle,2)+ovec_proj*dvec_proj
                                     sumxx360(icorr_angle,2)=sumxx360(icorr_angle,2)+ovec_proj**2
                                     sumyy360(icorr_angle,2)=sumyy360(icorr_angle,2)+dvec_proj**2
                                     sumx360(icorr_angle,2)=sumx360(icorr_angle,2)+ovec_proj
                                     sumy360(icorr_angle,2)=sumy360(icorr_angle,2)+dvec_proj
                                  end do

                                  !----------------------------------
                                  ! (3) The zx plane.
                                  !----------------------------------

                                  do icorr_angle = 1,360
                                     corr_angle=real(icorr_angle)*pi/180.0 
                                     unit_x = cos(corr_angle)
                                     unit_y = sin(corr_angle)
                                     ovec_proj = ovec(3)*unit_x+ovec(1)*unit_y
                                     dvec_proj = dvec(3)*unit_x+dvec(1)*unit_y
                                     sumxy360(icorr_angle,3)=sumxy360(icorr_angle,3)+ovec_proj*dvec_proj
                                     sumxx360(icorr_angle,3)=sumxx360(icorr_angle,3)+ovec_proj**2
                                     sumyy360(icorr_angle,3)=sumyy360(icorr_angle,3)+dvec_proj**2
                                     sumx360(icorr_angle,3)=sumx360(icorr_angle,3)+ovec_proj
                                     sumy360(icorr_angle,3)=sumy360(icorr_angle,3)+dvec_proj
                                  end do
                               end if
                            end if
                         end if
                      end if
                   end do
                end do
             end do
          end do
       end do

       if(otyp.lt.10)write(kip,'(a3,i1.1)')'000',otyp
       if(otyp.ge.10)write(kip,'(a2,i2.2)')'00',otyp
       if(otyp.ge.100)write(kip,'(a1,i3.3)')'0',otyp
       if(otyp.ge.1000)write(kip,'(i4.4)')otyp

       write(xyzfile(1),'(a2,a4,a1,a1,a4,a100)')'r_',kip,'_',oora,'_xy_',outname
       write(xyzfile(2),'(a2,a4,a1,a1,a4,a100)')'r_',kip,'_',oora,'_yz_',outname
       write(xyzfile(3),'(a2,a4,a1,a1,a4,a100)')'r_',kip,'_',oora,'_zx_',outname
       stdou1 = open(xyzfile(1))           
       stdou2 = open(xyzfile(2))           
       stdou3 = open(xyzfile(3))           

       !----------------------------------
       ! So now calc correlation function.
       !----------------------------------

       write(stdou1,*)'#---------------------------------------------------------'
       write(stdou1,*)'#xy plane, contact vector type: ',otyp
       write(stdou1,*)'# Angle to         correlation'
       write(stdou1,*)'# x axis           coefficient'
       do icorr_angle = 1,360
          corr1 = sumxy360(icorr_angle,1)-(sumx360(icorr_angle,1)*sumy360(icorr_angle,1))/n360
          corr2 = sumxx360(icorr_angle,1)-sumx360(icorr_angle,1)**2/n360
          corr3 = sumyy360(icorr_angle,1)-sumy360(icorr_angle,1)**2/n360
          corr = corr1/sqrt(corr2*corr3)
          write(stdou1,*)icorr_angle,corr
       end do

       write(stdou2,*)'#---------------------------------------------------------'
       write(stdou2,*)'#yz plane, contact vector type: ',otyp
       write(stdou2,*)'# Angle            correlation'
       write(stdou2,*)'# to y axis        coefficient'
       do icorr_angle = 1,360
          corr1 = n360*sumxy360(icorr_angle,2)-sumx360(icorr_angle,2)*sumy360(icorr_angle,2)
          corr2 = n360*sumxx360(icorr_angle,2)-sumx360(icorr_angle,2)*sumx360(icorr_angle,2)
          corr3 = n360*sumyy360(icorr_angle,2)-sumy360(icorr_angle,2)*sumy360(icorr_angle,2)
          corr = corr1/sqrt(corr2*corr3)
          write(stdou2,*)icorr_angle,corr
       end do

       write(stdou3,*)'#---------------------------------------------------------'
       write(stdou3,*)'#zx plane, contact vector type: ',otyp
       write(stdou3,*)'# Angle           correlation'
       write(stdou3,*)'# to z axis       coefficient'
       do icorr_angle = 1,360
          corr1 = n360*sumxy360(icorr_angle,3)-sumx360(icorr_angle,3)*sumy360(icorr_angle,3)
          corr2 = n360*sumxx360(icorr_angle,3)-sumx360(icorr_angle,3)*sumx360(icorr_angle,3)
          corr3 = n360*sumyy360(icorr_angle,3)-sumy360(icorr_angle,3)*sumy360(icorr_angle,3)
          corr = corr1/sqrt(corr2*corr3)
          write(stdou3,*)icorr_angle,corr
       end do

       close(stdou1)
       close(stdou2)
       close(stdou3)
    end do
    deallocate(globeave)
    deallocate(globeavecount)
    deallocate(working_array)

    return

  end subroutine peanut

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine pairs(outname)

    !--------------------------------------------------------------------------
    ! This simply puts out to files the variables for pairs of molecules separated by  
    ! each contact vector.  These can then be plotted, correlations looked for etc.
    !--------------------------------------------------------------------------

    character(len=4)                              :: kip
    character(len=4)                              :: lip
    character(len=116)                            :: pairfile
    character(len=100),intent(in)                 :: outname

    integer                                       :: ia2, ib2, ic2, typ,k,ddm,ddz,otyp
    integer                                       :: oa,ob,oc,om,oz,ol,dunit
    integer                                       :: i,iflag,iiflag

    real,dimension(3)                             :: x1,x2

    type (quaternion)                             :: q1,q2    

    if(.not.quiet) then
       write(stdout,*)'----------------------------------------------------------------------------'
       write(stdout,*)'Running subroutine to put out the molecular variables (x, q, i)'
       write(stdout,*)'for pairs of molecules connected by each type of contact vector.'     
       write(stdout,*)'The columns of numbers (a subset of the full crystal) can then be'       
       write(stdout,*)'imported into various programs (R, gnuplot, gnumeric, Excel).'
       write(stdout,*)'It assumes that a given contact vector always has the same types'
       write(stdout,*)'of z-matrices at its ends.  If this is not the case, interpret with care!'
       write(stdout,*)
       write(stdout,*)'Files will have names of the form pairs_X_NUM_outname.out'
       write(stdout,*)'where X is the location of the origin molecule and '
       write(stdout,*)'NUM is the contact vector type.'
       write(stdout,*)'----------------------------------------------------------------------------'
    end if

    do otyp = 1,ntyp
       if(otyp.lt.10)write(kip,'(a3,i1.1)')'000',otyp
       if(otyp.ge.10)write(kip,'(a2,i2.2)')'00',otyp
       if(otyp.ge.100)write(kip,'(a1,i3.3)')'0',otyp
       if(otyp.ge.1000)write(kip,'(i4.4)')otyp

       do ol = 1,num_loc
       if(ol.lt.10)write(lip,'(a3,i1.1)')'000',ol
       if(ol.ge.10)write(lip,'(a2,i2.2)')'00',ol
       if(ol.ge.100)write(lip,'(a1,i3.3)')'0',ol
       if(ol.ge.1000)write(lip,'(i4.4)')ol

          iiflag = 0
          write(pairfile,'(a6,a4,a1,a4,a1,a100)')'pairs_',lip,'_',kip,'_',outname
          dunit = open(pairfile)
          do oa = 2,(simsize(1)-2),2
             do ob = 2,(simsize(2)-2),2
                do oc = 2,(simsize(3)-2),2
                   q1 = qcrystal(oa,ob,oc,ol)
                   do i = 1,3
                      x1(i) = xyzcrystal(oa,ob,oc,ol,i)
                   end do
                   iflag = 0
                   oz = zocc(oa,ob,oc,ol)
                   om = mocc(oa,ob,oc,ol)

                   !------------------------------------------------------
                   ! So the origin molecule is at oa ob oc ol oz om
                   ! Now, get the first contact vector coming off this molecule
                   !------------------------------------------------------

                   do k = 1,conmol(ol,om,oz)
                      if (iflag.eq.0) then
                         typ = ty(ol,oz,om,k)
                         if(typ.eq.otyp) then

                            ! Go ahead and get the numbers.

                            ia2 = oa + da(ol,oz,om,k)
                            ib2 = ob + db(ol,oz,om,k)
                            ic2 = oc + dc(ol,oz,om,k)
                            ddz = zocc(ia2,ib2,ic2,dl(ol,oz,om,k))
                            if((oa.eq.ia2).and.(ob.eq.ib2).and.(oc.eq.ic2)    &   
                                 .and.(ol.eq.dl(ol,oz,om,k))) then

                               !do nothing

                            else

                               ! ddz is the zmatrix which is actually
                               ! there at the end of the contact vector

                               ddm = mocc(ia2,ib2,ic2,dl(ol,oz,om,k))
                               if( (ddz.eq.dz(ol,oz,om,k)).and. (ddm.eq.dm(ol,oz,om,k)) ) then

                                  ! Then this contact vector
                                  ! for this molecule us done.

                                  ! Write out the pair oa,ob,oc,ol and ia2,ib2,ic2,dl(ol,oz,om,k)

                                  q2 = qcrystal(ia2,ib2,ic2,dl(ol,oz,om,k))
                                  do i = 1,3
                                     x2(i) = xyzcrystal(ia2,ib2,ic2,dl(ol,oz,om,k),i)
                                  end do

                                  ! If not done already,
                                  ! write in a heading.

                                  if (iiflag.eq.0)then
                                     iiflag = 1
                                     write(dunit,'(7a8," ",a3," ")',advance='no')'x_1','y_1','z_1','q1_1','q1_2','q1_3','q1_4','im1'
                                     do i = 1,nintern(oz)
                                        write(dunit,'(a6,i2," ")', advance='no')'intl1_',i
                                     end do
                                     write(dunit,'(7a8," ",a3," ")',advance='no')'x_2','y_2','z_2','q2_1','q2_2','q2_3','q2_4','im2'
                                     do i = 1,nintern(ddz)
                                        write(dunit,'(a6,i2," ")', advance='no')'intl2_',i
                                     end do
                                     write(dunit,*)
                                  end if
                                  write(dunit,'(7f8.5," ",L3," ")',advance='no')(x1(i),i=1,3),(as_array(q1)),improper(q1)
                                  do i = 1,nintern(oz)
                                     write(dunit,'(f8.5," ")', advance='no')incrystal(oa,ob,oc,ol,i)
                                  end do
                                  write(dunit,'(7f8.5," ",L3," ")',advance='no')(x2(i),i=1,3),(as_array(q2)),improper(q2)
                                  do i = 1,nintern(ddz)
                                     write(dunit,'(f8.5," ")', advance='no')incrystal(ia2,ib2,ic2,dl(ol,oz,om,k),i)
                                  end do
                                  write(dunit,*)
                                  iflag = 1

                               end if
                            end if
                         end if
                      end if
                   end do
                end do
             end do
          end do
          close(dunit)
       end do
    end do

    return

  end subroutine pairs

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine calc_uij_new(xtal,outname)

    !------------------------------------------------------
    ! Calc Uij tensor for each atom that (might) exist in the unit cell.
    !------------------------------------------------------

    character(len=100)                            :: outname
    character(len=120)                            :: xyzfile
    character(len=6)                              :: lable,lablefinal
    character(len=2)                              :: lableshort

    integer                                       :: ia,ib,ic,im,iz,il,i,j,stdou2,inum

    real, dimension(3)                            :: ovec,fvec
    real, dimension(3)                            :: cellaxes, cellangles
    double precision,allocatable                  :: globeave(:,:,:,:,:)
    double precision,allocatable                  :: globeave2(:,:,:,:,:)
    double precision,allocatable                  :: globeavecount(:,:,:)
    double precision,allocatable                  :: uij(:,:,:,:,:)
    double precision                              :: Uiso

    type (crystallography_object)                 :: xtal

    allocate(globeave(3,num_loc,total_num_zm,maxval(mocc),maxval(num_at)),stat=ierr)
    if (ierr.ne.0) call allerr('globeave()                                          ')
    allocate(globeave2(3,num_loc,total_num_zm,maxval(mocc),maxval(num_at)),stat=ierr)
    if (ierr.ne.0) call allerr('globeave2()                                         ')
    allocate(globeavecount(num_loc,total_num_zm,maxval(mocc)),stat=ierr)
    if (ierr.ne.0) call allerr('globeavecount()                                     ')
    allocate(uij(num_loc,total_num_zm,maxval(mocc),maxval(num_at),6),stat=ierr)
    if (ierr.ne.0) call allerr('uij()                                               ')

    cellaxes(:)  = axes(xtal)
    cellangles(:)= angles(xtal)

    if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
    if (.not. quiet)write(stdout,*)'Routine to output a .cif file; the file should be readable by Mercury if '
    if (.not. quiet)write(stdout,*)'nothing else, since that was used to test the routine.      '

    !------------------------------------------------------
    ! To keep it simple, first loop over the molecules
    ! and calculate the mean position of each atom.
    !------------------------------------------------------

    globeavecount = 0.0
    globeave = 0.0
    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             do il = 1,num_loc
                iz = zocc(ia,ib,ic,il)
                im = mocc(ia,ib,ic,il)
                globeavecount(il,iz,im) = globeavecount(il,iz,im) + 1.0
                do i = 1,num_at(iz)
                   do j = 1,3
                      globeave(j,il,iz,im,i) = globeave(j,il,iz,im,i) + carts(j,ia,ib,ic,il,i) 
                   end do
                end do
             end do
          end do
       end do
    end do

    do il = 1, num_loc
       do iz = 1, total_num_zm
          do im = 1,maxval(mocc)
             if(globeavecount(il,iz,im).gt.0.5) then
                do i = 1,num_at(iz)
                   do j = 1,3
                      globeave(j,il,iz,im,i) = globeave(j,il,iz,im,i) / globeavecount(il,iz,im)
                   end do
                   ovec(:) = globeave(:,il,iz,im,i)
                   fvec(:) = as_fractional(xtal,ovec(:))
                   do j = 1,3
                      globeave2(j,il,iz,im,i) = fvec(j) * cellaxes(j)
                   end do
                end do
             end if
          end do
       end do
    end do

    !------------------------------------------------------
    ! Now the mean position of each atom in the unit cell is in globeave().
    ! and in the a*a coords used to define the Uij the means are in globeave2()
    ! (Trueblood et al, Acta. Crystallogr. A52 (1996) 770-781)
    !------------------------------------------------------

    uij = 0.0

    !------------------------------------------------------
    ! Now loop over the crystal again, calculating the various sums
    !------------------------------------------------------

    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             do il = 1,num_loc
                iz = zocc(ia,ib,ic,il)
                im = mocc(ia,ib,ic,il)
                do i = 1,num_at(iz)

                   !------------------------------------------------------
                   ! first, get the deviations from the average
                   !------------------------------------------------------

                   do j = 1,3
                      ovec(j) =  carts(j,ia,ib,ic,il,i) 
                   end do
                   fvec(:) = as_fractional(xtal,ovec(:))
                   do j = 1,3
                      ovec(j) = fvec(j) * cellaxes(j) - globeave2(j,il,iz,im,i) 
                   end do
                   uij(il,iz,im,i,1) = uij(il,iz,im,i,1) + ovec(1)**2
                   uij(il,iz,im,i,2) = uij(il,iz,im,i,2) + ovec(2)**2
                   uij(il,iz,im,i,3) = uij(il,iz,im,i,3) + ovec(3)**2
                   uij(il,iz,im,i,4) = uij(il,iz,im,i,4) + (ovec(2)*ovec(1)) 
                   uij(il,iz,im,i,5) = uij(il,iz,im,i,5) + (ovec(2)*ovec(3)) 
                   uij(il,iz,im,i,6) = uij(il,iz,im,i,6) + (ovec(3)*ovec(1)) 
                end do
             end do
          end do
       end do
    end do

    do il = 1, num_loc
       do iz = 1, total_num_zm
          do im = 1,maxval(mocc)
             if(globeavecount(il,iz,im).gt.0.5) then
                do i = 1,num_at(iz)
                   do j = 1,6
                      uij(il,iz,im,i,j) = uij(il,iz,im,i,j) / globeavecount(il,iz,im)
                   end do
                end do
             end if
          end do
       end do
    end do

    call  extension(outname,'cif',xyzfile)
    if (.not. quiet)write(stdout,*)'Spacegroup will be written as P1, .cif file is ',trim(xyzfile)
    stdou2 = open(xyzfile)
    write(stdou2,'(a)')'###########################################################################'
    write(stdou2,'(a)')'#                                                                         #'
    write(stdou2,'(a)')'#                 CIF from ZMC, version as of March 2026                  #'
    write(stdou2,'(a)')'#                                                                         #' 
    write(stdou2,'(a)')'###########################################################################'
    write(stdou2,'(a)')'#                                                                         #'
    write(stdou2,'(a)')'#  This CIF contains data generated by averaging across a model crystal   #'
    write(stdou2,'(a)')'#  as used in the program ZMC to model short-range order in crystalline   #'
    write(stdou2,'(a)')'#  materials.  The program assumes NO symmetry, so ALL sites in the       #'
    write(stdou2,'(a)')'#  unit cell are in the file, and space group is given as P1.             #'
    write(stdou2,'(a)')'#                                                                         #'
    write(stdou2,'(a)')'#  The routine keeps the first two characters of the atom label from      #'
    write(stdou2,'(a)')'#  the z-matrix and adds the ordinal number of the atom in the unit cell. #'
    write(stdou2,'(a)')'#                                                                         #'
    write(stdou2,'(a)')'#  This CIF generation routine has been tested by importing the output    #'
    write(stdou2,'(a)')'#  files into Mercury (http://www.ccdc.cam.ac.uk/products/mercury/)       #'
    write(stdou2,'(a)')'#  version 1.4.2 (Build 2).  No guarantees on the quality of the CIF      #'
    write(stdou2,'(a)')'#  files generated are given or implied.                                  #'
    write(stdou2,'(a)')'#                                                                         #'
    write(stdou2,'(a)')'#  For further information about ZMC, contact darren.goossens@gmail.com   #'
    write(stdou2,'(a)')'#                                                                         #'
    write(stdou2,'(a)')'###########################################################################'
    write(stdou2,'(a5,a52)')'data_',xyzfile
    write(stdou2,'(a)')"_audit_creation_method     'Hacked out from ZMC'"
    write(stdou2,'(a)')"_symmetry_space_group_name_H-M     'P 1'"
    write(stdou2,'(a)')'_symmetry_Int_Tables_number 1'
    write(stdou2,'(a21,f12.6)')'_cell_length_a       ',cellaxes(1)
    write(stdou2,'(a21,f12.6)')'_cell_length_b       ',cellaxes(2)
    write(stdou2,'(a21,f12.6)')'_cell_length_c       ',cellaxes(3)
    write(stdou2,'(a21,f12.6)')'_cell_angle_alpha    ',cellangles(1)
    write(stdou2,'(a21,f12.6)')'_cell_angle_beta     ',cellangles(2)
    write(stdou2,'(a21,f12.6)')'_cell_angle_gamma    ',cellangles(3)
    write(stdou2,'(a)')'loop_'
    write(stdou2,'(a)')'    _atom_site_label'
    write(stdou2,'(a)')'    _atom_site_fract_x'
    write(stdou2,'(a)')'    _atom_site_fract_y'
    write(stdou2,'(a)')'    _atom_site_fract_z'
    write(stdou2,'(a)')'    _atom_site_U_iso_or_equiv'
    write(stdou2,'(a)')'    _atom_site_adp_type'
    inum=0
    do il = 1, num_loc 
       do iz = 1, total_num_zm
          do im = 1,maxval(mocc) 
             if(globeavecount(il,iz,im).gt.0.5) then
                do i = 1, num_at(iz)
                   inum=inum+1
                   lable=labels(zmat(iz),i)
                   lableshort=lable(1:2)
                   Uiso = 0.0
                   do j = 1, 3
                      ovec(j) = globeave(j,il,iz,im,i)
                      Uiso = Uiso + uij(il,iz,im,i,j)
                   end do
                   fvec(:) = as_fractional(xtal,ovec(:))
                   Uiso = Uiso / 3.0 
                   if (len_trim(lable).ge.2) then
                      if (inum.le.9) write(lablefinal,'(a2,i1)')lableshort,inum
                      if ((inum.ge.10).and.(inum.le.99))write(lablefinal,'(a2,i2)')lableshort,inum
                      if ((inum.ge.100).and.(inum.le.999))write(lablefinal,'(a2,i3)')lableshort,inum
                      if ((inum.ge.1000).and.(inum.le.9999))write(lablefinal,'(a2,i4)')lableshort,inum
                      if ((inum.ge.10000).and.(inum.le.99999))write(lablefinal,'(a2,i5)')lableshort,inum
                      if ((inum.ge.100000).and.(inum.le.999999))write(lablefinal,'(a2,i6)')lableshort,inum
                      if ((inum.ge.1000000).and.(inum.le.9999999))write(lablefinal,'(a2,i7)')lableshort,inum
                   else
                      if (inum.le.9) write(lablefinal,'(a1,i1)')lableshort,inum
                      if ((inum.ge.10).and.(inum.le.99))write(lablefinal,'(a1,i2)')lableshort,inum
                      if ((inum.ge.100).and.(inum.le.999))write(lablefinal,'(a1,i3)')lableshort,inum
                      if ((inum.ge.1000).and.(inum.le.9999))write(lablefinal,'(a1,i4)')lableshort,inum
                      if ((inum.ge.10000).and.(inum.le.99999))write(lablefinal,'(a1,i5)')lableshort,inum
                      if ((inum.ge.100000).and.(inum.le.999999))write(lablefinal,'(a1,i6)')lableshort,inum
                      if ((inum.ge.1000000).and.(inum.le.9999999))write(lablefinal,'(a1,i7)')lableshort,inum
                   end if
                   write(stdou2,'(a10,4f9.5,a5)')lablefinal,(fvec(j),j=1,3),Uiso,' Uani'
                end do
             end if
          end do
       end do
    end do
    write(stdou2,'(a)')'loop_'
    write(stdou2,'(a)')'    _atom_site_aniso_label'
    write(stdou2,'(a)')'    _atom_site_aniso_U_11'
    write(stdou2,'(a)')'    _atom_site_aniso_U_22'
    write(stdou2,'(a)')'    _atom_site_aniso_U_33'
    write(stdou2,'(a)')'    _atom_site_aniso_U_12'
    write(stdou2,'(a)')'    _atom_site_aniso_U_23'
    write(stdou2,'(a)')'    _atom_site_aniso_U_13'
    inum=0
    do il = 1, num_loc
       do iz = 1, total_num_zm
          do im = 1,maxval(mocc) 
             if(globeavecount(il,iz,im).gt.0.5) then
                do i = 1, num_at(iz)
                   inum=inum+1
                   lable=labels(zmat(iz),i)
                   lableshort=lable(1:2)
                   if (len_trim(lable).ge.2) then
                      if (inum.le.9) write(lablefinal,'(a2,i1)')lableshort,inum
                      if ((inum.ge.10).and.(inum.le.99))write(lablefinal,'(a2,i2)')lableshort,inum
                      if ((inum.ge.100).and.(inum.le.999))write(lablefinal,'(a2,i3)')lableshort,inum
                      if ((inum.ge.1000).and.(inum.le.9999))write(lablefinal,'(a2,i4)')lableshort,inum
                      if ((inum.ge.10000).and.(inum.le.99999))write(lablefinal,'(a2,i5)')lableshort,inum
                      if ((inum.ge.100000).and.(inum.le.999999))write(lablefinal,'(a2,i6)')lableshort,inum
                      if ((inum.ge.1000000).and.(inum.le.9999999))write(lablefinal,'(a2,i7)')lableshort,inum
                   else
                      if (inum.le.9) write(lablefinal,'(a1,i1)')lableshort,inum
                      if ((inum.ge.10).and.(inum.le.99))write(lablefinal,'(a1,i2)')lableshort,inum
                      if ((inum.ge.100).and.(inum.le.999))write(lablefinal,'(a1,i3)')lableshort,inum
                      if ((inum.ge.1000).and.(inum.le.9999))write(lablefinal,'(a1,i4)')lableshort,inum
                      if ((inum.ge.10000).and.(inum.le.99999))write(lablefinal,'(a1,i5)')lableshort,inum
                      if ((inum.ge.100000).and.(inum.le.999999))write(lablefinal,'(a1,i6)')lableshort,inum
                      if ((inum.ge.1000000).and.(inum.le.9999999))write(lablefinal,'(a1,i7)')lableshort,inum
                   end if
                   write(stdou2,'(a10,6f10.6)')lablefinal, (uij(il,iz,im,i,j),j=1,6)
                end do
             end if
          end do
       end do
    end do

    close(stdou2)

    if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

    deallocate(globeave)
    deallocate(globeave2)
    deallocate(globeavecount)
    deallocate(uij)
    return

  end subroutine calc_uij_new

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine usage()  

    write(stdout,*)   
    write(stdout,*)'|---------------------------------------------------------------------|'
    write(stdout,*)'| Usage:                                                              |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| zmc [--option_1] [--option_2] ... [--option_n] infile [outfile]     |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| infile contains the parameters and additional filenames to run      |'
    write(stdout,*)'| the MC simulation.                                                  |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| outfile is the root name for most output; if not given the root     |'
    write(stdout,*)'| name will be infile (i.e., outfile = infile)                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| Options always begin with two dashes, and if the option can be      |'
    write(stdout,*)'| passed a value, the value must be indicated with an equals sign     |'
    write(stdout,*)'| and there must be no spaces.  E.g.: --summary=inline is right,      |'
    write(stdout,*)'| "--summary inline" is wrong, "--summary = inline" is wrong.         |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| Items enclosed in square brackets are optional. <text> denotes      |'
    write(stdout,*)'| that the user can specify any value.  For example                   |'
    write(stdout,*)'| --crystal[=<filename>] means you can choose "filename".             |'
    write(stdout,*)'| --summary[=inline] means you give "inline" or nothing.              |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| Options are:                                                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --crystal[=<filename>]                                              |'
    write(stdout,*)'|           Causes program to output a text file referred to as a     |'
    write(stdout,*)'|           "ZMC crystal file" which contains the variables for each  |'
    write(stdout,*)'|           molecule.  Can be analysed and/or read back in via        |'
    write(stdout,*)'|           --reread[=filename] to resume simulation.  If filename    |'
    write(stdout,*)'|           not given, will write to outfile.crystal.                 |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --diffuse[=<filename>]                                              |'
    write(stdout,*)'|           Causes program to output a file designed to be read into  |'
    write(stdout,*)'|           diffuse scattering calculation program DIFFUSE.           |'
    write(stdout,*)'|           Butler and  Welberry J.Appl.Cryst.(1992)25 391-399        |'
    write(stdout,*)'|           If filename not given, will write to outfile.diffuse.     |'
    write(stdout,*)'|           If outfile not given, will write to infile.diffuse        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --discus[=<filename>]                                               |'
    write(stdout,*)'|           Causes program to output a file designed to be read into  |'
    write(stdout,*)'|           DISCUS (http://discus.sourceforge.net/).  Also writes     |'
    write(stdout,*)'|           out a DISCUS macro file, .mac, and some info on usage.    |'
    write(stdout,*)'|           USE AT YOUR OWN RISK, and DO NOT TRUST THE DEFAULTS.      |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --summary[=inline]                                                  |'
    write(stdout,*)'|           Prints out some information about the model crystal at    |'
    write(stdout,*)'|           the end of the simulation.  If no argument passed, sends  |'
    write(stdout,*)'|           output to outfile.summary. inline writes output to stdout.|'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --pairs   Outputs variables for pairs of molecules connected by     |'
    write(stdout,*)'|           each contact vector to files pairs_L_C_outfile.out        |'
    write(stdout,*)'|           where L is the location number of the origin molecule     |'
    write(stdout,*)'|           C is contact number.  Many of these files will be         |'
    write(stdout,*)'|           duplicates.                                               |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --fracs   Outputs fractional coordinates of every atom in the model |'
    write(stdout,*)'|           crystal to outfile.fracs                                  |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --cartsn  Outputs Cartesian coordinates of every atom in the model  |'
    write(stdout,*)'|           to outfile.cartsn. Does not add on unit cell translations.|'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --cartst  Outputs Cartesian coordinates of every atom in the model  |'
    write(stdout,*)'|           to outfile.cartst. Adds on unit cell translations.        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --reread=<filename>                                                 |'
    write(stdout,*)'|           Reads in a crystal file and uses it as starting point for !'
    write(stdout,*)'|           the simulation.  Reads in from filename.                  !'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --corro   Calculates correlation ("peanut") diagrams and writes     |'
    write(stdout,*)'|           them to files r_C_o_XX_outname.out where C is the         |'
    write(stdout,*)'|           contact vector and XX is the plane (xy, yz, zx).          !'
    write(stdout,*)'|           Uses the position of the origin atom of each molecule     !'
    write(stdout,*)'|           in the calculation.  Many of these files will be          |'
    write(stdout,*)'|           duplicates.                                               |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --corra   As for --corro except file names are r_C_a_XX_outname.out |'
    write(stdout,*)'|           and it uses the average over all atoms in the molecule.   |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --energy  Outputs the MC energy of each molecule to outfile.energy  |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --cif     Outputs a crude CIF file to outfile.cif.                  |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --quiet   ZMC sends nothing to the screen except error messages and |'
    write(stdout,*)'|           summary information if --summary=inline is set.           |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --help    Prints this information and exits.                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --help2   Prints more help and exits.                               |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --version Prints version information and exits.                     |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --plot    Runs interactive plotting of the simulation; crude but    |'
    write(stdout,*)'|           sometimes useful. Exits after producing plots.            |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --getcontacts                                                       |'
    write(stdout,*)'|           Runs a simple routine to generate contact vectors         |'
    write(stdout,*)'|           based on parameters specified in keyword file.            |'
    write(stdout,*)'|           Then exits.                                               |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'|--modwave  runs a set of options to create modulations               |'
    write(stdout,*)'|           see README.modwave and modwave_tutorial_ZMC.ppt           |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --help, --help2 and --version can be run without infile             |'
    write(stdout,*)'| or outfile.out being specified. --quiet does not work with --help,  |'
    write(stdout,*)'| --help2 or --version. If --getcontacts is given, will not do MC     |'
    write(stdout,*)'| or plot.  If --plot is given, will not proceed to MC.               |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'|---------------------------------------------------------------------|'

  end subroutine usage

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  extension(filename,ext,newname)

  !----------------------------------------------------------------------------------------------
  ! Takes a filename, chops off whatever comes after the LAST stop, and adds the specified 
  ! extension. If no stop, just adds extension
  !----------------------------------------------------------------------------------------------

    implicit  none

    character*(*) :: filename,ext,newname
    integer iend

    do iend = len(filename),1,-1
       if(filename(iend:iend).eq.'.') goto 20
    end do
20  if(iend.lt.1) then
       newname = trim(filename) //'.'// ext
    else
       newname = filename(1:iend) // ext
    end if

    return

  end subroutine extension

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine plot(con_test)

    character(len=1),dimension(3)        :: view
    character(len=1)                     :: nmal,co1,co2,pcon,numatom,writenum,labcon,plmult
    character(len=12)                    :: psname
    character(len=80)                    :: ctext, f2name
    real                                 :: xxx,xxxx,yyy,yyyy,ca1,ca2,ca3
    real                                 :: dx,dy,dzz,x2,x3,y2,y3,z2,z3
    real                                 :: x2a,x3a,y2a,y3a
    real                                 :: xshift,fs, radius,dist
    integer                              :: icon, ii,mm,iil,iiz,im,ia,ib,ic,ia2,ib2,ic2,iiat,ranflag,slthick
    integer                              :: nummult,pltest, con_test
    type (zmatrix_object)                :: testzmat
    real,dimension(3)                    :: cvec,cvec2

    integer, dimension(3)                :: si

    integer, allocatable                 :: pltypes(:)
    real, allocatable                    :: conmin(:), conmax(:)
    real                                 :: c1(36),c2(36),c3(36)

    c1(1) = 0.0
    c2(1) = 0.0
    c3(1) = 1.0

    c1(2) = 0.0
    c2(2) = 1.0
    c3(2) = 0.0

    c1(3) = 1.0
    c2(3) = .0
    c3(3) = .0

    c1(4) = 1.0
    c2(4) = 0.0
    c3(4) = 1.0

    c1(5) = 1.0
    c2(5) = 1.0
    c3(5) = 0.0

    c1(6) = 0.0
    c2(6) = 1.0
    c3(6) = 1.0

    c1(7) = 0.5
    c2(7) = 0.0
    c3(7) = 1.0

    c1(8) = 0.5
    c2(8) = 1.0
    c3(8) = 0.0

    c1(9) = 1.0
    c2(9) = .5
    c3(9) = .0

    c1(10) = 1.0
    c2(10) = 0.5
    c3(10) = 1.0

    c1(11) = 1.0
    c2(11) = 1.0
    c3(11) = 0.5

    c1(12) = 0.5
    c2(12) = 1.0
    c3(12) = 0.5

    c1(13) = 0.5
    c2(13) = 0.5
    c3(13) = 1.0

    c1(14) = 1.0
    c2(14) = 0.5
    c3(14) = 0.5

    c1(15) = 0.5
    c2(15) = 1.0
    c3(15) = 1.0

    c1(16) = 0.0
    c2(16) = 0.5
    c3(16) = 1.0

    c1(17) =  .0
    c2(17) = 1.0
    c3(17) = 0.5

    c1(18) = 1.0
    c2(18) = 0.0
    c3(18) = 0.5

    do i = 1,18
       ca1  = 0.5*(c2(i)+c3(i))
       ca2  = 0.5*(c1(i)+c3(i))
       ca3  = 0.5*(c2(i)+c1(i))
       c1(i+18)=ca1
       c2(i+18)=ca2
       c3(i+18)=ca3
    end do

    ranflag = 0
    pltest = 0
    nummult = 0

    allocate(conmin(total_num_zm),stat=ierr)  
    if (ierr.ne.0) call allerr('conmin(total_num_zm)                                ')
    allocate(conmax(total_num_zm),stat=ierr)  
    if (ierr.ne.0) call allerr('conmax(total_num_zm)                                ')
    allocate(pltypes(ntyp),stat=ierr)
    if (ierr.ne.0) call allerr('pltypes(ntyp)                                        ')

    pltypes = 0 

112 write(stdout,*)'Number atoms  y or n'
    read(stdin,*)numatom
    if ((numatom.ne.'y').and.(numatom.ne.'n'))goto 112

113 write(stdout,*)'Write mol number on molecules y or n'
    read(stdin,*)writenum
    if ((writenum.ne.'y').and.(writenum.ne.'n'))goto 113

1191 write(stdout,*)'How thick (in unit cells) for slice down view dirn  '
    call read_buffer(stdin, buffer, finished, comment="!#", removeblanks=.TRUE.)
    read (buffer,*,iostat=ierr) slthick
    if ((slthick.lt.1).or.(ierr.ne.0)) then
       write(stderr,*)'Invalid number of cells!'
       goto 1191
    end if
    if(slthick.gt.6)slthick=6

111 write(stdout,*)'Plot contacts y or n'
    read(stdin,*)pcon
    if ((pcon.ne.'y').and.(pcon.ne.'n'))goto 111

    ! Here, test if we have the information we need -- is a contact vector
    ! list specified?

    if(con_test.eq.-1) then
       pcon = 'n'
       write(stdout,*)'No valid contact vector list specified. Not plotting contacts!'
    end if

    labcon = 'n'

    if (pcon.eq.'y') then
114    write(stdout,*)'Label contacts y or n'
       read(stdin,*)labcon
       if ((labcon.ne.'y').and.(labcon.ne.'n'))goto 114
    end if

    if (pcon.eq.'y') then
115    write(stdout,*)'Plot multiple contacts on same plot? y or n'
       read(stdin,*)plmult
       if ((plmult.ne.'y').and.(plmult.ne.'n'))goto 115
       if (plmult.eq.'y') then
116       write(stdout,*)'How many on the plot?'
          call read_buffer(stdin, buffer, finished, comment="!#", removeblanks=.TRUE.)
          read (buffer,*,iostat=ierr) nummult
          if ((nummult.lt.1).or.(ierr.ne.0).or.(nummult.gt.ntyp)) then
             write(stderr,*)'Invalid number of contacts!'
             goto 116
          end if
          pltypes = 0
117       write(stdout,*)'Which ones?'
          call read_buffer(stdin, buffer, finished, comment="!#", removeblanks=.TRUE.)
          read (buffer,*,iostat=ierr) (pltypes(i),i=1,nummult)
          if ((minval(pltypes(1:nummult)).lt.1).or.(ierr.ne.0).or.(maxval(pltypes).gt.ntyp)) then
             write(stderr,*)'Invalid contact(s) specified!'
             goto 117
          end if

       end if
    end if

    if (.not. quiet)write(stdout,*)'-------------------------------------------------------------------------------'
    if(.not.quiet)write(stdout,*)'Crystal will be viewed down x, y and z.'
    if (.not. quiet)write(stdout,*)'-------------------------------------------------------------------------------'

    !We can view the crystal down one of the three real-space Cartesian directions (NOT lattice vectors)

    view(1) = 'x'
    view(2) = 'y'
    view(3) = 'z'
    f2name = 'courier'
6546 FORMAT(4i2)
    do icon = 1,ntyp
       do ii = 1,3
          nmal = view(ii)
          !Name the output file
          if(icon.lt.10)then
             write(psname,'(a6,i1.1,a1,a1,a3)')'con_00',icon,'_',nmal,'.ps'
          else
             if(icon.lt.100)then
                write(psname,'(a5,i2.2,a1,a1,a3)')'con_0',icon,'_',nmal,'.ps'
             else
                write(psname,'(a4,i3.3,a1,a1,a3)')'con_',icon,'_',nmal,'.ps'
             end if  ! if icon < 100
          end if  !if icon < 10
          if(plmult.eq.'y') then
             write(psname,'(a8,a1,a3)')'con_mul_',nmal,'.ps'
          end if
          if(pcon.eq.'n') then
             write(psname,'(a8,a1,a3)')'no_cont_',nmal,'.ps'
          end if

          !-------------------------------------------------------------------------------------!
          ! End of naming of output postscript file (in psname)
          !-------------------------------------------------------------------------------------!
          ! Initialise the plot page.  In ps_init (ps_routines.f) I have changed the length of
          ! the psname to 12 chars when it used to be 32.
          ! Work out the min and max contact lengths for each z-matrix
          !-------------------------------------------------------------------------------------!

          DO iz = 1,total_num_zm
             mm = num_at(iz)
             testzmat = zmat(iz)
             conmax(iz) = maxval(parameter(testzmat,(/(i,i=2,mm)/),(/(1,i=2,mm)/)))+0.3
             if(conmax(iz).gt.1.8)conmax(iz)=1.8
             conmin(iz) = minval(parameter(testzmat,(/(i,i=2,mm)/),(/(1,i=2,mm)/)))
          end do

          call ps_init(psname)
          call ps_scale(15.,15.)
          call ps_translate(5.,5.)
          call ps_lw(0.05)
          call ps_rgbcolor(0.,0.,0.)
          radius = 0.1
          ! Draw a big circle to indicate the origin
          call ps_circ(0.,0.,1.5)


          si(1) = min(6,simsize(1))
          si(2) = min(6,simsize(2))
          si(3) = min(6,simsize(3))
          si(ii)=slthick

          DO ia = 1,si(1)
             DO ib = 1,si(2)
                DO ic = 1,si(3)
                   cvec(:) = as_cartesian(xtal,real((/(ia-1),(ib-1),(ic-1)/)))
                   DO iil = 1,num_loc
                      ! Now cycle over each molecule within a type (with the same z-matrix)
                      iiz = zocc(ia,ib,ic,iil)
                      im = mocc(ia,ib,ic,iil)
                      call ps_rgbcolor(c1(mod(iil,35)),c2(mod(iil,35)),c3(mod(iil,35)))
                      DO iat =1, num_at(iiz)
                         IF (nmal.eq.'z') then
                            xxx =carts(1,ia,ib,ic,iil,iat)+cvec(1)
                            yyy =carts(2,ia,ib,ic,iil,iat)+cvec(2)
                            co1 = 'x'
                            co2 = 'y'
                         END IF
                         IF (nmal.eq.'y') then
                            xxx =carts(3,ia,ib,ic,iil,iat) +cvec(3)
                            yyy = carts(1,ia,ib,ic,iil,iat)+cvec(1)
                            co1 = 'z'
                            co2 = 'x'
                         END IF
                         IF (nmal.eq.'x') then
                            xxx =carts(2,ia,ib,ic,iil,iat) +cvec(2)
                            yyy = carts(3,ia,ib,ic,iil,iat)+cvec(3)
                            co1 = 'y'
                            co2 = 'z'
                         END IF
                         call ps_circ(xxx,yyy,radius)
                         ! This bit numbers each atom.
                         if (numatom.eq.'y') then
                            fs = 0.39
                            call ps_setfont(f2name,fs)
                            write(ctext,'(i2)')iat
                            call ps_string(xxx+.05,yyy-.02,ctext,0.0,2)
                         end if
                         ! This bit draws in the intra molecular bonds -- not contacts
                         DO iiat = 1,num_at(iiz)
                            IF (nmal.eq.'z') then
                               xxxx =carts(1,ia,ib,ic,iil,iiat)+cvec(1)
                               yyyy =carts(2,ia,ib,ic,iil,iiat)+cvec(2)
                            END IF
                            IF (nmal.eq.'y') then
                               xxxx =carts(3,ia,ib,ic,iil,iiat)+cvec(3)
                               yyyy = carts(1,ia,ib,ic,iil,iiat)+cvec(1)
                            END IF
                            IF (nmal.eq.'x') then
                               xxxx =carts(2,ia,ib,ic,iil,iiat)+cvec(2)
                               yyyy = carts(3,ia,ib,ic,iil,iiat)+cvec(3)
                            END IF
                            dx = carts(1,ia,ib,ic,iil,iat)- carts(1,ia,ib,ic,iil,iiat)
                            dy = carts(2,ia,ib,ic,iil,iat)- carts(2,ia,ib,ic,iil,iiat)
                            dzz = carts(3,ia,ib,ic,iil,iat)- carts(3,ia,ib,ic,iil,iiat)
                            dist = sqrt(dx**2+dy**2+dzz**2)
                            ! If bond is in right range, draw it.
                            if ((dist.ge.conmin(iiz)-0.005).and.(dist.le.conmax(iiz)+0.005))   &
                                 call ps_line(xxx,yyy,xxxx,yyyy) 
                         END DO  !iiat
                      END DO   !iat
                      ! This bit writes unit cell and mol number next to each molecule.
                      call ps_rgbcolor(0.,0.,0.)
                      if (writenum.eq.'y') then
                         fs = 0.39
                         call ps_setfont(f2name,fs)
                         write(ctext,6546)im
                         xshift = float(ic)*.59
                         call ps_string(xxx-.09,yyy-xshift,ctext,0.0,12)
                      end if
                      ! This is where we are drawing contacts rather than just molecules
                      if((pcon.eq.'n').or. (pcon.eq.'N')) then
                         ! do nothing
                      else
                         ! plot contacts
                         do k = 1,conmol(iil,im,iiz)
                            pltest = 0
                            if((ty(iil,iiz,im,k).eq.icon).and.(plmult.eq.'n')) pltest = 2345
                            if((any(pltypes(1:nummult).eq.ty(iil,iiz,im,k))).and.(plmult.eq.'y')) pltest = 2345
                            if (pltest.eq.2345) then
                               pltest = 0
                               x2 = carts(1,ia,ib,ic,iil,oat(iil,iiz,im,k))+cvec(1)
                               y2 = carts(2,ia,ib,ic,iil,oat(iil,iiz,im,k))+cvec(2)
                               z2 = carts(3,ia,ib,ic,iil,oat(iil,iiz,im,k))+cvec(3)
                               ia2 = ia+da(iil,iiz,im,k)
                               IF((ia2.le.0).or.(ia2.gt.si(1))) goto 9087
                               ib2 = ib+db(iil,iiz,im,k)
                               IF((ib2.le.0).or.(ib2.gt.si(2))) goto 9087
                               ic2 = ic+dc(iil,iiz,im,k)
                               IF((ic2.le.0).or.(ic2.gt.si(3))) goto 9087
                               ! See what is actually there at the end of the contact we are plotting
                               if((zocc(ia2,ib2,ic2,dl(iil,iiz,im,k)).ne.dz(iil,iiz,im,k)).or.      &
                                    (mocc(ia2,ib2,ic2,dl(iil,iiz,im,k)).ne.dm(iil,iiz,im,k))) goto 9087
                               cvec2(:) = as_cartesian(xtal,real((/(ia2-1),(ib2-1),(ic2-1)/)))
                               x3 = carts(1,ia2,ib2,ic2,dl(iil,iiz,im,k),dat(iil,iiz,im,k))+cvec2(1)
                               y3 = carts(2,ia2,ib2,ic2,dl(iil,iiz,im,k),dat(iil,iiz,im,k))+cvec2(2)
                               z3 = carts(3,ia2,ib2,ic2,dl(iil,iiz,im,k),dat(iil,iiz,im,k))+cvec2(3)
                               IF (nmal.eq.'z') then
                                  x2a = x2
                                  y2a = y2
                                  x3a = x3
                                  y3a = y3
                               END IF
                               IF (nmal.eq.'y') then
                                  x2a = z2
                                  y2a = x2
                                  x3a = z3
                                  y3a = x3
                               END IF
                               IF (nmal.eq.'x') then
                                  x2a = y2
                                  y2a = z2
                                  x3a = y3
                                  y3a = z3
                               END IF
                               dist = sqrt((x3-x2)**2.0+(y3-y2)**2.0+(z3-z2)**2.0)
                               ! use this as a test compare with expected length
                               if(abs(dist-length(iil,iiz,im,k)).gt.0.005) ranflag = 99999
                               call ps_rgbcolor(0.,0.,0.)
                               !  This bit labels the contact 
                               if(labcon.eq.'y') then
                                  fs = 0.39
                                  call ps_setfont(f2name,fs)
                                  write(ctext,'(i3)')ty(iil,iiz,im,k)
                                  call ps_string((x2a+x3a)/2.0,(y2a+y3a)/2.0,ctext,0.0,12)
                               end if
                               call ps_rgbcolor(0.,0.,0.)
                               call ps_line(x2a,y2a,x3a,y3a)
                            end if ! if ty = icon
9087                        fs=fs 
                         end do ! k=1,conmol
                      end if  ! draw contacts?
                   END DO   !iil
                END DO   !ic
             END DO   !ib
          END DO   !ia
          ! Draw the axes and lable them
          call ps_rgbcolor(0.,0.,0.)
          call ps_lw(0.4)
          call ps_line(0.,0.,4.,0.)
          fs = 2.99
          call ps_setfont(f2name,fs)
          write(ctext,*)co1
          call ps_string(4.0,0.,ctext,0.0,2)
          call ps_line(0.,0.,0.,4.)
          write(ctext,*)co2
          call ps_string(0.0,4.0,ctext,0.0,2)
          call ps_page
          close(13)
          if (.not. quiet)write(stdout,*)
          if (.not. quiet)write(stdout,'(" Contact ",i3," viewed down ",a1," plotted to",a15)')icon,view(ii),psname
          if (.not. quiet)write(stdout,*)'-------------------------------------------------------------------------------'
          if (ranflag.eq.99999)then
             if (.not. quiet)write(stdout,*)' Calculated and stored contact lengths unequal.'
             if (.not. quiet)write(stdout,*)' If there is meant to be no randomness in the  '
             if (.not. quiet)write(stdout,*)' model then this is a problem...               '
             if (.not. quiet)write(stdout,*)' -----------------------------------------------------------'
             ranflag = 0
          end if
       end do !ii  
       if (pcon.eq.'n')exit
       if (plmult.eq.'y')exit
    end do ! ntyp

    deallocate(pltypes)
    deallocate(conmin)
    deallocate(conmax)

    stop

  end subroutine plot

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_w(init_w,flag)

    character(len=500)   :: buffer
    character(len=100)   :: keyword, keywords(3)
    character(len=100)   :: xyz, dummy
    character(len=2)     :: comp1
    character(len=1)     :: flag
    character(len=7)     :: prefix
    real, allocatable    :: width(:,:)
    integer              :: ierr,n,i,j,coord, zmat_int
    integer, allocatable :: fields_i(:)
    integer              :: numtimes
    real                 :: w, pi
    real                 :: init_w(:,:)

    pi = 3.141592654
    keywords(1) = 'XYZINITW'
    keywords(2) = 'QINITW'
    keywords(3) = 'ININITW'

    if (flag.eq.'w') then
       keywords(1) = 'XYZWIDTH'
       keywords(2) = 'QWIDTH'
       keywords(3) = 'INWIDTH'
    end if
    init_w = 0.0

    keyword = keywords(1)

    allocate(width(total_num_zm,3),stat=ierr)
    if (ierr.ne.0) call allerr('width(total_num_zm,3)                                ')

    width = 0.

    !--------------------------------------------------------------------------
    ! If first field is a letter, x, y, or z, then the second value is the
    ! init width for that direction.  If it is a number it applies to
    ! all of x, y and z                                               
    !
    ! First, scan for a case of setting default.
    !--------------------------------------------------------------------------

    n = num_value(key2, keyword)
    allocate(fields_i(n),stat=ierr)
    if (ierr.ne.0) call allerr('fields_i(n)                                          ')
    fields_i = 0
    numtimes = 0

    !--------------------------------------------------------------------------
    ! First, test all occurances of keyword to see if any do not have a character
    ! in the first field.  If there is one, it is used to set the default.
    ! BUT it gets more complicated -- what if there are multiple zmats?
    ! OK, there are three possibilities -- three fields follow the keyword, two fields or one.
    ! If only one field, it is must be real and is the initw for xy & z
    ! If two fields and the first is x, y or z then the second must be a real and
    ! must apply to all zmats
    !
    ! If two fields and the first is a integer, the second must be a real that applies to that zmatrix.
    ! If three fields, the first should be an integer, the second a latter and the third a real.
    ! So first, we want to work out which occurance of XYZINITWIDTH is followed by a single real.
    !--------------------------------------------------------------------------

    do i = 1,n
       buffer=' '
       buffer = get_value(key2,keyword,i)
       call get_num_fields(buffer, fields_i(i))
       ! If number of fields is one, use it!
       if (fields_i(i).eq.1) then
          numtimes = numtimes+ 1 
          ! It is a global one and use it
          read(buffer,*,iostat=ierr) w
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)
          width = w
       end if

    end do

    if(numtimes.gt.1) then
       if(.not.quiet) write(stdout,*)trim(keyword),' default set more than once.'
       if(.not.quiet) write(stdout,*)'Will use last occurance in keyword file.'
    end if

    !--------------------------------------------------------------------------
    ! OK, now that the default is established, we look for cases with two input fields.
    !--------------------------------------------------------------------------

    do i = 1,n
       if (fields_i(i).eq.2) then
          ! Process.  First, is the first field a number of a character?
          buffer=' '
          buffer = get_value(key2,keyword,i)
          read(buffer,*,iostat=ierr) xyz, dummy
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)

          xyz = .ucase.xyz
          if(.isletter.trim(xyz)) then      
             ! Then treat it as a global setting for a particular coordinate, x, y or z
             ! So read the value for w
             read(buffer,*,iostat=ierr) dummy, w
             if(ierr.ne.0) call fatal_keyword(keyword, buffer)

             coord = 0
             if (trim(xyz).eq.'X') coord = 1
             if (trim(xyz).eq.'Y') coord = 2
             if (trim(xyz).eq.'Z') coord = 3
             if(coord.eq.0)call fatal_keyword(keyword, buffer) 
             do j = 1, total_num_zm
                width(j,coord) = w
             end do
          else
             ! First field is a number in which case it must be an integer.
             ! See if we can read it into zmat_int
             read(buffer,*,iostat=ierr) zmat_int, w
             if(ierr.ne.0) call fatal_keyword(keyword, buffer)
             if ((zmat_int.lt.1).or.(zmat_int.gt.total_num_zm)) then
                call fatal_keyword(keyword,buffer)
             end if
             do j = 1,3
                width(zmat_int,j) = w
             end do
          end if ! if is letter
       end if ! if two fields
    end do

    !--------------------------------------------------------------------------
    ! OK, three input fields.
    !--------------------------------------------------------------------------

    do i = 1,n
       if (fields_i(i).eq.3) then
          ! Process.  First, fields should be int, char, real
          buffer=' '
          buffer = get_value(key2,keyword,i)
          read(buffer,*,iostat=ierr) zmat_int, xyz, w
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)
          ! Then, test what they've given us
          if ((zmat_int.lt.1).or.(zmat_int.gt.total_num_zm)) then
             call fatal_keyword(keyword,buffer)
          end if
          xyz = .ucase.xyz
          coord = 0
          if (trim(xyz).eq.'X') coord = 1
          if (trim(xyz).eq.'Y') coord = 2
          if (trim(xyz).eq.'Z') coord = 3
          if(coord.eq.0)call fatal_keyword(keyword, buffer)
          width(zmat_int,coord) = w
       end if
    end do

    !--------------------------------------------------------------------------
    ! now put width() into initw()
    !--------------------------------------------------------------------------

    do i = 1, total_num_zm
       do j = 1,3
          init_w(j,i) = width(i,j)
       end do
    end do

    deallocate(width)
    deallocate(fields_i)

    !--------------------------------------------------------------------------
    ! Now q widths.
    !--------------------------------------------------------------------------

    keyword = keywords(2)

    allocate(width(total_num_zm,4),stat=ierr)
    if (ierr.ne.0) call allerr('width(total_num_zm,4)                                ')

    width = 0.

    !--------------------------------------------------------------------------
    ! If first field is a string Q1...Q4 then the second value is the
    ! init width for that component.  If it is a number it applies to
    ! all of Q
    !--------------------------------------------------------------------------

    n = num_value(key2, keyword)
    allocate(fields_i(n),stat=ierr)
    if (ierr.ne.0) call allerr('fields_i(n)                                          ')
    fields_i = 0
    numtimes = 0
    do i = 1,n
       buffer=' '
       buffer = get_value(key2,keyword,i)
       call get_num_fields(buffer, fields_i(i))
       ! If number of fields is one, use it!
       if (fields_i(i).eq.1) then
          numtimes = numtimes+ 1 
          ! It is a global one and use it
          read(buffer,*,iostat=ierr) w
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)
          width = w
       end if

    end do

    if(numtimes.gt.1) then
       if(.not.quiet) write(stdout,*)trim(keyword),' default set more than once.'
       if(.not.quiet) write(stdout,*)'Will use last occurance in keyword file.'
    end if

    !--------------------------------------------------------------------------
    ! OK, now that the default is established, we look for cases with two input fields.
    !--------------------------------------------------------------------------

    do i = 1,n
       if (fields_i(i).eq.2) then
          ! Process.  First, is the first field a number of a character?
          buffer=' '
          buffer = get_value(key2,keyword,i)
          read(buffer,*,iostat=ierr) xyz, dummy
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)

          xyz = .ucase.xyz
          if(.isletter.xyz(1:1)) then      
             ! Then treat it as a global setting for a particular coordinate, Q1..Q4
             ! So read the value for w
             read(buffer,*,iostat=ierr) dummy, w
             if(ierr.ne.0) call fatal_keyword(keyword, buffer)

             coord = 0
             if (trim(xyz).eq.'Q1') coord = 1
             if (trim(xyz).eq.'Q2') coord = 2
             if (trim(xyz).eq.'Q3') coord = 3
             if (trim(xyz).eq.'Q4') coord = 4
             if(coord.eq.0)call fatal_keyword(keyword, buffer) 
             do j = 1, total_num_zm
                width(j,coord) = w
             end do
          else
             ! First field is a number in which case it must be an integer.
             ! See if we can read it into zmat_int
             read(buffer,*,iostat=ierr) zmat_int, w
             if(ierr.ne.0) call fatal_keyword(keyword, buffer)
             if ((zmat_int.lt.1).or.(zmat_int.gt.total_num_zm)) then
                call fatal_keyword(keyword, buffer)
             end if
             do j = 1,4
                width(zmat_int,j) = w
             end do
          end if ! if is letter
       end if ! if two fields
    end do

    !--------------------------------------------------------------------------
    ! OK, three input fields.
    !--------------------------------------------------------------------------

    do i = 1,n
       if (fields_i(i).eq.3) then
          ! Process.  First, fields should be int, char, real
          buffer=' '
          buffer = get_value(key2,keyword,i)
          read(buffer,*,iostat=ierr) zmat_int, xyz, w
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)
          ! Then, test what they've given us
          if ((zmat_int.lt.1).or.(zmat_int.gt.total_num_zm)) then
             call fatal_keyword(keyword, buffer) 
          end if
          xyz = .ucase.xyz
          coord = 0
          if (trim(xyz).eq.'Q1') coord = 1
          if (trim(xyz).eq.'Q2') coord = 2
          if (trim(xyz).eq.'Q3') coord = 3
          if (trim(xyz).eq.'Q4') coord = 4
          if(coord.eq.0)call fatal_keyword(keyword, buffer)
          width(zmat_int,coord) = w
       end if
    end do

    !--------------------------------------------------------------------------
    ! now put width() into initw()
    !--------------------------------------------------------------------------

    do i = 1, total_num_zm
       do j = 1,4
          init_w(j+3,i) = width(i,j)
       end do
    end do

    deallocate(width)
    deallocate(fields_i)

    !--------------------------------------------------------------------------
    ! Okay, internals now...
    !--------------------------------------------------------------------------

    keyword = keywords(3)

    !--------------------------------------------------------------------------
    ! Max possible number of internals is maxval(num_at), but
    ! we know how many there are - maxval(nintern)
    !--------------------------------------------------------------------------

    allocate(width(total_num_zm,maxval(nintern)),stat=ierr)
    if (ierr.ne.0) call allerr('width(total_num_zm,maxval(nintern))                  ')

    width = 0.

    n = num_value(key2, keyword)
    allocate(fields_i(n),stat=ierr)
    if (ierr.ne.0) call allerr('fields_i(n)                                          ')
    fields_i = 0
    numtimes = 0
    do i = 1,n
       buffer=' '
       buffer = get_value(key2,keyword,i)
       call get_num_fields(buffer, fields_i(i))
       ! If number of fields is one, use it!
       if (fields_i(i).eq.1) then
          numtimes = numtimes+ 1 
          ! It is a global one and use it
          read(buffer,*,iostat=ierr) w
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)
          width = w
       end if

    end do

    if(numtimes.gt.1) then
       if(.not.quiet) write(stdout,*)trim(keyword),' default set more than once.'
       if(.not.quiet) write(stdout,*)'Will use last occurance in keyword file.'
    end if

    !--------------------------------------------------------------------------
    ! OK, now that the default is established, we look for cases with two input fields.
    !--------------------------------------------------------------------------

    do i = 1,n
       if (fields_i(i).eq.2) then
          ! Process.  First, is the first field a number or a character?
          buffer=' '
          buffer = get_value(key2,keyword,i)
          read(buffer,*,iostat=ierr) xyz, dummy
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)

          xyz = .ucase.xyz
          if(.isletter.xyz(1:1)) then      
             ! Then treat it as a global setting for a particular coordinate, I1...IN
             ! So read the value for w
             read(buffer,*,iostat=ierr) dummy, w
             if(ierr.ne.0) call fatal_keyword(keyword, buffer)

             coord = 0
             ! Now here it gets tricky, because I don't know how many components there
             ! are in I.  So I have to do a loop from 1 to maxval(nintern)
             do k = 1, maxval(nintern)
                write(comp1,'(a1,i1)')'I',k
                if (trim(xyz).eq.comp1) coord = k
             end do

             if(coord.eq.0)call fatal_keyword(keyword, buffer) 
             do j = 1, total_num_zm
                width(j,coord) = w
             end do
          else
             ! First field is a number in which case it must be an integer.
             ! See if we can read it into zmat_int
             read(buffer,*,iostat=ierr) zmat_int, w
             if(ierr.ne.0) call fatal_keyword(keyword, buffer)
             if ((zmat_int.lt.1).or.(zmat_int.gt.total_num_zm)) then
                call fatal_keyword(keyword, buffer)
             end if
             do j = 1,maxval(nintern)
                width(zmat_int,j) = w
             end do
          end if ! if is letter
       end if ! if two fields
    end do

    !--------------------------------------------------------------------------
    ! OK, three input fields.
    !--------------------------------------------------------------------------

    do i = 1,n
       if (fields_i(i).eq.3) then
          ! Process.  First, fields should be int, char, real
          buffer=' '
          buffer = get_value(key2,keyword,i)
          read(buffer,*,iostat=ierr) zmat_int, xyz, w
          if(ierr.ne.0) call fatal_keyword(keyword, buffer)
          ! Then, test what they've given us
          if ((zmat_int.lt.1).or.(zmat_int.gt.total_num_zm)) then
             call fatal_keyword(keyword, buffer)
          end if
          xyz = .ucase.xyz
          coord = 0
          do k = 1, maxval(nintern)
             write(comp1,'(a1,i1)')'I',k
             if (trim(xyz).eq.comp1) coord = k
          end do
          if(coord.eq.0)call fatal_keyword(keyword, buffer)
          width(zmat_int,coord) = w
       end if
    end do

    !--------------------------------------------------------------------------
    ! now put width() into initw() 
    !--------------------------------------------------------------------------

    do i = 1, total_num_zm
       do j = 1,maxval(nintern)
          if(intern(j,1,i).eq.1) then
             ! It is bond length and leave it alone
             init_w(j+7,i) = width(i,j)
          else
             ! It is angle and convert to radians
             init_w(j+7,i) = width(i,j) * pi/180.0
          end if
       end do
    end do

    deallocate(width)
    deallocate(fields_i)

    if(.not.quiet) then
       prefix = 'Initial'
       if(flag.eq.'w') prefix = 'For  MC'
       do j = 1, total_num_zm

          write(stdout,'(1x,"Molecule type (z-matrix number): ",i3)')j
          write(stdout,'(1x,a7," widths for x,y,z are             ",3f9.5)') &
               prefix,(init_w(i,j),i=1,3)
          write(stdout,'(1x,a7," widths for q1,q2,q3,q4 are       ",4f9.5)') &
               prefix,(init_w(i,j),i=4,7)
          if(nintern(j).gt.0) then
             write(stdout,'(1x,a7," widths for internal d.f. are      ")',advance='no')prefix
             do k = 1,nintern(j)
                write(stdout,'(f8.5," ")', advance='no')init_w(k+7,j)
             end do
             write(stdout,*)
             write(stdout,*)'For internal d.f. that are angles, degrees have been converted to radians.'
          end if
          write(stdout,*)'----------------------------------------------------------------------------'
       end do
    end if

    return

  end  subroutine get_w

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_num_fields(buffer, fields)

    character(len=500) :: buffer
    character(len=100) :: test(100), comparator
    integer            :: ierr,j, fields

    comparator = '/.,mnbvcxz;lkjhgfdsaA][poiuytrewq=-0987654321?><MNBVCXZ:LKJHGFDSA}{POIUYTREWQ+_)(*^%$#@!mdeoityrnba9'
    do j = 1,100
       test(j) = comparator
    end do
    read(buffer,*,iostat=ierr)(test(j),j=1,100)
    do j =  1,100
       if(test(j).ne.comparator) fields = j
    end do

    return

  end subroutine get_num_fields

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine read_qxyz(quatern,translate)

    character(len=500)   :: buffer
    character(len=100)   :: dummy, qxyzname, keyword
    logical              :: improp
    integer              :: ierr,i,j,k,l,dunit
    integer              :: total_num_qxyz
    type (quaternion)    :: quatern(:,:,:)
    integer, allocatable :: zindex(:)
    real,dimension(4)    :: tempquat
    real                 :: translate(:,:,:,:)

    keyword = 'QXYZFILE'
    allocate(zindex(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('zindex(total_num_zm)                                 ')
    total_num_qxyz = num_value(key2,keyword)
    if (total_num_qxyz.eq.0) then
       if(.not.quiet)write(stdout,*)'No QXYZFILE given.  Program exiting.'
       stop
    end if
    if (total_num_qxyz.ne.total_num_zm) then
       if(.not.quiet)write(stdout,*)'Number of QXYZFILEs not same as number of z-matrices.'
       if(.not.quiet)write(stdout,*)'Exiting.'
       stop
    end if

    do i = 1,total_num_qxyz
       buffer = ' '
       buffer = get_value(key2,keyword,i)
       read(buffer,*,iostat=ierr)zindex(i), qxyzname
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)

       !---------------------------------------------------------------------
       ! OK, we have the name, let's see if we can open the file...
       !---------------------------------------------------------------------

       if (.not. quiet)write(stdout,'(" Quaternions and translations for z-matrix",i5)')zindex(i)
       dunit = open(qxyzname, status='old')

       if (.not. quiet)write(stdout,*)"are in file ",trim(qxyzname) 

       if(dunit.eq.-1) call file_does_not_exist(trim(qxyzname))

       !---------------------------------------------------------------------
       ! OK, the file is open.  Now, we have to read until we get an empty line.
       !---------------------------------------------------------------------

       do

          !---------------------------------------------------------------------
          ! j is 'which instance' and is read in here just to allow
          ! you to put that in explicitly and not rely on the ordering
          ! from the looping variables
          !---------------------------------------------------------------------

          call read_buffer(dunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
          if(finished) exit
          read(buffer,*,iostat=ierr)dummy,l,dummy,j,(tempquat(k),k=1,4),  &
               improp,(translate(l,zindex(i),j,k),k=1,3)
          if(ierr.ne.0) then
             write(stderr,*)'----------------------------------------------------------------------------'
             write(stderr,*)'Problem reading lines in .qxyz file: ',trim(qxyzname)
             write(stderr,*)'Offending line is: ',trim(buffer)
             write(stderr,*)'Program exiting.'
             write(stderr,*)'----------------------------------------------------------------------------'
             stop
          end if

          quatern(l,zindex(i),j) = tempquat

          improp=improper(quatern(l,zindex(i),j),improp)
          if(.not. quiet) then
             write(stdout,*)''
             write(stdout,'(1x,"                    Location): ",I4)') l
             write(stdout,'(1x," Sub structure (per location): ",I4)') j
             write(stdout,'(1x,"                   Quaternion: ",4F11.6)') as_array(quatern(l,zindex(i),j))
             write(stdout,'(1x,"                     Improper: ",L4)') improp
             write(stdout,'(1x,"              COM Translation: ",3F11.6)') translate(l,zindex(i),j,:)
          end if

       end do ! read through file

       close(dunit)

       if(.not. quiet) write(stdout,*)'----------------------------------------------------------------------------'

    end do ! loop over files

    deallocate(zindex)

    return

  end subroutine read_qxyz

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine read_in_contacts()

    character(len=500) :: buffer
    character(len=100) :: oname, keyword

    keyword = 'CONTACTFILE'
    buffer=''
    oname=''
    buffer = get_value(key2,keyword)
    read(buffer,*,iostat=ierr)oname
    if(ierr.ne.0) call fatal_keyword(keyword,buffer)

    call getcon(oname)

    return

  end subroutine read_in_contacts

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_maxval_mocc(maxval_mocc)

    character(len=100)     :: qxyzname
    character(len=500)     :: buffer
    integer                :: i,ierr, total_num_qxyz
    integer, intent(out)   :: maxval_mocc
    integer, allocatable   :: zindex(:)
    integer                :: l,j,dummy1, dummy2

    maxval_mocc = -9999
    allocate(zindex(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('zindex(total_num_zm)                                 ')

    total_num_qxyz = 0
    total_num_qxyz = num_value(key2,'QXYZFILE')
    if (total_num_qxyz.eq.0) then
       if(.not.quiet)write(stdout,*)'No QXYZFILE given, this may be a problem!'
    end if
    if (total_num_qxyz.ne.total_num_zm) then
       if(.not.quiet)write(stdout,*)'Number of QXYZFILEs not same as number of z-matrices.'
       if(.not.quiet)write(stdout,*)'This may be a problem!'
    end if
    do i = 1,total_num_qxyz
       buffer = get_value(key2,'QXYZFILE',i)
       read(buffer,*,iostat=ierr)zindex(i), qxyzname
       if(ierr.ne.0)then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'Illegal input in field ','QXYZFILE'
          write(stderr,*)'Offending input is: ',trim(buffer)
          write(stderr,*)'Program exiting.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if

       !----------------------------------------------------------------
       ! OK, we have the name, let's see if we can open the file...
       !----------------------------------------------------------------

       dunit = open(qxyzname, status='old')
       if(dunit.eq.-1) call file_does_not_exist(trim(qxyzname))

       !----------------------------------------------------------------
       ! OK, the file is open.  Now, we have to read until we get an empty line.
       !----------------------------------------------------------------

       do 
          call read_buffer(dunit, buffer, finished, comment="!#", removeblanks=.TRUE.)
          read(buffer,*,end=12,iostat=ierr)dummy1,l,dummy2,j
          if(ierr.ne.0) then
             write(stderr,*)'----------------------------------------------------------------------------'
             write(stderr,*)'Problem reading lines in .qxyz file: ',trim(qxyzname)
             write(stderr,*)'Offending line is: ',trim(buffer)
             write(stderr,*)'Program exiting.'
             write(stderr,*)'----------------------------------------------------------------------------'
             stop
          end if

          if (j.gt.maxval_mocc)maxval_mocc = j

       end do ! read through file

12     close(dunit)

    end do ! i=1,zm

    do i = 1,total_num_zm
       if(.not.(any(zindex.eq.i))) then
          write(stderr,*)'-----------------------------------------------'
          write(stderr,*)'QXYZ file zmatrix numbers illegal.  Check input file.'
          write(stderr,*)'Program exiting.'
          write(stderr,*)'-----------------------------------------------'
          stop
       end if

    end do
    if (maxval_mocc.ne.-9999) then
       if(.not. quiet)write(stdout,*)'Interrogation of qxyz files suggests there are'
       if(.not. quiet)write(stdout,'(" up to",i5," instances of a z-matrix on a location.")')maxval_mocc
       if(.not. quiet) write(stdout,*)'----------------------------------------------------------------------------'
    end if

    deallocate(zindex)

    return

  end subroutine get_maxval_mocc

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_new_contacts()

    character(len=100)                 :: oname, keyword
    character(len=500)                 :: buffer
    character(len=4)                   :: zmatext, text

    integer                            :: dunit
    integer                            :: typ
    integer                            :: iz,ol,iiz,iiat,om
    integer                            :: iat, oat,oom,ddm
    integer                            :: da,db,dc,dz,dm,dat,dl
    integer                            :: oz,j, ierr
    integer, allocatable               :: which_atoms_to_use(:,:)
    integer, allocatable               :: num_atoms_to_use(:)

    real                               :: vmin, vmax, separation
    real, dimension(3)                 :: dest, origin
    real, dimension(3)                 :: trans
    real, dimension(3,3)               :: celltrans
    real,allocatable                   :: cartscell(:,:,:,:,:)
    type (quaternion)                  :: q

    allocate(cartscell(num_loc,total_num_zm,maxval_mocc,maxval(num_at),3),stat=ierr)
    if (ierr.ne.0) call allerr('cartscell(etc etc etc)                              ')

    zmatext = 'ZMAT'

    if (exists(key2, 'NEWCONTACTS')) then
       buffer=''
       oname = ''
       buffer = get_value(key2,'NEWCONTACTS')
       read(buffer,*,iostat=ierr)oname
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       ! OK, we have the name, let's see if we can open the file...
       if (exists(oname)) call file_does_exist(trim(oname))
       dunit = open(oname)
    else
       write(stderr,*)'-----------------------------------------------'
       write(stderr,*)'Must supply a new contact file name using NEWCONTACTS'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'-----------------------------------------------'
       stop
    end if

    !----------------------------------------------------------------------------------------------
    ! OK, file is read in and opened.  now we need to find out which atoms of which zmats to use,
    ! and what the max and min lengths are.
    ! First, VMIN and VMAX
    !----------------------------------------------------------------------------------------------

    call get_required_real(vmin,'VMIN')
    call get_required_real(vmax,'VMAX')

    keyword = 'CONTATOMS'

    allocate(num_atoms_to_use(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('num_atoms_to_use(total_num_zm)                      ')
    num_atoms_to_use = 0
    if (exists(key2, keyword)) then
       do j = 1,num_value(key2,keyword)
          if (has_value(key2, keyword)) then
             buffer=' '
             buffer = get_value(key2,keyword,j)
             read(buffer,*,iostat=ierr)text,oz
             if (ierr.ne.0) call fatal_keyword(keyword, buffer)
             text = .ucase.text
             if(text.ne.zmatext)call fatal_keyword(keyword, buffer)
             call get_num_fields(buffer, ol)
             if((ol-2).gt.num_atoms_to_use(oz))num_atoms_to_use(oz)=(ol-2)
             if((ol-2).gt.num_at(oz))  call fatal_keyword(keyword, buffer)
          else
             write(stderr,*)keyword,' keyword given but it has no value!'
             write(stderr,*)'Program exiting.'
             write(stderr,*)'----------------------------------------------------------------------------'
             stop
          end if
       end do
    else
       write(stderr,*)'Value for ',keyword,' not specified explicitly.'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'----------------------------------------------------------------------------'
       stop
    end if

    allocate(which_atoms_to_use(total_num_zm,maxval(num_at)),stat=ierr)
    if (ierr.ne.0) call allerr('which_atoms_to_use(total_num_zm,maxval(num_at))     ')
    which_atoms_to_use = 0
    do j = 1,num_value(key2,keyword)
       buffer=' '
       buffer = get_value(key2,keyword,j)
       read(buffer,*,iostat=ierr)text,oz,(which_atoms_to_use(oz,i),i=1,num_atoms_to_use(oz))
       if (ierr.ne.0) call fatal_keyword(keyword, buffer)
       if (maxval(which_atoms_to_use(oz,:)).gt.num_at(oz)) call fatal_keyword(keyword, buffer)
    end do

    celltrans(1,:)=as_cartesian(xtal,real((/1,0,0/)))
    celltrans(2,:)=as_cartesian(xtal,real((/0,1,0/)))
    celltrans(3,:)=as_cartesian(xtal,real((/0,0,1/)))

    !----------------------------------------------------------------------------------------------
    ! The first thing is to fill up an array of Cartesian positions FOR A SINGLE
    ! UNIT CELL.
    !----------------------------------------------------------------------------------------------

    do ol = 1,num_loc
       do iz = 1,how_many_types(ol)
          oz = which_types(ol,iz)
          allocate(new_coords(3,(num_at(oz))),stat=ierr)
          if (ierr.ne.0) call allerr('new_coords(3,(num_at(oz)))    2         ')          
          call copy(z1,zmat(oz))
          do oom = 1,instances_types(ol,oz)
             om = which_mol(ol,oz,oom)
             ! So now we have defined the molecule of a certain z-matrix type on a certain location.
             q = quatern(ol,oz,om)
             trans = translate(ol,oz,om,1:3)
             new_coords = as_xyz(z1)
             rmat = as_rotmatrix(q)
             ! Pre-rotate the coordinates using this rotation matrix
             call rotate(rmat,new_coords)
             do j = 1,num_at(oz)
                do i = 1,3
                   cartscell(ol,iz,om,j,i)=new_coords(i,j)+trans(i)
                end do ! i = 1,3
             end do
             if(.not.quiet) then
                write(stdout,'(" Coordinates for location: ",i5)')ol
                write(stdout,'(" Occupied by zmatrix     : ",i5)')oz
                write(stdout,'(" Which is ordinal number : ",i5," on that location.")')iz
                write(stdout,'(" In orientation number   : ",i5," for that location.")')om
                do j = 1,num_at(oz)
                   write(stdout,'(i5,1x,3f7.3)')j,cartscell(ol,iz,om,j,1:3)
                end do
                write(stdout,*)'----------------------------------------------------'
             end if
          end do
         deallocate(new_coords)
       end do
    end do

    if(.not.quiet)write(stdout,*)'Initial Cartesian positions calculated.'
    if(.not.quiet)write(stdout,*)'----------------------------------------------------------------------------'

    !----------------------------------------------------------------------------------------------
    ! Now we loop over the unit cell and get the contact vectors.
    !
    ! We loop over all possible molecules in a single unit
    ! cell, and we look for atoms the appropriate distance
    ! away and satisfying the criteria put down in keyword file
    !----------------------------------------------------------------------------------------------

    typ = 0
    write(dunit,*)' ol  oz  om oat  da  db  dc  dl  dz  dm dat  length  type'
    do ol = 1,num_loc
       do iz = 1,how_many_types(ol)
          oz = which_types(ol,iz)

          do oom = 1,instances_types(ol,oz)
             om = which_mol(ol,oz,oom)
             do iat = 1, num_atoms_to_use(oz)
                oat = which_atoms_to_use(oz,iat)
                ! position of origin atom
                origin(1) = cartscell(ol,iz,om,oat,1)
                origin(2) = cartscell(ol,iz,om,oat,2)
                origin(3) = cartscell(ol,iz,om,oat,3)
                ! now we loop over neighbours looking for valid contact vectors
                do da=-3,3,1
                   do db=-3,3,1
                      do dc=-3,3,1
                         do dl = 1,num_loc
                            if((da.eq.0).and.(db.eq.0).and.(dc.eq.0).and.(dl.eq.ol))then
                               ! do nothing
                            else
                               do iiz = 1,how_many_types(dl)
                                  dz = which_types(dl,iiz)
                                  do ddm = 1,instances_types(dl,dz)
                                     dm = which_mol(dl,dz,ddm)
                                     do iiat = 1,num_atoms_to_use(dz)
                                        dat = which_atoms_to_use(dz,iiat)
                                        ! position of destination
                                        dest(1) =  cartscell(dl,iiz,dm,dat,1)
                                        dest(2) =  cartscell(dl,iiz,dm,dat,2)
                                        dest(3) =  cartscell(dl,iiz,dm,dat,3)
                                        ! separation, dding in unit cell shift
                                        dest = dest + da*celltrans(1,:) + db*celltrans(2,:) + dc*celltrans(3,:)
                                        separation = sqrt((dest(1)-origin(1))**2+      &
                                             (dest(2)-origin(2))**2+      &
                                             (dest(3)-origin(3))**2)
                                        if((separation.gt.vmin).and.(separation.lt.vmax)) then
                                           ! write out the vector
                                           typ = typ + 1
                                           if(separation.lt.1.5) then
                                              if(.not.quiet) then
                                                 write(stdout,*)'Contact is less than 1.5 A!!'
                                                 write(stdout,*)'Atom 1 ',origin
                                                 write(stdout,*)'Atom 2 ',dest
                                                 write(stdout,'(11i4,f8.5,i9)')ol,oz,om,oat,da,db,dc,dl,dz,dm,dat,separation,typ
                                              end if
                                           end if
                                           write(dunit,'(11i4,f8.5,i9)')ol,oz,om,oat,da,db,dc,dl,dz,dm,dat,separation,typ
                                        end if
                                     end do
                                  end do
                               end do
                            end if
                         end do
                      end do
                   end do
                end do ! oat loop
             end do
          end do
       end do
    end do

    close(dunit)

    if(.not.quiet) write(stdout,*)'Exiting normally.'
    if(.not.quiet) write(stdout,*)'----------------------------------------------------------------------------'
    stop

  end subroutine get_new_contacts

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  get_required_real(value_given, keywrd)

    character(len=500)          ::  buffer
    real, intent(out)           ::  value_given
    integer                     :: ierr
    character(len=*),intent(in) :: keywrd

    if (exists(key2, keywrd)) then
       if (has_value(key2, keywrd)) then
          buffer=' '
          buffer = get_value(key2,keywrd)
          read(buffer,*,iostat=ierr)value_given
          if(ierr.ne.0) call fatal_keyword(keywrd,buffer)
          if (.not. quiet)write(stdout,*)'Value for ',keywrd,' specified explicitly is: ',value_given
          if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
       else
          write(stderr,*)keywrd,' keyword given but it has no value!'
          write(stderr,*)'Program exiting.'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
    else
       write(stderr,*)'Value for ',keywrd,' not specified explicitly.'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'----------------------------------------------------------------------------'
       stop
    end if

    return

  end subroutine  get_required_real

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  get_required_integer(value_given, keywrd)

    character(len=500)          :: buffer    
    integer, intent(out)        :: value_given  
    integer                     :: ierr
    character(len=*),intent(in) :: keywrd 

    if (exists(key2, keywrd)) then   
       if (has_value(key2, keywrd)) then    
          buffer=' '    
          buffer = get_value(key2,keywrd)      
          read(buffer,*,iostat=ierr)value_given
          if(ierr.ne.0) call fatal_keyword(keywrd,buffer)    
          if (.not.quiet)write(stdout,*)'Value for ',keywrd,' specified explicitly is:',value_given    
          if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
       else
          write(stderr,*)keywrd,' keyword given but it has no value!'
          write(stderr,*)'Program exiting.'
          write(stderr,*)'----------------------------------------------------------------------------'
          stop
       end if
    else
       write(stderr,*)'Value for ',keywrd,' not specified explicitly.'
       write(stderr,*)'Program exiting.'
       write(stderr,*)'----------------------------------------------------------------------------'
       stop
    end if

    return

  end subroutine  get_required_integer

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  get_optional_string(opt_string,keyword)

    character(len=100)          :: opt_string
    character(len=*),intent(in) :: keyword

    if (exists(key2, keyword)) then
       if (has_value(key2, keyword)) then
          opt_string=' '
          opt_string = get_value(key2,keyword)
          if(.not.quiet)write(stdout,*)trim(opt_string)
          if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
       else
          if(.not.quiet)write(stdout,*)'Keyword ',trim(keyword),' is present but has no value'
          if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
       end if
    else
       if(.not.quiet)write(stdout,*)'Keyword ',trim(keyword),' not present.'
       if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if

    return

  end subroutine get_optional_string

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_optional_int_and_real(value_given,real_val,keyword)

    character(len=500)          :: buffer
    integer, intent(out)        :: value_given
    real, intent(out)           :: real_val
    integer                     :: ierr
    character(len=*),intent(in) :: keyword

    real_val = -9999.0
    value_given = -9999
    if (exists(key2, keyword)) then
       if (has_value(key2, keyword)) then
          ! Grab 1 integer and 1 real
          buffer=' '
          buffer = get_value(key2,keyword)
          read(buffer,*,iostat=ierr)value_given, real_val
          if(ierr.ne.0) then 
             call fatal_keyword(keyword,buffer)
          else
             if(value_given.lt.1)then
                if (.not. quiet)write(stdout,*)'A value for ',keyword,' is zero or negative.'
                if (.not. quiet)write(stdout,*)'Is this OK?'
                if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------' 
             end if
             if(real_val.lt.0.00001) call fatal_keyword(keyword,buffer)
             if (.not. quiet)  &  
                  write(stdout,'(" Values for ",a15," specified explicitly are: ",i5,f9.5)')trim(keyword),value_given, real_val
             if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
          end if
       else
          if (.not. quiet)write(stdout,*)keyword,' keyword given but it has no values!'
          if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
       end if
    else
       if (.not. quiet)write(stdout,*)'Values for ',keyword,' not specified explicitly.'
       if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if

    return

  end subroutine get_optional_int_and_real

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------
  !
  subroutine get_ncross(ncross_given)

    character(len=100)     :: keyword
    character(len=500)     :: buffer
    character(len=4)       :: zmat, text
    integer                :: i,ierr, n
    integer, intent(out)   :: ncross_given(:)
    integer, allocatable   :: zmat_num(:)

    allocate(zmat_num(total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('zmat_num(total_num_zm)                                 ')

    zmat = 'ZMAT'
    keyword = 'NUMCROSS'
    ncross_given = -9999
    if (exists(key2, keyword)) then
       n = num_value(key2, keyword)
       do i = 1,n
          buffer=' '
          buffer = get_value(key2,keyword,i)
          read(buffer,*,iostat=ierr)text,zmat_num(i),ncross_given(zmat_num(i))
          if(ierr.ne.0)then
             call fatal_keyword(keyword,buffer)
          else
             !              make sure they are acceptable numbers.
             if((zmat_num(i).lt.1).or.(nintern(zmat_num(i)).lt.0))call fatal_keyword(keyword,buffer)
             text = .ucase. text
             if(text.ne.zmat)call fatal_keyword(keyword,buffer)
             if (.not. quiet)write(stdout,'(" Number of cross terms for z-matrix ",i5)')zmat_num(i)
             if (.not. quiet)write(stdout,'(" specified explicitly is: ",i5)')ncross_given(zmat_num(i))
             if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
          end if
       end do ! i = 1,n

    end if

    if (minval(ncross_given).eq.-9999) then
       if (.not. quiet)write(stdout,*)'Value for ',trim(keyword),' not specified explicitly'
       if (.not. quiet)write(stdout,*)'for one or more z-matrices.'
       if (.not. quiet)write(stdout,*)'Will try to work it out from other inputs.'
       if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
    end if
    deallocate(zmat_num)

    return

  end subroutine get_ncross

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine  get_sprcon_size(sprcon,keyword, ntyp)

    character(len=500)           :: buffer
    character(len=*), intent(in) :: keyword
    integer                      :: ierr,n,i,ol,counter,j
    real, intent(inout)          :: sprcon(:)
    integer, intent(in)          :: ntyp
    real                         :: real_num
    integer, allocatable         :: fields(:)

    n = num_value(key2, keyword)
    counter = 0 
    do i = 1,n
       buffer=' '
       buffer = get_value(key2,keyword,i)
       call get_num_fields(buffer, ol)
       read(buffer,*,iostat=ierr) real_num
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       if(ol.eq.1) then
          sprcon = real_num
          counter = counter + 1
       end if
    end do ! i = 1,n
    if((.not.quiet).and.(counter.gt.1))    & 
         write(stdout,*)trim(keyword),' default set more than once.  Last instance will be used.'

    do i = 1,n
       buffer=' '
       buffer = get_value(key2,keyword,i)
       call get_num_fields(buffer, ol)
       if(ol.gt.1) then
          allocate(fields(ol-1),stat=ierr)
          if (ierr.ne.0) call allerr('fields(ol-1)                                           ')
          read(buffer,*,iostat=ierr) real_num,(fields(j),j=1,(ol-1))
          if(ierr.ne.0) call fatal_keyword(keyword,buffer)
          if(maxval(fields).gt.ntyp) call fatal_keyword(keyword,buffer)
          do j = 1,(ol-1)
             sprcon(fields(j)) = real_num
          end do
          deallocate(fields)
       end if
    end do ! i = 1,n

    return

  end subroutine  get_sprcon_size

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_inspr(insprcon, insprtype) 

    real, intent(out)              :: insprcon(:)
    integer, intent(out)           :: insprtype(:,:)
    character(len=100)             :: keyword
    character(len=500)             :: buffer
    character(len=10)              :: zmat, inte, text
    character(len=10), allocatable :: split_buffer(:)
    integer                        :: ierr,n,i,z,int_num, j, num_fields, jj

    keyword = 'INSPR'
    zmat = 'ZMAT'
    inte = 'INT'
    n = num_value(key2, keyword)
    do i = 1,n

       buffer=' '
       buffer = get_value(key2,keyword,i)
       call get_num_fields(buffer,num_fields)
       allocate(split_buffer(num_fields),stat=ierr)
       if (ierr.ne.0) call allerr('split_buffer(num_fields)                               ')
       split_buffer = ' '
       read(buffer,*,iostat=ierr)(split_buffer(j),j=1,num_fields)
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       read(split_buffer(1),*,iostat=ierr)insprcon(i)
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       read(split_buffer(2),*,iostat=ierr)text
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       text = .ucase.text
       if(trim(text).ne.trim(zmat)) call fatal_keyword(keyword,buffer)
       read(split_buffer(3),*,iostat=ierr)z 
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       if((z.lt.1).or.(z.gt.total_num_zm))call fatal_keyword(keyword,buffer)
       read(split_buffer(4),*,iostat=ierr) text
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       text = .ucase.text
       if(trim(text).ne.trim(inte)) call fatal_keyword(keyword,buffer)
       read(split_buffer(5),*,iostat=ierr)int_num
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       if((int_num.lt.1).or.(int_num.gt.nintern(z)))call fatal_keyword(keyword,buffer)
       insprtype(int_num,z) = i
       if (.not.quiet) write(stdout,'("     Internal d.f. #",i5," of z-matrix",i5)')int_num,z
       if (num_fields.gt.5) then
          j = 6
          do jj = 6, num_fields
             read(split_buffer(j),*,iostat=ierr) text
             if(ierr.ne.0) call fatal_keyword(keyword,buffer)
             if (.isnumber. text) then
                ! Must be another internal
                read(split_buffer(j),*,iostat=ierr) int_num
                insprtype(int_num,z) = i
                if (.not.quiet) write(stdout,'(" and internal d.f. #",i5," of z-matrix",i5)')int_num,z
                j = j + 1
                if(j.gt.num_fields)exit
             else
                ! It must be next ZMAT statement:
                text = .ucase.text
                if(trim(text).ne.trim(zmat)) call fatal_keyword(keyword,buffer)
                j = j + 1
                ! Sop must be zmat number next
                if(j.gt.num_fields)exit
                read(split_buffer(j),*,iostat=ierr)z
                j = j + 1
                if(ierr.ne.0) call fatal_keyword(keyword,buffer)
                if((z.lt.1).or.(z.gt.total_num_zm))call fatal_keyword(keyword,buffer)
                ! then 'int'
                if(j.gt.num_fields)exit
                read(split_buffer(j),*,iostat=ierr) text
                j = j + 1
                if(ierr.ne.0) call fatal_keyword(keyword,buffer)
                text = .ucase.text
                if(trim(text).ne.trim(inte)) call fatal_keyword(keyword,buffer)
                ! then int number
                if(j.gt.num_fields)exit
                read(split_buffer(j),*,iostat=ierr)int_num
                if(ierr.ne.0) call fatal_keyword(keyword,buffer)
                !                if((int_num.lt.1).or.(int_num.gt.nintern(z)))call fatal_keyword(keyword,buffer)
                insprtype(int_num,z) = i
                if (.not.quiet) write(stdout,'(" and internal d.f. #",i5," of z-matrix",i5)')int_num,z
                j = j + 1
                if(j.gt.num_fields)exit
             end if
          end do
       end if
       if(.not.quiet)write(stdout,'(" Has force constant of type",i5," and value",f8.3)')i,insprcon(i)

       deallocate(split_buffer)

    end do
    if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
    return

  end subroutine get_inspr

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine get_cross()
    !
    character(len=100)             :: keyword
    character(len=500)             :: buffer
    character(len=10)              :: zmat, inte, text, cross, text1,text2,text3
    integer                        :: n,i, num_fields,iz, zmat_num
    integer                        :: cross_num,int1,int2,z,j,jj
    character(len=10), allocatable :: split_buffer(:)

    !----------------------------------------------------------------------------------
    ! First get ncross from actual definitions in the input file and
    ! compare with ncross_given.

    ! this entails reading in CROSSDEF  lines (CROSSDEF ZMAT 1 CROSS 1 INT 2 3 means
    ! cross term 1 for ZMAT 1 is between internals 2 and 3.
    !----------------------------------------------------------------------------------

    zmat = 'ZMAT'
    inte = 'INT'
    cross = 'CROSS'
    keyword = 'CROSSDEF'
    allocate(ncross(total_num_zm),stat=ierr)   
    if (ierr.ne.0) call allerr('ncross(total_num_zm)                                ')
    ncross = 0
    n = num_value(key2,keyword)
    do i = 1,n
       buffer=' '
       buffer = get_value(key2,keyword,i)
       call get_num_fields(buffer,num_fields)
       if(num_fields.ne.7) call fatal_keyword(keyword,buffer)
       read(buffer,*,iostat=ierr)text,iz
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       text=.ucase.text
       if(trim(text).ne.trim(zmat)) call fatal_keyword(keyword,buffer)
       ncross(iz) = ncross(iz) +1    
    end do
    do i=1,total_num_zm
       if(ncross_given(i).ne.-9999) then
          if(ncross(i).ne.ncross_given(i)) then
             write(stderr,*)'----------------------------------------------------------------------------'
             write(stderr,*)'Number of cross terms specifed using CROSSDEF not same as number'
             write(stderr,*)'expected from NUMCROSS.  Program exiting.'
             write(stderr,*)'----------------------------------------------------------------------------'
             stop
          end if
       end if
    end do
    iz  = 0
    do i = 1,total_num_zm
       ncross(i) = ncross_given(i)
       iz = iz  + ncross(i)
    end do

    allocate(crosstype(maxval(ncross),total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('crosstype(maxval(ncross),total_num_zm)              ')
    allocate(cross_terms(total_num_zm,maxval(ncross),2),stat=ierr)
    if (ierr.ne.0) call allerr('cross_terms(total_num_zm,(maxval(ncross),2)         ')
    allocate(cross_spr(iz),stat=ierr)
    if (ierr.ne.0) call allerr('cross_spr(iz)                                       ')

    do i = 1,n
       buffer=' '
       buffer = get_value(key2,keyword,i)
       read(buffer,*,iostat=ierr)text1,zmat_num,text2,cross_num,text3,int1,int2
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       if((zmat_num.lt.1).or.(zmat_num.gt.total_num_zm))call fatal_keyword(keyword,buffer)
       if((min(int1,int2).lt.1).or.(max(int1,int2).gt.nintern(zmat_num)))call fatal_keyword(keyword,buffer)
       if((cross_num.lt.1).or.(cross_num.gt.ncross(zmat_num)))call fatal_keyword(keyword,buffer)
       text1 = .ucase. text1
       text2 = .ucase. text2
       text3 = .ucase. text3
       if(trim(text1).ne.trim(zmat))call fatal_keyword(keyword,buffer)
       if(trim(text2).ne.trim(cross))call fatal_keyword(keyword,buffer)
       if(trim(text3).ne.trim(inte))call fatal_keyword(keyword,buffer)

       cross_terms(zmat_num,cross_num,1) = int1
       cross_terms(zmat_num,cross_num,2) = int2

       if(.not.quiet)write(stdout,*)'Cross term',cross_num,' of z-matrix ',zmat_num,' is between internal d.f.',int1,' and ',int2
    end do

    !----------------------------------------------------------------------------------
    ! So we have cross_terms defined.  Now we need to fill crosstype and cross_spr
    ! CROSSSPR 123.0 ZMAT 1 CROSS 1 2 3 ZMAT 2 CROSS 3 4 5 etc
    !----------------------------------------------------------------------------------

    keyword = 'CROSSSPR'
    n = num_value(key2,keyword)
    do i = 1,n
       buffer=' '
       buffer = get_value(key2,keyword,i)
       call get_num_fields(buffer,num_fields)
       allocate(split_buffer(num_fields),stat=ierr)
       if (ierr.ne.0) call allerr('split_buffer(num_fields)                               ')
       split_buffer = ' '
       read(buffer,*,iostat=ierr)(split_buffer(j),j=1,num_fields)
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       read(split_buffer(1),*,iostat=ierr)cross_spr(i)
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       read(split_buffer(2),*,iostat=ierr)text
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       text = .ucase.text
       if(trim(text).ne.trim(zmat)) call fatal_keyword(keyword,buffer)
       read(split_buffer(3),*,iostat=ierr)z 
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       if((z.lt.1).or.(z.gt.total_num_zm))call fatal_keyword(keyword,buffer)
       read(split_buffer(4),*,iostat=ierr) text
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       text = .ucase.text
       if(trim(text).ne.trim(cross)) call fatal_keyword(keyword,buffer)
       read(split_buffer(5),*,iostat=ierr)cross_num
       if(ierr.ne.0) call fatal_keyword(keyword,buffer)
       if((cross_num.lt.1).or.(cross_num.gt.ncross(z)))call fatal_keyword(keyword,buffer)
       crosstype(cross_num,z) = i
       if (.not.quiet)write(stdout,*)'     Cross term    #',cross_num, 'of z-matrix',z
       if (num_fields.gt.5) then
          j = 6
          do jj = 6, num_fields
             read(split_buffer(j),*,iostat=ierr) text
             if(ierr.ne.0) call fatal_keyword(keyword,buffer)
             if (.isnumber. text) then
                read(split_buffer(j),*,iostat=ierr) cross_num
                crosstype(cross_num,z) = i
                if (.not.quiet) write(stdout,'(" and cross term    #",i5," of z-matrix",i5)')cross_num,z
                j = j + 1
                if(j.gt.num_fields)exit
             else
                ! It must be next ZMAT statement:
                text = .ucase.text
                if(trim(text).ne.trim(zmat)) call fatal_keyword(keyword,buffer)
                j = j + 1
                ! Sop must be zmat number next
                if(j.gt.num_fields)exit
                read(split_buffer(j),*,iostat=ierr)z
                j = j + 1
                if(ierr.ne.0) call fatal_keyword(keyword,buffer)
                if((z.lt.1).or.(z.gt.total_num_zm))call fatal_keyword(keyword,buffer)
                ! then 'cross'
                if(j.gt.num_fields)exit
                read(split_buffer(j),*,iostat=ierr) text
                j = j + 1
                if(ierr.ne.0) call fatal_keyword(keyword,buffer)
                text = .ucase.text
                if(trim(text).ne.trim(cross)) call fatal_keyword(keyword,buffer)
                ! then int number
                if(j.gt.num_fields)exit
                read(split_buffer(j),*,iostat=ierr)cross_num
                if(ierr.ne.0) call fatal_keyword(keyword,buffer)
                crosstype(cross_num,z) = i
                if (.not.quiet) write(stdout,'(" and cross term    #",i5," of z-matrix",i5)')cross_num,z
                j = j + 1
                if(j.gt.num_fields)exit
             end if
          end do
       end if
       if(.not.quiet)write(stdout,'(" Has force constant of type",i5," and value",f8.3)')i,cross_spr(i)

       deallocate(split_buffer)
    end do
    if((n.gt.0).and.(.not. quiet))write(stdout,*)'----------------------------------------------------------------------------'

    return

  end subroutine get_cross

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine longusage()

    write(stdout,*)'|---------------------------------------------------------------------|'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| The main input file is a series of keywords and values.  Some       |'
    write(stdout,*)'| are mandatory, some are not.  The order in the keyword file does    |'
    write(stdout,*)'| not matter.  The details are outlined below.                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| It is best to begin all filenames called within the keyword file    |'
    write(stdout,*)'| with alphaetical characters, not with numbers or other characters.  |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| ZMC uses a z-matrix to describe a molecule.  If carefully           |'
    write(stdout,*)'| constructed, this allows segmented motion of the molecules in a     |'
    write(stdout,*)'| simple way.  The convention for naming these files is to give       |'
    write(stdout,*)'| them the extension ".zmat".  While a z-matrix defines the           |'
    write(stdout,*)'| molecular geometry, the molecule must be oriented and positioned    |'
    write(stdout,*)'| within the unit cell.  This is done using a 3-vector (x, y, z)      |'
    write(stdout,*)'| to specify the origin of the molecule (the position of the first    |'
    write(stdout,*)'| atom) and a quaternion (a normalised 4-vector (q1,q2,q3,q4)) to     |'
    write(stdout,*)'| give the orientation.  This information is in files with            |'
    write(stdout,*)'| extension ".qxyz". If one uses a third vector to hold the values    |'
    write(stdout,*)'| of the internal degrees of freedom for the molecule (this will      |'
    write(stdout,*)'| be an n-vector if there are n internal d.f.) then the molecule      |'
    write(stdout,*)'| is completely specified by the 3-vector, the 4-vector and the n-    |'
    write(stdout,*)'| vector.  Given that a substantial molecule may have 50 atoms in     |'
    write(stdout,*)'| it, requiring 150 coordinates, this is a great economy of           |'
    write(stdout,*)'| variables.                                                          |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| A molecule is said to occupy a location rather than a site          |'
    write(stdout,*)'| simply because site is a word commonly attached to atomic           |'
    write(stdout,*)'| position, and there is a desire to avoid confusion by using a       |'
    write(stdout,*)'| word (hopefully) without connotations.                              |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| A single Monte Carlo (MC) step consists of choosing a molecule      |'
    write(stdout,*)'| at random, calculating its energy by summing over any               |'
    write(stdout,*)'| interactions pertaining to it (whether internal to the molecule     |'
    write(stdout,*)'| or between the molecule and its environment), then randomly         |'
    write(stdout,*)'| altering the molecules configuration (putting random shifts on      |'
    write(stdout,*)'| the three vectors noted above), calculating the new energy and      |'
    write(stdout,*)'| then accepting or rejecting the new configuration based on a MC     |'
    write(stdout,*)'| algorithm.                                                          |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| Hence the keywords allow the user to specify what the degrees of    |'
    write(stdout,*)'| freedom are, what the force constants (spring constants) acting     |'
    write(stdout,*)'| in the system are, and how wide the distribution of random          |'
    write(stdout,*)'| shifts can be (referred to as various kinds of widths).             |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| It is my intention that a fuller manual for ZMC will become         |'
    write(stdout,*)'| available with time, and that it will be distributed with           |'
    write(stdout,*)'| example input files and simulations that work.  Also, some          |'
    write(stdout,*)'| documentation exists in the reviewed literature, as noted in        |'
    write(stdout,*)'| the messages on running the program.                                |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| --------------------------------------------------------------      |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| CRYSTAL      Specifies the size of the model crystal in unit        |'
    write(stdout,*)'|              cells. Must be followed by three integers.             |'
    write(stdout,*)'|              e.g. CRYSTAL 32 16 5                                   |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| HEADER       Specifies (if given) a header (up to 100               |'
    write(stdout,*)'|              characters long) for the keyword file and which        |'
    write(stdout,*)'|              is then written into the output file.                  |'
    write(stdout,*)'|              e.g. HEADER This is the header.                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| TEMPERATURE  Specifies the temperature of the Monte Carlo           |'
    write(stdout,*)'|              simulation. This affects how likely it is that a       |'
    write(stdout,*)'|              MC move which increases the system energy will be      |'
    write(stdout,*)'|              accepted.  Must be followed by one real.               |'
    write(stdout,*)'|              e.g. TEMPERATURE 1.0                                   |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| MCCYCLES     Specifies number of Monte Carlo cycles to do. One      |'
    write(stdout,*)'|              cycle consists of a number of MC steps sufficient      |'
    write(stdout,*)'|              to visit every location in the crystal once.           |'
    write(stdout,*)'|              Must be followed by one integer.                       |'
    write(stdout,*)'|              e.g. MCCYCLES 500                                      |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| INCUPDATE    Specifies how many MC cycles between adjusting         |'
    write(stdout,*)'|              the sizes of the increments (widths) that govern       |'
    write(stdout,*)'|              the shifts put on the variables in the                 |'
    write(stdout,*)'|              simulation.  The initial randomness is governed        |'
    write(stdout,*)'|              by XYZINITW, QINITW and ININITW, while the sizes       |'
    write(stdout,*)'|              of the random shifts applied to the variables          |'
    write(stdout,*)'|              during the MC are governed by XYZWIDTH, QWIDTH         |'
    write(stdout,*)'|              and INWIDTH.  These latter parameters can be           |'
    write(stdout,*)'|              dynamically adjusted to give a rejection ratio of      |'
    write(stdout,*)'|              about 50% (half moves rejected) which is a rule        |'
    write(stdout,*)'|              of thumb often applied to MC simulation.               |'
    write(stdout,*)'|              Further, the subroutine adjusts the shifts on the      |'
    write(stdout,*)'|              individual variables to ensure that all are            |'
    write(stdout,*)'|              having approximately equal influence on the            |'
    write(stdout,*)'|              acceptance/rejection ratio. If it turns out that       |'
    write(stdout,*)'|              one variable (for example the x position of the        |'
    write(stdout,*)'|              origin) has a very big width, this indicates that      |'
    write(stdout,*)'|              the simulation is only weakly sensitive to shifts      |'
    write(stdout,*)'|              in this direction.  This may reflect reality, or       |'
    write(stdout,*)'|              may imply more or stronger constraints                 |'
    write(stdout,*)'|              (interactions) are needed.  Must be followed by a      |'
    write(stdout,*)'|              single integer.  e.g. INCUPDATE 1                      |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| BADJUST      Specifies how many MC cycles between adjusting         |'
    write(stdout,*)'|              the strengths of the interactions in the system        |'
    write(stdout,*)'|              to achieve a specified average B-factor across         |'
    write(stdout,*)'|              ALL atoms.  Interactions are scaled globally, so       |'
    write(stdout,*)'|              this is analogous to changing the simulation           |'
    write(stdout,*)'|              temperature.  It is helpful because it gives a         |'
    write(stdout,*)'|              good indication of the scale of the interactions.      |'
    write(stdout,*)'|              Since it scales them all at once their relative        |'
    write(stdout,*)'|              values do not change.  The B-factor specified          |'
    write(stdout,*)'|              (B = 8pi^2U) is an average across all atoms in         |'
    write(stdout,*)'|              the z-matrix, so is a very crude number, and is        |'
    write(stdout,*)'|              calculated as isotropic, which is also a crude         |'
    write(stdout,*)'|              approximation.  Must be followed by one integer        |'
    write(stdout,*)'|              and one real, where the real is the desired            |'
    write(stdout,*)'|              B-factor.  e.g. BADJUST 1 2.5                          |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| NUMSPRCON    Specifies the number of different spring               |'
    write(stdout,*)'|              constants (interaction constants) acting on            |'
    write(stdout,*)'|              contact vectors in the model.  Will be deduced         |'
    write(stdout,*)'|              from other inputs if not given.  Must be followed      |'
    write(stdout,*)'|              by one integer.  e.g. NUMSPRCON 23                     |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| NUMINSPRCON  Specifies the number of different spring               |'
    write(stdout,*)'|              constants (interaction constants) acting on            |'
    write(stdout,*)'|              internal degrees of freedom in the model.              |'
    write(stdout,*)'|              Will be deduced from other inputs if not given.        |'
    write(stdout,*)'|              Must be followed by one integer.                       |'
    write(stdout,*)'|              e.g. NUMINSPRCON 23                                    |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| NUMLOCS      Specifies the number of locations in the unit          |'
    write(stdout,*)'|              cell.  Location is used instead of site because        |'
    write(stdout,*)'|              site tends to be used to refer to atoms; a             |'
    write(stdout,*)'|              location is a position in the unit cell which is       |'
    write(stdout,*)'|              occupied by a molecule.  The occupation of a           |'
    write(stdout,*)'|              given location may vary from cell to cell -- the       |'
    write(stdout,*)'|              type or orientation of z-matrix on the location        |'
    write(stdout,*)'|              may vary.  Will be deduced from other inputs if        |'
    write(stdout,*)'|              not given.  Must be  followed by one integer.          |'
    write(stdout,*)'|              e.g. NUMLOCS 23                                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| NUMZMATS     Specifies the number of different z-matrices           |'
    write(stdout,*)'|              present in the simulation.  Will be deduced from       |'
    write(stdout,*)'|              other inputs if not given.  Must be followed by        |'
    write(stdout,*)'|              one integer.  e.g. NUMZMATS 2                          |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| NUMINTERNAL  Specifies the number of internal degrees of            |'
    write(stdout,*)'|              freedom a given z-matrix is given.  Must be            |'
    write(stdout,*)'|              followed by the sub-keyword ZMAT then two              |'
    write(stdout,*)'|              integers, the ZMAT number (as given in the             |'
    write(stdout,*)'|              ZMATRIX lines of the keyword file) and how many        |'
    write(stdout,*)'|              internal degrees of freedom that z-matrix has          |'
    write(stdout,*)'|              been given.  Will be deduced from other inputs if      |'
    write(stdout,*)'|              not given.  To specify that z-matrix 2 has 3           |'
    write(stdout,*)'|              internal degrees of freedom, write:                    |'
    write(stdout,*)'|              NUMINTERNAL ZMAT 2 3                                   |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| NUMCROSS     Specifies the number of cross-terms a given            |'
    write(stdout,*)'|              z-matrix is given. A cross-term is an interaction      |'
    write(stdout,*)'|              between internal degrees of freedom.  Must be          |'
    write(stdout,*)'|              followed by the sub-keyword ZMAT then two              |'
    write(stdout,*)'|              integers, the ZMAT number (as given in the             |'
    write(stdout,*)'|              ZMATFILE lines of the keyword file) and how many       |'
    write(stdout,*)'|              cross terms that z-matrix has been given.  Will        |'
    write(stdout,*)'|              be deduced from other inputs if not given.  To         |'
    write(stdout,*)'|              specify that z-matrix 2 has 3 internal cross           |'
    write(stdout,*)'|              terms, write: NUMCROSS ZMAT 2 3                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| OCCFILE      Specifies the name of the file the occupancies         |'
    write(stdout,*)'|              are to be read from.  If not given it is assumed       |'
    write(stdout,*)'|              that it is not needed (a model in which there is       |'
    write(stdout,*)'|              no occupancy disorder, generally).                     |'
    write(stdout,*)'|              e.g. OCCFILE occupancies.txt                           |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| ZMATFILE     Specifies the file containing a z-matrix.  Must        |'
    write(stdout,*)'|              be followed by an integer, which is the z-matrix       |'
    write(stdout,*)'|              number, and a filename.                                |'
    write(stdout,*)'|              e.g. ZMATFILE 1 paraterphenyl.zmat                     |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| QXYZFILE     Specifies the file containing the variables            |'
    write(stdout,*)'|              describing the average origin of the z-matrix          |'
    write(stdout,*)'|              (xyz) and its orientation (quaternion, q).  Must       |'
    write(stdout,*)'|              be followed by an integer, which is the z-matrix       |'
    write(stdout,*)'|              number, and a filename.                                |'
    write(stdout,*)'|              e.g.  QXYZFILE 1 paraterphenyl.qxyz                    |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| CELL         Specifies unit cell parameters (angles in              |'
    write(stdout,*)'|              degrees). Must be followed by 6 reals.                 |'
    write(stdout,*)'|              e.g.  CELL 12.34 5.126 8.902 90.0  109.35 90.0         |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| CONTACTFILE  Specifies the file containing the contact              |'
    write(stdout,*)'|              vectors.  e.g. CONTACTFILE Contacts.txt                |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| SPRCON       Specifies value of spring constants (interaction       |'
    write(stdout,*)'|              constants).  If given with a single real               |'
    write(stdout,*)'|              (e.g. SPRCON 12.0) sets a default value that           |'
    write(stdout,*)'|              applies to all springs.  Subsequent invocations,       |'
    write(stdout,*)'|              with a real followed by integers,                      |'
    write(stdout,*)'|              (e.g. SPRCON 50.2 1 2 3 4 5 6 7) sets the spring       |'
    write(stdout,*)'|              constants for those contact vectors specified.         |'
    write(stdout,*)'|              The two invocations here would set all springs to      |'
    write(stdout,*)'|              be 12.0, then set the spring constants on vectors      |'
    write(stdout,*)'|              1 to 7 to be 50.2, leaving the others, if any, at      |'
    write(stdout,*)'|              12.0.                                                  |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| SIZE         Specifies value of size-effects on contact             |'
    write(stdout,*)'|              vectors.  Invocation rules are the same as for         |'
    write(stdout,*)'|              SPRCON.                                                |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| INSPR        Specifies spring constants on internal degrees of      |'
    write(stdout,*)'|              freedom.  An example would be:                         |'
    write(stdout,*)'|              INSPR  50.0  ZMAT 1 INT 1 ZMAT 2 INT 2 3               |'
    write(stdout,*)'|              This would set the spring constant (interaction        |'
    write(stdout,*)'|              constant) for internal degree of freedom 1 on          |'
    write(stdout,*)'|              z-matrix 1 to be 50.0 and internal spring              |'
    write(stdout,*)'|              constants 2 and 3 of z-matrix 2 would also be          |'
    write(stdout,*)'|              50.0.  These must all be set explicitly; there is      |'
    write(stdout,*)'|              no defaulting invocation as for SPRCON.                |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| INTERNAL     Specifies the nature of an internal degree of          |'
    write(stdout,*)'|              freedom.  To specify internal degree of freedom        |'
    write(stdout,*)'|              number 1 for z-matrix number 2 to be the               |'
    write(stdout,*)'|              dihedral angle of atom 17, type:                       |'
    write(stdout,*)'|              INTERNAL ZMAT 2 INT 1 dihedral 17                      |'
    write(stdout,*)'|              The types of internal d.f. are bond lengths (use       |'
    write(stdout,*)'|              "length" instead of "dihedral"), bond angle            |'
    write(stdout,*)'|              (sub-keyword "angle") or dihedral angle as shown.      |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| CROSSSPR     Specifies spring constants on cross terms.  An         |'
    write(stdout,*)'|              example would be:                                      |'
    write(stdout,*)'|              CROSSSPR  50.0  ZMAT 1 CROSS 1 ZMAT 2 CROSS 2 3        |'
    write(stdout,*)'|              Implemented much like INSPR                            |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| CROSSDEF     Specifies cross-terms.  For example:                   |'
    write(stdout,*)'|              CROSSDEF ZMAT 1 CROSS 1 INT 1 2                        |'
    write(stdout,*)'|              would define the first cross-term of z-matrix 1        |'
    write(stdout,*)'|              is an interaction between internal d.f. 1 and 2.       |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| XYZWIDTH     Sets the sizes of the random increments on x, y        |'
    write(stdout,*)'|              and z coordinates of z-matrices.  Can be invoked       |'
    write(stdout,*)'|              in several ways.  "XYZWIDTH 0.2" would set widths      |'
    write(stdout,*)'|              to be 0.2 for all x,y and z on all z-matrices.         |'
    write(stdout,*)'|              "XYZWIDTH X 0.1" would set the widths on the           |'
    write(stdout,*)'|              x-shifts on all z-matrices to be 0.1.                  |'
    write(stdout,*)'|              "XYZWIDTH 1 0.15" would set widths on x, y and z       |'
    write(stdout,*)'|              for all of x, y and z on z-matrix 1 to be 0.15         |'
    write(stdout,*)'|              and (lastly) "XYZWIDTH 1 Y 0.4" would set the          |'
    write(stdout,*)'|              width of the y coordinate on z-matrix 1  to            |'
    write(stdout,*)'|              be 0.4.  "XTZWIDTH 1 X Y 0.15" would not work;         |'
    write(stdout,*)'|              use separate invocations for X and Y in this case.     |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| XYZINITW     Sets the initial randomness on the variables;          |'
    write(stdout,*)'|              invocation as for XYZWIDTH.                            |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| QWIDTH       Sets the sizes of the random increments on             |'
    write(stdout,*)'|              quaternion components.  As for XYZWIDTH but refer      |'
    write(stdout,*)'|              to Q1, Q2, Q3 and Q4 rather than X, Y and Z.           |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| QINITW       See QWIDTH and XYZINITW                                |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| INWIDTH      Sets the sizes of the random increments on             |'
    write(stdout,*)'|              internal degrees of freedom.  As for XYZWIDTH but      |'
    write(stdout,*)'|              refer to I1, I2, ..., IN rather than X, Y and Z.       |'
    write(stdout,*)'|              NOTE that it works best with one line per d.f..        |'
    write(stdout,*)'|              For example "INWIDTH 1 I2 15.0" would set the width    |'
    write(stdout,*)'|              for internal 2 on ZMAT 1 to be 15.0 degrees.           |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| ININITW      See previous couple of entries.                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| VMIN         When creating some contact vectors, this is the        |'
    write(stdout,*)'|              minimum length (Angstrom). e.g. VMIN 1.9               |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| VMAX         When creating some contact vectors, this is the        |'
    write(stdout,*)'|              maximum length.                                        |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| NEWCONTACTS  Specifies the new contact vector file name.            |'
    write(stdout,*)'|              e.g. NEWCONTACTS newfilename.out                       |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'| CONTATOMS    Specified which atoms to look at when generating       |'
    write(stdout,*)'|              contact vectors.  "CONTATOMS ZMAT 1 17 25" would       |'
    write(stdout,*)'|              look for contacts between atoms 17 and 25 on           |'
    write(stdout,*)'|              ZMAT 1.  To look for contacts involving more than      |'
    write(stdout,*)'|              one z-matrix type, use multiple instances of           |'
    write(stdout,*)'|              CONTATOMS.                                             |'
    write(stdout,*)'|                                                                     |'
    write(stdout,*)'|---------------------------------------------------------------------|'
    stop

  end subroutine longusage

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine write_readat(filename)

    !----------------------------------------------------------------------------------------------
    ! This suggests that what we need to do is (1) Work out number of possible sites in a cell.    
    ! (2) Work out the maximum number of sites that can be occupied in a cell.                   
    ! (3) Write these things and simsize() to the output file.                                  
    ! (4) Loop over the POSSIBLE sites in the cell (by looping over locations then possible     
    !     z-matrix then possible orientations) and increment a counter n for each non-dummy atom.
    !     Before doing this, set correspond() to zeroes; then, when on a site which IS occupied  
    !     increment a second counter, m, and set correspond(i,j,k,n) = m.  Then set            
    !     x(i,j,k,m), y, z and label() as well, so that they can be accessed by using         
    !     x(i,j,k,correspond(i,j,k,n)) for example.  Then, when doing readat() it is simply  
    !     that if correspond(i,j,k,isite)=0, it is unoccupied.  If .ne.0, then get label and
    !     compare to type passed in.  If agreement, set the flag, set x,y,z and get out. Easy!
    !----------------------------------------------------------------------------------------------

    character(len=107)                   :: filename
    integer :: possible_sites,il,iz,oz,im,iat,occupied_sites,max_occ,num_nondummy,dunit
    integer :: z_there, m_there,ia,ib,ic,in2,ind
    type (zmatrix_object)                :: z1
    character(len=6)                     :: lable
    real,dimension(3)                    :: cvec,fvec

    allocate(how_many_types(num_loc),stat=ierr)
    if (ierr.ne.0) call allerr('how_many_types(num_loc)                             ')
    allocate(which_types(num_loc,total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('which_types(num_loc,total_num_zm)                   ')
    allocate(instances_types(num_loc,total_num_zm),stat=ierr)
    if (ierr.ne.0) call allerr('instances_types(num_loc,total_num_zm)               ')
    allocate(which_mol(num_loc,total_num_zm,maxval(mocc)),stat=ierr)
    if (ierr.ne.0) call allerr('which_mol(num_loc,total_num_zm,maxval(mocc))         ')
    call get_what_goes_where(how_many_types,which_types,instances_types,which_mol,total_num_zm,num_loc,quiet,key2)

    !----------------------------------------------------------------------------------------------
    ! (1) Work out number of possible sites in a cell.
    !----------------------------------------------------------------------------------------------

    possible_sites = 0
    do il = 1,num_loc
       do iz = 1,how_many_types(il)
          oz = which_types(il,iz)
          call copy(z1,zmat(oz))
          do im = 1,instances_types(il,oz)
             ! now, we have to loop over the atoms, because we have to factor out dummies
             do iat = 1,num_at(oz)
                lable = ''
                lable=labels(z1,iat)
                if((lable(1:1)).ne.'x') possible_sites = possible_sites + 1
             end do
          end do
       end do
    end do

    !----------------------------------------------------------------------------------------------
    ! (2) Work out the maximum number of sites that can be occupied in a cell.                    !
    ! For each location, what is the biggest z-matrix that can go on it?
    ! Should we worry about dummy atoms at this stage?  It will mean our arrays
    ! come out a little bigger than needed...do we need a subroutine to work out
    ! how many dummy atoms?
    !----------------------------------------------------------------------------------------------

    occupied_sites = 0
    do il = 1,num_loc
       max_occ = 0
       do iz = 1,how_many_types(il)
          oz = which_types(il,iz)
          call copy(z1,zmat(oz))
          num_nondummy = 0
          do iat = 1,num_at(oz)
             lable=labels(z1,iat)
             if((lable(1:1)).ne.'x') num_nondummy = num_nondummy + 1
          end do
          if(max_occ.lt.num_nondummy) max_occ = num_nondummy 
       end do
       occupied_sites = occupied_sites + max_occ 
    end do

    write(stdout,*)'Possible atom sites:', possible_sites
    write(stdout,*)'Of which at most',occupied_sites,' can be occupied.'

    !----------------------------------------------------------------------------------------------
    ! Open the output file and write these to it.
    !----------------------------------------------------------------------------------------------

    dunit = open(filename)
    write(dunit,*)simsize,possible_sites,occupied_sites

    !----------------------------------------------------------------------------------------------
    ! Now the guts of it.
    ! First write out lookup table array; readat() can interrogate that to find out
    ! what it needs to know,and then read in the actual atoms.
    !----------------------------------------------------------------------------------------------

    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             write(dunit,*)ia,ib,ic
             occupied_sites = 0
             possible_sites = 0
             do il = 1,num_loc
                z_there = zocc(ia,ib,ic,il)
                m_there = mocc(ia,ib,ic,il)
                do iz = 1,how_many_types(il)
                   oz = which_types(il,iz)
                   call copy(z1,zmat(oz))
                   do im = 1,instances_types(il,oz)
                      om = which_mol(il,oz,im)
                      do iat = 1,num_at(oz)
                         lable = ''
                         lable=labels(z1,iat)
                         if((lable(1:1)).ne.'x') then 
                            possible_sites = possible_sites + 1
                            write(dunit,'(i7)',advance='no')possible_sites
                            if((om.eq.m_there).and.(oz.eq.z_there)) then
                               ! then we actually have an occupied site!
                               occupied_sites = occupied_sites + 1
                               write(dunit,'(i7)')occupied_sites
                            else
                               write(dunit,*) 0  
                            end if
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    ! Now write out the actual atoms
    do ia = 1,simsize(1)
       do ib = 1,simsize(2)
          do ic = 1,simsize(3)
             write(dunit,*)ia,ib,ic
             do il = 1,num_loc
                z_there = zocc(ia,ib,ic,il)
                call copy(z1,zmat(z_there))
                do iat = 1,num_at(z_there)
                   lable = ''
                   lable=labels(z1,iat)
                   if((lable(1:1)).ne.'x') then 
                      do i = 1,3
                         cvec(i)=carts(i,ia,ib,ic,il,iat)
                      end do
                      fvec(:)=as_fractional(xtal,cvec(:))
                      in2 = max((scan(lable,'0123456789')-1),1)
                      if (in2.lt.1) then
                         write(stdout,*)'All atoms in zmats need to be numbered!!'
                         write(stdout,*)'E.g. Nb1 or O23'
                         write(stdout,*)'Not halting, but check your output files!'
                      end if
                      ind = min(2,in2)
                      write(dunit,'(3f10.5,a4)')(fvec(i),i=1,3),lable(1:ind)
                   end if
                end do
             end do
          end do
       end do
    end do

    close(dunit)

    return

  end subroutine write_readat

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  subroutine write_discus(xyzfile,xtal)

    !----------------------------------------------------------------
    ! Write out xyz coordinates of all atoms to some vast ASCII file. 
    !----------------------------------------------------------------

    character(len=107),intent(in)       :: xyzfile
    character(len=107)                  :: pdffile,macfile

    integer                             :: ia,ib,ic,il,dunit,i,ind

    real,dimension(3)                   :: cvec,fvec,cellaxes,cellangles

    type (crystallography_object)       :: xtal
    integer                             :: iat,in2
    integer                             :: z_there
    type (zmatrix_object)               :: z1
    character(len=6)                    :: lable
    character(len=4)                    :: text

    if (.not.  quiet)write(stdout,*)'----------------------------------------------------------------------------'
    if (.not. quiet)write(stdout,*)'Writing DISCUS output to ',trim(xyzfile)
    dunit = open(xyzfile)

    write(dunit,'(a21)')'title Output from ZMC'
    write(dunit,'(a8)')'spcgr P1'
    cellaxes(:)  = axes(xtal)
    cellangles(:)= angles(xtal)
    write(dunit,'(a5,6f10.4)')'cell ',cellaxes,cellangles
    write(dunit,'(a22)')'PUT ncells LINE HERE!!'
    write(dunit,'(a5)')'atoms'


          do ic = 1,simsize(3)
       do ib = 1,simsize(2)
    do ia = 1,simsize(1)
             do il = 1,num_loc
                z_there = zocc(ia,ib,ic,il)
                call copy(z1,zmat(z_there))
                do iat = 1,num_at(z_there)
                   lable = ''
                   lable=labels(z1,iat)
                   if((lable(1:1)).ne.'x') then
                      do i = 1,3
                         cvec(i)=carts(i,ia,ib,ic,il,iat)
                      end do
                      fvec(:)=as_fractional(xtal,cvec(:))
                      fvec(1) = fvec(1) + ia - 1
                      fvec(2) = fvec(2) + ib  - 1
                      fvec(3) = fvec(3) + ic  - 1
                      in2 = max((scan(lable,'0123456789')-1),1)
                      if (in2.lt.1) then
                         write(stdout,*)in2,'All atoms in zmats need to be numbered!!'
                         write(stdout,*)'E.g. Nb1 or O23'
                         write(stdout,*)'Not halting, but check your output files!'
                      end if
                      ind = min(2,in2)
                      if(ind.eq.1) then
                         write(text,'(a1,a3)')lable(1:1),'   '
                      else
                         write(text,'(a2,a2)')lable(1:2),'  '
                      end if
                      write(dunit,'(a4,f10.5,a2,f10.5,a2,f10.5,a2,f5.2,a2,f5.2)')   &
                              text,fvec(1),', ',fvec(2),', ',fvec(3),', ',1.0,', ',1.0
                   end if
                end do
             end do
          end do
       end do
    end do
    close(dunit)


    call extension(xyzfile,'pdf',pdffile)
    call extension(xyzfile,'mac',macfile)
    if (.not. quiet)write(stdout,*)'Writing DISCUS pdf macro to ',trim(macfile)
    dunit = open(macfile)
    write(dunit,'(a63)')'# This comes from a very simple routine in ZMC                 '
    write(dunit,'(a63)')'# to output a DISCUS-compatible file. Also, outputs a          '
    write(dunit,'(a63)')'# DISCUS macro file that ought to give a first pass at an input'
    write(dunit,'(a63)')'# file for pdf and powder calculations. NOT THOROUGHLY TESTED! '
    write(dunit,'(a1)')'#'
    write(dunit,'(a63)')'# See http://discus.sourceforge.net/                           '
    write(dunit,'(a1)')'#'
    write(dunit,'(a63)')'# Note the spacegroup is always given as P1 by ZMC             ' 
    write(dunit,'(a1)')'#'
    write(dunit,'(a63)')'# Note the Biso is always given as unity and is then           '
    write(dunit,'(a63)')'# disabled in the macro.                                       '
    write(dunit,'(a1)')'#'
    write(dunit,'(a63)')'# Note also that DISCUS may not be able to read in all the     '
    write(dunit,'(a63)')'# atoms in the simulation, so you may have to reduce the       '
    write(dunit,'(a63)')'# simulation size or recompile DISCUS with bigger limits.      '
    write(dunit,'(a63)')'# Just truncating the .discus file may not work.               '
    write(dunit,'(a1)')'#'
    write(dunit,'(a63)')'# To get a quick look at the outputs, plot in a spreadsheet   '
    write(dunit,'(a63)')'# or just run gnuplot and type plot "filename.pdf" (in quotes).'
    write(dunit,'(a1)')'#'
    write(dunit,'(a4)')'read'
    write(dunit,'(a6)',advance='no')'stru  '
    write(dunit,*)trim(xyzfile)
    write(dunit,'(a1)')'#'
    write(dunit,'(a3)')'pdf'
    write(dunit,'(a1)')'#'
    write(dunit,'(a13)')'set therm,off'
    write(dunit,'(a1)')'#'
    write(dunit,'(a20)')'set rang, 20.0, 0.02'
    write(dunit,'(a14)')'set qmax, 26.0'
    write(dunit,'(a13)')'set qsig, 0.0'
    write(dunit,'(a13)')'set rad, xray'
    write(dunit,'(a1)')'#'
    write(dunit,'(a4)')'calc'
    write(dunit,'(a10)',advance='no')'save pdf, '
    write(dunit,*)trim(pdffile)
    write(dunit,'(a4)')'exit'
    write(dunit,'(a1)')'#'
    write(dunit,'(a6)')'powder'
    write(dunit,'(a10)')'set dh,0.2'
    write(dunit,'(a10)')'set dk,0.2'
    write(dunit,'(a10)')'set dl,0.2'
    write(dunit,'(a13)')'set dtth,0.05'
    write(dunit,'(a15)')'set tthmin, 0.1'
    write(dunit,'(a15)')'set tthmax,95.0'
    write(dunit,'(a15)')'set axis,tth'
    write(dunit,'(a1)')'#'
    write(dunit,'(a4)')'xray'
    write(dunit,'(a12)')'set wvle,1.514'
    write(dunit,'(a15)')'set temp,ignore'
    write(dunit,'(a1)')'#'
    write(dunit,'(a3)')'run'
    write(dunit,'(a1)')'#'
    write(dunit,'(a4)')'exit'
    write(dunit,'(a1)')'#'
    write(dunit,'(a6)')'output'
    write(dunit,'(a13)')'format powder'
    write(dunit,'(a5)',advance='no')'outf '
    call extension(xyzfile,'powder',pdffile)
    write(dunit,*)trim(pdffile)
    write(dunit,'(a3)')'run'
    write(dunit,'(a1)')'#'
    write(dunit,'(a4)')'exit'
    write(dunit,'(a1)')'#'
    write(dunit,'(a4)')'exit'
    write(dunit,'(a1)')'#'

    close(dunit)

    return

  end subroutine write_discus

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------


  !----------------------------------------------------------------------------------------------
  ! The following subroutine was provided by Eric Chan, <echanj@gmail.com>
  ! and merged into the main ZMC code by Darren Goossens, August 2014.  It has not been as heavily
  ! tested as some of the other tools in this program.  Please refer to README.modwave and/or
  ! modwave_tutorial_ZMC.ppt for help. Contact Eric about it, Darren does not know!
  !----------------------------------------------------------------------------------------------
 
   subroutine modwave(key,xtal)

!    use myparams
!    use zmatrix_class
!    Use file_functions
!    use quaternion_class
!    use crystallography_class
     use vector_class !, only: operator(.cross.)
!    use varmod

    implicit none

    ! character(len=100), intent(in) :: outname
    ! character(len=100), intent(in) :: header
    ! character(len=500)             :: buffer
      character(len=30)              :: outfile


    !Integers!

   ! integer, intent(in),dimension(3)                    :: simsize
   ! integer                                             :: ia,ib,ic,il,j
    integer                                             :: ia,ib,ic,j
   ! integer                                             :: iia,iib,iic,iil,iz,im,i,dunit
    integer                                             :: iia,iib,iic,iil,im
   ! integer, intent(in)                                 :: num_loc
   ! integer,intent(in),dimension(amax,bmax,cmax,loc_max):: zocc,mocc
   ! integer, intent(in),dimension(zm_max)               :: nintern
    integer                                              :: qzocc,modtype
    integer            :: zmat_new, zmat_old, mocc_new


    !Reals!

   ! real,intent(out), dimension(amax,bmax,cmax,loc_max,3) :: xyzcrystal
   ! real, intent(out), dimension(amax,bmax,cmax,loc_max,in_max) :: incrystal
    real, dimension(3) :: trans
    real, dimension(4) :: ddd
    real, dimension(6)                                   :: cellpar
   ! real, dimension(in_max)                             ::  interns

    logical :: improp

    !Custom Types!

   ! type (quaternion)                                                   :: qcrystal
    type (quaternion) :: q_shift_factor,current_q_shift
   ! type (crystallography_object),intent(in)                           :: xtal
     type (crystallography_object)                                      :: xtal
     type (keyword_object)                                              :: key

!       real,dimension(3)  :: corr_vect

!   1st q vector parameters
       real,dimension(3)  :: qvect,qdir
       real               :: qamp,qconc
       real,dimension(3)  :: qpol
       real,dimension(3)  :: q_vect_1
       real,dimension(3)  :: Axyz_1
       real,dimension(3)  :: Axyz_norm_1
       real               :: amp_1
       real               :: Ax1,Ay1,Az1
       real,dimension(3)  :: q_carts_1, q_carts_norm_1
!   2nd q vector parameters
       real,dimension(3)  :: q_vect_2
       real,dimension(3)  :: Axyz_2
       real,dimension(3)  :: Axyz_norm_2
       real               :: amp_2
       real               :: Ax2,Ay2,Az2
       real,dimension(3)  :: q_carts_2, q_carts_norm_2

! parameters for the jana input
!       real,dimension(3)  :: cshift_real_a
!       real,dimension(3)  :: cshift_real_b
!       real,dimension(3)  :: cshift_real_c
!       real,dimension(3)  :: total_shift

       real               :: qw_shift
       real               :: qx_shift
       real               :: qy_shift
       real               :: qz_shift

       real,dimension(3)  :: basis1,basis1c
       real,dimension(3)  :: basis2,basis2c
       real,dimension(3)  :: basis3,basis3c

       real,dimension(3)  :: Astar_carts, Astar_carts_norm, astar_vect
       real,dimension(3)  :: Bstar_carts, Bstar_carts_norm, bstar_vect
       real,dimension(3)  :: Cstar_carts, Cstar_carts_norm, cstar_vect
       real,dimension(3)  :: cvec,cvec2
       real               :: current_shift_a,current_shift_b,current_shift_c
       real               :: current_shift
       real               :: astar,bstar,cstar
       real               :: stol
       real               :: dspace1,dspace2,dspace3
       real               :: phase_shift,conc

  real,parameter  ::  pi = 3.141592654

! test you can punch a hole in
          buffer = get_value(key,'QMODTYPE')
          read(buffer,*,iostat=ierr)modtype
          buffer = get_value(key,'QVECTOR')
          read(buffer,*,iostat=ierr)qvect(1),qvect(2),qvect(3)
          buffer = get_value(key,'QAMP')
          read(buffer,*,iostat=ierr)qamp
          buffer = get_value(key,'QPOL')
          read(buffer,*,iostat=ierr)qpol(1),qpol(2),qpol(3)
          buffer = get_value(key,'QZOCC')
          read(buffer,*,iostat=ierr)qzocc
          buffer = get_value(key,'QCONC')
          read(buffer,*,iostat=ierr)qconc
          buffer = get_value(key,'QDIR')
          read(buffer,*,iostat=ierr)qdir(1),qdir(2),qdir(3)

    write(stdout,*)'punched a hole in'
    write(stdout,*)'qvect' , qvect
    write(stdout,*)'check variables are passing ok '
    write(stdout,*)'simsize' , simsize
    write(stdout,*)'qamp' , qamp
    write(stdout,*)'qpol' , qpol
    write(stdout,*)'qpol' , qzocc

    cellpar(1:3) = axes(xtal)
    cellpar(4:6) = angles(xtal)

    write(stdout,*)'cellpar' , cellpar

! setup cyclic boundry conditions

!       integer bnd1(-(simsize(1)-1):simsize(1)*2)
!       integer bnd2(-(simsize(2)-1):simsize(2)*2)
!       integer bnd3(-(simsize(3)-1):simsize(3)*2)

!       bnd1 = (/ (i,i=1,simsize(1)), (i,i=1,simsize(1) ), (i,i=1,simsize(1) ) /)
!       bnd2 = (/ (i,i=1,simsize(2)), (i,i=1,simsize(2) ), (i,i=1,simsize(2) ) /)
!       bnd3 = (/ (i,i=1,simsize(3)), (i,i=1,simsize(3) ), (i,i=1,simsize(3) ) /)

     write (*,*)
     write (*,*)
     write (*,*)  'in subroutine xyzcorr'
     write (*,*)  'modified creating modulations wave vectors in the crystal'
     write (*,*)
     write (*,*)

     !  q_vect_1 = (/0.21,0.0,0.21/) !
       q_vect_1 = qvect
  !     q_vect_2 = (/-0.45,0.0,0.95/) !

 write (*,*) 'q_vect_1(:) = ' , q_vect_1(:)
! write (*,*) 'q_vect_2(:) = ' , q_vect_2(:)

basis1 = (/1.0,0.0,0.0 /)
basis2 = (/0.0,1.0,0.0 /)
basis3 = (/0.0,0.0,1.0 /)

call find_stol(basis1,stol,dspace1,astar,bstar,cstar,cellpar)
call find_stol(basis2,stol,dspace2,astar,bstar,cstar,cellpar)
call find_stol(basis3,stol,dspace3,astar,bstar,cstar,cellpar)

write (*,*) 'dspace1' , dspace1
write (*,*) 'dspace2' , dspace2
write (*,*) 'dspace3' , dspace3
write (*,*) 'astar' , astar
write (*,*) 'bstar' , bstar
write (*,*) 'cstar' , cstar

basis1c = as_cartesian(xtal,basis1)
basis2c = as_cartesian(xtal,basis2)
basis3c = as_cartesian(xtal,basis3)

astar_carts(:) = basis2c .cross. basis3c
bstar_carts(:) = basis3c .cross. basis1c
cstar_carts(:) = basis1c .cross. basis2c

! write (*,*) astar_carts(:)
! write (*,*) astar_carts/sqrt(sum(astar_carts**2))   !  normalised astar as cartesian coordinates
! write (*,*) bstar_carts(:)
! write (*,*) bstar_carts/sqrt(sum(bstar_carts**2))   !  normalised bstar as cartesian coordinates
! write (*,*) cstar_carts(:)
! write (*,*) cstar_carts/sqrt(sum(cstar_carts**2))   !  normalised cstar as cartesian coordinates

astar_carts_norm =  astar_carts/sqrt(sum(astar_carts**2))   !  normalised astar as cartesian coordinates
bstar_carts_norm =  bstar_carts/sqrt(sum(bstar_carts**2))   !  normalised astar as cartesian coordinates
cstar_carts_norm =  cstar_carts/sqrt(sum(cstar_carts**2))   !  normalised astar as cartesian coordinates

! now express  the A*, B*, C*  magnitudes in terms of cartesian coordinate vectors
! eg. the components of the unit vector Astar_carts_norm act as direction cosines for the scalar magnitude of A*

astar_vect(:) = (/astar*astar_carts_norm(1),astar*astar_carts_norm(2),astar*astar_carts_norm(3)/)   !  normalised astar as cartesian coordinates
bstar_vect(:) = (/bstar*bstar_carts_norm(1),bstar*bstar_carts_norm(2),bstar*bstar_carts_norm(3)/)
cstar_vect(:) = (/cstar*cstar_carts_norm(1),cstar*cstar_carts_norm(2),cstar*cstar_carts_norm(3)/)

write (*,*) astar_vect(:)
write (*,*) bstar_vect(:)
write (*,*) cstar_vect(:)

 ! now sort out the q components of each of these vectors and place in one final cartesian vector

  q_carts_1(:) =   (/q_vect_1(1)*astar_vect(1),q_vect_1(1)*astar_vect(2),q_vect_1(1)*astar_vect(3)/)  &
                 + (/q_vect_1(2)*bstar_vect(1),q_vect_1(2)*bstar_vect(2),q_vect_1(2)*bstar_vect(3)/)  &
                 + (/q_vect_1(3)*cstar_vect(1),q_vect_1(3)*cstar_vect(2),q_vect_1(3)*cstar_vect(3)/)

 write (*,*) 'q_carts_1 = ' , q_carts_1(:)

 ! normmalise the direction into a unit vector
        q_carts_norm_1(:) = q_carts_1(:)/sqrt(sum(q_carts_1(:)**2))   !  normalised qvector

 write (*,*) 'q_carts_norm_1 = ' , q_carts_norm_1(:)

 if ((qpol(1)+qpol(2)+qpol(3)) .eq. 0.0) then
 ! this makes the displacements in the direction q  i.e. longitudinal
       Axyz_norm_1(:) = q_carts_norm_1(:)
 else

 ! this makes displacments transverse  but is ad hoc depending on the q vector
     !   Axyz_norm_1(:) = (/-q_carts_norm_1(2),q_carts_norm_1(1),q_carts_norm_1(3)/)

 ! this makes displacments some component of q-carts norm
    !    Axyz_1(:) = (/ qpol(1)*q_carts_norm_1(2),qpol(2)*q_carts_norm_1(1),qpol(3)*q_carts_norm_1(3)/)
    !this results in errors

 ! you could specify some other arbitary direction if you like
 ! this you specify near the q-vector  (next two lines)
    !    Axyz_1(:)= (/1.0,0.0,0.0 /)    ! polarisation in the c* direction
        Axyz_1(:)= qpol(:)    ! polarisation in the reletive to k-space setting
        Axyz_norm_1(:) = Axyz_1/sqrt(sum(Axyz_1(:)**2))   !
 endif

       amp_1 = qamp
       Ax1=Axyz_norm_1(1)*amp_1
       Ay1=Axyz_norm_1(2)*amp_1
       Az1=Axyz_norm_1(3)*amp_1

 ! beginging of parameter code for second q_vector

  !    q_carts_2(:) = (/q_vect_2(1)*astar_vect(1),q_vect_2(1)*astar_vect(2),q_vect_2(1)*astar_vect(3)/)  &
  !                 + (/q_vect_2(2)*bstar_vect(1),q_vect_2(2)*bstar_vect(2),q_vect_2(2)*bstar_vect(3)/)  &
  !                 + (/q_vect_2(3)*cstar_vect(1),q_vect_2(3)*cstar_vect(2),q_vect_2(3)*cstar_vect(3)/)

 !write (*,*) 'q_carts_2 = ' , q_carts_2(:)

 !      q_carts_norm_2(:) = q_carts_2(:)/sqrt(sum(q_carts_2(:)**2))   !  normalised qvector

 !write (*,*) 'q_carts_norm_2 = ' , q_carts_norm_2(:)

 ! this makes displacments transverse  but is ad hoc depending on the q vector
 !      Axyz_norm_2(:) = q_carts_norm_2(:)


  !     amp_2 = 1.0
  !     Ax2=Axyz_norm_2(1)*amp_2
  !     Ay2=Axyz_norm_2(2)*amp_2
  !     Az2=Axyz_norm_2(3)*amp_2

 !  seperate inputs modified for the Jana parameters
 !  which i have not checked hold true 100%
 ! but its worth keeping some of the code structure in case the problem arizes
  !     Ax = 0.1
  !     Ay = 0.065
  !     Az = 0.00


if (modtype.eq.1) then

 ! test the direction

    write(*,*) 'dot_product(q_carts_norm_1(:),Axyz_1(:))' , dot_product(q_carts_norm_1(:),Axyz_norm_1(:))

        current_shift_a = 0.0   ! reinitialize current shift
        current_shift_b = 0.0   ! reinitialize current shift
        current_shift_c = 0.0   ! reinitialize current shift

   write (*,*)  cellpar

      do ia = 1,simsize(1)
         do ib = 1,simsize(2)
           do ic = 1,simsize(3)
               do il = 1,num_loc

   if (zocc(ia,ib,ic,il).eq.qzocc) then
  !  if (il.eq.1.or.il.eq.2.or.il.eq.3.or.il.eq.4.or.il.eq.9.or.il.eq.10.or.il.eq.11.or.il.eq.12) then

             cvec(:) = as_cartesian(xtal,real((/(ia-1),(ib-1),(ic-1)/)))
             cvec2(:) = cvec(:)+xyzcrystal(ia,ib,ic,il,:)

      current_shift_a = Ax1*cos((2*pi)*dot_product(q_carts_1,cvec2)) !
      current_shift_b = Ay1*cos((2*pi)*dot_product(q_carts_1,cvec2)) !
      current_shift_c = Az1*cos((2*pi)*dot_product(q_carts_1,cvec2)) !

    !     current_shift_a = current_shift_a + Ax2*cos((2*pi)*dot_product(q_carts_2,cvec2))  ! cosine modulation
    !     current_shift_b = current_shift_b + Ay2*cos((2*pi)*dot_product(q_carts_2,cvec2))  ! cosine modulation
    !     current_shift_c = current_shift_c + Az2*cos((2*pi)*dot_product(q_carts_2,cvec2))  ! cosine modulation

          xyzcrystal(ia,ib,ic,il,1) = xyzcrystal(ia,ib,ic,il,1) + current_shift_a
          xyzcrystal(ia,ib,ic,il,2) = xyzcrystal(ia,ib,ic,il,2) + current_shift_b
          xyzcrystal(ia,ib,ic,il,3) = xyzcrystal(ia,ib,ic,il,3) + current_shift_c

 !  Variation of the shift parameters if needed for the Jana style input

      !    cshift_real_a(:) = as_cartesian(xtal,(/current_shift_a,0.0,0.0/))
      !    cshift_real_b(:) = as_cartesian(xtal,(/0.0,current_shift_b,0.0/))
      !    cshift_real_c(:) = as_cartesian(xtal,(/0.0,0.0,current_shift_c/))

      !    total_shift(:) = cshift_real_a + cshift_real_b + cshift_real_c
      !    xyzcrystal(ia,ib,ic,il,:) = xyzcrystal(ia,ib,ic,il,:) + total_shift(:)

                   endif ! if iz


                enddo
              enddo
            enddo
          enddo
 !
 !
 ! ! now just have to try write those coords to a new crystal
 !
 !    no longer have to write crystal but need to make sure parameter changes are made to the crystal object
 !      call writecrystal(incrystal,xyzcrystal,qcrystal,simsize,num_loc                        &
 !          ,            outname,header,nintern,zocc)

elseif (modtype.eq.2) then

! begin occupancy modulation

! test the direction

   write(*,*) 'dot_product(q_carts_norm_1(:),Axyz_1(:))' , dot_product(q_carts_norm_1(:),Axyz_norm_1(:))

       current_shift_a = 0.0   ! reinitialize current shift
       current_shift_b = 0.0   ! reinitialize current shift
       current_shift_c = 0.0   ! reinitialize current shift

  write (*,*)  cellpar

        outfile = "modulate.occ"

 dunit = open(outfile, status='new')

     do ia = 1,simsize(1)
        do ib = 1,simsize(2)
          do ic = 1,simsize(3)
              do il = 1,num_loc
            !     elseif (il .gt. 108) then

 ! initialise variables  to prevent memory leaks
        zmat_new = zocc(ia,ib,ic,il)
        mocc_new = mocc(ia,ib,ic,il)

  if (zocc(ia,ib,ic,il).eq.qzocc) then

           ! zmat_old = zocc(ia,ib,ic,il)  ! only need this for a different type of modulation

            cvec(:) = as_cartesian(xtal,real((/(ia-1),(ib-1),(ic-1)/)))
            cvec2(:) = cvec(:)+xyzcrystal(ia,ib,ic,il,:)

! in the case of occupancy I only need to worry about one modulation

        phase_shift = 0.0 ! this will complicate things so leave it out for now
      !   conc   = 0.10 ! current setting that works
        conc = qconc
        mocc_new = mocc(ia,ib,ic,il)

        current_shift = cos((2*pi)*dot_product(q_carts_1,cvec2) + phase_shift )   ! cosine modulation
        current_shift = (current_shift+1)/2 ! calibrate values between 0.0 to 1.0

         if (current_shift.lt.conc) zmat_new =  qzocc   !zmat_old
         if (current_shift.ge.conc) zmat_new =  qzocc+1

                  endif ! if iz = 2


 write (dunit,*) ia, ib, ic, il, zmat_new, mocc_new

               enddo
             enddo
           enddo
         enddo
!
 close(dunit)

!  begin optional quaternion modulation

elseif (modtype.eq.3) then

 ! initialise variables
         current_shift = amp_1
         current_q_shift = (/0.5, 0.5,  0.5,   0.5/)

!        outfile = "test.occ"
           write (*,*)   'current shift is', current_shift  , 'degrees - use qamp to adjust magnitude '

! dunit = open(outfile, status='new')

     do ia = 1,simsize(1)
        do ib = 1,simsize(2)
          do ic = 1,simsize(3)
              do il = 1,num_loc

  if (zocc(ia,ib,ic,il).eq.qzocc) then

            cvec(:) = as_cartesian(xtal,real((/(ia-1),(ib-1),(ic-1)/)))
            cvec2(:) = cvec(:)+xyzcrystal(ia,ib,ic,il,:)

! in the case of occupancy I only need to worry about one modulation

        phase_shift = 0.0

         current_shift = amp_1*cos((2*pi)*dot_product(q_carts_1,cvec2) + phase_shift )   ! cosine modulation
      !  current_shift = sin(current_shift) ! sin of current shift

         ! convert to degrees to radians
              current_shift = (current_shift/180)*pi

      qw_shift = cos((current_shift)/2)
      qx_shift = qdir(1)*sin((current_shift)/2)
      qy_shift = qdir(2)*sin((current_shift)/2)
      qz_shift = qdir(3)*sin((current_shift)/2)

      current_q_shift = (/qw_shift,  qx_shift,  qy_shift,   qz_shift/)

!      q1=(/-0.315403,0.677179,0.658814,0.088957/) ! test quaternion
!      q2 = current_q_shift*q1

!         if (current_shift.lt.conc) zmat_new = 1
!         if (current_shift.ge.conc) zmat_new = 2

! write (dunit,*) ia, ib, ic, il, zmat_new, mocc_new

!  write (*,*) ia, ib, ic, il, as_array(q2)

!   how to manipulate a quaternion
                    qcrystal(ia,ib,ic,il) = current_q_shift*qcrystal(ia,ib,ic,il)

        endif ! zocc eq qzocc

!  in order to print the final quaternion
!           write (*,*)   as_array(q1)

               enddo
             enddo
           enddo
         enddo
!
! close(dunit)



! end occupancy modulation
endif


    return

      end subroutine modwave

  !----------------------------------------------------------------------------------------------
  ! The following subroutine was provided by Eric Chan, <echanj@gmail.com>
  ! and merged into the main ZMC code by Darren Goossens, August 2014.  It has not been as heavily
  ! tested as some of the other tools in this program.  Please refer to README.modwave and/or
  ! modwave_tutorial_ZMC.ppt for help. Contact Eric about it. Darren does not know!
  !----------------------------------------------------------------------------------------------


SUBROUTINE  find_stol(hkl,stol,d_space,astar,bstar,cstar,cellpar)

! use myparams
! use fundamental_constants
implicit none

!  integer,intent(IN)  :: numu, numv, numw

 integer                                          :: i,j,k
 real ,intent(OUT)                                :: stol
 real ,intent(OUT)                                :: d_space
! real, dimension(3),intent(IN)                    :: hkl
 real, dimension(3)                               :: hkl
 real                                             :: vol
 real                                             :: alphar,betar,gammr
 real ,intent(OUT)                                :: astar,bstar,cstar
 real                                             :: alphastar,betastar,gammastar
 real                                             :: h_in,k_in,l_in
 real                                             :: ssqtheta,stheta
! real,intent(in), dimension(6)                    :: cellpar
 real, dimension(6)                    :: cellpar

! ad hoc input of cell params really messy so fix this for later

  real  :: cell_a
  real  :: cell_b
  real  :: cell_c
  real  :: alpha
  real  :: beta
  real  :: gamm
  real  :: lambda

  real,parameter  ::  pi = 3.141592654

! use the hkl values to find the stol value given the standard cell dimensions
! based on equations for sstheta for triclinic cell

! print *, 'in find stol'

  cell_a =  cellpar(1)
  cell_b =  cellpar(2)
  cell_c =  cellpar(3)
  alpha  =  cellpar(4)
  beta   = cellpar(5)
  gamm   = cellpar(6)

alphar = alpha*pi/180.0
betar  = beta*pi/180.0
gammr  = gamm*pi/180.0

 Vol = cell_a*cell_b*cell_c* sqrt(1+2*cos(alphar)*cos(betar)*cos(gammr) &
       -cos(alphar)**2  &
       -cos(betar)**2  &
       -cos(gammr)**2  )

 astar =  (cell_b*cell_c/Vol)*sin(alphar)
 bstar =  (cell_c*cell_a/Vol)*sin(betar)
 cstar =  (cell_a*cell_b/Vol)*sin(gammr)

! in radians now
 alphastar = acos((cos(betar)*cos(gammr)-cos(alphar))/(sin(betar)*sin(gammr)))
 betastar = acos((cos(gammr)*cos(alphar)-cos(betar))/(sin(gammr)*sin(alphar)))
 gammastar = acos((cos(alphar)*cos(betar)-cos(gammr))/(sin(alphar)*sin(betar)))

! print *, 'alphar' , alphar
! print *, 'betar' , betar
! print *, 'gammr' , gammr
! print *, 'vol' , vol
!
! print *, 'astar' , astar
! print *, 'bstar' , bstar
! print *, 'cstar' , cstar

! print *, 'alphastar' , alphastar
! print *, 'betastar' , betastar
! print *, 'gammastar' , gammastar

! print *, 'h' , hkl(1)
! print *, 'k' , hkl(2)
! print *, 'l' , hkl(3)

! print  *, 'calculating stol and 2theta values'

h_in = hkl(1)
k_in = hkl(2)
l_in = hkl(3)

     ssqtheta = ((lambda**2)/4)*                             &
               (                                             &
               ((h_in**2) * (astar**2))+                     &
               ((k_in**2) * (bstar**2))+                     &
               ((l_in**2) * (cstar**2))+                     &
               (2*k_in*l_in * bstar*cstar*cos(alphastar))+   &
               (2*l_in*h_in * cstar*astar*cos(betastar))+    &
               (2*h_in*k_in * astar*bstar*cos(gammastar))    &
               )

  stheta=sqrt(ssqtheta)

     d_space=lambda/(2.*stheta)
!       d_star=1/d

!   print *, stheta/lambda

  stol = stheta/lambda

   print  *, stol
   print  *, '.'

! print *, stol(i,j,k)

 return

end SUBROUTINE  find_stol

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

end module mainmod


program ZMC

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------

  use cmdline_arguments
  use mainmod

  implicit none

  ! Version string.  Update when change program, maintain its length

  cver = '06 March 2026                                     '

  !------------------------------------------------------------------------------!
  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  !------------------------------------------------------------------------------!

  myoptions(1)  = 'reread'
  myoptions(2)  = 'help'
  myoptions(3)  = 'pairs'
  myoptions(4)  = 'fracs'
  myoptions(5)  = 'cartsn'
  myoptions(6)  = 'cartst'
  myoptions(7)  = 'corro'
  myoptions(8)  = 'corra'
  myoptions(9)  = 'energy'
  myoptions(10) = 'quiet' 
  myoptions(11) = 'crystal'
  myoptions(12) = 'version'
  myoptions(13) = 'cif'
  myoptions(14) = 'summary'
  myoptions(15) = 'plot'
  myoptions(16) = 'getcontacts'
  myoptions(17) = 'help2'
  myoptions(18) = 'diffuse'
  myoptions(19) = 'discus'
  myoptions(20) = 'modwave'

  !------------------------------------------------------------------------------!
  !     Keywords for use in input files.  See subroutine longusage() for details
  !------------------------------------------------------------------------------!

  mykeywords(1) = 'CRYSTAL'      
  mykeywords(2) = 'HEADER'      
  mykeywords(3) = 'TEMPERATURE'
  mykeywords(4) = 'MCCYCLES'  
  mykeywords(5) = 'INCUPDATE'
  mykeywords(6) = 'BADJUST' 
  mykeywords(7) = 'NUMSPRCON'      
  mykeywords(8) = 'NUMINSPRCON'   
  mykeywords(9) = 'NUMLOCS'      
  mykeywords(10) = 'NUMZMATS'   
  mykeywords(11) = 'NUMINTERNAL'
  mykeywords(12) = 'NUMCROSS'  
  mykeywords(13) = 'OCCFILE'  
  mykeywords(14) = 'ZMATFILE'
  mykeywords(15) = 'QXYZFILE'     
  mykeywords(16) = 'CELL'        
  mykeywords(17) = 'CONTACTFILE'
  mykeywords(18) = 'SPRCON'    
  mykeywords(19) = 'SIZE'     
  mykeywords(20) = 'INSPR'   
  mykeywords(21) = 'INTERNAL'     
  mykeywords(22) = 'CROSSSPR'    
  mykeywords(23) = 'CROSSDEF'   
  mykeywords(24) = 'XYZWIDTH'  
  mykeywords(25) = 'XYZINITW' 
  mykeywords(26) = 'QWIDTH'  
  mykeywords(27) = 'QINITW' 
  mykeywords(28) = 'INWIDTH'      
  mykeywords(29) = 'ININITW'  
  mykeywords(30) = 'POTENTIAL'
  mykeywords(31) = 'VMIN' 
  mykeywords(32) = 'VMAX'
  mykeywords(33) = 'NEWCONTACTS'
  mykeywords(34) = 'CONTATOMS'
  mykeywords(35) = 'QVECTOR'
  mykeywords(36) = 'QAMP'
  mykeywords(37) = 'QPOL'
  mykeywords(38) = 'QZOCC'
  mykeywords(39) = 'QMODTYPE'
  mykeywords(40) = 'QCONC'
  mykeywords(41) = 'QDIR'

  ! By using a seeded random number generator with known seeds, you can compare
  ! exactly. BUT it means there is little point in averaging the results of 
  ! identical simulations, because  they will be the same.
  ! Hence, TO BE modified to red in seeds from the input file IF they are
  ! given in the file. NOT YET DONE.
  ! Needs a new keyword (RANDOMSEEDS) which is attached to 2 integers.

  i1 = 12375
  i2 = 9875

  call rseed(i1,i2)

  !------------------------------------------------------------------------------!
  ! This call parses the command line arguments for command line options
  !------------------------------------------------------------------------------!

  call get_options(myoptions, ierr)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (ierr > 0) then
     write(stderr,*) 'ERROR! Unknown option(s): ',join(bad_options()," ")
     call usage()
     stop
  end if

  if((option_exists('version')).or.(option_exists('help')).or.(option_exists('help2')))then
     if (option_exists('version')) then
        write(stdout,*)'----------------------------------------------------------------------------'
        write(stdout,*)'ZMC version: ',cver
        write(stdout,*)'----------------------------------------------------------------------------'
     end if
     if (option_exists('help')) call usage()
     if (option_exists('help2')) call longusage()
     stop
  end if

  quiet = .FALSE.
  if (option_exists('quiet')) quiet = .TRUE.

  !------------------------------------------------------------------------------!
  !                      Initialise keyword list.                                !
  !------------------------------------------------------------------------------!

  call initialise(key2, mykeywords)

  !------------------------------------------------------------------------------!
  !                      Open keyword file (command line argument)
  !------------------------------------------------------------------------------!

     write(stdout,*)'----------------------------------------------------------------------------'
     write(stdout,*)'This is ZMC version dated ',trim(cver),', a program for implementing '
     write(stdout,*)'Monte Carlo simulations of crystal structures, primarily for the'
     write(stdout,*)'analysis of diffuse scattering.'
     write(stdout,*)'                                 darren.goossens@gmail.com'
     write(stdout,*)' '
     write(stdout,*)'The modwave() option is provided by Eric Chan, echanj@gmail.com.'
     write(stdout,*)' '
     write(stdout,*)'All use should reference the paper:'
     write(stdout,*)"D. J. Goossens, A. P. Heerdegen, E. J. Chan & T. R. Welberry, 'Monte"
     write(stdout,*)"Carlo Modelling of Diffuse Scattering from Single Crystals: The "
     write(stdout,*)"Program ZMC' Metallurgical and Materials Transactions A, vol 42A, "
     write(stdout,*)"(2011) 23-31.  DOI:10.1007/s11661-010-0199-1"
     write(stdout,*)""
     write(stdout,*)"ZMC is licenced under the Academic Free License version 3.0; details"
     write(stdout,*)"in the files DISCLAIMER.txt and COPYRIGHT.txt and at"
     write(stdout,*)"http://opensource.org/licenses/AFL-3.0"
     write(stdout,*)' '
     write(stdout,*)'Invoke it with ZMC --help for some help, and ZMC --help2 for more help.'
  if ( .not. have_args() ) then
     write(stdout,*)""
     write(stderr,*) "Must provide a keyword file filename or a command line option."
     write(stderr,*) 'Program exiting.'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if   ! arg test
     write(stdout,*)'----------------------------------------------------------------------------'

  if(.not. quiet) then       
  end if
  keyfile = next_arg()  

  if (.not. quiet)write(stdout,*)'Input file is ',trim(keyfile)
  ! Read values from text file
  if (.not. read(key2, keyfile)) then
     write(stderr,*)'----------------------------------------------------------------------------'
     write(stderr,*)'Problem with keyword file ',trim(keyfile)
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if
  if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

  !------------------------------------------------------------------------------!
  ! Now get root of outname. If not given, use keyfile
  !------------------------------------------------------------------------------!

  if ( .not. have_args() ) then
     if(.not.quiet) write(stdout,*) "Keyword filename will be used as root name for output files."
     outname = keyfile
  else
     outname  = next_arg()  
  end if   ! arg test
  call extension(outname,' ',fname2) 
  if(.not.quiet) write(stdout,*) "Output files have root name ",trim(fname2)
  if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

  !------------------------------------------------------------------------------!
  ! Initialisation common to all three subsections:
  !------------------------------------------------------------------------------!

  call get_cell(xtal,key2,quiet)
  call get_optional_integer(num_loc_given,'NUMLOCS',key2,quiet)
  call get_optional_integer(given_total_num_zm,'NUMZMATS',key2,quiet)
  call get_total_zm_files(total_num_zm,given_total_num_zm,key2,quiet)
  allocate(zmat(total_num_zm),stat=ierr)
  if (ierr.ne.0) call allerr('zmat(total_num_zm)                                  ')
  allocate(num_at(total_num_zm),stat=ierr)
  if (ierr.ne.0) call allerr('num_at(total_num_zm)                                ')
  call get_zmats(zmat,total_num_zm,num_at,key2,quiet)

  !------------------------------------------------------------------------------!
  ! If --getcontacts is specified, generate some contact vectors and exit.
  !------------------------------------------------------------------------------!

  if (option_exists('getcontacts'))  then

     i1 = num_value(key2, 'QXYZFILE')
     if((i1.lt.1).and.(num_loc_given.eq.-9999)) then
        write(stderr,*)'Must give NUMLOCS and/or QXYZFILE and/or OCCFILE to specify number of locations.'
        write(stderr,*)'----------------------------------------------------------------------------'
        stop
     end if

     call get_num_loc_qxyz(num_loc_given, num_loc_qxyz, total_num_zm,key2,quiet)

     if((num_loc_given.eq.-9999).and.(num_loc_qxyz.eq.-9999)) then
        write(stderr,*)'----------------------------------------------------------------------------' 
        write(stderr,*)'Number of locations (molecular sites) not given!'
        write(stderr,*)'Program exiting.'             
        write(stderr,*)'----------------------------------------------------------------------------' 
        stop
     else
        num_loc = max(num_loc_given,num_loc_qxyz)
     end if

     call get_maxval_mocc(maxval_mocc)
     allocate(quatern(num_loc,total_num_zm,maxval_mocc),stat=ierr)
     if (ierr.ne.0) call allerr('quatern(num_loc,total_num_zm,maxval_mocc)           ')
     allocate(translate(num_loc,total_num_zm,maxval_mocc ,3),stat=ierr)
     if (ierr.ne.0) call allerr('translate(num_loc,total_num_zm,maxval_mocc, 3)      ')
     call read_qxyz(quatern,translate)
     allocate(how_many_types(num_loc),stat=ierr)
     if (ierr.ne.0) call allerr('how_many_types(num_loc)                             ')
     allocate(which_types(num_loc,total_num_zm),stat=ierr)
     if (ierr.ne.0) call allerr('which_types(num_loc,total_num_zm)                   ')
     allocate(instances_types(num_loc,total_num_zm),stat=ierr)
     if (ierr.ne.0) call allerr('instances_types(num_loc,total_num_zm)               ')
     allocate(which_mol(num_loc,total_num_zm,maxval_mocc),stat=ierr)
     if (ierr.ne.0) call allerr('which_mol(num_loc,total_num_zm,maxval_mocc)         ')
     call get_what_goes_where(how_many_types,which_types,instances_types,which_mol,total_num_zm,num_loc,quiet,key2)
     call get_new_contacts()

     stop

  end if ! getcontacts

  !------------------------------------------------------------------------------!
  ! Initialisation common to plot and MC
  !------------------------------------------------------------------------------!

  i1 = (num_value(key2, 'QXYZFILE')+num_value(key2, 'OCCFILE'))
  if((i1.lt.1).and.(num_loc_given.eq.-9999)) then
     write(stderr,*)'Must give NUMLOCS and/or QXYZFILE and/or OCCFILE to specify number of locations.'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if

  call get_num_loc_qxyz(num_loc_given, num_loc_qxyz, total_num_zm,key2,quiet)
  allocate(nintern(total_num_zm),stat=ierr)
  if (ierr.ne.0) call allerr('nintern(total_num_zm)                               ')
  call get_ninternal(nintern,total_num_zm,key2,quiet)
  allocate(intern(maxval(num_at),2,total_num_zm),stat=ierr)
  if (ierr.ne.0) call allerr('intern(maxval(num_at),2,total_num_zm)               ')
  call get_internals(nintern,total_num_zm,intern,key2,quiet,num_at) 

  call get_simsize(simsize,key2,quiet)
  num_loc_occ = -9999
  i1 = num_value(key2, 'OCCFILE')
  if(i1.lt.1) then
     if(.not.quiet)write(stdout,*)'No occupancy file given.  Will try to work out'
     if(.not.quiet)write(stdout,*)'whether one is needed from QXYZ files.'
     if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
     i2 = num_value(key2, 'QXYZFILE')
     if(i2.lt.1) then
        write(stderr,*)'No occupancy or QXYZ files given.  Exiting.'
        stop
     end if
  else
     if(i1.gt.1) then
        write(stderr,*)'Too many occupancy files given! Exiting.'
        stop
     end if
     call get_num_loc_occ(num_loc_given,num_loc_qxyz,num_loc_occ,total_num_zm,mol_max,key2,quiet,simsize)
  end if

  if((num_loc_occ.eq.-9999).and.(num_loc_given.eq.-9999).and.(num_loc_qxyz.eq.-9999)) then
     write(stderr,*)'----------------------------------------------------------------------------' 
     write(stderr,*)'Number of locations (molecular sites) not given!'
     write(stderr,*)'Program exiting.'             
     write(stderr,*)'----------------------------------------------------------------------------' 
     stop
  else
     num_loc = max(num_loc_occ,num_loc_given,num_loc_qxyz)
  end if

  allocate(mocc(simsize(1),simsize(2),simsize(3),num_loc),stat=ierr)
  if (ierr.ne.0) call allerr('mocc(simsize(1),simsize(2),simsize(3),num_loc)      ')
  allocate(zocc(simsize(1),simsize(2),simsize(3),num_loc),stat=ierr)
  if (ierr.ne.0) call allerr('zocc(simsize(1),simsize(2),simsize(3),num_loc)      ')
  allocate(xyzcrystal(simsize(1),simsize(2),simsize(3),num_loc,3),stat=ierr)
  if (ierr.ne.0) call allerr('xyzcrystal(simsize(1),simsize(2),simsize(3),numloc,3')
  allocate(incrystal(simsize(1),simsize(2),simsize(3),num_loc,maxval(num_at)),stat=ierr)
  if (ierr.ne.0) call allerr('incrystal()                                         ')
  allocate(qcrystal(simsize(1),simsize(2),simsize(3),num_loc),stat=ierr)
  if (ierr.ne.0) call allerr('qcrystal(simsize(1),simsize(2),simsize(3),num_loc)  ')
  allocate(carts(3,simsize(1),simsize(2),simsize(3),num_loc,maxval(num_at)),stat=ierr)
  if (ierr.ne.0) call allerr('carts()                                             ')

  call get_zocc_mocc(zocc,mocc,simsize,num_loc,total_num_zm,key2,quiet)

  if (option_exists('reread')) then
     if(has_value('reread'))then
        inname = get_value('reread')
        call  readcrystal(inname,outname)
     else
        write(stderr,*)'Please specify file to reread the simulation from.'
        write(stderr,*)'(--reread needs a filename, --reread=filename.)'
        write(stderr,*) 'Program exiting.'
        write(stderr,*)'----------------------------------------------------------------------------'
        stop
     end if
  else
     ! Construct a starting state from the given parameters
     allocate(initw((7+maxval(num_at)),total_num_zm),stat=ierr)
     if (ierr.ne.0) call allerr('initw((7+maxval(num_at)),total_num_zm)              ')
     if(.not.quiet)write(stdout,*)'Constructing crystal'
     initw = 0.0
     call get_w(initw,'i')
     allocate(quatern(num_loc,total_num_zm,maxval(mocc)),stat=ierr)
     if (ierr.ne.0) call allerr('quatern(num_loc,total_num_zm,maxval(mocc))          ')
     allocate(translate(num_loc,total_num_zm,maxval(mocc),3),stat=ierr)
     if (ierr.ne.0) call allerr('translate(num_loc,total_num_zm,maxval(mocc),3)      ')
     call read_qxyz(quatern,translate)
     call get_maxval_mocc(maxval_mocc)
     call replicator()
     call replinternal()
  end if

  if (.not. quiet)write(stdout,*)'Model crystal populated.'
  if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

  call getcarts()

  !------------------------------------------------------------------------------!
  ! If option --plot is given, plot some stuff interactively and exit.
  !------------------------------------------------------------------------------!

  if (option_exists('plot'))  then

     ntyp = 1
     con_test=-1

     if (exists(key2,'CONTACTFILE')) then
        allocate(conmol(num_loc,maxval(mocc),total_num_zm),stat=ierr)
        if (ierr.ne.0) call allerr('conmol(num_loc,maxval(mocc),total_num_zm)           ')
        con_test = 1
        call read_in_contacts()
     end if
     call plot(con_test)
     stop
  end if !! if option exists plot

  !------------------------------------------------------------------------------!
  ! OK, if not getting contacts or plotting, proceed to MC.                                
  ! Must specify at least one of --crystal or --diffuse (and later --discus)
  !------------------------------------------------------------------------------!

  if ((.not.option_exists('crystal')).and.(.not.option_exists('diffuse')).and.(.not.option_exists('discus')))  then
     write(stderr,*)'----------------------------------------------------------------------------'
     write(stderr,*) 'Must specify one or more of --crystal, --diffuse, --discus'
     write(stderr,*) '--getcontacts or --plot or there will be no output!'
     write(stderr,*) 'Program exiting.'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if

  if (exists(key2,'CONTACTFILE')) then
     allocate(conmol(num_loc,maxval(mocc),total_num_zm),stat=ierr)
     if (ierr.ne.0) call allerr('conmol(num_loc,maxval(mocc),total_num_zm)           ')
     call read_in_contacts()
  else
     write(stderr,*)'----------------------------------------------------------------------------'
     write(stderr,*) 'Must provide a list of contact vectors via CONTACTFILE.'
     write(stderr,*) 'Program exiting.'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if

  if(.not.quiet) write(stdout,'(" Contact vector file contains",i4," contact vector types.")')ntyp
  if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

  allocate(mcw((7+maxval(num_at)),total_num_zm),stat=ierr)
  if (ierr.ne.0) call allerr('  mcw((7+maxval(num_at)),total_num_zm)              ')
  mcw   = 0.0
  call get_w(mcw,'w')
  allocate(widths(7,total_num_zm),stat=ierr)   
  if (ierr.ne.0) call allerr('widths(7,total_num_zm)                              ')
  allocate(inwidths(maxval(num_at),total_num_zm),stat=ierr)
  if (ierr.ne.0) call allerr('inwidths(maxval(num_at),total_num_zm)               ')

  inwidths = 0.0
  do i1 = 1,total_num_zm
     do i = 1,7
        widths(i,i1)= mcw(i,i1)
     end do
     do i = 1,nintern(i1)
        inwidths(i,i1) = mcw(i+7,i1)
     end do
  end do

  call get_required_integer(ncy,'MCCYCLES')
  call get_required_real(temp,'TEMPERATURE') 
  adjcy = 0
  call get_optional_integer(adjcy,'INCUPDATE',key2,quiet)
  if(adjcy.eq.-9999) adjcy = 0
  radjcy  = real(adjcy)
  badjcy = 0 
  call get_optional_int_and_real(badjcy,B_aim,'BADJUST') 
  if(badjcy.eq.-9999) badjcy = 0
  bradjcy = real(badjcy)
  call get_optional_integer(ninspr_given,'NUMINSPRCON',key2,quiet)
  call get_optional_integer(ntyp_given,'NUMSPRCON',key2,quiet)

  if ((ntyp_given.ne.ntyp).and.(ntyp_given.ne.-9999)) then
     write(stderr,*)'----------------------------------------------------------------------------'
     write(stderr,*)'Number of contact types in list is: ',ntyp
     write(stderr,*)'This is different to number specified by NUMSPRCON: ',ntyp_given
     write(stderr,*)'Program exiting.'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if

  allocate(sprcon(ntyp),stat=ierr)
  if (ierr.ne.0) call allerr('sprcon(ntyp)                                        ')
  allocate(size_effect(ntyp),stat=ierr)
  if (ierr.ne.0) call allerr('size_effect(ntyp)                                   ')

  sprcon = -9999.0
  call get_sprcon_size(sprcon,'SPRCON', ntyp)
  if (minval(sprcon).eq.-9999.0) then
     write(stderr,*)'----------------------------------------------------------------------------'
     write(stderr,*)'Some contact vector force constants not set!'
     write(stderr,*)'Program exiting.'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if

  size_effect = 0.0
  call get_sprcon_size(size_effect,'SIZE', ntyp)

  if(.not.quiet) then
     do i = 1,ntyp
        write(stdout,'(" Contact type: ",i4,"  Sprcon: ",f7.3,"  Size effect: ",f7.3)')i,sprcon(i),size_effect(i)
     end do
     write(stdout,*)'----------------------------------------------------------------------------'
  end if

  ninspr = num_value(key2,'INSPR')
  if((ninspr.ne.ninspr_given).and.(ninspr_given.ne.-9999)) then
     write(stderr,*)'----------------------------------------------------------------------------'
     write(stderr,*)'Number of  internal force constants appears to be: ',ninspr
     write(stderr,*)'This is different to number specified by NUMINSPRCON: ',ninspr_given
     write(stderr,*)'Program exiting.'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if
  if((ninspr.lt.1).and.(num_value(key2,'INTERNAL').gt.0)) then
     write(stderr,*)'----------------------------------------------------------------------------'
     write(stderr,*)' ',num_value(key2,'INTERNAL'),' internal degrees of freedom defined.'
     write(stderr,*)'but no internal force constants!'
     write(stderr,*)'Program exiting.'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if

  allocate(insprcon(ninspr),stat=ierr)
  if (ierr.ne.0) call allerr('insprcon(ninspr)                                    ')
  insprcon = 0.
  allocate(insprtype(maxval(num_at),total_num_zm),stat=ierr)
  if (ierr.ne.0) call allerr('insprtype(maxval(num_at),total_num_zm)              ')
  call get_inspr(insprcon, insprtype) 
  allocate(ncross_given(total_num_zm),stat=ierr)
  if (ierr.ne.0) call allerr('ncross_given(total_num_zm)                          ')
  call get_ncross(ncross_given) 
  call get_cross()

  !------------------------------------------------------------------------------!
  ! At last we have all the information we need to proceed.
  !------------------------------------------------------------------------------!

  allocate(irej(total_num_zm),stat=ierr)  
  if (ierr.ne.0) call allerr('irej(total_num_zm)                                  ')
  allocate(iacc(total_num_zm),stat=ierr)  
  if (ierr.ne.0) call allerr('iacc(total_num_zm)                                  ')
  allocate(iup(total_num_zm),stat=ierr)   
  if (ierr.ne.0) call allerr('iup(total_num_zm)                                   ')

  iacc = 0
  iup = 0
  irej = 0
  rdiv_max = real(div_max)

  allocate(tempin(maxval(nintern)),stat=ierr)
  if (ierr.ne.0) call allerr('tempin(maxval(nintern))                             ')
  allocate(tempcarts(3,maxval(num_at)),stat=ierr)
  if (ierr.ne.0) call allerr('tempcarts(3,maxval(num_at))                         ')
  allocate(tweak(7+maxval(nintern)),stat=ierr)
  if (ierr.ne.0) call allerr('tweak(7+maxval(nintern))                            ')
  allocate(invals_ge(maxval(nintern)),stat=ierr)
  if (ierr.ne.0) call allerr('invals_ge(maxval(nintern))                          ')

  mtyp = 0
  do il = 1,num_loc
     do iz = 1,total_num_zm
        do m = 1,maxval(mocc)
           do k = 1,conmol(il,m,iz)
              typ =  ty(il,iz,m,k)
              if(typ.gt.mtyp)mtyp=typ
           end do
        end do
     end do
  end do

  if(mtyp.ne.ntyp) then
     write(stderr,*)'----------------------------------------------------------------------------'
     write(stderr,*)'Number of contact vector types in contact vector file does not'
     write(stderr,*)'match the number of types specified in the MC parameters file.'
     write(stderr,*)'Program exiting. Code NVT1'
     write(stderr,*)'----------------------------------------------------------------------------'
     stop
  end if

  allocate(wrap(-1*maxval(simsize):2*maxval(simsize),3),stat=ierr)
  if (ierr.ne.0) call allerr('wrap(-1*maxval(simsize):2*maxval(simsize),3)        ')
  do i = 1,3
     rsimsize(i)=real(simsize(i))
     do  m=-simsize(i),simsize(i)+simsize(i)
        wrap(m,i)=m
        if (wrap(m,i).le.0) wrap(m,i)=simsize(i)+wrap(m,i)
        if (wrap(m,i).gt.simsize(i)) wrap(m,i)=wrap(m,i)-simsize(i)
     end do
  end do

  mcell   = num_loc*simsize(1)*simsize(2)*simsize(3)
  radjcy  = real(adjcy)
  bradjcy = real(badjcy)
  allocate(reject(total_num_zm,7+maxval(num_at),div_max),stat=ierr)
  if (ierr.ne.0) call allerr('reject(total_num_zm,7+maxval(num_at),div_max)       ')
  reject = 0

  celltrans(1,:) = as_cartesian(xtal,real((/1,0,0/)))
  celltrans(2,:) = as_cartesian(xtal,real((/0,1,0/)))
  celltrans(3,:) = as_cartesian(xtal,real((/0,0,1/)))
  one_on_temp = 1.0/temp

  if (.not. quiet)write(stdout,*)'Initialisation completed.'
  if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

  if (.not. quiet)write(stdout,*)'Header of input file is:'
  header = 'ZMC simulation'
  call get_optional_string(header,'HEADER')

  !------------------------------------------------------------------------------!
  ! MC loop starts here, outer (icy) and inner (ic2)
  !------------------------------------------------------------------------------!

  do icy = 1,ncy

     if(icy.gt.1)then
        if (.not. quiet)write(stdout,'(1x," acc:",i8,", up:",i8,", rej:",i8,", R/T:",f5.2)')  &
             sum(iacc),sum(iup),sum(irej),real(sum(irej)) / real(sum(iacc)+sum(iup)+sum(irej))
        if ((adjcy.gt.0).and.(icy.gt.1)) then
           if (real(icy)/radjcy.eq.real(icy/adjcy)) call incupdate() 
        end if
     end if
     iacc = 0
     iup = 0
     irej = 0
     reject = 0

     if ((badjcy.gt.0).and.(icy.gt.1)) then
        if (real(icy)/bradjcy.eq.real(icy/badjcy)) then
           call  B_adjust(B_aim)
        end if
     end if

     do ic2 = 1, mcell

        !------------------------------------------------------------------------------!
        !        Now we have to identify the molecule to modify
        !------------------------------------------------------------------------------!

        call rannum(rvar,1)
        oa = wrap(nint(rvar(1)*rsimsize(1)),1)
        call rannum(rvar,1)
        ob = wrap(nint(rvar(1)*rsimsize(2)),2)
        call rannum(rvar,1)
        oc = wrap(nint(rvar(1)*rsimsize(3)),3)
        call rannum(rvar,1)
        ol = int(rvar(1)*real(num_loc))+1
        om = mocc(oa,ob,oc,ol)
        oz = zocc(oa,ob,oc,ol)
        allocate(new_coords(3,(num_at(oz))),stat=ierr)
        if (ierr.ne.0) call allerr('new_coords(3,(num_at(oz)))     3   ')

        !------------------------------------------------------------------------------!
        !         Store its variables for temporary backup 
        !------------------------------------------------------------------------------!

        nin = nintern(oz)
        tempin(1:nin)             = incrystal(oa,ob,oc,ol,1:nin)
        tempxyz                   = xyzcrystal(oa,ob,oc,ol,:)
        tempq                     = qcrystal(oa,ob,oc,ol)
        tempcarts(:,1:num_at(oz)) = carts(:,oa,ob,oc,ol,1:num_at(oz))

        !------------------------------------------------------------------------------!
        !         Get its energy
        !------------------------------------------------------------------------------!

        call get_energy(ol,om,oz,oa,ob,oc,old_e)

        !------------------------------------------------------------------------------!
        !         Generate the modification
        !------------------------------------------------------------------------------!

        call gettweak(oz)

        !------------------------------------------------------------------------------!
        !         Apply the modification
        !------------------------------------------------------------------------------!

        qcrystal(oa,ob,oc,ol)        = qcrystal(oa,ob,oc,ol)        + tweak(4:7)
        xyzcrystal(oa,ob,oc,ol,:)    = xyzcrystal(oa,ob,oc,ol,:)    + tweak(1:3)
        incrystal(oa,ob,oc,ol,1:nin) = incrystal(oa,ob,oc,ol,1:nin) + tweak(7+1:7+nin)

        !------------------------------------------------------------------------------!
        !         Get the coordinates of the modified  molecule and
        !         Update (a temporary copy of) the z-matrix
        !------------------------------------------------------------------------------!

        call copy(z1,zmat(oz))
        incrystal(oa,ob,oc,ol,1:nin)=                                             &
             parameter(z1,intern(1:nin,1,oz),intern(1:nin,2,oz),incrystal(oa,ob,oc,ol,1:nin))
        new_coords = as_xyz(z1)
        rmat = as_rotmatrix(qcrystal(oa,ob,oc,ol))
        call rotate(rmat,new_coords)

        do iat = 1, num_at(oz)
           carts(:,oa,ob,oc,ol,iat)= new_coords(:,iat)+xyzcrystal(oa,ob,oc,ol,:)
        end do

        !------------------------------------------------------------------------------!
        !         Get its energy 
        !------------------------------------------------------------------------------!

        call get_energy(ol,om,oz,oa,ob,oc,new_e)

        !------------------------------------------------------------------------------!
        !         Do the MC step
        !------------------------------------------------------------------------------!

        edif = new_e-old_e

        if(edif.lt.0.0) then
           ! Increment the acceptance counter
           iacc(oz) = iacc(oz)+1
        else
           call rannum(rvar,1)
           !!           if(rvar(1).le.exp(-edif*one_on_temp)) then
           if(rvar(1).le.exp(-edif/temp)) then
              ! Conditional acceptance
              iup(oz) = iup(oz)+1
           else
              ! Put it all back the way it was.
              incrystal(oa,ob,oc,ol,1:nin)     = tempin(1:nin) 
              xyzcrystal(oa,ob,oc,ol,:)        = tempxyz
              qcrystal(oa,ob,oc,ol)            = tempq
              carts(:,oa,ob,oc,ol,1:num_at(oz))= tempcarts(:,1:num_at(oz))
              irej(oz) = irej(oz)+1
              do iivar = 1,7
                 ibin  = min(int(abs(tweak(iivar))*2.0*rdiv_max/widths(iivar,oz))+1,div_max)
                 ibin = max(ibin,1)
                 reject(oz,iivar,ibin) = reject(oz,iivar,ibin)+1
              end do
              do iivar = 1,nin
                 ibin  = min(int(abs(tweak(iivar+7))*2.0*rdiv_max/inwidths(iivar,oz))+1,div_max)
                 ibin = max(ibin,1)
                 reject(oz,iivar+7,ibin) = reject(oz,iivar+7,ibin)+1
              end do
           end if
        end if
        deallocate(new_coords)
     end do
     if (.not. quiet)write(stdout,'(1x," %done:",f5.1)')  real(icy)/real(ncy)*100.0
  end do

  !------------------------------------------------------------------------------!
  !        Write out final variable increments ('widths')
  !------------------------------------------------------------------------------!

  if(ncy.gt.0) then
     if (adjcy.gt.0) then
        do iz = 1,total_num_zm
           if (.not. quiet) then
              write(stdout,*)'----------------------------------------------------------------------------'
              write(stdout,'(1x,"Molecule type (z-matrix number): ",i3)')iz
              write(stdout,'(1x,"Widths for increments on x,y,z are       ",3f9.5)') &
                   (widths(i,iz),i=1,3)
              write(stdout,'(1x,"Widths for increments on q1,q2,q3,q4 are ",4f9.5)') &
                   (widths(i,iz),i=4,7)
           end if
           if(nintern(iz).gt.0) then
              do i = 1,nintern(iz)
                 if (.not. quiet) write(stdout,*)'Widths for internal degrees of freedom:'
                 if(intern(i,1,iz).eq.1) then
                    ! It is bond length and leave it alone
                    if (.not. quiet) write(stdout,*)'                                        ',inwidths(i,iz)
                 else
                    ! It is angle and convert to degrees
                    if (.not. quiet) write(stdout,*)'                                        ',inwidths(i,iz)*180.0/3.14159
                 end if
              end do
              write(stdout,*)'For internal d.f. that are angles, widths are expressed in degrees.'
           else
              if (.not. quiet) write(stdout,*)'No internal d.f.'
           end if
        end do
     end if

     if (badjcy.gt.0) then
        if (.not. quiet) write(stdout,*)'----------------------------------------------------------------------------'
        if (.not. quiet)  write(stdout,*)'Globally scaled intermolecular spring constants are:'
        do iz = 1,ntyp
           if (.not. quiet)  write(stdout,*)iz,sprcon(iz)
        end do

        if(ninspr.gt.0) then
           if (.not. quiet) write(stdout,*)'----------------------------------------------------------------------------'
           if (.not. quiet) write(stdout,*)'Globally scaled intramolecular spring constants are:'
           do iz = 1,ninspr
              if (.not. quiet)write(stdout,*)iz,insprcon(iz)
           end do
        end if
        if(maxval(ncross).gt.0) then
           if (.not. quiet) write(stdout,*)'----------------------------------------------------------------------------'
           if (.not. quiet) write(stdout,*)'Globally scaled cross term spring constants are:'
           do iz = 1,maxval(crosstype)
              if (.not. quiet)write(stdout,*)iz,cross_spr(iz)
           end do
        end if
     end if

  !---------------------------------------------------------------------------!
  ! Place the call to modwave subroutine here  - Eric Chan
  ! just before the program decides to exit
  !
  !---------------------------------------------------------------------------!

     if (option_exists('modwave')) then
         call modwave(key2,xtal)
     endif


     !------------------------------------------------------------------------------!
     !        Output the final coordinates
     !------------------------------------------------------------------------------!

     if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
     ! Output to a ZMC crystal file...
     if (option_exists('crystal'))  then
        ! get the name from the value of crystal; if it does not have a value, derive it from outname.
        if(has_value('crystal')) then
           fname2 = ''
           fname2 = get_value('crystal')
           call extension(fname2,'crystal',fname2)
        else
           call extension(outname,'crystal',fname2) 
        end if
        if (.not. quiet)write(stdout,*)'Writing ZMC crystal file to ',trim(fname2)
        call writecrystal(fname2 ,header)
     end if
     ! Output to file for DIFFUSE
     if (option_exists('diffuse'))  then
        if(has_value('diffuse')) then
           fname2 = ''
           fname2 = get_value('diffuse')
           call extension(fname2,'diffuse',fname2)
        else
           call extension(outname,'diffuse',fname2)
        end if
        if (.not. quiet)write(stdout,*)'Writing file for DIFFUSE calculation ',trim(fname2)
        call write_readat(fname2)
     end if
     ! And for discus 
     if (option_exists('discus'))  then
        if(has_value('discus')) then
           fname2 = ''
           fname2 = get_value('discus')
           call extension(fname2,'discus',fname2)
        else
           call extension(outname,'discus',fname2)
        end if
        if (.not. quiet)write(stdout,*)'Writing files for DISCUS calculation '
        call write_discus(fname2,xtal)
     end if
  else
     if (.not. quiet)write(stdout,*)'No output since no MC cycles done.'
  end if

  !------------------------------------------------------------------------------!
  ! When the MC has ended, we need to dip into a routine to do some analysis of it.
  ! There are various subroutines which can be used                                   
  !------------------------------------------------------------------------------!

  if (option_exists('summary'))  then
     if(has_value('summary')) then
        inline  = get_value('summary')
        if(inline.ne.'inline') then
           write(stderr,*)'Invalid value passed into summary.'
           write(stderr,*)'Writing summary to file.'
           call get_info(outname,'y')
        else
           call get_info(outname,'n')
        end if
     else
        call get_info(outname,'y')
     end if
  end if

  if (option_exists('cartsn')) then
     call extension(outname,'cartsn',fname2) 
     call write_xyz(fname2,xtal,'n')
  end if

  if (option_exists('cartst')) then
     call extension(outname,'cartst',fname2) 
     call write_xyz(fname2,xtal,'y')
  end if

  if (option_exists('fracs')) then
     call extension(outname,'fracs',fname2) 
     call write_frac(fname2,xtal)
  end if

  if (option_exists('corro')) call peanut(outname,'o')
  if (option_exists('corra')) call peanut(outname,'a')
  if (option_exists('pairs')) call pairs(outname)
  if (option_exists('cif'))  call calc_uij_new(xtal,outname)
  if (option_exists('energy')) then
     !        Write out molecular energies
     call extension(outname,'energy',fname2)
     if (.not. quiet)write(stdout,*)'Writing out molecular energies to ',fname2
     dunit = open(fname2)
     do oa=1,simsize(1)
        do ob=1,simsize(2)
           do oc=1,simsize(3)
              do ol = 1,num_loc
                 om = mocc(oa,ob,oc,ol)
                 oz = zocc(oa,ob,oc,ol)
                 !         Get its energy                                                      !
                 call get_energy(ol,om,oz,oa,ob,oc,old_e)
                 write(dunit,'(6i6,f11.4)')oa,ob,oc,ol,om,oz,old_e
              end do
           end do
        end do
     end do
     close(dunit)
     if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'
  end if

  if (.not. quiet)write(stdout,*)'Exiting normally.'
  if (.not. quiet)write(stdout,*)'----------------------------------------------------------------------------'

  stop

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------


end program ZMC
