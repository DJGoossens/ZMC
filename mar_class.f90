module mar_class

  use binary_io
  use file_functions, only: stderr, exists, freeunit
  use variable_array, only: push, splice

  implicit none

  private

  character(len=*), parameter :: version = "$Id: mar_class.f90,v 1.5 2008/05/12 00:22:05 aidan Exp $"

  !! Converts a mar345 file into a pgm file

  !! $Log: mar_class.f90,v $
  !! Revision 1.5  2008/05/12 00:22:05  aidan
  !! Added back the "new" routine to allow mar objects to be made from scratch.
  !!
  !! Revision 1.4  2008/05/08 06:35:00  aidan
  !! Massive update. Changed C routines for reading in mar files (now only supports
  !! mar345, packed (pck) format). Can now also write these files.
  !!
  !! Revision 1.3  2008/04/16 05:09:32  aidan
  !! Deleted unecessary gumpf dealing with other formats supported by the XDS
  !! library.
  !!
  !! Added in a masking routine.
  !!
  !! Revision 1.2  2006/12/06 03:11:21  aidan
  !! Added info routine to output information about te mar object.
  !!
  !! Revision 1.1  2005/07/28 03:19:42  aidan
  !! Initial revision
  !!

  type mar_object 
     private
     character(len=200) :: filename
     integer :: info(15)
     integer, dimension(:,:), pointer :: data
     character(len=3968) :: asciiheader
     logical :: littleendian
  end type mar_object

  interface new
     module procedure new_mar_object
  end interface

  interface read
     module procedure read_mar_object
  end interface

  interface write
     module procedure write_mar_object
  end interface

  interface destroy
     module procedure destroy_mar_object
  end interface

  interface get_data
     module procedure get_mar_data
  end interface

  interface mask_data
     module procedure mask_mar_data
  end interface

  interface size
     module procedure get_mar_size
  end interface

  interface info
     module procedure print_mar_info
  end interface

  interface operator (==)
     module procedure mar_eq_mar
  end interface

  interface operator (/=)
     module procedure mar_neq_mar
  end interface

  public :: mar_object
  public :: new, read, write, get_data, mask_data, size, destroy, info
  public :: operator(==), operator(/=)

contains

  subroutine destroy_mar_object (mar, status)

    type (mar_object), intent(inout) :: mar

    integer, optional, intent(inout) :: status

    integer :: error

    if (associated(mar%data)) deallocate(mar%data, stat=error)

    mar%filename = ' '
    mar%info = 0
    mar%littleendian = .FALSE.

    if (present(status)) then
       status = 0
       if (error > 0) status = error
    else
       if (error > 0) stop "MAR_CLASS :: Error deallocating memory for mar object"
    end if

  end subroutine destroy_mar_object

  subroutine new_mar_object (mar, data, pixlength, pixheight, wavelength, distance, phistart, phiend, header, status)

    ! Interface variables
    type (mar_object), intent(inout) :: mar
    integer, intent(in)              :: data(:,:)
    real, optional, intent(in)       :: pixlength, pixheight, wavelength, distance, phistart, phiend
    integer, optional, intent(inout) :: status
    character(len=*), optional, intent(in) :: header

    ! Local variables
    integer :: error, nx, numhigh
    logical :: make_header
    real :: lambda, d, phis, phie
    character(len=64) :: buf

    nx = size(data,1)
    numhigh = count(data>65535)

    ! if (associated(mar%data)) deallocate(mar%data, stat=error)

    call destroy(mar)
    allocate(mar%data(nx,nx), stat=error)

    mar%data = data

    ! Initialise all the info fields to a bogus value
    mar%info = 9999999

    mar%info(1) = nx
    mar%info(2) = numhigh
    mar%info(3) = 0         ! Compressed (pck) is the only format we support
    mar%info(4) = 1         ! Collection mode is set to TIME, we don't really care anyway
    mar%info(5) = nx**2   ! Total pixels
    if (present(pixlength))  mar%info(6) = pixlength*1000
    if (present(pixheight))  mar%info(7) = pixheight*1000
    if (present(wavelength)) mar%info(8) = wavelength*1.0e6
    if (present(distance))   mar%info(9) = distance*1000
    if (present(phistart))   mar%info(10) = phistart*1000
    if (present(phiend))     mar%info(11) = phiend*1000

    ! mar%filename = ' '
    mar%littleendian = .TRUE.

    if (present(header)) then
       mar%asciiheader = header
    else
       ! Write the ascii part of the header. 
       mar%asciiheader = add_to_header('mar research')
       mar%asciiheader = trim(mar%asciiheader) // add_to_header('PROGRAM MAR_CLASS')
       write(buf,'(A,I0,A,I0)') 'FORMAT ',nx,' PCK ',mar%info(5)
       mar%asciiheader = trim(mar%asciiheader) // add_to_header(buf)
       write(buf,'(A,I0)') 'HIGH ',numhigh
       mar%asciiheader = trim(mar%asciiheader) // add_to_header(buf)
       write(buf,'(A,F0.3,A,F0.3,A)') 'PHI START ',real(mar%info(10))/1000.,' END ',real(mar%info(11))/1000.,' OSC ?'
       mar%asciiheader = trim(mar%asciiheader) // add_to_header(buf)
       write(buf,'(A,F0.3,A,F0.3)') 'PIXEL LENGTH ',real(mar%info(6))/1000.,' HEIGHT ',real(mar%info(7))/1000.
       mar%asciiheader = trim(mar%asciiheader) // add_to_header(buf)
       write(buf,'(A,F0.6)') 'WAVELENGTH ',real(mar%info(8))/1.0e6
       mar%asciiheader = trim(mar%asciiheader) // add_to_header(buf)
       write(buf,'(A,F0.3)') 'DISTANCE ',real(mar%info(9))/1000.
       mar%asciiheader = trim(mar%asciiheader) // add_to_header(buf)
       ! The mar programs segfaulted if I didn't have this 'END OF HEADER' line
       mar%asciiheader = trim(mar%asciiheader) // add_to_header('END OF HEADER')
       
       ! Pad the rest of the header with space
       mar%asciiheader(len_trim(mar%asciiheader)+1:) = ' '
    end if

    if (present(status)) then
       status = 0
       if (error > 0) status = error
    else
       if (error > 0) stop "MAR_CLASS :: Error deallocating memory for mar object"
    end if

  end subroutine new_mar_object

  integer function get_mar_size (mar) result(size)

    type (mar_object), intent(in) :: mar

    size = mar%info(1)

  end function get_mar_size

  function get_mar_data (mar) result(data)

    type (mar_object), intent(in) :: mar

    integer, dimension(mar%info(1), mar%info(1)) :: data

    data = mar%data

  end function get_mar_data

  subroutine read_header (file, mar)

    character(len=*), intent(in)     :: file
    type (mar_object), intent(inout) :: mar

    ! local variables
    integer :: flag, i
    type (binary_filehandle) :: fh

    mar%littleendian = .TRUE.

    ! Open the mar file
    call open(fh, file)

    ! Read 4 bytes into integer variable flag
    call read(fh, flag, 4, littleendian=mar%littleendian)

    ! This should be the number 1234, if not, we'll swap our endianess
    ! and try again
    if (flag /= 1234) then

       call close(fh)

       mar%littleendian = .NOT. mar%littleendian

       call open(fh, file)

       ! Read 4 bytes into integer variable flag
       call read(fh, flag, 4, littleendian=mar%littleendian)

       ! Can't find 1234 reading big or little endian -- time to give up
       if (flag /= 1234) then
          write(stderr,*) 'MAR_CLASS :: Not a mar345 file: ',file
          stop
       end if

    end if

    ! Read in header information 
    do i = 1, size(mar%info)
       call read(fh, mar%info(i), 4, littleendian=mar%littleendian)
    end do

    ! Read in 16 more bytes of nothing
    do i = 1, 16
       call read(fh, flag, 4, littleendian=mar%littleendian)
    end do

    call read(fh, mar%asciiheader)

    call close(fh)

  end subroutine read_header

  subroutine write_header (fh, mar)

    ! Write the mar-header to the file. The description of the mar345
    ! header from the man page is:

    !  The image header consists of "lines" of 64 bytes each with exception of
    !  the first line, that contains 16 4-byte words (32  bit  integers).  The
    !  very  first  value  always  is 1234 for all scanning modes, i.e.  it is
    !  just a marker for the necessity of byte swapping and does NOT give  the
    !  size of the image (in contrast to the first value of the old image for-
    !  mats which is 1200 or 2000 or its byte swapped  equivalent).   Programs
    !  of the mar345 suite rely on this value to decide wether the image is in
    !  old format or in the new one and/or if byte-swapping is required.   The
    !  remaining 15 32-bit integers are:
    !  2) Size of image in one dimension
    !  3) Number of high intensity pixels
    !  4) Image format (1=COMPRESSED, 2=SPIRAL)
    !  5) Collection mode (0=DOSE, 1=TIME)
    !  6) Total number of pixels in image
    !  7) Pixel length (in mm * 1.000)
    !  8) Pixel height (in mm * 1.000)
    !  9) Used wavelength (in Ang * 1.000.000 )
    !  10) Used distance (in mm * 1.000 )
    !  11) Used starting PHI (in deg. * 1.000 )
    !  12) Used ending   PHI (in deg. * 1.000 )
    !  13) Used starting OMEGA (in deg. * 1.000 )
    !  14) Used ending   OMEGA (in deg. * 1.000 )
    !  15) Used CHI (in deg. * 1.000 )
    !  16) Used TWOTHETA (in deg. * 1.000 )
    !
    !  If  the first value is not 1234, the bytes of the following 15 integers
    !  must be swapped.
    !
    !  The next 64 character line contains a  general  identifier  string  for
    !  this  type  of file: "mar research" (bytes 65 to 76 in the image file).
    !  This is for possible use with the "file" command under Unix.
    !
    !  All following lines contain keyworded  information.  The  last  keyword
    !  should  be "END OF HEADER". All keywords are in capital letters and all
    !  "lines" are pure ASCII, so they are not affected by the byte  order  of
    !  different  computer  platforms.  Processing  of  the  keywords  is  not
    !  required. For using the formats correctly, the most important  informa-
    !  tion  is  contained  in  the  second (bytes 5-8) and third (bytes 9-12)
    !  header value: the size of the image and the number  of  high  intensity
    !  pixels!

    ! Interface variables
    type (binary_filehandle), intent(inout) :: fh
    type (mar_object), intent(inout)        :: mar

    ! Local variables
    integer :: unit_num, flag, i
    character(len=3968) :: asciiheader

    ! Write integer variable flag
    call write(fh, 1234, 4, littleendian=mar%littleendian)

    ! Write binary file info
    do i = 1, size(mar%info)
       call write(fh, mar%info(i), 4, littleendian=mar%littleendian)
    end do

    ! Despite the docs saying that the first line of a mar file
    ! contains 16 4-byte words I am padding it by 16 more bytes
    ! as this is what the mar files from the detector look like.
    do i = 1, 16
       call write(fh, 0, 4, littleendian=mar%littleendian)
    end do

    ! Write the minimum into the ascii part of the header. The
    ! mar programs segfaulted if I didn't have the 'END OF HEADER'
    ! line
    ! asciiheader = add_to_header('mar research')
    ! asciiheader = trim(asciiheader) // add_to_header('PROGRAM MAR_CLASS')
    ! asciiheader = trim(asciiheader) // add_to_header('END OF HEADER')

    ! Pad the rest of the header with space
    ! asciiheader(len_trim(asciiheader)+1:) = ' '

    call write(fh, mar%asciiheader)

  end subroutine write_header

  function add_to_header(text) result(line)

    ! A simple function to pad out or truncate text to fit in
    ! the 64 byte limit of a mar header file

    ! Interface variables
    character(len=*), intent(in)    :: text

    ! Return value
    character(len=64) :: line

    line = ' '
    line(1:min(63,len(text))) = text
    line(64:64) = achar(10)

  end function add_to_header
  
  subroutine read_mar_object (mar, file)

    ! Read a mar object from a file

    ! Interface variables
    type (mar_object), intent(inout) :: mar
    character(len=*), intent(in)     :: file

    ! Local variables
    integer, dimension(:), allocatable :: tmpdata
    type (binary_filehandle) :: fh
    integer :: i, header(1024), location, intensity, status

    ! Make sure the mar object is empty
    call destroy(mar, status)

    if (status > 0) write(stderr,*) 'MAR_CLASS :: Possible memory leak!'

    ! Read in the header so we know how big it is
    call read_header(file, mar)

    ! Allocate the space required for the data
    allocate(tmpdata(mar%info(1)*mar%info(1)))
    allocate(mar%data(mar%info(1),mar%info(1)))

    mar%filename = file

    ! Call the C routine to read in the packed data first
    call readpack( tmpdata, trim(file), len_trim(file))

    ! Check to see if we have any overflow records
    if (mar%info(2) > 0) then

       ! Now read in overflow records -- these replace the
       ! values we read in previously

       ! Open the mar file
       call open(fh, file)
       
       ! Position data file at overflow records = read 
       ! 4096 bytes into integer variable flag
       call read(fh, header, 4, littleendian=mar%littleendian)
       
       do i = 1, mar%info(2)
          ! Read location 
          call read(fh, location, 4, littleendian=mar%littleendian)
          ! Read intensity
          call read(fh, intensity, 4, littleendian=mar%littleendian)
          ! Set this location to this value (note the +1 accounts for
          ! fortran 1 based indexing, but the mar file assumed c-style
          ! 0 indexing)
          tmpdata(location+1) = intensity
       end do

       call close(fh)

    end if

    ! The C routine returns the data as a one dimensional array,
    ! so we make it into a 2-D matrix at this point
    mar%data = reshape(tmpdata, (/ mar%info(1), mar%info(1) /))

    deallocate(tmpdata)

  end subroutine read_mar_object

  subroutine write_mar_object (mar, file)

    type (mar_object), intent(inout) :: mar
    character(len=*), intent(in), optional :: file

    integer, dimension(size(mar%data)) :: tmpdata, indices
    type (binary_filehandle) :: fh
    character(len=1000) :: filename
    integer :: i, nov, header(1024)
    integer, pointer :: locations(:)

    if (present(file)) then
       filename = file
    else
       filename = mar%filename
    end if

    ! Open the mar file
    call open(fh, trim(filename), action='write')

    ! Write the header
    call write_header(fh, mar)

    ! Make a one-dimensional version of the data
    ! tmpdata = pack(mar%data, mask=.true.)
    tmpdata = reshape(mar%data, (/ size(mar%data) /))

    ! Find all the pixels greater than 65535
    ! nov = push(locations, pack((/(i,i=1,(tmpdata))/),mask=tmpdata>65535)) -> TOO SLOW!
    do i = 1, size(indices)
       indices(i) = i
    end do
    nov = push(locations, pack(indices,mask=tmpdata>65535))

    if (nov /= mar%info(2)) then
       write(stderr,'(A,I0,A,I0,A)') 'WARNING: number of overflow records (',nov,') does not agree with mar345 header (',mar%info(2),')'
    end if

    do i = 1, nov 
       ! Write location (offset by -1 to account for C based indexing)
       call write(fh, locations(i)-1, 4, littleendian=mar%littleendian)
       ! Write intensity
       call write(fh, tmpdata(locations(i)), 4, littleendian=mar%littleendian)
    end do
    if (mod(nov,8) /= 0) then
       call write(fh, (/(0,i=1,mod(nov,8))/), 4, littleendian=mar%littleendian)
    end if

    ! Close now, as the C routine below opens the file in append mode
    call close(fh)

    call writepack( tmpdata, size(mar%data,1), size(mar%data,2), trim(file), len_trim(file))

    nov = splice(locations, 0)

  end subroutine write_mar_object

  subroutine print_mar_info (mar)

    type (mar_object), intent(in) :: mar

    character(len=*), parameter :: image_format(2) = (/ 'COMPRESSED', 'SPIRAL    ' /)
    character(len=*), parameter :: collection_mode(2) = (/ 'DOSE', 'TIME' /)

    ! print *,'marinfo: ',mar%info

    print *,'Size of image:         ',mar%info(1),' x ',mar%info(1)
    print *,'Number of overflows:   ',mar%info(2)
    print *,'Image format:          ',image_format(mar%info(3))
    print *,'Collection mode:       ',collection_mode(mar%info(4))
    print *,'Total pixels:          ',mar%info(5)
    print *,'Pixels size lxh (micometres):  ',mar%info(6)/1000.,'x',mar%info(7)/1000.
    print *,'Wavelength (Angstroem):',mar%info(8)/1.0e6
    print *,'Sample/detector distance (mm):',mar%info(9)/1000.
    print *,'Phi range (deg):              ',mar%info(10)/1000.,'-',mar%info(11)/1000.

  end subroutine print_mar_info

  subroutine mask_mar_data (mar, maskfile, maskvalue)

    ! Interface variables
    type (mar_object), intent(inout) :: mar
    character(len=*), intent(in), optional :: maskfile
    integer, intent(in), optional :: maskvalue

    ! Local variables
    character(len=1000) :: fname
    integer :: i, ix, iy, nx, unitnum, num_masked_pixels, value
    integer, allocatable :: indices(:)

    if (present(maskfile)) then
       fname = maskfile
    else
       fname = mar%filename
       i = index(fname,".mar",back=.true.)
       if (i == 0) i = len_trim(fname) + 1
       fname(i:) = ".mask"
    end if

    if (exists(fname)) then

       unitnum = freeunit()

       open (unit=unitnum, file=trim(fname), form='unformatted', status='old')
       read(unitnum) num_masked_pixels
       allocate(indices(num_masked_pixels))
       read(unitnum) indices
       close (unitnum)

       value = 0
       if (present(maskvalue)) value = maskvalue

       nx = size(mar)

       ! Now we mask off the data ...
       do i = 1, num_masked_pixels
          ix = mod(indices(i),nx)
          iy = (indices(i)/nx) + 1
          mar%data(ix, iy) = value
       end do

       deallocate(indices)
    end if

  end subroutine mask_mar_data

  logical function mar_eq_mar (mar1, mar2) 

    ! Interface variables
    type (mar_object), intent(in) :: mar1, mar2

    mar_eq_mar = .false.

    if (any(mar1%info /= mar2%info)) then
       !print *,(mar1%info /= mar2%info)
       return
    end if
    if (any(mar1%data /= mar2%data)) then
       print *,count(mar1%data /= mar2%data)
       print *,count(mar1%data == mar2%data)
       return
    end if

    mar_eq_mar = .true.
    
  end function mar_eq_mar

  logical function mar_neq_mar (mar1, mar2) 

    ! Interface variables
    type (mar_object), intent(in) :: mar1, mar2

    mar_neq_mar = .not. (mar1 == mar2)
    
  end function mar_neq_mar

end module mar_class
