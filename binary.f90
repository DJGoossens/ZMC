module binary

  use precision
  use iso_varying_string
  
  implicit none

  private

  character(len=*), parameter :: version = "$Id: binary.f90,v 1.2 2003/10/14 05:07:44 aidan Exp aidan $"

  !! $Log: binary.f90,v $
  !! Revision 1.2  2003/10/14 05:07:44  aidan
  !! Cleaned up the internal logic -- conversion between the byte format and
  !! integer and raw bytes (integer(kind=1)) is elemental, i.e. input an array,
  !! get an array out. Conversion between byte and character is not elemental.
  !! We convert a 'byte vector' (array of byte type) into a character string
  !! of the same length as the size of the array of bytes and vice versa.
  !!
  !! Revision 1.1  2003/10/08 04:40:56  aidan
  !! Initial revision
  !!
  !!

  !-- For most purposes the components of byte_type should be private.
  !-- The exact components can vary from implementation to implementation,
  !-- so users of this module should avoid referencing the components.
  !-- However, I/O (even unformatted) requires public components.
  type, public :: byte_type
     integer(i1_kind) :: data
  end type
  type(byte_type), public, parameter :: zero_byte = byte_type(0)

  interface assignment(=)
     module procedure assign_byte_from_rawbytes, assign_rawbytes_from_byte, &
          assign_byte_from_int, assign_int_from_byte, &
          assign_bytevector_from_string, assign_string_from_bytevector, & 
          assign_bytevector_from_varstr, assign_varstr_from_bytevector 
  end interface

  interface char
     module procedure byte_to_char, bytevector_to_string 
  end interface

  interface int
     module procedure byte_to_int
  end interface

  interface byte
     module procedure int_to_byte, rawbytes_to_byte, string_to_bytevector !, char_to_byte
  end interface

  interface operator(==)
     module procedure byte_eq_byte
  end interface

  interface operator(/=)
     module procedure byte_neq_byte
  end interface

  ! Public routines
  public :: byte
  
  ! Over-loaded routines
  public :: char, int

  ! Over-loaded operators
  public :: assignment(=), operator(==), operator(/=)

contains

  ! Equivalence

  elemental logical function byte_eq_byte (bytel, byter) result(equal)

    ! Tests to see if the two bytes are equal. Elemental so that it
    ! can be used in by any and all on arrays of bytes

    type(byte_type), intent(in) :: bytel, byter

    equal = (bytel%data == byter%data)

  end function byte_eq_byte

  elemental logical function byte_neq_byte (bytel, byter) result(notequal)

    ! Tests to see if the two bytes are not equal. Elemental so that it
    ! can be used in by any and all on arrays of bytes

    type(byte_type), intent(in) :: bytel, byter

    notequal = (bytel%data /= byter%data)

  end function byte_neq_byte


  ! Conversion functions
  
  elemental function rawbytes_to_byte (rawbytes) result(byteval)

    ! Convert a raw byte to a byte type--just copies it into
    ! the byte_type

    ! Input variable
    integer(kind=1), intent(in) :: rawbytes

    ! Return value
    type(byte_type) :: byteval
    
    byteval%data = rawbytes
    
  end function rawbytes_to_byte
  
  elemental function byte_to_rawbytes (byteval) result(rawbytes)
    
    ! Not used. Will result in an ambiguous function reference if included
    ! in the 'int' interface.

    ! Input variable
    type(byte_type), intent(in) :: byteval

    ! Return value
    integer(kind=1) :: rawbytes
    
    rawbytes = byteval%data
    
  end function byte_to_rawbytes


  pure function byte_to_char (byteval) result(charval)

    ! Convert from a byte to a character. Not elemental because we
    ! want a vector of bytes to transform to a single string, not
    ! an array of single characters

    type(byte_type), intent(in) :: byteval
    character(len=1) :: charval

    charval = transfer(byteval%data,charval)

  end function byte_to_char

  elemental function char_to_byte (charval) result(byteval)

    ! Convert from a character to a byte. Not used. Results in an 
    ! ambiguous function reference. It would be nice if this 
    ! elemental function could zip along an entire string, but a
    ! character string is not an array .. so must define the
    ! function below for this

    character(len=1), intent(in) :: charval

    type(byte_type) :: byteval
    
    byteval%data = transfer(charval,byteval%data)

  end function char_to_byte

  pure function bytevector_to_string (bytevector) result(string)

    ! Convert from an array of bytes to a string

    type(byte_type), dimension(:), intent(in) :: bytevector
    character(len=size(bytevector)) :: string

    string = transfer(bytevector(:)%data,string)

  end function bytevector_to_string

  pure function string_to_bytevector (string) result(bytevector)

    ! Convert from a string to a byte

    character(len=*), intent(in) :: string

    type(byte_type), dimension(len(string)) :: bytevector
    
    integer :: i

    do i = 1,len(string)
       bytevector(i)%data = transfer(string(i:i),bytevector(i)%data)
    end do

  end function string_to_bytevector

  elemental function byte_to_int (byteval) result(intval)

    ! Convert from byte to integer

    type(byte_type), intent(in) :: byteval

    ! Return type
    integer :: intval
    
    intval = int(byteval%data)
    if (intval < 0) intval = intval + 256
    
  end function byte_to_int

  elemental function int_to_byte (intval) result(byteval)

    ! Convert from integer to byte
    integer, intent(in)          :: intval

    ! Return type
    type(byte_type) :: byteval

    ! Local copy
    integer :: intcopy

    intcopy = intval

    if (intcopy > 255) then
       intcopy = -1
    else if (intcopy < 0) then
       intcopy = 0
    else if (intcopy > 128) then
       intcopy = intcopy - 256
    end if
    byteval%data = int(intcopy,1)
    
  end function int_to_byte

  
  ! Assignment routines

  elemental subroutine assign_rawbytes_from_byte (rawbytes, byteval)

    ! rawbytes = byte

    type(byte_type), intent(in)  :: byteval
    integer(kind=1), intent(out) :: rawbytes
    
    rawbytes = byteval%data
    
  end subroutine assign_rawbytes_from_byte

  elemental subroutine assign_byte_from_rawbytes (byteval, rawbytes)

    ! byte = rawbytes

    integer(kind=1), intent(in)  :: rawbytes
    type(byte_type), intent(out) :: byteval
    
    byteval%data = byte(rawbytes)
    
  end subroutine assign_byte_from_rawbytes

  elemental subroutine assign_int_from_byte (intval, byteval)

    ! int = byte

    type(byte_type), intent(in) :: byteval
    integer, intent(out)        :: intval
    
    intval = int(byteval)
    
  end subroutine assign_int_from_byte

  elemental subroutine assign_byte_from_int (byteval, intval)

    ! byte = int

    integer, intent(in)          :: intval
    type(byte_type), intent(out) :: byteval

    byteval = byte(intval)
    
  end subroutine assign_byte_from_int

  subroutine assign_string_from_bytevector (string, bytevector)

    ! string = bytevector

    type(byte_type), dimension(:), intent(in) :: bytevector
    character(len=size(bytevector)), intent(out) :: string

    ! character(len=1), dimension(size(bytevector), intent(out) :: string

    string = char(bytevector)

  end subroutine assign_string_from_bytevector

  subroutine assign_bytevector_from_string (bytevector, string)

    ! bytevector = string

    character(len=*), intent(in)               :: string
    type(byte_type), dimension(len(string)), intent(out) :: bytevector

    bytevector = byte(string)

  end subroutine assign_bytevector_from_string

  subroutine assign_varstr_from_bytevector (varstring, bytevector)

    ! varstring = bytevector

    type(byte_type), dimension(:), intent(in) :: bytevector
    type(varying_string), intent(out)         :: varstring

    varstring = char(bytevector)

  end subroutine assign_varstr_from_bytevector

  subroutine assign_bytevector_from_varstr (bytevector, varstring)

    ! bytevector = varstring

    type(varying_string), intent(in)           :: varstring
    type(byte_type), dimension(:), intent(out) :: bytevector

    bytevector = byte(char(varstring))

  end subroutine assign_bytevector_from_varstr

end module binary
