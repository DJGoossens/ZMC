module crystallography_class

  use fundamental_constants, only: radian
  use mol2_class, only: mol2_object, cell, group, setting
  use file_functions, only: stderr

  implicit none

  private

  character(len=*), parameter :: version = "$Id: crystallography_class.f90,v 1.1 2003/08/21 06:48:59 aidan Exp $" 

  !! $Log: crystallography_class.f90,v $
  !! Revision 1.1  2003/08/21 06:48:59  aidan
  !! Initial revision
  !!

  type cell_object
     private
     ! The cell parameters a, b, c, alpha, beta, gamma and their 
     ! reciprocal counterparts, a*, b*, c*, alpha*, beta*, gamma*
     real, dimension(3) :: length, angle, reciprocal_length, reciprocal_angle
     ! Cell volume
     real :: volume
  end type cell_object

  type crystallography_object
     private
     ! The space group and setting (integer values as defined in the Intl Tables).
     integer            :: space_group, setting
     type (cell_object) :: cell
     ! The transformation matrices for frac -> orthogonal and vice versa
     real, dimension(3,3) :: orth_matrix, frac_matrix
  end type crystallography_object

  ! Define a tolerance for equivalence and "effectively zero"
  real :: tolerance = 1e-5 

  ! Procedure to make a new crystallography object
  interface new
     module procedure new_cryst, new_cryst_from_array
  end interface

  ! Procedure to make a crystallography object from another object
  interface as_crystallography
     module procedure mol2_to_crystallography
  end interface

  ! Procedure to return cell volume
  interface volume
     module procedure cellvolume_from_xtal, cellvolume
  end interface

  interface as_cartesian
     module procedure fractional_to_cartesian, fractional_array_to_cartesian 
  end interface

  interface as_fractional
     module procedure cartesian_to_fractional, cartesian_array_to_fractional 
  end interface

  interface axes
     module procedure get_cell_axes
  end interface

  interface reciprocal_axes
     module procedure get_reciprocal_cell_axes
  end interface

  interface angles
     module procedure get_cell_angles
  end interface

  interface reciprocal_angles
     module procedure get_reciprocal_cell_angles
  end interface

  interface group
     module procedure get_space_group
  end interface

  interface setting
     module procedure get_setting
  end interface

  ! Define equivalence operator (also defines the .eq. operator by default)
  interface operator (==)
     module procedure xtal_eq_xtal
  end interface

  ! Define non-equivalence operator (also defines the .neq. operator by default)
  interface operator (/=)
     module procedure xtal_neq_xtal
  end interface

  ! Define assignment
  interface assignment (=)
     module procedure assign_array, assign_mol2
  end interface

  ! Public data types
  public :: crystallography_object, cell_object

  ! Public routines
  public :: new, as_crystallography, volume, as_cartesian, as_fractional
  public :: group, setting, axes, reciprocal_axes, angles, reciprocal_angles

  ! Overloaded operators
  public :: assignment(=), operator(==), operator(/=)

contains

  subroutine new_cryst(xtal_obj, a, b, c, alpha, beta, gamma, group, setting)

    ! Input variables
    type(crystallography_object), intent(out) :: xtal_obj

    ! Cell parameters, in angstroms and degrees respectively
    real, intent(in) :: a, b, c, alpha, beta, gamma

    ! Space group and setting as per the International tables
    integer, optional, intent(in) :: group, setting

    xtal_obj%cell%length = (/ a, b ,c /)
    xtal_obj%cell%angle  = (/ alpha, beta, gamma /) * (1.0d0/radian)

    ! Initialise space group and setting to a 'null' value
    xtal_obj%space_group = -1 ; xtal_obj%setting = -1
    if (present(group)) then
       xtal_obj%space_group = group
    end if
    if (present(setting)) then
       xtal_obj%setting = setting
    end if

    call calculate_internal_parameters(xtal_obj)

  end subroutine new_cryst

  subroutine new_cryst_from_array(xtal_obj, cell_array, group, setting)

    ! This initialises a crystallography object but uses an array
    ! of length 6 instead of indvidually naming them. This is not
    ! a wrapper to new_cryst, but a duplicate, as it would have
    ! been painful to deal with the optional arguments

    ! Input variables
    type(crystallography_object), intent(out) :: xtal_obj

    ! Cell parameters, in angstroms and degrees respectively
    real, intent(in) :: cell_array(6)

    ! Space group and setting as per the International tables
    integer, optional, intent(in) :: group, setting

    xtal_obj%cell%length = cell_array(1:3)
    xtal_obj%cell%angle  = cell_array(4:6) * (1.0d0/radian)

    ! Initialise space group and setting to a 'null' value
    xtal_obj%space_group = -1 ; xtal_obj%setting = -1
    if (present(group)) then
       xtal_obj%space_group = group
    end if
    if (present(setting)) then
       xtal_obj%setting = setting
    end if

    call calculate_internal_parameters(xtal_obj)

  end subroutine new_cryst_from_array

  !
  ! ASSIGNMENT
  !

  subroutine assign_array(xtal_obj,array)

    ! Assign to crystallography object the elements in a cell array

    type(crystallography_object), intent(out) :: xtal_obj
    real, intent(in) :: array(6)

    ! Make a new crystallography object
    call new(xtal_obj,array)

  end subroutine assign_array

  subroutine assign_mol2(xtal_obj,mol2)

    ! Assign to crystallography object the elements in a cell array

    type(crystallography_object), intent(out) :: xtal_obj
    type (mol2_object), intent(in) :: mol2

    ! Coerce a mol2 object into a new crystallography object
    xtal_obj = as_crystallography(mol2)

  end subroutine assign_mol2



  subroutine calculate_internal_parameters(xtal_obj)

    type(crystallography_object), intent(inout) :: xtal_obj

    ! Calculate the cell volume
    xtal_obj%cell%volume = volume(xtal_obj%cell%length, xtal_obj%cell%angle)

    ! Calculate the reciprocal lattice parameters
    xtal_obj%cell%reciprocal_length = xtal_obj%cell%length 
    xtal_obj%cell%reciprocal_angle  = xtal_obj%cell%angle
    call make_reciprocal(xtal_obj%cell%reciprocal_length,xtal_obj%cell%reciprocal_angle,xtal_obj%cell%volume)

    ! Now calculate the transformation matrices for 
    ! fractional -> orthogonal and vice versa
    call make_trans_matrices(xtal_obj)

  end subroutine calculate_internal_parameters

  function mol2_to_crystallography(mol2) result(xtal_obj)

    ! Extracts all the necessaries out of the mol2 object and makes a 
    ! crystallography object.
    
    ! Return type of function
    type(crystallography_object) :: xtal_obj

    ! Input variable
    type (mol2_object) :: mol2

    call new(xtal_obj, cell(mol2), group = group(mol2), setting = setting(mol2))

  end function mol2_to_crystallography

  real function cellvolume_from_xtal(xtal_obj) result(vol)

    ! Return the cell volume from a crystallography object.

    ! Input type
    type(crystallography_object), intent(in) :: xtal_obj

    vol = xtal_obj%cell%volume 

  end function cellvolume_from_xtal

  function cellvolume(length,angle)

    ! Calculate the cell volume, V, given by
    !
    ! V = a * b * c * sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2
    !                           + 2 * cos(alpha) * cos(beta) * cos(gamma))

    ! Return type
    real :: cellvolume

    ! Axis lengths and cell angles input as two length 3 arrays
    real, dimension(3) :: length, angle

    real(kind=8) :: check

    check = 1.0d0 + 2.0d0 * product(cos(real(angle,8))) - sum((cos(real(angle,8)))**2)

    ! Check for wonky angles
    if (check < 0.0d0) then
       write(stderr,*) 'CRYSTALLOGRAPHY_CLASS :: cellvolume : Error in unit cell angles ',angle
       stop
    end if

    cellvolume = real(product(real(length,8)) * sqrt(check)) 

    ! Check for sensible answer
    if (cellvolume < 0.0d0) then
       write(stderr,*) 'CRYSTALLOGRAPHY_CLASS :: cellvolume : Volume < 0 ',cellvolume
       stop
    end if
    
  end function cellvolume

  subroutine make_reciprocal(length,angle,volume)

    ! Calculate the reciprocal cell axes (see p64 of 'Fundamentals of
    ! Crystallography' by Giacovazzo et al for details):
    !
    !
    !      a* = b c sin(gamma)
    !           -------------
    !                V
    !
    !      b* = a c sin(beta)
    !           -------------
    !                V
    !
    !      c* = a b sin(alpha)
    !           -------------
    !                V
    !
    !
    !      cos(alpha*) = cos(beta) cos(gamma) - cos(alpha)
    !                    --------------------------------
    !                          sin(beta) sin(gamma)
    !
    !       cos(beta*) = cos(alpha) cos(gamma) - cos(beta)
    !                    --------------------------------
    !                          sin(alpha) sin(gamma)
    !
    !      cos(gamma*) = cos(alpha) cos(beta) - cos(gamma)
    !                    --------------------------------
    !                          sin(alpha) sin(beta)
    !

    ! Axis lengths and cell angles input as two length 3 arrays
    real, dimension(3), intent(inout) :: length, angle
    real, intent(in) :: volume

    ! Local variables

    ! Need some local versions of the axis lengths and angles
    real, dimension(3) :: oldlength, oldangle

    integer :: i

    oldlength = length
    oldangle = angle

    do i=1,3
       length(i) = oldlength(mod(i,3)+1) * oldlength(mod(i+1,3)+1) * sin(oldangle(i))
    end do
    length = length / volume

    do i=1,3
       angle(i) = (cos(oldangle(mod(i,3)+1)) * cos(oldangle(mod(i+1,3)+1)) - cos(oldangle(i))) / &
            (sin(oldangle(mod(i,3)+1)) * sin(oldangle(mod(i+1,3)+1)))
    end do
    angle = acos(angle)
    ! print *
    ! print *,oldlength, length
    ! print *,oldangle*radian, angle*radian

  end subroutine make_reciprocal

  subroutine make_trans_matrices(xtal_obj)

    ! Calculate the fractionalisation and orthogonalisation matrices for a 
    ! crystallography object.  
    ! 
    ! We use the convention that, given the crystallographic basis, A = {a,b,c},
    ! and the cartesian basis, E = {i,j,k}, then
    !
    !                       i || a
    !                       j is in (a,b) plane
    !                       k = i x j
    !
    ! The orthogonalisation matrix, M, is given on p68 of Giacovazzo et al:
    !
    ! 
    !         /                                                       \
    !        |           1/a                      0               0    |
    !        |                                                         |
    !    M = | -cos(gamma)/(a sin(gamma))  1/(b sin(gamma))       0    |
    !        |                                                         |
    !        |      a* cos(beta*)           b* cos(alpha*)        c*   |
    !         \                                                       /
    !
    ! and it's inverse (fractionalisation matrix):
    ! 
    !         /                                                       \
    !        |            a                       0               0    |
    !        |                                                         |
    ! M^-1 = |       b cos(gamma)           b sin(gamma))         0    |
    !        |                                                         |
    !        |       c cos(beta)      -c sin(beta) cos(alpha*)   1/c*  |
    !         \                                                       /
    !

    ! Input type
    type(crystallography_object), intent(inout) :: xtal_obj

    integer :: i

    real(kind=8) :: RR(3,3,6)

    real (kind=8) :: A,ALPH,AS,B,BET,BS,C,COSA,COSAS,COSB,COSBS,COSG,COSGS,CS,GAMM
    real (kind=8) :: SINA,SINAS,SINB,SINBS,SING,SINGS
    INTEGER J,K,N,NCODE

    ! xtal_obj%cell%orth_matrix(1,1) = 1. / xtal_obj%cell%length(1)  
    ! xtal_obj%cell%orth_matrix(2,1) = -cos(gamma)/(a sin(gamma))  
    ! xtal_obj%cell%orth_matrix(3,1) = 

!!$    xtal_obj%orth_matrix = reshape( &
!!$         (/ 1/xtal_obj%cell%length(1), 0., 0., &
!!$          -cos(xtal_obj%cell%angle(3))/(xtal_obj%cell%length(1)*sin(xtal_obj%cell%angle(3) )), &
!!$          1/( xtal_obj%cell%length(2)*sin(xtal_obj%cell%angle(3))), 0., &
!!$          xtal_obj%cell%reciprocal_length(1)*cos(xtal_obj%cell%reciprocal_angle(2)), &
!!$          xtal_obj%cell%reciprocal_length(2)*cos(xtal_obj%cell%reciprocal_angle(1)), &
!!$          xtal_obj%cell%reciprocal_length(3) /), (/3,3/), order = (/2,1/))

    ALPH = xtal_obj%cell%angle(1)
    BET  = xtal_obj%cell%angle(2)
    GAMM = xtal_obj%cell%angle(3)

    SINA = SIN(ALPH)
    COSA = COS(ALPH)
    SINB = SIN(BET)
    COSB = COS(BET)
    SING = SIN(GAMM)
    COSG = COS(GAMM)
    COSAS = cos(xtal_obj%cell%reciprocal_angle(1))   !(COSG*COSB-COSA)/ (SINB*SING)
    SINAS = sin(xtal_obj%cell%reciprocal_angle(1))   !SQRT(1.0-COSAS*COSAS)
    COSBS = cos(xtal_obj%cell%reciprocal_angle(2))   !(COSA*COSG-COSB)/ (SINA*SING)
    SINBS = sin(xtal_obj%cell%reciprocal_angle(2))   !SQRT(1.0-COSBS*COSBS)
    COSGS = cos(xtal_obj%cell%reciprocal_angle(3))   !(COSA*COSB-COSG)/ (SINA*SINB)
    SINGS = sin(xtal_obj%cell%reciprocal_angle(3))   !SQRT(1.0-COSGS*COSGS)

    A = xtal_obj%cell%length(1)
    B = xtal_obj%cell%length(2)
    C = xtal_obj%cell%length(3)

    AS = xtal_obj%cell%reciprocal_length(1)
    BS = xtal_obj%cell%reciprocal_length(2)
    CS = xtal_obj%cell%reciprocal_length(3)

    RR = 0.0
!
!---- Calculate matrices
!
!---- XO along a  Zo along c*
!
    NCODE = 1
    RR(1,1,NCODE) = A
    RR(1,2,NCODE) = B*COSG
    RR(1,3,NCODE) = C*COSB
    RR(2,2,NCODE) = B*SING
    RR(2,3,NCODE) = -C*SINB*COSAS
    RR(3,3,NCODE) = C*SINB*SINAS
!
!---- XO along b  Zo along a*
!
    NCODE = 2
    RR(1,1,NCODE) = A*COSG
    RR(1,2,NCODE) = B
    RR(1,3,NCODE) = C*COSA
    RR(2,1,NCODE) = -A*SING*COSBS
    RR(2,3,NCODE) = C*SINA
    RR(3,1,NCODE) = A*SING*SINBS
!
!---- XO along c  Zo along b*
!
    NCODE = 3
    RR(1,1,NCODE) = A*COSB
    RR(1,2,NCODE) = B*COSA
    RR(1,3,NCODE) = C
    RR(2,1,NCODE) = A*SINB
    RR(2,2,NCODE) = -B*SINA*COSGS
    RR(3,2,NCODE) = B*SINA*SINGS
!
!---- trigonal only - XO along a+b  YO alon a-b  Zo along c*
!
    NCODE = 4
    RR(1,1,NCODE) = A/2.0
    RR(1,2,NCODE) = A/2.0
    RR(2,1,NCODE) = -A*SING
    RR(2,2,NCODE) = A*SING
    RR(3,3,NCODE) = C
!
!---- XO along a*   ZO along c
!
    NCODE = 5
    RR(1,1,NCODE) = A*SINB*SINGS
    RR(2,1,NCODE) = -A*SINB*COSGS
    RR(2,2,NCODE) = B*SINA
    RR(3,1,NCODE) = A*COSB
    RR(3,2,NCODE) = B*COSA
    RR(3,3,NCODE) = C
!
!---- Grr*! to  Gerard Bricogne - his setting for P1 in SKEW.
!     XO along a  Yo along b*
!
    NCODE = 6
    RR(1,1,NCODE) = A
    RR(1,2,NCODE) = B*COSG
    RR(1,3,NCODE) = C*COSB
    RR(2,2,NCODE) = B*SING*SINAS
    RR(3,2,NCODE) = -B*SING*COSAS
    RR(3,3,NCODE) = C*SINB
!
    xtal_obj%orth_matrix(:,:) = RR(:,:,1)
    ! do i=1,3
    !    print *,'orth',xtal_obj%orth_matrix(i,:)
    ! end do

!!$    xtal_obj%frac_matrix = reshape( &
!!$         (/ xtal_obj%cell%length(1), 0., 0., &
!!$          xtal_obj%cell%length(2)*cos(xtal_obj%cell%angle(3)), &
!!$          xtal_obj%cell%length(2)*sin(xtal_obj%cell%angle(3)), 0., &
!!$          xtal_obj%cell%length(3)*cos(xtal_obj%cell%angle(2)), &
!!$          -xtal_obj%cell%length(3)*sin(xtal_obj%cell%angle(2))*cos(xtal_obj%cell%reciprocal_angle(1)), &
!!$          1/xtal_obj%cell%reciprocal_length(3) /), (/3,3/), order = (/2,1/))

    xtal_obj%frac_matrix = invert_matrix(xtal_obj%orth_matrix)
    ! xtal_obj%frac_matrix = xtal_obj%orth_matrix
    ! call gauss(xtal_obj%frac_matrix,3)
    ! do i=1,3
    !    print *,'frac',xtal_obj%frac_matrix(i,:)
    ! end do

  end subroutine make_trans_matrices

  function invert_matrix(M1)

    ! Nicked from spicelib
 
    ! Generate the inverse of a 3x3 matrix.

    real, dimension(3,3) :: M1, invert_matrix
 
    !      M1         I   Matrix to be inverted.
    !      MOUT       O   Inverted matrix (M1)**-1.  If M1 is singular, then
    !                     MOUT will be the zero matrix.   MOUT can
    !                     overwrite M1.

    !      M1    An arbitrary 3x3 matrix.  The limits on the size of
    !            elements of M1 are determined by the process of calculating
    !            the cofactors of each element of the matrix.  For a 3x3
    !            matrix this amounts to the differencing of two terms, each
    !            of which consists of the multiplication of two matrix
    !            elements.  This multiplication must not exceed the range
    !            of double precision numbers or else an overflow error will
    !            occur.
    !
    !      MOUT  is the inverse of M1 and is calculated explicitly using
    !            the matrix of cofactors.  MOUT is set to be the zero matrix
    !            if M1 is singular.
    !
    !      First the determinant is explicitly calculated using the
    !      fundamental definition of the determinant.  If this value is less
    !      that 10**-16 then the matrix is deemed to be singular and the
    !      output value is filled with zeros.  Otherwise, the output matrix
    !      is calculated an element at a time by generating the cofactor of
    !      each element.  Finally, each element in the matrix of cofactors
    !      is multiplied by the reciprocal of the determinant and the result
    !      is the inverse of the original matrix.  Since a temporary matrix
    !      is used, the output matrix may overwrite the input matrix.  NO
    !      INTERNAL CHECKING ON THE INPUT MATRIX M1 IS PERFORMED EXCEPT ON
    !      THE SIZE OF ITS DETERMINANT.  THUS IT IS POSSIBLE TO GENERATE A
    !      FLOATING POINT OVERFLOW OR UNDERFLOW IN THE PROCESS OF
    !      CALCULATING THE MATRIX OF COFACTORS.
    !
    !$ Examples
    !
    !      Suppose that M1 is given by the following matrix equation:
    !
    !           | 0   -1    0 |
    !      M1 = | 0.5  0    0 |  then if INVERT is called according to the
    !           | 0    0    1 |    FORTRAN code:
    !
    !      CALL INVERT (M1, M1)
    !
    !      then M1 will be set to be:
    !
    !           | 0    2    0 |
    !      M1 = |-1    0    0 |
    !           | 0    0    1 |
    !
    !$ Restrictions
    !
    !      The input matrix must be such that generating the cofactors will
    !      not cause a floating point overflow or underflow.  The
    !      strictness of this condition depends, of course, on the computer
    !      installation and the resultant maximum and minimum values of
    !      double precision numbers.

    real(kind=8) :: MTEMP(3,3)
    real :: MDET

    !  Find the determinant of M1 and check for singularity
    mdet =  det(M1)
    IF ( ABS(MDET) .LT. 1.D-16 ) THEN
       invert_matrix = 0.
       RETURN
    END IF
    !
    !  Get the cofactors of each element of M1
    !
    MTEMP(1,1) =  (M1(2,2)*M1(3,3) - M1(3,2)*M1(2,3))
    MTEMP(1,2) = -(M1(1,2)*M1(3,3) - M1(3,2)*M1(1,3))
    MTEMP(1,3) =  (M1(1,2)*M1(2,3) - M1(2,2)*M1(1,3))
    MTEMP(2,1) = -(M1(2,1)*M1(3,3) - M1(3,1)*M1(2,3))
    MTEMP(2,2) =  (M1(1,1)*M1(3,3) - M1(3,1)*M1(1,3))
    MTEMP(2,3) = -(M1(1,1)*M1(2,3) - M1(2,1)*M1(1,3))
    MTEMP(3,1) =  (M1(2,1)*M1(3,2) - M1(3,1)*M1(2,2))
    MTEMP(3,2) = -(M1(1,1)*M1(3,2) - M1(3,1)*M1(1,2))
    MTEMP(3,3) =  (M1(1,1)*M1(2,2) - M1(2,1)*M1(1,2))

    ! Divide the cofactor matrix by MDET to obtain the inverse
    invert_matrix = mtemp / real(mdet,8)

  end function invert_matrix

  function fractional_to_cartesian(xtal_obj,frac_coords)

    ! Return cartesian coordinates given fractional and a transformation matrix

    ! Input type
    type(crystallography_object), intent(in) :: xtal_obj
    real, dimension(:,:) :: frac_coords

    ! Return type
    real, dimension(size(frac_coords,1),size(frac_coords,2)) :: fractional_to_cartesian

    fractional_to_cartesian = matmul(xtal_obj%orth_matrix,frac_coords)

  end function fractional_to_cartesian

  function fractional_array_to_cartesian(xtal_obj,frac_coords)

    ! Return cartesian coordinates given fractional and a transformation matrix

    ! Input type
    type(crystallography_object), intent(in) :: xtal_obj
    real :: frac_coords(3)

    ! Return type
    real :: fractional_array_to_cartesian(3)

    fractional_array_to_cartesian = matmul(xtal_obj%orth_matrix,frac_coords)

  end function fractional_array_to_cartesian

  function cartesian_to_fractional(xtal_obj,cart_coords)

    ! Return cartesian coordinates given fractional and a transformation matrix

    ! Input type
    type(crystallography_object), intent(in) :: xtal_obj
    real, dimension(:,:) :: cart_coords

    ! Return type
    real, dimension(size(cart_coords,1),size(cart_coords,2)) :: cartesian_to_fractional

    cartesian_to_fractional = matmul(xtal_obj%frac_matrix,cart_coords)

  end function cartesian_to_fractional

  function cartesian_array_to_fractional(xtal_obj,cart_coords)

    ! Return cartesian coordinates given fractional and a transformation matrix

    ! Input type
    type(crystallography_object), intent(in) :: xtal_obj
    real :: cart_coords(3)

    ! Return type
    real :: cartesian_array_to_fractional(3)

    cartesian_array_to_fractional = matmul(xtal_obj%frac_matrix,cart_coords)

  end function cartesian_array_to_fractional

  function det(matrix)

    ! Determinant of a 3x3 matrix
 
    real :: det
    real, intent(in) :: matrix(3,3)

    real(kind=8) :: tmp

    tmp = matrix(1,1) * ( matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2) ) & 
         - matrix(1,2) * ( matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1) ) & 
         + matrix(1,3) * ( matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1) )

    det = real(tmp)

  end function det

  ! Access routines

  function get_cell_axes(xtal_obj, array) result(axes)

    ! Return type
    real, dimension(3) :: axes

    ! Input variables
    type(crystallography_object), intent(inout) :: xtal_obj

    real, dimension(3), intent(in), optional :: array

    if (present(array)) then
       xtal_obj%cell%length = array
       call calculate_internal_parameters(xtal_obj)
    end if

    axes = xtal_obj%cell%length

  end function get_cell_axes

  function get_reciprocal_cell_axes(xtal_obj) result(axes)

    ! Return type
    real, dimension(3) :: axes

    ! Input variables
    type(crystallography_object), intent(inout) :: xtal_obj

    axes = xtal_obj%cell%reciprocal_length

  end function get_reciprocal_cell_axes

  function get_cell_angles(xtal_obj, array) result(angles)

    ! Return type
    real, dimension(3) :: angles

    ! Input variables
    type(crystallography_object), intent(inout) :: xtal_obj

    real, dimension(3), intent(in), optional :: array

    if (present(array)) then
       xtal_obj%cell%angle = array/radian
       call calculate_internal_parameters(xtal_obj)
    end if

    angles = xtal_obj%cell%angle*radian

  end function get_cell_angles

  function get_reciprocal_cell_angles(xtal_obj) result(angles)

    ! Return type
    real, dimension(3) :: angles

    ! Input variables
    type(crystallography_object), intent(inout) :: xtal_obj

    angles = xtal_obj%cell%reciprocal_angle*radian

  end function get_reciprocal_cell_angles

  function get_space_group(xtal_obj, groupnum)

    ! Return type
    integer :: get_space_group

    ! Input variables
    type(crystallography_object), intent(inout) :: xtal_obj

    integer, intent(in), optional :: groupnum

    if (present(groupnum)) then
       xtal_obj%space_group = groupnum
    end if

    get_space_group = xtal_obj%space_group

  end function get_space_group

  function get_setting(xtal_obj, settingnum)

    ! Return type
    integer :: get_setting

    ! Input variables
    type(crystallography_object), intent(inout) :: xtal_obj

    integer, intent(in), optional :: settingnum

    if (present(settingnum)) then
       xtal_obj%setting = settingnum
    end if

    get_setting = xtal_obj%setting

  end function get_setting



  !
  ! EQUIVALENCE
  !

  function xtal_eq_xtal(xtal_objL, xtal_objR)

    ! Determines if two quaternions are equivalent (within
    ! an arbitrary tolerance)

    logical :: xtal_eq_xtal

    ! Input variables
    type(crystallography_object), intent(in) :: xtal_objL, xtal_objR

    ! Initialise the return to false, so that any return below will
    ! indicate failure
    xtal_eq_xtal = .FALSE.

    if (xtal_objL%space_group /= xtal_objR%space_group) return
    if (xtal_objL%setting /= xtal_objR%setting) return

    if (sum(abs(xtal_objL%cell%length - xtal_objR%cell%length)) > tolerance) return
    if (sum(abs(xtal_objL%cell%angle - xtal_objR%cell%angle)) > tolerance) return

    ! If we checked all the entries in the zmatrix and had no problems
    xtal_eq_xtal = .TRUE.

  end function xtal_eq_xtal

  function xtal_neq_xtal(xtal_objL, xtal_objR)

    ! Determines if two quaternions are equivalent (within
    ! an arbitrary tolerance)

    logical :: xtal_neq_xtal

    ! Input variables
    type(crystallography_object), intent(in) :: xtal_objL, xtal_objR

    ! Initialise the return to false
    xtal_neq_xtal = .TRUE.

    ! If we checked all the entries in the zmatrix and had no problems
    if (xtal_objL == xtal_objR) xtal_neq_xtal = .FALSE.

  end function xtal_neq_xtal

  ! --------------------------------------------------------------------
  SUBROUTINE Gauss (a,n)       ! Invert matrix by Gauss method
    ! --------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: n
    REAL    :: a(n,n)

    ! - - - Local Variables - - -
    REAL    :: b(n,n), c, d, temp(n)
    INTEGER :: i, j, k, m, imax(1), ipvt(n)
    ! - - - - - - - - - - - - - -
    
    b = a
    ipvt = (/ (i, i = 1, n) /)

    DO k = 1,n

       imax = MAXLOC(ABS(b(k:n,k)))
       m = k-1+imax(1)

       IF (m /= k) THEN
          ipvt( (/m,k/) ) = ipvt( (/k,m/) )
          b((/m,k/),:) = b((/k,m/),:)
       END IF
       d = 1/b(k,k)

       temp = b(:,k)
       DO j = 1, n
          c = b(k,j)*d
          b(:,j) = b(:,j)-temp*c
          b(k,j) = c
       END DO
       b(:,k) = temp*(-d)
       b(k,k) = d
    END DO

    a(:,ipvt) = b

  END SUBROUTINE Gauss

end module crystallography_class
