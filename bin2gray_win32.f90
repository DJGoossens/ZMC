program bin2gray

  use precision
  use fundamental_constants
  use pnm_class
  use cmdline_arguments
  use iso_varying_string
  use string_functions, only: join
  use file_functions, only: stderr, stdin, stdout
  use image_transforms


  implicit none

  character(len=*), parameter :: version = "$Id: bin2gray.f90,v 1.6 2007/06/04 05:23:11 aidan Exp $"

  ! This program will take the file 'intensity.bin' which was
  ! created by the program DIFFUSE and will convert the computed
  ! diffuse intensities to a series of 8-bit gray-scale images
  ! that can be viewed by a range of software packages.  The output
  ! will be to a series of files named img###.pgm where ### is the
  ! number of the reciprocal plane it contains.  The output format
  ! is Jef Poskanzer's portable graymap (pgm) file format which can
  ! be read by many programs and converted to just about any other
  ! image format.

  !! $Log: bin2gray.f90,v $
  !! Revision 1.6  2007/06/04 05:23:11  aidan
  !! Added option '--square' to keep the output square shape (ignore the
  !! sin(theta)/lambda) limits in the intensity file output from diffuse).
  !!
  !! Revision 1.5  2007/04/04 05:45:03  aidan
  !! Added option autoscale to force automatic scaling of the output to
  !! average + 4*stddev. Previously this scaling was only done if the
  !! intensity overflowed the dynamic range of the pgm (be it 8-bit or
  !! 16-bit).
  !!
  !! Revision 1.4  2005/05/30 05:38:32  aidan
  !! Added an option to n-fold average the data. This uses the image_transforms
  !! module to rotate the image and then adds up the contributions. Seems to work
  !! ok but does not produce an *identical* picture to the --twofold option (there
  !! are small differences around the edge, the middle is identical).
  !!
  !! Revision 1.3  2004/09/09 23:44:09  aidan
  !! New working version of bin2gray. Should be fully backwards compatible,
  !! though there are issues with the 'addbottom' and 'addleft' functionality --
  !! Brent assumed the central pixels were an origin of sorts, so only replicated
  !! numv-1 and numu-1 pixels to the bottom and left respectively. We have just
  !! doubled the size of the image in this version.
  !!
  !! Also added horizontal and vertical mirroring (--hmirror, --vmirror) and two
  !! fold averaging (--twofold).
  !!
  !! If no command line options are specified we assume "compatibility mode" and
  !! ask for input from stdin and default to 8bit output. If any command line
  !! options are used, even '--' to signify no more options, or intensity files
  !! are specified as command line arguments then we have different behaviour.
  !! The default depth is then 16bit, we only normalise the data if it will
  !! overflow (i.e. > 65535 or > 255 if the --8bit options is used) and then we
  !! default to mean + 4*stddev, with the --norm=value option available if th
  !! user wishes to normalise to another value. When in compatibility mode the
  !! program will only open 'intensity.bin' files and outputs img00x.pgm. When
  !! we specify an intensity file on the command line the program will attempt
  !! to name the output pgm based on the intensity file, replacing '.bin' with
  !! '_xxx.pgm'. There is an exception to this rule however. Specifying an
  !! intensity file with the name 'intensity.bin' means we fall back to the
  !! img00x.pgm naming scheme. This is so the user can take advantage of some of
  !! the command line options but still get the same output files as they did
  !! before.
  !!
  !! Revision 1.2  2004/09/06 00:11:14  aidan
  !! A *just* working version of bin2gray which uses the pnm_class to output the
  !! pgm file. Does not support adding horizontal or vertical mirrors. Just a bug
  !! fix release for Darren.
  !!

  integer :: numu, numv, numw, nu, nv, nw, nsites, ntypes, nlots
  integer :: i, j, irec, irlen, iseed, jseed, error

  real, dimension(:,:), allocatable :: dsi, stl, dsitmp

  real ah_o(3),ah_u(3),ah_v(3),ah_w(3)
  real cell(6), stlmax, anorm, xmax, xmin, xnorm
  real uin(3),vin(3),win(3)

  integer :: maxvalue = 65535

  double precision xmean,xstd
  integer csize(3),ls_xyz(3), nfold
  character(len=1) qalat,qhorz,qvert,qnorm,qval,qbnd
  character(len=70) descrip
  character(len=1000) fname 

  character(len=1000) intensity_file, tmpprefix

  logical :: verbose = .TRUE., compatibility_mode = .TRUE.
  logical :: addleft = .FALSE., addbottom = .FALSE.
  logical :: fixed_max = .FALSE., fixed_name = .FALSE.
  logical :: twofold, hmirror, vmirror, rotaverage, got_prefix = .FALSE.
  logical :: auto_scale = .FALSE., layer_scale_separate = .FALSE.
  logical :: make_round = .TRUE.

  type (pnm_object) :: pnm
  type (varying_string) :: myoptions(14), prefix

  real, parameter :: no_data_value = -9999

!  integer, external :: nb

  common/sizes/ah_o,uin,vin,win,numu,numv,numw,stlmax

  ! Before we process command line arguments check if we have ANYTHING AT ALL
  ! on the command line. If we do then we assume we are not using bin2gray in
  ! 'compatibility mode', which means we will not look for input on stdin and
  ! we might generate different output files 
  if (num_args() /= 0) compatibility_mode = .FALSE.

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'version'
  myoptions(3) = 'norm'
  myoptions(4) = 'quiet'
  myoptions(5) = 'addbottom'
  myoptions(6) = 'addleft'
  myoptions(7) = 'hmirror'
  myoptions(8) = 'vmirror'
  myoptions(9) = '8bit'
  myoptions(10) = 'twofold'
  myoptions(11) = 'prefix'
  myoptions(12) = 'rotave'
  myoptions(13) = 'autoscale'
  myoptions(14) = 'square'

  ! This call parses the command line arguments for command line options
  call get_options(myoptions, error)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (error > 0) then
     write(stderr,*) 'ERROR! Unknown option(s): ',join(bad_options()," ")
     call usage
     STOP
  end if

  ! Check if we want to print the version number
  if (option_exists('version')) then
     call print_version
     STOP
  end if

  ! Check if we just want to print the usage
  if (option_exists('help')) then
     call usage
     STOP
  end if

  ! We will switch to 8bit output if we asked for it via an option
  ! or we are running in compatibility mode
  if (option_exists('8bit') .OR. compatibility_mode) then
     maxvalue = 255
  end if

  ! Check if we don't want verbose output
  verbose = (.NOT. option_exists('quiet')) 

  ! Do we want to add a left half?
  addleft = option_exists('addleft')
  ! Do we want to add a bottom half?
  addbottom = option_exists('addbottom')

  hmirror = option_exists('hmirror')
  vmirror = option_exists('vmirror')
  twofold = option_exists('twofold')
  auto_scale = option_exists('autoscale')
  make_round = .not. option_exists('square')

  if (option_exists('rotave')) then
     if (twofold) then 
        write(stderr,*) 'Requested rotave and twofold simultaneously. Ignoring twofold option.'
        twofold = .FALSE.
     end if
     if (.NOT. has_value('rotave')) then
        write(stderr,*) 'Option rotave must have a value!'
        call usage
        stop
     end if
     ! Get the n-fold rotation ...
     nfold = get_value('rotave')
     rotaverage = .TRUE.
  end if

  ! See if we have specified a maximum we wish to scale the data to
  if (option_exists('norm')) then
     if (auto_scale) write(stderr,*) 'Option autoscale ignored when option norm also specified'
     ! Make sure we have a value
     if (.NOT. has_value('norm')) then
        write(stderr,*) 'Option norm must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     xnorm = get_value('norm')
     fixed_max = .TRUE.
  else
     ! Default to normalising by ave + 4*stddev if we are not in
     ! compatibility mode -- in this latter case we will read in 
     ! a value for qnorm
     if (.NOT. compatibility_mode) qnorm = 'y'
  end if

  ! See if we have specified a prefix for our output file name
  if (option_exists('prefix')) then
     ! Make sure we have a value
     if (.NOT. has_value('prefix')) then
        write(stderr,*) 'Option prefix must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     prefix = get_value('prefix')
     got_prefix = .TRUE.
  end if


  ! Check if we have specified an intensity file as a command line argument
  if (num_args() == 0) then
     intensity_file = 'intensity.bin'
  else
     intensity_file = next_arg()
  end if

  do

     ! Open the input file and find out what record size it is ...
     open(unit=1,file=intensity_file,status='old',               &
          &          form='unformatted',access='direct',recl=74)
     read(1,rec=1)irlen

     close (unit=1)

     ! Re-open the file (record size now known) and read header ... 
     open(unit=1,file=intensity_file,status='old',               &
          &          form='unformatted',access='direct',recl=irlen)
     read(1,rec=1)irlen,descrip
     read(1,rec=2)ah_o,ah_u,ah_v,ah_w,numu,numv,numw
     read(1,rec=3)cell,stlmax
     read(1,rec=4)iseed,jseed,csize,qbnd,nsites,ntypes
     read(1,rec=5)ls_xyz,nlots,qalat
     close (unit=1)

     if (verbose) then
        ! Now echo to standard output what was in the header ...
        ! call print_version(stdout)
        write(stdout,'(A,A,A,I0)') ' File ',trim(intensity_file),' opened and irlen = ',irlen
        write(stdout,*)
        write(stdout,*) '---------------------------------------------------'
        write(stdout,*) descrip
        write(stdout,*) '---------------------------------------------------'
        write(stdout,*)
        write(stdout,*) '  The computation volume is defined by:'
        print 101,(ah_o(i), i=1,3),(ah_u(i), i=1,3)
        print 101,(ah_o(i), i=1,3),(ah_v(i), i=1,3)
        print 101,(ah_o(i), i=1,3),(ah_w(i), i=1,3)
        write(stdout,'(A,I0,A,I0,A,I0)') '   Image size is ',numu,' X ',numv,' X ',numw
        write(stdout,*) '  sin(theta)/lambda maximum = ',stlmax
        write(stdout,'(A,1x,3(I0,2X))') '   Crystal size is: ',(csize(i), i=1,3)
        write(stdout,'(A,1x,3(I0,2X))') '   Lot Size is: ',(ls_xyz(i), i=1,3)
        write(stdout,'(A,I0)') '   Number of Lots Computed: ',nlots
        write(stdout,*) '  Subtract Average Lattice? ',qalat
        write(stdout,*) ' '
101     format(5x,'(',2(f7.2,','),f7.2,') => (',2(f7.2,','),f7.2,')' )
     end if
     
     ! Get the increments along each of the three axes ...
     !vocl loop,nopreex
     do i=1,3
        uin(i)=0.d0
        vin(i)=0.d0
        win(i)=0.d0
        if (numu.gt.1) uin(i)= (ah_u(i)-ah_o(i))/real(numu-1)
        if (numv.gt.1) vin(i)= (ah_v(i)-ah_o(i))/real(numv-1)
        if (numw.gt.1) win(i)= (ah_w(i)-ah_o(i))/real(numw-1)
     end do

     ! Check standard input if we are in compatibility mode
     ! if (nb() > 0) then
     if (compatibility_mode) then

        ! Ask some pertinent questions about how to do this ...
        write(stdout,*) '  Put a bottom half on using a mirror?'
        read(stdin,'(a)')qhorz
        write(stdout,*) '  Put a left half on using a mirror?'
        read(stdin,'(a)')qvert
        write(stdout,*) '  Use the default ave+4*std normalization?' 
        read(stdin,'(a)')qnorm
        if ( (qnorm.ne.'Y').and.(qnorm.ne.'y') ) then
           write(stdout,*) '  Normalize all planes by the same value?'
           read(stdin,'(a)')qval
           if ( (qval.eq.'Y').or.(qval.eq.'y') ) then
              write(stdout,*) '  O.K., what value should I normalize to?'
              read(stdin,*)xnorm
              fixed_max = .TRUE.
           else
              layer_scale_separate = .TRUE.
           end if
        else
           auto_scale = .TRUE.
        end if

        addbottom = ( (qhorz.eq.'Y').or.(qhorz.eq.'y') )
        addleft   = ( (qvert.eq.'Y').or.(qvert.eq.'y') )
        
     end if
  
     ! Allocate our arrays
     allocate(dsi(numu,numv),stl(numu,numv))

     ! Check if we want to append mirrors to the left and bottom
     nu = numu
     if (addleft)   nu = 2 * numu
     nv = numv
     if (addbottom) nv = 2 * numv
     
     ! Allocate space for a temporary version which may or may not
     ! be bigger than the original data
     allocate(dsitmp(nu,nv))
        
     ! Loop over all reciprocal planes that are in 'intensity.bin'
     do nw=1,numw
        if (verbose) write(stdout,'(A,I0,A,I0,A)')'  Working on plane # ',nw,'/',numw,' ...'
        ! Extract from file the values associated with plane nw ...
        open(unit=1,file=intensity_file,status='old',             &
             &            form='unformatted',access='direct',recl=irlen)
        do j=1,Numv
           irec=numv*(nw-1)+j+5
           read(1,rec=irec)(dsi(i,j), i=1,Numu)
        end do
        close(unit=1)

        ! Find the Max, Min, Ave, and Std Dev of this image ...
        xmax=-1.e32
        xmin=1.e32
        xmean=0.d0
        xstd=0.d0
        do j=1,numv
           do i=1,numu
              if (dsi(i,j).gt.xmax) xmax=dsi(i,j)
              if (dsi(i,j).lt.xmin) xmin=dsi(i,j)
              xmean=xmean+dsi(i,j)
              xstd=xstd+dsi(i,j)**2
           end do
        end do
        xmean=xmean/(real(numu)*real(numv))
        xstd=xstd/(real(numu)*real(numv)) - xmean**2
        xstd=sqrt(xstd)
        if (verbose) then
           write(stdout,*) '    Range = ',xmin,' => ',xmax
           write(stdout,*) '    Mean , STD  = ', xmean,xstd
        end if
        
        ! What do we normalize this image to? ...
        ! If we haven't specified a value to normalise to then we
        ! check out the possibilities
        if (.NOT. fixed_max) then

           if ( compatibility_mode .and. layer_scale_separate ) then
              write(stdout,*) '    Normalize to what? (0 gives ave + 4*std) '
              read(stdin,*)xnorm
              if (xnorm.lt.1.e-6) xnorm=xmean+4.*xstd
           else
              ! Only normalise data if it overflows our accuracy or
              ! we specified auto scaling in the options and 
              ! then we normalise to mean + 4*std.dev
              if ((xmax > real(maxvalue)) .or. auto_scale) then
                 xnorm=xmean+4.*xstd
              else
                 xnorm=maxvalue
              end if
           end if
        end if
        if (verbose) write(stdout,*) '    Normalizing to ',xnorm

        if (make_round) then
        
           ! Normalize all values in the image but set to white anything
           ! outside of the maximum sin(theta)/lambda ...
           call stltab(stl,cell,nw)
           do j=1,numv
              do i=1,numu
                 ! dsi(i,j)=dsi(i,j)*255./xnorm		!altered to avoid zero's may-98
                 dsi(i,j)=(dsi(i,j)*(real(maxvalue)-1.)/xnorm)+1.
                 ! if (dsi(i,j).gt.real(maxvalue)) dsi(i,j)=real(maxvalue)
                 ! if (stl(i,j).gt.stlmax) dsi(i,j)=real(maxvalue)
                 if (stl(i,j).gt.stlmax) dsi(i,j)=no_data_value
                 ! inc(i,j)=char( nint(dsi(i,j)) )
              end do
           end do

        end if
        
        ! Fill in the original data
        dsitmp((nu-numu)+1:,(nv-numv)+1:) = dsi
        
        ! Left mirror (top left quadrant only)
        if (addleft) dsitmp(1:numu,(nv-numv)+1:nv) = dsi(numu:1:-1,:)
        ! dsitmp((nu-numu)+1:nu,1:numv) = dsi(:,numv:1:-1)

        ! Bottom mirror (entire bottom half)
        if (addbottom) dsitmp(:,1:numv) = dsitmp(:,nv:(nv-numv)+1:-1)
        ! dsitmp(1:numu,:) = dsitmp(nu:(nu-numu)+1:-1,:)

        ! 2-fold average the data ? (This is hokey two-fold rotation averaging
        ! where we run the x and y directions in reverse order, i.e. a flip
        ! and a flop)
        if (twofold) call average(dsitmp,dsitmp(nu:1:-1,nv:1:-1))

        ! n-fold average the data ?
        if (rotaverage) call rotation_average(dsitmp, nfold)
        
        ! Apply horizontal mirror?
        if (hmirror) call average(dsitmp,dsitmp(nu:1:-1,:))

        ! Apply vertical mirror?
        if (vmirror) call average(dsitmp,dsitmp(:,nv:1:-1))
        ! if (vmirror) dsi = (dsi + dsi(:,numv:1:-1))/2.
        
        ! Map missing values to zero
        where(dsitmp == no_data_value) dsitmp = 0.0

        ! Truncate to maximum value
        where(dsitmp > real(maxvalue)) dsitmp = real(maxvalue)
           
        ! Stuff the data into a pnm object
        ! pnm = nint(dsitmp)
        call new(pnm, nint(dsitmp), maximum=maxvalue)
        
        ! If the program was invoked the 'old style' way or we have
        ! an intensity file named the old way (intensity.bin) then we'll
        ! generate an old style pgm file. Otherwise we'll generate
        ! the pgm file name from the intensity file name
        if (got_prefix) then
           ! Build a file name from the specified prefix
           call getname(char(prefix),nw,fname)
        else if (compatibility_mode .OR. (intensity_file == 'intensity.bin')) then
           ! Build a file name with prefix 'img'
           call getname('img',nw,fname)
        else
           ! Extract a prefix from the intensity file name
           i = index(intensity_file,'.bin',back=.TRUE.) - 1
           if (i == -1) i = len_trim(intensity_file)
           call getname(intensity_file(1:i),nw,fname)
        end if
        
        if (verbose) write(stdout,*) 'Writing pgm ',trim(fname)
        call write(pnm, fname)
        
     end do
     ! We have finished with all of the reciprocal planes now.
     
     deallocate(dsi,dsitmp,stl)
     close(unit=1)
        
     if (.NOT. have_args()) exit
     intensity_file = next_arg()
     
  end do

contains
    
  subroutine usage
    
    write(stderr,*)
    write(stderr,*) 'Usage: bin2gray --help --quiet --addleft --addbottom --8bit --norm'
    write(stderr,*) '                --twofold --rotave=n --hmirror --vmirror binfile(s)'
    write(stderr,*)
    write(stderr,*) '  --help      - print this message'
    write(stderr,*) '  --version   - print the version number'
    write(stderr,*) '  --quiet     - non-verbose output'
    write(stderr,*) '  --norm      - scale the data to this value'
    write(stderr,*) '  --autoscale - scale the data average + 4*stddev'
    write(stderr,*) '  --twofold   - two-fold average the data.'
    write(stderr,*) '  --rotave=n  - n-fold average the data (n is any integer).'
    write(stderr,*) '  --hmirror   - apply horizonal mirror averaging'
    write(stderr,*) '  --vmirror   - apply vertical mirror averaging'
    write(stderr,*) '  --addleft   - add a left half using a vertical mirror'
    write(stderr,*) '  --addbottom - add a bottom half using a horizonal mirror'
    write(stderr,*) '  --8bit      - specify an 8 bit output image (defaults to 16 bit)'
  ! write(stderr,*) '  --outfile   - specify an output file (default is replace .bin with _xxx.pgm).'
    write(stderr,*) '  --prefix    - specify a prefix for your image files (default is "img").'
    write(stderr,*) '  --square    - leave image square (do not clip to sin(theta)/lambda limits)'
    write(stderr,*)
    
  end subroutine usage
    
  subroutine print_version(unit)

    integer, optional :: unit

    integer :: myunit

    myunit = stderr
    if (present(unit)) myunit = unit

    write(myunit,*) version(6:index(version,":",back=.true.)-1)
    
  end subroutine print_version
    
  subroutine stltab(stl,cell,nw)
    
    ! Calculates the value of sin(theta)/lambda for each element 
    ! of the diffuse scattering array.  The cell parameters (a,b,c,
    ! cos(bc),cos(ac),cos(ab)), computation spacings and dimmensions,
    ! and the reciprocal section must be provided.  This table is 
    ! used simply to provide a white mask over values that have 
    ! exceeded the parameter stlmax.
    
    real, dimension(:,:), intent(inout) :: stl
    real, intent(in)                    :: cell(6)
    integer, intent(in)                 :: nw
    
    integer numu, numv, numw, ii, jj
    real ah_o(3),uin(3),vin(3),win(3), stlmax
    real a1, a2, a3, c1, c2, c3, s1, s2, s3
    real V, S11, S22, S33, S12, S23, S13, XH1, XH2, XH3
    
    common/sizes/ah_o,uin,vin,win,numu,numv,numw,stlmax
    
    ! The cell parameters ...
    a1=cell(1)
    a2=cell(2)
    a3=cell(3)
    c1=cell(4)
    c2=cell(5)
    c3=cell(6)
    s1=sin(acos(c1))
    s2=sin(acos(c2))
    s3=sin(acos(c3))
    
    ! Some constants that I need (reference Cullity) ...
    V=a1*a2*a3*sqrt(1.-c1**2-c2**2-c3**2+2.*c1*c2*c3)
    S11=(a2**2)*(a3**2)*(s1**2)
    S22=(a1**2)*(a3**2)*(s2**2)
    S33=(a1**2)*(a2**2)*(s3**2)
    S12=(a1)*(a2)*(a3**2)*(c1*c2-c3)
    S23=(a2)*(a3)*(a1**2)*(c2*c3-c1)
    S13=(a1)*(a3)*(a2**2)*(c1*c3-c2)
    
    ! Loop over all pixels in the image ...
    do jj=1,Numv
       do ii=1,Numu
          xh1=ah_o(1)+real(ii-1)*uin(1)+real(jj-1)*vin(1)+nw*win(1)
          xh2=ah_o(2)+real(ii-1)*uin(2)+real(jj-1)*vin(2)+nw*win(2)
          xh3=ah_o(3)+real(ii-1)*uin(3)+real(jj-1)*vin(3)+nw*win(3)
          stl(ii,jj)=sqrt(S11*xh1**2+S22*xh2**2+S33*xh3**2            &
               &                  +2.*S12*xh1*xh2+2.*S13*xh1*xh3+2.*S23*xh2*xh3)
          stl(ii,jj)=0.5*stl(ii,jj)/V
       end do
    end do
    return
  end subroutine stltab
  
  subroutine getname(prefix,nw,fname)

    ! Builds the character constant 'img###.pgm' where ###=nw
    character(len=*), intent(in)  :: prefix
    integer, intent(in)           :: nw
    character(len=*), intent(out) :: fname

!!$    ihun=nw/100
!!$    iten=nw/10 - ihun*10
!!$    ione=nw - 100*ihun - 10*iten
!!$    ihun=ihun+48
!!$    iten=iten+48
!!$    ione=ione+48
!!$    fname='img'//char(ihun)//char(iten)//char(ione)//'.pgm'

    if ((numw > 1) .OR. compatibility_mode) then
       write(fname,'(A,I3.3,".pgm")') prefix,nw
    else
       write(fname,'(A,".pgm")') prefix
    end if
    return
  end subroutine getname

  subroutine average (original, transform)

    ! Interface variables
    real, dimension(:,:), intent(inout) :: original
    real, dimension(:,:), intent(in)    :: transform

    ! Local variables
    real, dimension(size(original,1),size(original,2)) :: flipflopcopy

    flipflopcopy = transform

    where (flipflopcopy == no_data_value)
       ! we have no data in our copy, get data from the original
       flipflopcopy = original
    end where
           
    where (original == no_data_value)
       ! we have no data in our original, get data from the copy
       original = flipflopcopy
    elsewhere
       ! we have data in both now -- take average
       original = (original + flipflopcopy) / 2.
    end where

  end subroutine average
    
  subroutine rotation_average(input, nfold)

    ! Interface variables
    real, dimension(:,:), intent(inout) :: input
    integer, intent(in)                 :: nfold

    ! Local variables
    real, dimension(size(input,1),size(input,2))    :: dsitot, dsibuffer
    integer, dimension(size(input,1),size(input,2)) :: bufcount, totcount
    real(kind=rd_kind) :: step
    integer :: i, j, k

    dsitot = input
    totcount = 0
    where(dsitot /= no_data_value) totcount = 1
    
    step = (360.d0/real(nfold,rd_kind))/radian
    
    do i = 1, nfold-1
       ! print *,'number ',i,real(i,rd_kind)*step*radian
       dsibuffer = input
       call rotate_image(dsibuffer,real(i,rd_kind)*step,no_data_value)
       do j = 1, size(input,1)
          do k = 1, size(input,2)
             if (dsibuffer(j,k) /= no_data_value) then
                if (dsitot(j,k) /= no_data_value) then
                   dsitot(j,k) = dsitot(j,k) + dsibuffer(j,k)
                else
                   dsitot(j,k) = dsibuffer(j,k)
                end if
                totcount(j,k) = totcount(j,k) + 1
             end if
          end do
       end do
    end do

    where(totcount > 0) dsitot = dsitot/real(totcount)
    input = dsitot

  end subroutine rotation_average

end program bin2gray
