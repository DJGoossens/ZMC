PROGRAM DIFFUSE

  use file_functions

  implicit none

  !----------------------------------------------------------------
  ! August 2009: DJGoossens.  Making all arrays allocatable to
  ! avoid need to edit code.  Also, I don't like common blocks,
  ! so I will probably 'contain' the subroutines.
  !
  ! Modify use of cur.count (comment out)
  !
  ! Also, remove use of implicit typing.
  !
  ! Have not dealt with MAT
  !----------------------------------------------------------------


  !c      Computes the diffuse diffracted intensity from a simulated
  !c      disordered crystal.  Input parameters are given in the file
  !c      'diffuse.in' - output is to the file 'intensity.bin'.  An
  !c      interface routine that provides access to information in the
  !c      simulation must be linked to this program.

  !c       ----------------------------------------
  !c       Version 2.1a  May 11, 1995   *BDB*
  !c       Version 2.1 adds option to subtract exact ave. lattice
  !c       version 2.0b computes Biso with orthoganal axes.
  !c       ----------------------------------------



  !        parameter (MI=501,MJ=501,MT=10,MS=155,MAT=6000)
  !parameter (MAT=6000)
  integer, parameter :: MAT=6000

  !parameter (MASK=2**14-1)
  integer, parameter  :: MASK=2**14-1

  !        complex :: csf(MI,MJ),acsf(MI,MJ),tcsf(MI,MJ)
  complex, allocatable :: csf(:,:),acsf(:,:),tcsf(:,:)

  ! Moved this to here from scatf
  complex, allocatable :: cftab(:,:)
  integer, allocatable :: itest(:)

  complex :: cex(0:MASK),cfact(0:1999)

  !        real :: dsi(MI,MJ),xat(MAT,3),cell(6)
  real, allocatable :: dsi(:,:)
  real :: xat(MAT,3),cell(6)

  real :: ah_o(3),ah_u(3),ah_v(3),ah_w(3)
  double precision :: uin(3),vin(3),win(3)

  !        real a1(MT),b1(MT),a2(MT),b2(MT),a3(MT),b3(MT)
  real,allocatable :: a1(:),b1(:),a2(:),b2(:),a3(:),b3(:)

  !        real a4(MT),b4(MT),c1(MT),fp(MT),fpp(MT)
  real, allocatable :: a4(:),b4(:),c1(:),fp(:),fpp(:)

  !real :: Biso(MS,MT),wt(MS,MT),ax(MS,MT),ay(MS,MT),az(MS,MT)
  real, allocatable :: Biso(:,:),wt(:,:),ax(:,:),ay(:,:),az(:,:)

  !integer :: istl(MI,MJ),csize(3),ls_xyz(3),lbeg(3),l(3)
  integer, allocatable :: istl(:,:)
  integer :: csize(3),ls_xyz(3),lbeg(3),l(3)

  character*70 :: descrip

  character*4 ,allocatable :: c_at(:)
  character*4 :: cat

  character*1 :: qalat,qbnd

  character*32 :: fname

  integer ::  numv, numu, numw, iseed, jseed, nsites, ntypes, nlots, ierr, isite
  integer ::  inum, irlen,itype
  integer :: i,j, ii, jj, natm,nw,mp,l1,l2,l3,ir,kk,nlot,ncell,irec

  real :: stlmax,xx,yy,zz,wx,B,x,y,z,denom,xnorm

  !        common/sizes/ah_o,uin,vin,win,numu,numv,numw,stlmax
  !        common/tabels/cex,istl,cfact
  !        common/scatter/a1,b1,a2,b2,a3,b3,a4,b4,c1,fp,fpp
  !        common/crystal/csize,cell
  !        common/info1/ls_xyz,iseed,jseed,nlots,nsites,ntypes
  !        common/info2/ah_u,ah_v,ah_w
  !        common/info3/c_at,descrip,qalat,qbnd

  !c       Find out everything we need to know to run the program ...


  write(6,*)'-----------------------------------------------------------------'
  write(6,*)'This is DZMC, top calculte diffuse sections from output from ZMC'
  write(6,*)'Distributed under the Academic Free License version 3.0; details in'
  write(6,*)'the files DISCLAIMER.txt and COPYRIGHT.txt'
  write(6,*)'(http://opensource.org/licenses/AFL-3.0)'
  write(6,*)'-----------------------------------------------------------------'

  call READINF

  allocate(cftab(0:1999,ntypes),stat=ierr)
  if (ierr.ne.0) call allerr('cftab(0:1999,ntypes)              ')
  allocate(itest(ntypes),stat=ierr)
  if (ierr.ne.0) call allerr('itest(ntypes)              ')
  itest = 0


  !c       Write the input info to standard output ...

  call WRITEINF

  !c       Initialize pseudo-random number generator ...

  call RSEED(iseed,jseed)

  !c       Initialize the complex exponent table ...

  call cexpt(cex)

  !c       Open the output file and write the appropriate header ...
  !c       For 4byte record size (i.e. DEC) use these lines
  !c       irlen=numu
  !c       if (irlen.lt.20) irlen=20
  !c       Else use these lines

  irlen=numu*4
  if (irlen.lt.74)irlen=74

  write(6,*) 'intensity file?'
  read(5,1122) fname
  write(6,1122) fname
1122 format(a32)

  !c         open(unit=1,file='intensity.bin',status='unknown',

  open(unit=1,file=fname,status='unknown',     &
       form='unformatted',access='direct',recl=irlen)
  write(1,rec=1)irlen,descrip
  write(1,rec=2)ah_o,ah_u,ah_v,ah_w,numu,numv,numw
  write(1,rec=3)cell,stlmax
  write(1,rec=4)iseed,jseed,csize,qbnd,nsites,ntypes
  write(1,rec=5)ls_xyz,nlots,qalat
  close(unit=1)

  !c       Get average atom positions, occupancies, and Debye factors ...

  print *,'  Average Structure ...'
           do 10, isite=1,nsites
            print *,'   Site # ',isite,' :'
            do 10, itype=1,ntypes
             cat=c_at(itype)
             call AVINFO(isite,cat,xx,yy,zz,wx,B)
             ax(isite,itype)=xx
             ay(isite,itype)=yy
             az(isite,itype)=zz
             wt(isite,itype)=wx
             Biso(isite,itype)=B
             if (nint(wx*100).gt.0) then
              print 101,cat,xx,yy,zz,wt(isite,itype),Biso(isite,itype)
101           format (6x,A,'atom at (',3f6.3,'), xocc = ',     &
                      f5.2,', Biso = ',f6.3)
             end if
10         continue



!c       If there is no periodic boundary limit the crystal size ...

         if ( (qbnd.ne.'y').and.(qbnd.ne.'Y') ) then
           csize(1)=csize(1)-ls_xyz(1)
           csize(2)=csize(2)-ls_xyz(2)
           csize(3)=csize(3)-ls_xyz(3)
         end if

!c       -------------------------------------------------------------
!c       For each reciprocal plane along the w-axis compute scattering
!c       -------------------------------------------------------------

         do 100, nw=1,numw
           print *,'  '
           print *,'  --------------------------'
           print *,'  Reciprocal plane # ',nw
           print *,'  --------------------------'

!c          Initialize table of sin(theta)/lambda ...

            call STLTAB(istl,cell,real(nw))

!c          Zero some arrays ...

            do 110, j=1,Numv
             do 110, i=1,Numu
              csf(i,j)=cmplx(0.d0,0.d0)
              acsf(i,j)=cmplx(0.d0,0.d0)
              dsi(i,j)=0.d0
110        continue

!c         Calculate the average structure factor if asked to do so...
           if( (qalat.eq.'y').or.(qalat.eq.'Y') ) then
             print *,'  Computing Average Scattering using Biso ...'

!c            SCATTERING FROM AVERAGE CELL
             do 120, isite=1,nsites
!cthp             print *,'   Site # ',isite,' :'
              do 120, itype=1,ntypes
               cat=c_at(itype)
               xat(1,1)=ax(isite,itype)
               xat(1,2)=ay(isite,itype)
               xat(1,3)=az(isite,itype)
               wx=wt(isite,itype)
               B=Biso(isite,itype)
	       natm=1
	       if (nint(wx*1000).gt.1) then
                call SCATF(cfact,itype)
                call DEBYE(cfact,B,wx)
                call STRUCF(tcsf,xat,natm,real(nw))
                do 125, j=1,Numv
                 do 125, i=1,Numu
                  acsf(i,j)=acsf(i,j)+tcsf(i,j)
125             continue
               end if
120          continue

!c            INTERFERENCE FUNCTION OF LOT SHAPE
             xx=0.
             yy=0.
             zz=0.
             call GETAV(xat,natm,ls_xyz,xx,yy,zz)
             do 130, i=0,1999
	       cfact(i)=cmplx(1.,0.)
130	     continue
	     call STRUCF(tcsf,xat,natm,real(nw))
             do 131, j=1,Numv
              do 131, i=1,Numu
               acsf(i,j)=acsf(i,j)*tcsf(i,j)
131	     continue
             print *,' '
           end if

!c         A more exact method if asked ...
           if( (qalat.eq.'e').or.(qalat.eq.'E') ) then
	     if(qalat.eq.'e') then
	       mp=( csize(1)*csize(2)*csize(3) )/20
               print *,'  Average Scattering from 5% of Crystal ...'
!cthp               print *,'  This could take some time ...'
!c              SCATTERING FROM AVERAGE CELL
               do 140, isite=1,nsites
!cthp                print *,'   Site # ',isite,' :'
                do 140, itype=1,ntypes
                 cat=c_at(itype)
                 call SCATF(cfact,itype)
		 natm=0
		 do 145, ir=1,mp
		   call RANLOC(csize,l)
		   call READAT(l,isite,cat,inum,x,y,z)
		   if(inum.eq.1) then
		     natm=natm+1
		     xat(natm,1)=x
		     xat(natm,2)=y
		     xat(natm,3)=z
		     if(natm.eq.MAT) then
                      call STRUCF(tcsf,xat,natm,real(nw))
                      do 146,  j=1,Numv
                       do 146, i=1,Numu
                        acsf(i,j)=acsf(i,j)+tcsf(i,j)
146                   continue
		      natm=0
		     end if
		   end if
145              continue
                 call STRUCF(tcsf,xat,natm,real(nw))
                 do 147,  j=1,Numv
                  do 147, i=1,Numu
                   acsf(i,j)=acsf(i,j)+tcsf(i,j)
147              continue
140            continue
	      else
	       l1=csize(1)
	       l2=csize(2)
	       l3=csize(3)
	       mp=l1*l2*l3
               print *,'  Computing EXACT Average Scattering ...'
!cthp               print *,'  This could take a loooooooong time ...'
!c              SCATTERING FROM AVERAGE CELL
               do 150, isite=1,nsites
!cthp                print *,'   Site # ',isite,' :'
                do 150, itype=1,ntypes
                  cat=c_at(itype)
                  call SCATF(cfact,itype)
		  natm=0
		  do 155, kk=1,l3
		   l(3)=kk
		   do 155, jj=1,l2
		    l(2)=jj
		    do 155, ii=1,l1
		     l(1)=ii
		     call READAT(l,isite,cat,inum,x,y,z)
		     if(inum.eq.1) then
		       natm=natm+1
		       xat(natm,1)=x
		       xat(natm,2)=y
		       xat(natm,3)=z
		       if(natm.eq.MAT) then
                        call STRUCF(tcsf,xat,natm,real(nw))
                        do 156,  j=1,Numv
                         do 156, i=1,Numu
                          acsf(i,j)=acsf(i,j)+tcsf(i,j)
156                     continue
		        natm=0
		       end if
		     end if
155		    continue
                    call STRUCF(tcsf,xat,natm,real(nw))
                    do 157,  j=1,Numv
                     do 157, i=1,Numu
                       acsf(i,j)=acsf(i,j)+tcsf(i,j)
157                 continue
150            continue
             end if

!c            INTERFERENCE FUNCTION OF LOT SHAPE
             xx=0.
             yy=0.
             zz=0.
             call GETAV(xat,natm,ls_xyz,xx,yy,zz)
             do 160, i=0,1999
	       cfact(i)=cmplx(1.,0.)
160	     continue
	     call STRUCF(tcsf,xat,natm,real(nw))
	     denom=1./mp
             do 161, j=1,Numv
              do 161, i=1,Numu
               acsf(i,j)=acsf(i,j)*tcsf(i,j)*cmplx(denom,0.)
161	     continue
             print *,' '
           end if

!c         Now compute the diffuse scattering by averaging the total
!c         intensities from 'nlots' regions of the simulated crystal.
!c         Loop over all atom types and then over all atom sites.

           do 200 nlot=1,nlots
             call RANLOC(csize,lbeg)
!c	             print *,'    Lot number = ',nlot,' (l,m,n) = (',
!c     *               (lbeg(i),i=1,3),')'
             do 210, itype=1,ntypes
                cat=c_at(itype)
                call SCATF(cfact,itype)
                do 210, isite=1,nsites
	         if (nint(wt(isite,itype)*1000).gt.1) then
                 call GETATM(xat,natm,lbeg,csize,ls_xyz,isite,cat,ncell)
!c                 print *,'        # of ',c_at(itype),' atoms on lattice'
!c     *                 ,' site ',isite,' = ',natm
                 call STRUCF(tcsf,xat,natm,real(nw))
!c                Add this part of the structure factor to the total ...
                  do 215, j=1,Numv
                   do 215, i=1,Numu
                    csf(i,j)=csf(i,j)+tcsf(i,j)
215               continue
                 end if
210          continue

!c            Subtract average (Bragg) scattering amplitude ...
              do 220, j=1,Numv
               do 220, i=1,Numu
                csf(i,j)=csf(i,j)-acsf(i,j)
220           continue

!c            Convert to intensity, add to total, and zero csf() ) ...
              do 230, j=1,Numv
               do 230, i=1,Numu                           
                 dsi(i,j)=dsi(i,j)+real( csf(i,j)*conjg(csf(i,j)) )
                 csf(i,j)=cmplx(0.d0,0.d0)
230           continue

!c          Save the diffuse intensity to disk ...
!c          Comment out if you don't want to save this often.
!c          There is a pretty big performance penalty on VP!!
!c            open(unit=1,file='intensity.bin',status='old',
!c     *           form='unformatted',access='direct',recl=irlen)
!c              xnorm=1./(real(ncell*nlot))
!c              do 240 j=1,Numv
!c               irec=numv*(nw-1)+j+5
!c               write(1,rec=irec)(dsi(i,j)*xnorm, i=1,Numu)
!c240           continue
!c            close(unit=1)

!c            Tell the file 'cur.txt' where we are ...
	if((nlot/50)*50.eq.nlot.or.nlots.lt.200) then
              open(unit=1,file='cur.txt',status='unknown')
               write(1,*)' Now at lot # ',nlot,'/',nlots
               write(1,*)' Now at plane # ',nw,'/',numw
              close(unit=1)
	endif
200        continue

!c          Finished doing 'nlots' regions on this one plane.

!c          Check to see if we did zero lots and do average I
	   if (nlots.eq.0) then
              write(*,*)'  This will only have the average in it!'
              do 300, j=1,Numv
               do 300, i=1,Numu
                dsi(i,j)=real( acsf(i,j)*conjg(acsf(i,j)) )
300           continue

              nlot=1
	      ncell=1
	   end if

!c          Save the diffuse intensity to disk ...
!c            open(unit=1,file='intensity.bin',status='old',

            open(unit=1,file=fname,status='old',        &
                form='unformatted',access='direct',recl=irlen)
              xnorm=1./(real(ncell*nlots))
              do 170 j=1,Numv
               irec=numv*(nw-1)+j+5
               write(1,rec=irec)(dsi(i,j)*xnorm, i=1,Numu)
170           continue
            close(unit=1)

100      continue

!c        Finished.

         print *,'  All done!'
         print *,' '
!c

!	errcnt=0
!556            open(unit=1,file='cur.count',status='unknown',err=555)
!               read(1,*) icount
!	       rewind 1
!	       icount=icount+1
!c
!               write(1,*) icount
!              close(unit=1)
!		goto 558
!555	continue
!	errcnt=errcnt+1
!	if(errcnt.lt.100) goto 556
!	print *, 'unable to write cur.count'
!c
!558	continue

      stop
  
      contains

        
        SUBROUTINE READINF

          ! DJG, AUgust 2009; Modified by CONTAINing within DIFFUSE


          !c      This routine reads in all the information necessary to run 
          !c      this program from the file 'diffuse.in' and passes all of this
          !c      information, through common blocks, back to the main routine
          !c      for use as it sees fit.

          !        parameter (MI=501,MJ=501)
          !        parameter (MT=10,MS=155)
          !        real ah_o(3),ah_u(3),ah_v(3),ah_w(3)
          !        real cell(6),uin(3),vin(3),win(3)
          !        integer csize(3),ls_xyz(3)
          !        real a1(MT),b1(MT),a2(MT),b2(MT),a3(MT),b3(MT)
          !        real a4(MT),b4(MT),c1(MT),fp(MT),fpp(MT)

          character*70 :: descrip

          !        character*4 c_at(MT)
          !        character*1 qalat,qbnd

          character*32 :: diffin
          integer :: i,k

          !common/sizes/ah_o,uin,vin,win,numu,numv,numw,stlmax
          !        !common/scatter/a1,b1,a2,b2,a3,b3,a4,b4,c1,fp,fpp
          !common/crystal/csize,cell
          !common/info1/ls_xyz,iseed,jseed,nlots,nsites,ntypes
          !common/info2/ah_u,ah_v,ah_w
          !common/info3/c_at,descrip,qalat,qbnd

          !c       Find out everything we need to know to run the program ...

          print *,' '
          print *,'DIFFUSE -  Ver. 2.1a  - May 11, 1995'
          print *,' '

          !c         open(unit=1,file='diffuse.in',status='old')
          write(6,*)'diffuse.in filename?'
          read(5,*)diffin
          open(unit=1,file=diffin,status='old')
          !c         Get the run description

          read(1,*)descrip
          !c         Random Number Seeds [ < 31328 and 30081 respectively ]

          read(1,*)iseed,jseed

          !c         Get the Lattice Parameter Info (Axial lengths, Cosines)

          read(1,*)(cell(i), i=1,6)        

          !c         Get the Crystal Size in Unit Cells

          read(1,*)(csize(i), i=1,3)

          !c         Does this crystal contain a periodic boundary?

          read(1,'(A1)')qbnd

          !c         Get the origin of the image computation

          read(1,*)(ah_o(i),i=1,3)

          !c         Get the maximum horizontal point and number of divisions

          read(1,*)(ah_u(i),i=1,3),numu

          !c         Get the maximum vertical point and Number of divisions

          read(1,*)(ah_v(i),i=1,3),numv

          !c         Get the maximum w-axis point and Number of divisions

          read(1,*)(ah_w(i),i=1,3),numw

          !           if ( (Numu.gt.MI).or.(Numv.gt.MJ) ) then
          !
          !             stop '  Your image is too big!'
          !
          !           end if

          !c         Get maximum sin(theta)/lambda 

          ! Allocate arrays that used to be dimensioned by MI, MJ


          allocate(csf(numu,numv),stat=ierr)
          if (ierr.ne.0) call allerr('csf(numu,numv)                    ')
          allocate(acsf(numu,numv),stat=ierr)
          if (ierr.ne.0) call allerr('acsf(numu,numv)                    ')
          allocate(tcsf(numu,numv),stat=ierr)
          if (ierr.ne.0) call allerr('tcsf(numu,numv)                    ')
          allocate(dsi(numu,numv),stat=ierr)
          if (ierr.ne.0) call allerr('dsi(numu,numv)                    ')
          allocate(istl(numu,numv),stat=ierr)
          if (ierr.ne.0) call allerr('istl(numu,numv)                    ')

          read(1,*)stlmax

          !c         Get the Lot Size (Size of Lots << csize() )

          read(1,*)(ls_xyz(i), i=1,3)

          !c         Get the Number of Lots

          read(1,*)nlots

          !c         How many Atom Site Positions are there per Cell? 

          read(1,*)nsites


          !           if (nsites.gt.MS) stop ' Too many sites!'

          !c         How many Atom Types are we to Deal With?  (<10)

          read(1,*)ntypes

          !           if (ntypes.gt.MT) stop ' Too many atom types!'
          !          allocate the needed arrays
          allocate(a1(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('a1(ntypes)                        ')
          allocate(a2(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('a2(ntypes)                        ')
          allocate(a3(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('a3(ntypes)                        ')
          allocate(a4(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('a4(ntypes)                        ')
          allocate(b1(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('b1(ntypes)                        ')
          allocate(b2(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('b2(ntypes)                        ')
          allocate(b3(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('b3(ntypes)                        ')
          allocate(b4(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('b4(ntypes)                        ')
          allocate(c1(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('c1(ntypes)                        ')
          allocate(fp(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('fp(ntypes)                        ')
          allocate(fpp(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('fpp(ntypes)                        ')
          allocate(c_at(ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('c_at(ntypes)                        ')
          allocate(Biso(nsites,ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('Biso(nsites,ntypes)                 ')
          allocate(wt(nsites,ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('wt(nsites,ntypes)                 ')
          allocate(ax(nsites,ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('ax(nsites,ntypes)                 ')
          allocate(ay(nsites,ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('ay(nsites,ntypes)                 ')
          allocate(az(nsites,ntypes),stat=ierr)
          if (ierr.ne.0) call allerr('az(nsites,ntypes)                 ')

          !c         Do we Want to Subtract an Average Lattice (Y/N) ?

          read(1,'(A1)')qalat
          !c         Read in the Scattering Factor Info ...

           do 10, k=1,ntypes

!c            Get these from Table 2.2B, p. 99, Vol. 3, Int. Tables

              read(1,*)c_at(k)
              read(1,*)a1(k),b1(k),a2(k),b2(k)
              read(1,*)a3(k),b3(k),a4(k),b4(k),c1(k)
!c            And the Dispersion Corrections

              read(1,*)fp(k),fpp(k)
10         continue
         close(unit=1)



!c        Compute the increments along the u,v,w directions

!c        This statement prevents a divide by zero on the vp2200

!*vocl loop,nopreex

          do 20, i=1,3

           uin(i)=0.d0

           vin(i)=0.d0

           win(i)=0.d0

           if (numu.gt.1) uin(i)= (ah_u(i)-ah_o(i))/real(numu-1)

           if (numv.gt.1) vin(i)= (ah_v(i)-ah_o(i))/real(numv-1)

           if (numw.gt.1) win(i)= (ah_w(i)-ah_o(i))/real(numw-1)

20        continue
      RETURN

      END subroutine READINF



      SUBROUTINE STRUCF(csf,xat,n,nw)

!c      Computes the complex stucture factor 'csf()' of n identical 
!c      atoms at the positions 'xat()' on the reciprocal plane 'nw'.
!c      This is the work-horse routine of the program DIFFUSE. Any
!c      real speed improvement will come from improving the inner loop
!c      of this subroutine.



!        parameter (MI=501,MJ=501,MAT=6000)
        integer, parameter :: MAT=6000
        !parameter (MAT=6000)

        !parameter (i2pi=2**14,MASK=2**14-1)
        integer, parameter  :: i2pi=2**14
        integer, parameter  :: MASK=2**14-1

!        complex csf(MI,MJ)

        complex ::  csf(:,:)

!        complex cex(0:MASK),cfact(0:1999)

!        real ah_o(3),uin(3),vin(3),win(3)

        real :: xm(3)
        real :: xat(MAT,3)

        double precision :: xarg0,xincu,xincv

!        integer istl(MI,MJ)

         real :: nw
         integer :: n,i,j,k,iarg0,iincu,iincv,iarg,iadd

        !common/sizes/ah_o,uin,vin,win,numu,numv,numw,stlmax

        !common/tabels/cex,istl,cfact

        

!c       Compute the origin of this reciprocal plane ...

         do 10, i=1,3

          xm(i)=ah_o(i)+(nw-1)*win(i)

10       continue



!c       Zero the complex scattering factor array 'csf()' ...

         do 20, j=1,Numv

          do 20, i=1,Numu

           csf(i,j)=cmplx(0.d0,0.d0)

20       continue



         if(n.eq.0) RETURN



!c       Loop over all of the atoms we are handling now ...

         do 100 k=1,n

!c          Get initial argument to the exponent and increments along

!c          the two axies 'u' and 'v'

            xarg0=  xm(1)*xat(k,1) +  xm(2)*xat(k,2) +  xm(3)*xat(k,3)

            xincu= uin(1)*xat(k,1) + uin(2)*xat(k,2) + uin(3)*xat(k,3)

            xincv= vin(1)*xat(k,1) + vin(2)*xat(k,2) + vin(3)*xat(k,3)



!c          Convert to high precision integers (64*i2pi=2^20) ...

            iarg0=nint( 64*i2pi*( xarg0-aint(xarg0)+1.d0 ) )

            iincu=nint( 64*i2pi*( xincu-aint(xincu)+1.d0 ) )

            iincv=nint( 64*i2pi*( xincv-aint(xincv)+1.d0 ) )

            iarg=iarg0



!c          Loop over all image pixels.  'iadd' is the address of the

!c          argument to the complex exponent (in the table 'cex()').

!c          The ISHFT operation divides out the 64 and the IAND

!c          is equivelent to a MOD and is used so that the argument to

!c          the complex exponent is inside our table which has range

!c          0=>2pi. 99.5% of the time the CPU will be busy with one of

!c          the four statements inside loop 210.

            do 200, j=1,Numv

              do 210, i=1,Numu

               iadd=ISHFT(iarg,-6)

               iadd=IAND(iadd,mask)

               csf(i,j)=csf(i,j)+cex(iadd)

               iarg=iarg+iincu

210           continue

              iarg=iarg0 + j*iincv

200         continue

100      continue

!c        We are through with all N atoms ...



!c        Multiply the complex scattering factor by the atomic

!c        scattering factor for this atom type and return ...

          do 40, j=1,Numv

            do 40, i=1,Numu

             csf(i,j)=csf(i,j)*cfact(istl(i,j))

40        continue

      RETURN

      END subroutine STRUCF



      SUBROUTINE WRITEINF

!c      Takes the information that was read from the file 'diffuse.in'

!c      and writes it back out to standard output so that a record of

!c      the run can be kept.



!        parameter (MT=10)

!        real ah_o(3),ah_u(3),ah_v(3),ah_w(3)

!        real cell(6),uin(3),vin(3),win(3)

!        integer csize(3),ls_xyz(3)

!        real a1(MT),b1(MT),a2(MT),b2(MT),a3(MT),b3(MT)

!        real a4(MT),b4(MT),c1(MT),fp(MT),fpp(MT)

!        character*70 descrip

!       character*4 c_at(MT)

!        character*1 qalat,qbnd



        !common/sizes/ah_o,uin,vin,win,numu,numv,numw,stlmax

        !common/scatter/a1,b1,a2,b2,a3,b3,a4,b4,c1,fp,fpp

        !common/crystal/csize,cell

        !common/info1/ls_xyz,iseed,jseed,nlots,nsites,ntypes

        !common/info2/ah_u,ah_v,ah_w

        !common/info3/c_at,descrip,qalat,qbnd



        print *,descrip

        print *,' '

        print *,'  The computation volume is defined by:'

        print 101,(ah_o(i), i=1,3),(ah_u(i), i=1,3)

        print 101,(ah_o(i), i=1,3),(ah_v(i), i=1,3)

        print 101,(ah_o(i), i=1,3),(ah_w(i), i=1,3)

        print *,'  Image size is ',numu,' X ',numv,' X ',numw

        print *,'  sin(theta)/lambda maximum = ',stlmax

        print *,'  Random number seeds are: ',iseed,jseed

        print *,'  Crystal size is: ',(csize(i), i=1,3)

        print *,'  Periodic boundary? ',qbnd

        print *,'  Number of atom sites per cell : ',nsites

        print *,'  Number of Atom types: ',ntypes

        print *,'  Lot Size is: ',(ls_xyz(i), i=1,3)

        print *,'  Number of Lots to Compute: ',nlots

        print *,'  Subtract Average Lattice? ',qalat

        print *,' '

101     format(5x,'(',2(f6.2,','),f6.2,') => (',2(f6.2,','),f6.2,')' )



      RETURN

      END subroutine WRITEINF


      SUBROUTINE RANLOC(csize,lbeg)

!c      Returns a pseudo-random cell from within the simulated crystal

!c      which is 'csize()' cells on edge.



        integer csize(3),lbeg(3)

        real xran(3)



        call rannum(xran,3)

        lbeg(1)=int(xran(1)*csize(1))+1

        lbeg(2)=int(xran(2)*csize(2))+1

        lbeg(3)=int(xran(3)*csize(3))+1



      RETURN

      END subroutine RANLOC



      SUBROUTINE DEBYE(cfact,Biso,wx)

!c      Multiplies a table of scattering lengths 'cfact()' by a 

!c      static displacement Debye factor and by a weighting factor

!c      'wx' which is usually a site occupancy fraction.



        !parameter (xinc=0.001)
        real, parameter :: xinc=0.001
        complex :: cfact(0:1999)
     
        real :: wx, B, Biso, stl, stl2

        integer :: i


!c       Loop over the whole tabel cfact() ...

         B=-1.*Biso

         do 10, i=0,1999

          stl=float(i)*xinc

          stl2=stl**2

          cfact(i)=wx*cfact(i)*exp(B*stl2)

10       continue



      RETURN

      END subroutine DEBYE




      SUBROUTINE SCATF(cfact,itype)

!c      This routine computes the complex atomic scattering factor as
!c      a function of sin(theta)/lambda for the element itype given
!c      the parameters in the common block.  If this has already been
!c      computed for this atom type then this routine just returns the
!c      previously computed values.  (That is what the extra array
!c      cftab() holds and what the test is all about.)


        !parameter (xinc=0.001)
         real, parameter :: xinc=0.001
!        parameter (MT=10)

        complex ::  cfact(0:1999)

        real :: stl, stl2, sf, sfp, sfpp

        integer :: i, itype

!        real a1(MT),b1(MT),a2(MT),b2(MT),a3(MT),b3(MT)
!        real a4(MT),b4(MT),c1(MT),fp(MT),fpp(MT)

!        complex cftab(0:1999,MT)
!        complex, allocatable :: cftab(:.:)




        !common/scatter/a1,b1,a2,b2,a3,b3,a4,b4,c1,fp,fpp

      !  data itest/10*0/

      !  save cftab,itest ! just make cftab global so no need to save.
      !  save itest

!c       Test if I have already computed it for this element 

         if (itest(itype).eq.1234) then

           do 5, i=0,1999

            cfact(i)=cftab(i,itype)

5          continue

          else

!c          I haven't done this element yet so compute it ...

            do 10, i=0,1999

             stl=float(i)*xinc

             stl2=stl**2

             sf =      a1(itype)*exp(-1.*b1(itype)*stl2)

             sf = sf + a2(itype)*exp(-1.*b2(itype)*stl2)

             sf = sf + a3(itype)*exp(-1.*b3(itype)*stl2)

             sf = sf + a4(itype)*exp(-1.*b4(itype)*stl2) + c1(itype)

             sfp=fp(itype)

             sfpp=fpp(itype)

             cfact(i)=cmplx(sf+sfp,sfpp)

             cftab(i,itype)=cfact(i)

10          continue

!c           Now I have computed it so let future call know this ...

             itest(itype)=1234

         end if



      RETURN

      END subroutine SCATF


                                  

      SUBROUTINE STLTAB(istl,cell,nw)

!c      Calculates the value of sin(theta)/lambda for each element
!c      of the diffuse scattering array.  The cell parameters (a,b,c,
!c      cos(bc),cos(ac),cos(ab)), computation spacings and dimmensions,
!c      and the reciprocal section must be provided.  This table is
!c      used for quick lookup of the scattering factor curves.
!c      The output table is integer with istl=nint(stl*1000).



!        parameter (MI=501,MJ=501)
!        integer istl(MI,MJ)
!        real ah_o(3),cell(6)
!        real uin(3),vin(3),win(3)

         integer :: istl(:,:)
         real :: cell(6)

        !common/sizes/ah_o,uin,vin,win,numu,numv,numw,stlmax

         real :: a1,a2,a3,c1,c2,c3,s1,s2,s3,V,S11,S22,S33,S12,S23,S13
         real :: xh1,xh2,xh3,stl,nw

!c       The cell parameters ...

         a1=cell(1) 

         a2=cell(2) 

         a3=cell(3)

         c1=cell(4)

         c2=cell(5)

         c3=cell(6)

         s1=sin(acos(c1))

         s2=sin(acos(c2))

         s3=sin(acos(c3))



!c       I got this out of Culity ...

         V=a1*a2*a3*sqrt(1.-c1**2-c2**2-c3**2+2.*c1*c2*c3)

         S11=(a2**2)*(a3**2)*(s1**2)

         S22=(a1**2)*(a3**2)*(s2**2)

         S33=(a1**2)*(a2**2)*(s3**2)

         S12=(a1)*(a2)*(a3**2)*(c1*c2-c3)

         S23=(a2)*(a3)*(a1**2)*(c2*c3-c1)

         S13=(a1)*(a3)*(a2**2)*(c1*c3-c2)



!c      Loop over the all elements in the output image ...

        do 10, jj=1,Numv

          do 10, ii=1,Numu

           xh1=ah_o(1)+real(ii-1)*uin(1)+real(jj-1)*vin(1)+nw*win(1)

           xh2=ah_o(2)+real(ii-1)*uin(2)+real(jj-1)*vin(2)+nw*win(2)

           xh3=ah_o(3)+real(ii-1)*uin(3)+real(jj-1)*vin(3)+nw*win(3)

           stl=sqrt(S11*xh1**2+S22*xh2**2+S33*xh3**2       &
                   +2.*S12*xh1*xh2+2.*S13*xh1*xh3+2.*S23*xh2*xh3)

           stl=0.5*stl/V

           istl(ii,jj)=nint(stl*1000.)

           if (istl(ii,jj).gt.1999) then

             write(*,*)ii,jj,istl(ii,jj)

             stop '  sin(theta)/lambda is greater than 2!'

           end if

10      continue



      RETURN

      END subroutine STLTAB


      SUBROUTINE CEXPT(cex)

!c      This routine computes a table of exp(i*2pi*delta) that is used
!c      to save time in computing trig functions in the routine STRUCF.
!c      The spacing of the values in the table is 1/2**N where N is 
!c      chosen to compromise between precision and the size of the 
!c      table.  Too large of a table size will slow down the computation
!c      considerably.



        !parameter (i2pi=2**14,MASK=2**14-1)
        integer, parameter  :: i2pi=2**14
        integer, parameter  :: MASK=2**14-1

        complex cex(0:MASK)

        double precision twopi,xmult,xarg,xt



         xt=1.d0/i2pi

         twopi=8.d0*datan(1.d0)

         do 10, i=0,MASK

           xmult=real(i)*xt

           xarg=twopi*xmult

           cex(i)=cmplx( cos(xarg),sin(xarg) )

10       continue



      RETURN

      END subroutine CEXPT


      SUBROUTINE GETATM(xat,natm,lbeg,csize,ls_xyz,isite,cat,ncell)

!c       Composes a list of 'natm' atom positions 'xat()' for a site 
!c       position 'isite' and atom label 'cat' that lie inside inside 
!c       unit cells with centers inside an elipsoid with major axes 
!c       (along the crysallographic (a,b,c) directions) of ls_xyz unit 
!c       cells.  The ellipsoid center is placed using the variable 
!c       'lbeg()'.  The crystal size 'csize()' is used to avoid over-
!c       running the array and the number of unit cells 'ncell' 
!c       contained in the elipsoid is also returned so that we know
!c       what to normalize to later.



        !parameter (MAT=6000)
        integer, parameter :: MAT=6000

        real :: xat(MAT,3)

        integer :: lbeg(3),csize(3),ls_xyz(3),natm

        integer :: l(3)

        character*4 :: cat

         real :: x01, x02, x03, xtest1, xtest2, xtest3, x, y, z, xtest
         integer :: ncell, ii, jj, kk, isite

!c       Initialize a few variables ...

         x01=real(ls_xyz(1))/2.

         x02=real(ls_xyz(2))/2.

         x03=real(ls_xyz(3))/2.

         natm=0

         ncell=0



!c       Loop over all cells that might be inside the ellipsoid ...

         do 100, kk=0,ls_xyz(3)-1

           l(3)=kk+lbeg(3)

           if(l(3).gt.csize(3)) l(3)=l(3)-csize(3)

           xtest3=( real(kk)-x03+0.5 )**2/x03**2

           do 100, jj=0,ls_xyz(2)-1

            l(2)=jj+lbeg(2)

            if(l(2).gt.csize(2)) l(2)=l(2)-csize(2)

            xtest2=( real(jj)-x02+0.5 )**2/x02**2

            do 100, ii=0,ls_xyz(1)-1

             l(1)=ii+lbeg(1)

             if(l(1).gt.csize(1)) l(1)=l(1)-csize(1)

             xtest1=( real(ii)-x01+0.5 )**2/x01**2



!c            Is it inside the ellipsoid?

              xtest=xtest1+xtest2+xtest3

              if (xtest.le.1.) then

                ncell=ncell+1

                call READAT(l,isite,cat,inum,x,y,z)

                if(inum.eq.1) then

                  natm=natm+1

                  xat(natm,1)=real(ii)+x

                  xat(natm,2)=real(jj)+y

                  xat(natm,3)=real(kk)+z

                 else if (inum.ne.0) then

                  stop 'The variable inum from readat must be 0 or 1'

                end if

              end if

100      continue



      RETURN

      END subroutine GETATM


      SUBROUTINE GETAV(xat,natm,ls_xyz,xx,yy,zz)

!c       Composes a list of atom positions 'xat()' for a site position
!c       (xx,yy,zz).  There are 'natm' positions returned and these lie
!c       inside unit cells with centers inside an elipsoid with major
!c       axes (along the crysallographic (a,b,c) directions) of ls_xyz
!c       unit cells.



        ! parameter (MAT=6000)
        integer, parameter :: MAT=6000

         real :: xat(MAT,3)

         integer :: ls_xyz(3)

         real :: x01, x02, x03, xtest1, xtest2, xtest3, xtest, xx, yy, zz

         integer :: natm, ii, jj, kk

         x01=real(ls_xyz(1))/2.

         x02=real(ls_xyz(2))/2.

         x03=real(ls_xyz(3))/2.

         natm=0



!c        Loop over all cells that might be in the elipsoid but only
!c        keep those cells that have a center inside the elipsoid.

          do 20, kk=0,ls_xyz(3)-1

            xtest3=( real(kk)-x03+0.5 )**2/x03**2

            do 20, jj=0,ls_xyz(2)-1

             xtest2=( real(jj)-x02+0.5 )**2/x02**2

             do 20, ii=0,ls_xyz(1)-1

              xtest1=( real(ii)-x01+0.5 )**2/x01**2

              xtest=xtest1+xtest2+xtest3

              if (xtest.le.1.) then

               natm=natm+1

               xat(natm,1)=real(ii)+xx

               xat(natm,2)=real(jj)+yy

               xat(natm,3)=real(kk)+zz

              end if

20        continue

      RETURN

      END subroutine GETAV


      SUBROUTINE AVINFO(isite,cat,xx,yy,zz,wx,Biso)

!c       Takes a look at all of the unit cells in the simulated
!c       crystal and computes the average position (xx,yy,zz) of
!c       an atom labeled 'cat' on the lattice site 'isite'.  It
!c       also returns the occupation fraction 'wx' and an
!c       isotropic static displacement Debye factor 'Biso'


         integer :: l(3)

         double precision :: sumx,sumy,sumz,sumua,sumub,sumuc
         double precision :: sumaa,sumab,sumac

         real :: a1,a2,a3,c1,c2,c3,s3,ugh,xpi2,x,y,z,wx,denom,xx,yy,zz 
         real :: Ba, Bb, Bc, Biso

         character*4 :: cat

         integer :: ntot, kk, jj, ii, isite, inum


         !common/crystal/csize,cell

!c        These are the cell dimensions and cosines

          a1=cell(1)
          a2=cell(2)
          a3=cell(3)
          c1=cell(4)
          c2=cell(5)
          c3=cell(6)

	  s3=sqrt(1.-c3**2)

	  ugh=(1.-c2**2+c1**2*s3**2)



!c        Define some constants (using double precision for variables
!c        that loop real sums) ...

          xpi2=(4.*atan(1.))**2

          ntot=0

          sumx=0.d0

          sumy=0.d0

          sumz=0.d0

          sumua=0.d0

          sumub=0.d0

          sumuc=0.d0


!c         Loop over all cells and do the sums ...

           do 10, kk=1,csize(3)

            l(3)=kk

            do 10, jj=1,csize(2)

             l(2)=jj

             do 10, ii=1,csize(1)

              l(1)=ii

!c             This is a call to the user supplied subroutine.  We pass
!c             it a cell 'l()', site 'isite', and atom description 'cat'
!c             and it tells us if there is such an atom there (i.e.
!c             inum=1) and if so where (x,y,z).
               call READAT(l,isite,cat,inum,x,y,z)

               if (inum.eq.1) then

                 ntot=ntot+1

                 sumx=sumx+x

                 sumy=sumy+y

                 sumz=sumz+z

                 sumua=sumua+( a1*x +a2*y*c3 + a3*z*c2 )**2

                 sumub=sumub+( a2*y*s3 + a3*z*c1*s3 )**2

                 sumuc=sumuc+(a3*z)**2*ugh

                else if (inum.ne.0) then

                 stop 'The variable inum from readat must be 0 or 1'

               end if

10         continue



!c         Normalize, compute Debye factor, and get outa here ...

           wx=real(ntot)/real(csize(1)*csize(2)*csize(3))

	   if (wx*1000.gt.1) then

	    denom=1./real(ntot)

            xx=sumx*denom

            yy=sumy*denom

            zz=sumz*denom

            sumua=sumua*denom

            sumub=sumub*denom

            sumuc=sumuc*denom

            sumaa=( a1*xx +a2*yy*c3 + a3*zz*c2 )**2

            sumab=( a2*yy*s3 + a3*zz*c1*s3 )**2

            sumac=(a3*zz)**2*ugh

            Ba=8.*xpi2*(sumua-sumaa)

            Bb=8.*xpi2*(sumub-sumab)

            Bc=8.*xpi2*(sumuc-sumac)

            Biso=(Ba+Bb+Bc)/3.

	   end if

      RETURN

      END subroutine AVINFO



  subroutine allerr(xyzfile)          

    character(len=*)                             :: xyzfile

    write(stderr,*)'----------------------------------------------------------------------------'
    write(stderr,*) 'Problem allocating memory to ',trim(xyzfile)
    write(stderr,*) 'Do you have enough for the simulation size you want?'
    write(stderr,*) 'Exiting.'
    write(stderr,*)'----------------------------------------------------------------------------'

    stop

  end subroutine allerr


      END program DIFFUSE 
