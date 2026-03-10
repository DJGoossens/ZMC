! make_random_occ.f90
! program to do just that for a specified set of locations, zmats and so on.
!------------------------------------------------------------------------------!
module myparams
  !Integer Parameters!
  integer, parameter:: mol_max = 4     !max number of molecules of each z-mat
  integer, parameter:: amax = 32       !max number of unit cells along a dirn
  integer, parameter:: bmax = 32       !max number of unit cells along b dirn
  integer, parameter:: cmax = 32       !max number of unit cells along c dirn
  integer, parameter:: zm_max = 6  
  integer, parameter:: loc_max = 20 
end module myparams


!------------------------------------------------------------------------------!
!          Program proper starts here                                          !
!------------------------------------------------------------------------------!
program make_random_occ

  use cmdline_arguments, only: have_args, next_arg
  Use file_functions
  use myparams

  implicit none

  !------------------------------------------------------------------------------!
  !          Variables declared by Darren                                        !
  !------------------------------------------------------------------------------!


  !Characters!
  character(len=30)                   ::  outfile
  character(len=1)                    ::  response

  !Integers!
  integer                              :: dunit, asize, bsize, csize, num_loc, num_zm, il,iz,ia,ib,ic,iiz
  integer,dimension(loc_max,zm_max)    :: how_many_instances


  !Reals!
  real                                 :: total_prob
  real                                 :: rvar
  real, dimension(loc_max,zm_max,mol_max) :: prob

  !Custom Types!

  call rseed(8739,7635)

  write(stdout,*)'How big is the model crystal? (asize bsize csize)'
  read(stdin,*)asize,bsize,csize
  write(stdout,*)'How many locations in the unit cell?'
  read(stdin,*)num_loc
  write(stdout,*)'How many types of zmatrix in the structure?'
  read(stdin,*)num_zm

  do il = 1,num_loc
     do iz = 1,num_zm
        write(stdout,'("How many instances of zmatrix ",i3," on location ",i3)')iz,il
        read(stdin,*)how_many_instances(il,iz)
     end do
  end do

  do il = 1, num_loc
     total_prob = 0.0
     do iz = 1, num_zm
        do iiz = 1,how_many_instances(il,iz)
           write(stdout,'("What is the probability of getting instance ",i3)')iiz
           write(stdout,'("of zmatrix ",i3," on location ",i3)')iz,il
           read(stdin,*)prob(il,iz,iiz)
           total_prob = total_prob + prob(il,iz,iiz)
        end do
     end do
     if(abs(1-total_prob).gt.0.000001) then
        write(stderr,*)'Probabilities do not add to 1 for location ',il
        write(stderr,*)'Exiting.'
        stop
     end if
  end do
  
  write(stdout,*)'What is output filename?'
  read(stdin,*)outfile
  
  dunit = open(outfile, status='new')
  if(dunit.lt.1) then
     write(stdout,*)'Ouput file alrady exists.  Overwrite?'
11   write(stdout,*)'y or n'
     read(stdin,*)response
     if (response.eq.'n') stop
     if (response.eq.'y') then
        dunit = open(outfile)
     else
        goto 11
     end if
  end if
     
  do ia = 1,asize
     do ib = 1,bsize
        do ic = 1,csize
           do il = 1,num_loc
              ! So we are on location l.  Let's generate a random number
              call rannum(rvar,1)
              ! The total probability on il is 1
              total_prob = 0.0
              zmloop:  do iz = 1, num_zm
                 do iiz = 1,how_many_instances(il,iz)
                    total_prob = total_prob + prob(il,iz,iiz)
                    if(rvar.lt.total_prob) then
                       write(dunit,*)ia,ib,ic,il,iz,iiz
                       exit zmloop   
                    end if
                 end do
              end do zmloop
           end do
        end do
     end do
  end do
  
  close(dunit)
  write(stdout,*)'Done.'

end program make_random_occ

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
