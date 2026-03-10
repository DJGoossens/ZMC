c
        subroutine rseed(ij,kl)

c this is the initialization routine for the random number 
c generator rannum() and must be called once prior to rannum().
c
c note: the seed variables must have values between: 0 <= ij <= 31328
c                                                    0 <= kl <= 30081
      real u(97), c, cd, cm
      integer i97, j97
      
      common /raset1/ u, c, cd, cm, i97, j97
      save
      
      if( ij .lt. 0  .or.  ij .gt. 31328  .or.
     *    kl .lt. 0  .or.  kl .gt. 30081 ) then
          print '(a)', ' the first random number seed must have a value 
     *between 0 and 31328'
          print '(a)',' the second seed must have a value between 
     *0 and 30081'
            stop
      endif

      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169) 

      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
3        continue
         u(ii) = s
2     continue

      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0

      i97 = 97
      j97 = 33

      return
      end

      subroutine rannum(rvec, len)
      
c this is a random number generator proposed by george marsaglia in 
c florida state university report: fsu-scri-87-50
c it was slightly modified by f. james to produce an array of 
c pseudorandom numbers.  

      real rvec(*)
      real u(97), c, cd, cm
      integer i97, j97
      integer ivec
      
      common /raset1/ u, c, cd, cm, i97, j97  
      save
 
      do 100 ivec = 1, len
         uni = u(i97) - u(j97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         u(i97) = uni
         i97 = i97 - 1
         if(i97 .eq. 0) i97 = 97
         j97 = j97 - 1
         if(j97 .eq. 0) j97 = 97
         c = c - cd
         if( c .lt. 0.0 ) c = c + cm
         uni = uni - c
         if( uni .lt. 0.0 ) uni = uni + 1.0
         rvec(ivec) = uni
100   continue
      return
      end
c----------------------------------
c-----------------------gauss routine-------------------------
c---------------------------------------------------------
        subroutine gaussm(ix,iy,s,am,v)
        real :: x(1),y(1)
        a=0.
        do 50 i=1,12
        call rannum(y(1),1)
50        a=a+y(1)
        v=(a-6.)*s+am
        return
        end
c
        subroutine random(i1,i2,x)
        dimension a(100)
        real :: x(1),y(1)
        data init/0/
        if(init.ne.0) goto 10
        do 5 n=1,100
        call rannum(x(1),1)
5        a(n)=x(1)
        init=-1
10        call rannum(x(1),1)
        call rannum(y(1),1)
        n=int(99.99*x(1))+1
        x(1)=a(n)
        a(n)=y(1)
        return
        end
