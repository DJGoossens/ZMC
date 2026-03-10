module simpose

  use quaternion_class
  use rotmatrix_class

  implicit none

contains

  real function superimpose(fixed, moving ,quat, trans)                                                    

    ! Based on the program pdbsup by B.Rupp and S.Parkin (1996), itself based
    ! an published method by S.K.Kearsley, Acta Cryst. A45, 208 (1989)       
    ! ----------------------------------------------------------------------
    ! Determines rotation matrix and translation vector for best fit 
    ! superimposition of two coordinate lists by solving the quaternion 
    ! eigenvalue problem.
    ! ----------------------------------------------------------------------

    ! Coordinates of fixed (xf) and moving (xm) coordinate lists
    ! real, dimension(:,:) :: xf,xm
    real, dimension(:,:), intent(in) :: fixed, moving

    type(quaternion), intent(out) :: quat
    real, intent(out) :: trans(3)

    real(8) :: xtemp(3)

    ! real, allocatable, dimension(:,:) :: dxp, dxm ,dxpsq ,dxmsq
    real, dimension(size(fixed,1),size(fixed,2)) :: xf, xm
    real(8), dimension(size(fixed,1),size(fixed,2)) :: dxp, dxm ,dxpsq ,dxmsq
    ! real, allocatable, dimension(:) :: rmsd
    real(8), dimension(size(fixed,2)) :: rmsd

    real(8) :: dm(4), vm(4,4), cm(3), cf(3)
    real(8) :: t(3,3), q(4,4)
    integer i, j, n, ns, nmrot
    
    real :: rotcm(3)

    type(rotmatrix) :: rmat

    ! Make local copies of the coordinates -- make explicit cast
    ! to double precision reals
    xf = fixed
    xm = moving

    ! Initialise the quaternion matrix
    q = 0.0
    
    ! Sum up all coordinates (in dble precision) to find centres
    cm = sum(xm,DIM=2)/size(moving,DIM=2)
    cf = sum(xf,DIM=2)/size(fixed,DIM=2)

    ! write(*,'(/a,3f8.3)')' Centre of target molecule  =',(cf(i),i=1,3)
    ! write(*,'(a,3f8.3)') ' Centre of moving molecule  =',(cm(i),i=1,3)
    ! write(*,'(a,3f8.3/)')' T - vector probe -> target =',(trans(i),i=1,3)
    
    ! Make coordinates relative to centre of mass
    do i=1,3                                                       
       xf(i,:)=xf(i,:)-real(cf(i))
       xm(i,:)=xm(i,:)-real(cm(i))
    end do

    ! Create coordinate differences delta x plus (dxp) and minus (dxm)  
    dxm=xm-xf
    dxp=xm+xf
    dxmsq=(xm-xf)**2
    dxpsq=(xm+xf)**2
    
    ! --- fill upper triangle of (symmetric) quaternion matrix --           

    q = 0.d0

    ! ---    diags are sums of squared cyclic coordinate differences        
    q(1,1)=sum(dxmsq)
    q(2,2)=sum(dxpsq(2:3,:))+sum(dxmsq(1,:))
    q(3,3)=sum(dxpsq(1:3:2,:))+sum(dxmsq(2,:))
    q(4,4)=sum(dxpsq(1:2,:))+sum(dxmsq(3,:))

    ! ---    cross differences                                              
    q(1,2)=sum(dxp(2,:)*dxm(3,:))-sum(dxm(2,:)*dxp(3,:))             
    q(1,3)=sum(dxm(1,:)*dxp(3,:))-sum(dxp(1,:)*dxm(3,:))             
    q(1,4)=sum(dxp(1,:)*dxm(2,:))-sum(dxm(1,:)*dxp(2,:))             
    q(2,3)=sum(dxm(1,:)*dxm(2,:))-sum(dxp(1,:)*dxp(2,:))             
    q(2,4)=sum(dxm(1,:)*dxm(3,:))-sum(dxp(1,:)*dxp(3,:))             
    q(3,4)=sum(dxm(2,:)*dxm(3,:))-sum(dxp(2,:)*dxp(3,:))              

    ! --- fill the rest by transposing it onto itself                       
    ! call trpmat(4,q,q)
    where(q == 0) q = q + transpose(q)

    ! write (*,'(/a)')                                                  &
    !      &     '       q(1)         q(2)         q(3)        q(4)'          
    ! do i=1,4                                                          
    !    write(*,'(4e13.5)') (q(i,j),j=1,4)                             
    ! end do
    !

    n = 4
    ns = 4

    ! --- orthogonalization by jacobi rotation = solution of EV -problem -- 
    call jacobi(q,n,ns,dm,vm,nmrot)                                   

    ! --- sort eigenvectors after eigenvalues, descending --                
    ! write (*,'(a/)')' Sorting eigenvalues/vectors .......'            
    call eigsrt(dm,vm,n,ns)                                           

    ! write (*,'(a,i2,a)')' Eigenvalues and Eigenvectors (',            &
    !      & nmrot,' Jacobi rotations)'                                       
    ! write (*,'(a)') '      e(1)        e(2)        e(4)        e(4)'  
    ! write (*,'(4e12.5,i5)') (dm(i),i=1,4)                             
    ! write (*,'(a)') '      ev(1)       ev(2)       ev(3)       ev(4)' 
    ! do i=1,4                                                          
    !    write(*,'(4f12.6)') (vm(i,j),j=1,4)                            
    ! end do
    
    ! --- the smallest eigenvector contains best fit srs                    
    ! rmsd=sqrt(abs(dm(4)/imov))                                        
    ! write(*,'(/a/)')                                                  &
    !      & ' The smallest eigenvalue represents s.r.s. of best fit'         
    ! write(*,'(a)')                                                    &
    !      & ' Constructing the best fit rotation matrix from associated'     
    ! write(*,'(a/)') ' eigenvector elements (last column).....'        
    
    !
    ! --- fill the rotation matrix which is made of elements from 4th EV    
    t(1,1)=vm(1,4)**2+vm(2,4)**2-vm(3,4)**2-vm(4,4)**2                
    t(2,1)=2*(vm(2,4)*vm(3,4)+vm(1,4)*vm(4,4))                        
    t(3,1)=2*(vm(2,4)*vm(4,4)-vm(1,4)*vm(3,4))                        
    t(1,2)=2*(vm(2,4)*vm(3,4)-vm(1,4)*vm(4,4))                        
    t(2,2)=vm(1,4)**2+vm(3,4)**2-vm(2,4)**2-vm(4,4)**2                
    t(3,2)=2*(vm(3,4)*vm(4,4)+vm(1,4)*vm(2,4))                        
    t(1,3)=2*(vm(2,4)*vm(4,4)+vm(1,4)*vm(3,4))                        
    t(2,3)=2*(vm(3,4)*vm(4,4)-vm(1,4)*vm(2,4))                        
    t(3,3)=vm(1,4)**2+vm(4,4)**2-vm(2,4)**2-vm(3,4)**2                
    
    ! print *,'Rotation matrix'
    ! do i=1,3                                                          
    !    write(*,'(3f11.5)') (t(i,j),j=1,3)                             
    ! end do
    ! Test we get the same result from quaternion module
    quat = real(vm(:,4))
    ! print *,as_array(quat)
    rmat = as_rotmatrix(quat)
    ! print *,as_matrix(rmat)
    ! --- xm and xf are translated                                      
    ! print *,'rmsds'

    ! Rotate the centre of mass of the moving molecule to put it
    ! in the same 'frame' as the fixed molecule
    rotcm = real(cm)
    call rotate(rmat,rotcm)

    ! Find net translation vector required to put moving coordinates
    ! at the same position as the fixed coordinates
    trans = real(cf) - rotcm ! real(cf - cm)
    
    do i=1,size(moving,2)
       call rotate(rmat,xm(:,i))
       rmsd(i)=sqrt(sum((xf(:,i)-xm(:,i))**2))
       ! write(*,'(I4F11.5)') i,rmsd(i)
       ! write(*,'(2(3F11.5,3X))') xf(:,i),xm(:,i)
    end do

    ! superimpose = sum(rmsd)/size(rmsd)
    superimpose = sum(rmsd)

  end function superimpose

  subroutine trpmat(n,t,tr)                                         
    ! --- transpose matrix -------------------------------------------------
    real(8) :: t(n,n), tr(n,n)                                              
    integer i,j,n
    do i=1,n                                                          
       do j=1,n                                                       
          tr(j,i)=t(i,j)                                              
       end do
    end do
    return                                                            
  end subroutine trpmat
  !
  subroutine rotvec (n,v,t)                                         
    ! --- multiply vector with matrix --------------------------------------
    real t(n,n), v(n),s(n)                                            
    integer i,j,n
    !
    do i=1,n                                                          
       s(i)=v(i)                                                      
       v(i)=0.0                                                       
    end do
    do i=1,n                                                          
       do j=1,n                                                       
          v(i)=v(i)+s(j)*t(i,j)                                       
       end do
    end do
    return                                                            
  end subroutine rotvec
  !
  SUBROUTINE eigsrt(d,v,n,np)                                       
    ! ----------------------------------------------------------------------
    INTEGER n,np                                                      
    REAL(8) :: d(np),v(np,np)                                               
    INTEGER i,j,k                                                     
    REAL(8) :: p                                                            
    do i=1,n-1                                                     
       k=i                                                             
       p=d(i)                                                          
       do j=i+1,n                                                   
          if(d(j).ge.p)then                                             
             k=j                                                         
             p=d(j)                                                      
          endif
       end do
       if(k.ne.i)then                                                  
          d(k)=d(i)                                                     
          d(i)=p                                                        
          do j=1,n                                                   
             p=v(j,i)                                                    
             v(j,i)=v(j,k)                                               
             v(j,k)=p                                                    
          end do
       endif
    end do
    return                                                            
  end subroutine eigsrt
  !  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             
  !
  SUBROUTINE jacobi(a,n,np,d,v,nrot)                                
    ! ----------------------------------------------------------------------
    !     modified from numerical recipes book                              
    !     one needs to set the threshold for sm from sm.eq.0 to sm.lt.10E-30
    !     (anything in this range would be ok) due to underflow errors on   
    !     some computers/compilers.                                         
    ! ----------------------------------------------------------------------
    INTEGER nmax
    PARAMETER (nmax=500)                                              
    !
    INTEGER n,np,nrot                                                 
    REAL(8) ::  a(np,np),d(np),v(np,np)                                      
    INTEGER i,ip,iq,j,maxrot                                          
    REAL(kind=8) :: c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax),zero          
    !
    ! --- zero set and iteration maximum                                    
    zero=10E-30                                                       
    maxrot=50                                                         
    !
    do ip=1,n                                                      
       do iq=1,n                                                    
          v(ip,iq)=0.                                                   
       end do
       v(ip,ip)=1.                                                     
    end do
    do ip=1,n                                                      
       b(ip)=a(ip,ip)                                                  
       d(ip)=b(ip)                                                     
       z(ip)=0.                                                        
    end do
    nrot=0                                                            
    do i=1,maxrot                                                  
       sm=0.                                                           
       do ip=1,n-1                                                  
          do iq=ip+1,n                                               
             sm=sm+abs(a(ip,iq))                                         
          end do
       end do
       ! ---   modified convergence threshold ---                              
       if(sm.lt.zero)return                                            
       if(i.lt.4)then                                                  
          tresh=0.2*sm/n**2                                             
       else                                                            
          tresh=0.                                                      
       endif
       do ip=1,n-1                                                  
          do iq=ip+1,n                                               
             g=100.*abs(a(ip,iq))                                        
             if((i.gt.4).and.(abs(d(ip))+                                &
                  &g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then            
                a(ip,iq)=0.                                               
             else if(abs(a(ip,iq)).gt.tresh)then                         
                h=d(iq)-d(ip)                                             
                if(abs(h)+g.eq.abs(h))then                                
                   t=a(ip,iq)/h                                            
                else                                                      
                   theta=0.5*h/a(ip,iq)                                    
                   t=1./(abs(theta)+sqrt(1.+theta**2))                     
                   if(theta.lt.0.)t=-t                                     
                endif
                c=1./sqrt(1+t**2)                                         
                s=t*c                                                     
                tau=s/(1.+c)                                              
                h=t*a(ip,iq)                                              
                z(ip)=z(ip)-h                                             
                z(iq)=z(iq)+h                                             
                d(ip)=d(ip)-h                                             
                d(iq)=d(iq)+h                                             
                a(ip,iq)=0.                                               
                do j=1,ip-1                                            
                   g=a(j,ip)                                               
                   h=a(j,iq)                                               
                   a(j,ip)=g-s*(h+g*tau)                                   
                   a(j,iq)=h+s*(g-h*tau)                                   
                end do
                do j=ip+1,iq-1                                         
                   g=a(ip,j)                                               
                   h=a(j,iq)                                               
                   a(ip,j)=g-s*(h+g*tau)                                   
                   a(j,iq)=h+s*(g-h*tau)                                   
                end do
                do j=iq+1,n                                            
                   g=a(ip,j)                                               
                   h=a(iq,j)                                               
                   a(ip,j)=g-s*(h+g*tau)                                   
                   a(iq,j)=h+s*(g-h*tau)                                   
                end do
                do j=1,n                                               
                   g=v(j,ip)                                               
                   h=v(j,iq)                                               
                   v(j,ip)=g-s*(h+g*tau)                                   
                   v(j,iq)=h+s*(g-h*tau)                                   
                end do
                nrot=nrot+1                                               
             endif
          end do
       end do
       do ip=1,n                                                    
          b(ip)=b(ip)+z(ip)                                             
          d(ip)=b(ip)                                                   
          z(ip)=0.                                                      
       end do
    end do
    stop 'too many iterations in jacobi'                             
    return                                                            
  END SUBROUTINE jacobi
  !  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             
end module simpose
