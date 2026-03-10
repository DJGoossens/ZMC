c
      subroutine ps_init(fname)
c      this is the initialization routine for the postscript output.
c      it opens the file 'fname' as unit 13, and sets the current
c      pen position to (0.,0.).  you must call this routine before any
c      other ps_ routines.

       character*12 fname
       character*72 ct

        open(unit=13,file=fname,status='unknown')
        write(13,'(a)')'%!ps-adobe-3.0 epsf-3.0'
        write(13,'(a)')'%%boundingbox: 0 0 612 792'
        write(13,'(a)')'%%creator ps_subs.f, b.d. butler'
        write(13,'(a)')'%%creationdate: mm/dd/yy'
        write(13,'(a)')'%%title: whatever'
        write(13,'(a)')'%%for: whoever'
        write(13,'(a)')'%%endcomments'
        write(13,'(a)')' '
        ct='/psline {newpath moveto lineto closepath stroke} def'
        write(13,'(a)')ct
        ct='/pscirc {newpath 0 360 arc closepath stroke} def'
        write(13,'(a)')ct
        ct='/pscircf {newpath 0 360 arc closepath fill} def'
        write(13,'(a)')ct
        write(13,'(a)')'/psrect {newpath moveto lineto lineto lineto'
        write(13,'(a)')'         closepath stroke} def'
        write(13,'(a)')'/psrectf{newpath moveto lineto lineto lineto'
        write(13,'(a)')'         closepath fill} def'
        write(13,'(a)')'/psbox {gsave translate rotate scale newpath '
        write(13,'(a)')'        -0.5 -0.5 moveto 0.5 -0.5 lineto '
        write(13,'(a)')'        0.5 0.5 lineto -0.5 0.5 lineto'
        write(13,'(a)')'        closepath stroke grestore} def'
        write(13,'(a)')'/psboxf {gsave translate rotate scale newpath '
        write(13,'(a)')'         -0.5 -0.5 moveto 0.5 -0.5 lineto '
        write(13,'(a)')'         0.5 0.5 lineto -0.5 0.5 lineto'
        write(13,'(a)')'         closepath stroke grestore} def'
        write(13,'(a)')'/pselip {gsave translate rotate scale '
        write(13,'(a)')'         newpath 0 0 1.0 0 360 arc '
        write(13,'(a)')'         closepath stroke grestore} def'
        write(13,'(a)')'/pselipf{gsave translate rotate scale '
        write(13,'(a)')'         newpath 0 0 1.0 0 360 arc '
        write(13,'(a)')'         closepath fill grestore} def'
        write(13,'(a)')'/psstring{gsave translate rotate 0 0 moveto'
        write(13,'(a)')'          show grestore} def'
        write(13,'(a)')'0.0 setgray'
        write(13,'(a)')' '
        write(13,'(a)')'/times-roman findfont 12. scalefont setfont'
        write(13,'(a)')' '
       return
      end

      subroutine ps_page
c      issues the newpage "showpage" command. this must be included
c      at the end of the file or the page won't necessarily print!!!!!
        write(13,'(a)')'showpage'
       return
      end

      subroutine ps_scale(xs,ys)
c      sets the x and y scales of the plot. it starts in units
c      of 1/72 inches (points) so to change to mm use
c      call ps_scale(2.835,2.835) etc.
        write(13,101)xs,ys
101     format(2f9.4,' scale')
       return
      end

      subroutine ps_translate(x,y)
c      translates the origin to (x,y).
        write(13,101)x,y
101     format(2f9.4,' translate')
       return
      end

      subroutine ps_rotate(xphi)
c      rotates the coordinate system by xphi degrees ccw
        write(13,101)xphi
101     format(f9.4,' rotate')
       return
      end

      subroutine ps_rgbcolor(r,g,b)
c      sets new color for drawing
        write(13,101)r,g,b
101     format(3f9.4,' setrgbcolor')
       return
      end

      subroutine ps_setgray(g)
c      sets new color for drawing
        write(13,101)g
101     format(f9.4,' setgray')
       return
      end

      subroutine ps_lw(xw)
c      sets a new default line width in current user scale.  0.0
c      is allowed and gives the narrowest line the output device
c      can produce.
        write(13,101)xw
101     format(f9.4,' setlinewidth')
       return
      end

      subroutine ps_line(x0,y0,xf,yf)
c      draws a line from (x0,y0) to (xf,yf)
        write(13,101)x0,y0,xf,yf
101     format(4f9.4,' psline')
       return
      end

      subroutine ps_circf(x,y,r)
c      draws a filled circle centered at x,y with radius r (mm)
        write(13,101)x,y,r
101     format(3f9.4,' pscircf')
       return
      end

      subroutine ps_circ(x,y,r)
c      draws a circle centered at x,y with radius r (mm)
        write(13,101)x,y,r
101     format(3f9.4,' pscirc')
       return
      end

      subroutine ps_rect(xbl,ybl,xtr,ytr)
c      draws a rectangle defined by bottom-left and top-right coords
c      (xbl,ybl), (xtr,ytr)
        write(13,101)xbl,ybl,xtr,ybl
        write(13,102)xtr,ytr,xbl,ytr
101     format(4f9.4)
102     format(4f9.4,' psrect')
       return
      end

      subroutine ps_rectf(xbl,ybl,xtr,ytr)
c      fills a rectangle defined by bottom-left and top-right coords
        write(13,101)xbl,ybl,xtr,ybl
        write(13,102)xtr,ytr,xbl,ytr
101     format(4f9.4)
102     format(4f9.4,' psrectf')
       return
      end

      subroutine ps_box(x,y,xs,ys,xphi)
c      draws a rectangle centered at x,y with sides xs,ys big
c      at an angle of xphi.
        write(13,101) xs,ys,xphi,x,y
101     format(4f9.4,' psbox')
       return
      end

      subroutine ps_boxf(x,y,xs,ys,xphi)
c      draws a filled rectangle just like ps_box.
        write(13,101)xs,ys,xphi,x,y
101     format(4f9.4,' psboxf')
       return
      end

      subroutine ps_elip(x,y,xs,ys,xphi)
c      draws an elipse centered at x,y with axes xs,ys big at
c      an angle of xphi.
        write(13,101)xs,ys,xphi,x,y
101     format(4f9.4,' pselip')
       return
      end

      subroutine ps_elipf(x,y,xs,ys,xphi)
c      fills an elipse centered at x,y with axes xs,ys big at
c      an angle of xphi.
        write(13,101)xs,ys,xphi,x,y
101     format(4f9.4,' pselipf')
       return
      end

      subroutine ps_rgbimage(iw,ih,image,x,y,xw,yw)
c      writes rgb image out 3hex bytes per pixel (use ps_rgb85 to
c      instead to save space unless your printer does not have
c      an ascii85 decoder.
       character*1 image(*)
       character*72 cline
       character*6  cpix
       character*2  ccolor

       ibpc=8
       cline=' '
       nlines=(iw*ih)/12        !each line will be 12 pixels (72 chars)
       nlast=(iw*ih) - ( (iw*ih)/12 )*12  !number of pixels in last line

       write(13,*)' '
       write(13,100)
100    format(' gsave')
       write(13,101)iw*3
101    format('/picstr ',i4,' string def')
       write(13,102)x,y,xw,yw
102    format(2f9.4,' translate',2f9.4,' scale')
       write(13,103)iw,ih,ibpc,iw,ih
103    format(3i5,' [',i4,' 0 0 ',i4,' 0 0 ]')
       write(13,104)
104    format('{currentfile picstr readhexstring pop} false 3')
       write(13,105)
105    format('colorimage')

       do j=0,nlines-1                !do all of the complete lines
        do i=0,11
            do ic=1,3
              icolor=ichar( image(j*36+i*3+ic) )
              write(ccolor,'(z2.2)')icolor
              cpix( ic*2-1:ic*2 ) = ccolor
            end do
            cline( i*6+1:i*6+6 ) = cpix
          end do
          write(13,'(a)')cline
       end do
       do i=0,nlast-1                !and then the last incomplete line
        do ic=1,3
         icolor=ichar( image(nlines*36+i*3+ic) )
           write(ccolor,'(z2.2)')icolor
         cpix( ic*2-1:ic*2 ) = ccolor
          end do
          cline( i*6+1:i*6+6 ) = cpix
       end do
       write(13,'(a)')cline(1:nlast*6)

       write(13,106)
106    format(' grestore')
       write(13,*)' '

      return
      end

      subroutine ps_rgb85(iw,ih,image,x,y,xw,yw)
c      routine that allows you to write an rgb image to the page.
c      iw,ih = width and hight of image in pixels
c      image = row-major order r-g-b image bytes
c              from bottom left to top right
c      x,y = position of bottom left of picture
c      xw,yw = width and hight of box to fill with the image
c      this routine requires the ps interpreter to have an
c      ascii85 decoder (officially level 2 stuff).  if yours
c      does not use ps_rgbimage instead.
        character*1 image(*)
        integer*4 ic(5)
        character*72 cline

        ibpc=8
        cline=' '

        write(13,*)' '
        write(13,100)
100     format(' gsave')
        write(13,102)x,y,xw,yw
102     format(2f9.4,' translate',2f9.4,' scale')
        write(13,103)iw,ih,ibpc,iw,ih
103     format(3i5,' [',i4,' 0 0 ',i4,' 0 0 ]')
        write(13,104)
104     format('currentfile /ascii85decode filter false 3')
        write(13,105)
105     format('colorimage')

        num_4t=(iw*ih*3)/4      !number of clean 4-tuples
        nrem_4t=mod(iw*ih*3,4)  !number of excess bytes.

        nchar=1
        do ihex=1,num_4t    !write all but the last few bytes
        i1=ichar( image(4*(ihex-1)+1) )
        i2=ichar( image(4*(ihex-1)+2) )
        i3=ichar( image(4*(ihex-1)+3) )
        i4=ichar( image(4*(ihex-1)+4) )
        i256=ishft(i1,24)
        i256=ior(ishft(i2,16),i256)
        i256=ior(ishft(i3,8),i256)
        i256=ior(i4,i256)
        if (i256.eq.0) then
          cline(nchar:nchar)='z'
          nchar=nchar+1
          if(nchar.eq.73) then
            write(13,'(a)')cline(1:72)
          cline=' '
          nchar=1
          end if
         else
          continue
          if (i256.gt.0) then
            ic(5) = mod(i256,85)
            i256 = i256/85
           else if (i256.lt.0) then
            ic(5)=mod(i256,85)+85+1     !mod(i256,85)
            ic(5)=mod(ic(5),85)
            itmp=i256/85+50529027       !i256/85 ; don't ask :)
            if( (ic(5)-mod(i256,85)).eq.86 ) itmp=itmp-1
            i256 = itmp
          end if
          do i=4,1,-1
           ic(i) = mod(i256,85)
           i256=i256/85
          end do
          do i=1,5
            cline(nchar:nchar)=char(33+ic(i))
            nchar=nchar+1
            if(nchar.eq.73) then
              write(13,'(a)')cline(1:72)
            cline=' '
            nchar=1
            end if
            end do
        end if
        end do

        i2=0
        i3=0
        if (nrem_4t.gt.0) then
          i1=ichar( image(4*(ihex-1)+1) )   !do remaining bytes
          if(nrem_4t.gt.1) i2=ichar( image(4*(ihex-1)+2) )
          if(nrem_4t.gt.2) i3=ichar( image(4*(ihex-1)+3) )
          i256=ishft(i1,24)
          i256=ior(ishft(i2,16),i256)
          i256=ior(ishft(i3,8),i256)
          if (i256.ge.0) then        !get first char
            ic(5) = mod(i256,85)
            i256 = i256/85
           else
            ic(5)=mod(i256,85)+85+1     !mod(i256,85)
            ic(5)=mod(ic(5),85)
            itmp=i256/85+50529027       !i256/85 ; don't ask :)
            if( (ic(5)-mod(i256,85)).eq.86 ) itmp=itmp-1
            i256 = itmp
          end if
          do i=4,1,-1
           ic(i) = mod(i256,85)
           i256=i256/85
          end do
          do i=1,nrem_4t+1
            cline(nchar:nchar)=char(33+ic(i))
            nchar=nchar+1
            if(nchar.eq.73) then
             write(13,'(a)')cline(1:72)
             cline=' '
             nchar=1
            end if
          end do
        end if
        write(13,'(a)')cline(1:nchar)

        cline(1:2)='~>'      !put on eof marker
        write(13,'(a)')cline(1:2)

        write(13,106)
106     format(' grestore')
        write(13,*)' '

       return
      end

      subroutine ps_image(iw,ih,image,x,y,xw,yw)
c      routine that allows you to write a gray-scale image to the page.
c      iw,ih = width and hight of image in pixels
c      image = row-major order image bytes
c              from bottom left to top right
c      x,y = position of bottom left of picture
c      xw,yw = width and hight of box to fill with the image
        character*1 image(*)
        character*72 cline
        character*2  cpix

        ibpc=8
        cline=' '
        nlines=(iw*ih)/36    !each line will be 36 pixels (72 chars)
        nlast=(iw*ih) - ( (iw*ih)/36 )*36 !number of pixels in last line

        write(13,*)' '
        write(13,100)
100     format(' gsave')
        write(13,101)iw
101     format('/picstr ',i4,' string def')
        write(13,102)x,y,xw,yw
102     format(2f9.4,' translate',2f9.4,' scale')
        write(13,103)iw,ih,ibpc,iw,ih
103     format(3i5,' [',i4,' 0 0 ',i4,' 0 0 ]')
        write(13,104)
104     format('{currentfile picstr readhexstring pop}')
        write(13,105)
105     format('image')

        do j=0,nlines-1    !do all of the complete lines
         do i=1,36
          icolor=ichar( image(j*36+i) )
          write(cpix,'(z2.2)')icolor
          cline( i*2-1:i*2 ) = cpix
         end do
         write(13,'(a)')cline(1:72)
        end do
        do i=1,nlast    !and then the last incomplete line
          icolor=ichar( image(nlines*36+i) )
          write(cpix,'(z2.2)')icolor
          cline( i*2-1:i*2 ) = cpix
        end do
        write(13,'(a)')cline(1:nlast*2)

        write(13,106)
106     format(' grestore')
        write(13,*)' '

       return
      end

      subroutine ps_string(x,y,ctext,xphi,n)
c      puts a character string 'ctext' of length 'n', at the location
c      (x,y) and rotated counter clockwise by the angle xphi.
        character*80 ctext
        write(13,101)ctext(:n),xphi,x,y
101     format(' (',a,') ',3f9.4,' psstring')
       return
      end

      subroutine ps_setfont(fname,fs)
c      loads a new font with name fname and size fs in current user
c      coordinates.  unless you are sure the font is
c      currently in the printer only use the times, helvetica, courier,
c      and symbol fonts.

c      times-roman  times-italic  times-bold  times-bolditalic

c      helvetica,helvetica-oblique,helvetica-bold,helvetica-boldoblique

c      courier  courier-oblique  courier-bold  courier-boldoblique

c      symbol

       character*80 fname
        n=0
        do i=1,80
          if(fname(i:i).eq.' ') goto 10
        n=n+1
        end do
10      continue
        write(13,101)fname(1:n),fs
101     format(' /',a,' findfont ',f9.4,' scalefont setfont')
       return
      end



