      Program MATCOV1
      include 'MAXSZ.INC'
      dimension r(2),ichans(MAXCH),istuff(64),a(MAXSZ),b(MAXSZ)
      logical bo(2),errflg
      character*80 strs(2)
      character*71 cstr,dstr
      character*2 bi(2)
      integer avgs, updt,nseg
      data bi/'-a','-c'/, dstr/' '/
C      include 'FASTFPA.INC'
      cstr='Usage:  matcov2 [-ac] avgs update inputfile outputfile'
      call prompts(cstr,2,bi,bo,2,2,r,strs,dstr,errflg)
      
      if (errflg)stop
      avgs= r(1)
      updt = r(2)
      call ncsetup(np ,nc ,ichans,strs(1))

      nseg = (np-avgs)/updt + 1
      npc = np*nc
      npout = nc*nseg*nc
c     if complex then half the number of channels
      if ( bo(2) ) npout = npout/2
      
c     confirm arrays will fit
      if (npc .gt. MAXSZ/2) then
          errflg=.true.
          call errsio(' Error - Inputfile too large ',npc)
      end if

C  If -c option is set, make sure Nc is even; 
      if (bo(2) .and. (mod(nc,2).ne.0)) then
          errflg=.true.
          call errsio(' Error - #channels not divisible by 2 ',nc)
      end if

      if (.not. errflg) then

        if (.not. bo(2)) then
          call  rdsio(a,np,nc,1,ichans,irl,strs(1),istuff)
          call  matcovr(a,np,nc,b,avgs,updt,nseg) 
        else
          call  rdsio(a,np,nc,1,ichans,irl,strs(1),istuff)
          call  matcovc(a,np,nc,b,avgs,updt,nseg) 
        end if

        if (.not. bo(2)) then
                   irl=irlget(nc*avgs,4)
           call  mksio(b,nc*nseg,nc,irl,strs(2),dstr,istuff)
        else    
                   irl=irlget(nc/2*avgs,4)
           call  mksio(b,nc/2*nseg,nc,irl,strs(2),dstr,istuff)
        end if

      end if
      stop
      end

      subroutine matcovr(a,np,nc,b,avgs,updt,nseg)
      integer avgs,updt
      dimension a(np,nc),b(1)
      double precision c
      nseg = (np-avgs)/updt + 1

        do 45 i=1,nc
          do 45 j=1,nc
            do 45 k=0,nseg-1
              c=0.
              do 40 l=1,avgs
40            c = c + a(l+(k*updt),i)*a(l+(k*updt),j)
              b(i+((j-1)*nseg*nc)+(k*nc)) = c/avgs
45          continue

      return
      end

      subroutine matcovc(a,np,nc,b,avgs,updt,nseg)
      integer avgs,updt
      dimension a(np,nc)
      dimension b(1)
      double precision c1,c2
      nseg = (np-avgs)/updt + 1
      nc2=nc/2
      np1 = nc2*nseg
        do 45 i=1,nc2
                i2=i*2
                i1=i2-1
          do 45 j=1,nc2
                   j2=j*2
                   j1=j2-1
            do 45 k=0,nseg-1
                     c1=0.
                     c2=0.
              do 40 l=1,avgs
                c1 = c1 + a(l+k*updt,i1)*a(l+k*updt,j1)+
     *                    a(l+k*updt,i2)*a(l+k*updt,j2)
                c2 = c2 + a(l+k*updt,i2)*a(l+k*updt,j1)-
     *                    a(l+k*updt,i1)*a(l+k*updt,j2)
40            continue
             ij = i + ((j-1)*2*nseg + k)*nc2
              ik = i + ((j-1)*2*nseg + k+nseg)*nc2
              b(ij) = c1/avgs
              b(ik) = c2/avgs
45            continue
      return
      end
