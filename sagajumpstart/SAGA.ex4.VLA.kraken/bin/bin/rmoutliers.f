      program rmoutliers

C  compare estimated slant ranges to sfile slant ranges
C  either removing the whole point if differences in 
C  slant ranges are too large or changing sfile to not
C  use that point.
C  will read in transponder coordinates from 'pfile.1' 
C  will read in array info from 'pfile.1'

      parameter( maxnp   = 1000 )
      parameter( maxxpdr = 31 )
 
      dimension pf(maxnp*3),lf(maxnp*3),xpdr(maxxpdr*3)
      dimension sf(maxnp*maxxpdr),se(maxnp*maxxpdr)
      dimension Opf((maxnp+maxxpdr)*3),Olf(maxnp*3)
      dimension Osf(maxnp*maxxpdr),Ose(maxnp*maxxpdr)
      dimension r(1),ic(1000)
      dimension ic7(1000),slant(1000),diff(1000)
      logical bo(3),errflg
      character*80 strs(4)
      character*71 cstr,dstr
      character*2 bi(3)
      logical iflg, lflg
      data bi/'-a','-i','-l'/, dstr/' '/
      cstr='Usage:  rmoutliers [-ail] thres pfile lfile sfile sefile'
      call prompts(cstr,3,bi,bo,1,4,r,strs,dstr,errflg)
      if (errflg)stop

      call ncsetup(np ,nc ,ic,strs(1))
      call ncsetup(np1,nc1,ic,strs(3))

      nele = np1
      nxpdr= np - np1
      thres = r(1)
C     -i option = remove slant range instead of whole point if it 
C                 exceeds limit. (default is to remove whole point)
C     -l option = loquacious mode
      iflg = bo(2)
      lflg = bo(3)

      if ( nele .GT. maxnp) then
           print *,'Error: maxnp',maxnp,' < ',nele,' # of obs'
           stop
      endif
      if ( iflg ) then
           call rmslant(nxpdr,nele,thres,lflg,strs,pf,lf,sf,se,xpdr,
     *          slant,diff,ic7)
      else
           call rmpoint(nxpdr,nele,thres,lflg,strs,pf,lf,sf,se,xpdr,
     *          slant,diff,ic7,Opf,Olf,Osf,Ose)
      endif

      end

      subroutine rmpoint(nxpdr,nele,thres,lflg,strs,pf,lf,sf,se,xpdr,
     *           slant,diff,ic7,Opf,Olf,Osf,Ose)
      integer nxpdr,nele
      real    thres
      logical lflg
      character*80 strs(4)
      real xpdr(nxpdr,3),pf(nele,3),lf(nele*3)
      real sf(nele,nxpdr),se(nele,nxpdr)
      real Opf(nele+nxpdr,3),Olf(nele*3)
      real Osf(nele,nxpdr),Ose(nele,nxpdr)
      real diff(nxpdr),slant(nxpdr)
      integer ic7(nxpdr)

      real maxdif,sum
      integer i,xy,ic3(3),badpt,head(64),rl,Onp
      integer np1,nc1,ic1(1),cnt,totslt
      data ic1/1/,ic3/1,2,3/

      do 10 i = 1,nxpdr
10       ic7(i)=i

C  read in lfile,sefile,sfile,xpdrs,pfile
      call ncsetup(np1,nc1,ic,strs(2))
      call rdsio( lf , np1 ,  1  ,   1   ,ic1,rl,strs(2),head)
      call rdsio( se ,nele ,nxpdr,   1   ,ic7,rl,strs(4),head)
      call rdsio( sf ,nele ,nxpdr,   1   ,ic7,rl,strs(3),head)
      call rdsio(xpdr,nxpdr,  3  ,   1   ,ic3,rl,strs(1),head)
      call rdsio( pf ,nele ,  3  ,nxpdr+1,ic3,rl,strs(1),head)
C  fill the output pfile array with the transponder positions
      call rdsiom(Opf,nele+nxpdr,3,  1   ,ic3,rl,strs(1),head,nxpdr)

C  determine if lfile allows XY or XYZ positions to vary
C    xy = 2 means XY only, xy = 3 means XYZ can vary
      xy = np1/nele

C  write out xpdr numbers if desired
      if (lflg) write(6,91) ' ',1,2,3,4,5,6,7
91    format(A1,7I9)
      Onp = 0
C  loop over all observations 
      do 20 i = 1,nele
         sum    = 0.
         badpt  = 0
         totslt = nxpdr
         maxdif = 0.
C  loop over all transponders computing estimated slant range from
C  current transponder (j) to current observation location (i)
         do 25 j = 1,nxpdr
            slant(j) = sqrt( ( pf(i,1) - xpdr(j,1) )**2 +
     *                       ( pf(i,2) - xpdr(j,2) )**2 +
     *                       ( pf(i,3) - xpdr(j,3) )**2 )
            diff(j)  = abs(slant(j) - sf(i,j))
C           throw away result if sefile says point not being used
            if ( se(i,j) .NE. 1. ) then
                  diff(j) = 0.
                  totslt= totslt-1
            endif
            sum = sum + diff(j)
C           compute maximum difference
            if (diff(j) .GT. maxdif) then
                  maxdif = diff(j)
                  cnt = j
            endif
25       continue
C  avgerr=sum of (estimated slant range - true slant range)/# of slant ranges
         avgerr = sum/totslt
C  if average error is < threshold then add the point to the output nav 
C  files (pfile,sfile,sefile, and lfile). Note this will alter the 
C  input files.
         if (avgerr .LT. thres) then
                  Onp = Onp + 1
                  do 110 j = 1,3
110                  Opf(Onp+nxpdr,j) = pf(i,j)
C                  do 120 j = 1,7
                  do 120 j = 1,nxpdr
                     Osf(Onp,j) = sf(i,j)
120                  Ose(Onp,j) = se(i,j)
                  do 130 j = 1,xy
130                  Olf((Onp-1)*xy+j) = lf((i-1)*xy+j)
         endif
C  if a point was changed and printing is desired then output point number
C  the differences and the number of the point that was changed
         if ((lflg).AND.(avgerr.GT.thres)) write(6,90) i,diff,avgerr,badpt
90       format(I4,32F9.2,I4)
20    continue
C
      np1 = Onp*xy
      call mksio(lf,np1,1,rl,strs(2),' ',head)
      call mksiom(Ose,nele,nxpdr,rl,strs(4),' ',head,1,Onp,nxpdr,ic7)
      call mksiom(Osf,nele,nxpdr,rl,strs(3),' ',head,1,Onp,nxpdr,ic7)
      call mksiom(Opf,nele+nxpdr,3,rl,strs(1),' ',head,1,Onp+nxpdr,
     *               3,ic3)
      write(6,*) 'removed',nele-Onp,' points,',Onp,' points remaining'
C
      return
      end

      subroutine rmslant(nxpdr,nele,thres,lflg,strs,pf,lf,sf,se,xpdr,
     *           slant,diff,ic7)
      integer nxpdr,nele
      real    thres
      logical lflg
      character*80 strs(4)
      real xpdr(nxpdr,3),pf(nele,3),lf(nele*3)
      real sf(nele,nxpdr),se(nele,nxpdr)
      real diff(nxpdr),slant(nxpdr)
      integer ic7(nxpdr)

      real maxdif
      integer i,ic3(3),badpt,head(64),rl,totfix
      integer np1,nc1,ic1(1),cnt,totslt
      data ic1/1/,ic3/1,2,3/,totfix/0/

      do 11 i = 1,nxpdr
11       ic7(i)=i
C  read in lfile,sefile,sfile,xpdrs,pfile
      call ncsetup(np1,nc1,ic,strs(2))
      call rdsio( lf , np1 ,  1  ,   1   ,ic1,rl,strs(2),head)
      call rdsio( se ,nele ,nxpdr,   1   ,ic7,rl,strs(4),head)
      call rdsio( sf ,nele ,nxpdr,   1   ,ic7,rl,strs(3),head)
      call rdsio(xpdr,nxpdr,  3  ,   1   ,ic3,rl,strs(1),head)
      call rdsio( pf ,nele ,  3  ,nxpdr+1,ic3,rl,strs(1),head)

C  write out xpdr numbers if desired
      if (lflg) write(6,91) ' ',1,2,3,4,5,6,7
91    format(A1,7I9)
C  loop over all observations 
      do 20 i = 1,nele
         badpt  = 0
         totslt = nxpdr
         maxdif = 0.
C  loop over all transponders computing estimated slant range from
C  current transponder (j) to current observation location (i)
         do 25 j = 1,nxpdr
            slant(j) = sqrt( ( pf(i,1) - xpdr(j,1) )**2 +
     *                       ( pf(i,2) - xpdr(j,2) )**2 +
     *                       ( pf(i,3) - xpdr(j,3) )**2 )
            diff(j)  = abs(slant(j) - sf(i,j))
C           throw away result if sefile says point not being used
            if ( se(i,j) .NE. 1. ) then
                  diff(j) = 0.
                  totslt= totslt-1
            endif
C           compute maximum difference
            if (diff(j) .GT. maxdif) then
                  maxdif = diff(j)
                  cnt = j
            endif
25       continue
C  if maximum difference is > threshold and > 3 slant ranges present
C  then fix sefile so that nav will not use bad slant range
         if ((maxdif .GT. thres) .AND. (totslt .GT. 3))then
                  badpt = cnt
                  se(i,badpt)=10000.
                  totfix = totfix+1
         endif
C  if a point was changed and printing is desired then output point number
C  the differences and the number of the point that was changed
         if ((lflg).AND.(badpt.NE.0)) write(6,90) i,diff,badpt
90       format(I4,31F9.2,I4)
20    continue
C
      call mksio(se,nele,nxpdr,rl,strs(4),' ',head)
      write(6,*) 'fixed',totfix,' slant ranges out of',nele,' possible'
C
      return
      end
