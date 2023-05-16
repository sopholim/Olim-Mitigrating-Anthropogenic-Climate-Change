      subroutine OTEC (ts, i, j, joff, nd, kmt)

#if defined O_mom && defined O_OTEC
      implicit none

      integer i, ie, is, j, je, joff, jrow, js, k, kb, kbo
      integer kt, ktot, kven, la, lb, l, m, n, nd

      logical chk_la, chk_lb

      real rl, ru, tmx(3), tsm(3), zsm, dens, tq, sq, drodt, drods
      real drhodt, drhods, ddensdtdt, ddensdtds, ddensdsds

      include "size.h"
      include "param.h"
      include "scalar.h"
      include "coord.h"
      include "grdvar.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "state.h"
      include "dens.h"
      include "OTEC.h"

      integer kmt(imt,jmt)

      integer nti
      real ts(km,nt), dt, tk, tkm1

      jrow = j + joff
      kbo = kmt(i,jrow)

      nti = 1 + c2dtts*disvel(nd)/dzt(1)
      dt = c2dtts/nti

      do m=1,nti

!-----------------------------------------------------------------------
!       upsteam advection over the dischage area
!
!       draw water up from the bottom (uptake) of the pipe then advect
!       it down over the discharge area before doing lateral mixing
!-----------------------------------------------------------------------

        do n=1,nt
          k = levdis
          tk = ts(k,n)
          ts(k,n) = ts(k,n) + dt*disvel(nd)*(ts(levupt,n) 
     &            - ts(levdis,n))/dzt(k)
          tkm1 = tk
          do k=levdis+1,levupt-1
            tk = ts(k,n)
            ts(k,n) = ts(k,n) + dt*disvel(nd)*(tkm1 - tk)/dzt(k)
            tkm1 = tk
          enddo
          k = levupt
          ts(k,n) = ts(k,n) + dt*disvel(nd)*(tkm1 - ts(k,n))/dzt(k)
        enddo

!-----------------------------------------------------------------------
!       convect over the dischage area after advection
!-----------------------------------------------------------------------

!       search for unstable regions starting from the top
        kt = 1
        kb = 2
        do while (kt .lt. kbo)
          ru = dens (ts(kt,1)-to(kb), ts(kt,2)-so(kb), kb)
          rl = dens (ts(kb,1)-to(kb), ts(kb,2)-so(kb), kb)

!         sum the first pair found in an unstable region
          if (ru .gt. rl) then
            chk_la = .true.
            chk_lb = .true.
            zsm = dztxcl(kt) + dztxcl(kb)
            tsm(1) = ts(kt,1)*dztxcl(kt) + ts(kb,1)*dztxcl(kb)
            tmx(1) = tsm(1)/zsm
            tsm(2) = ts(kt,2)*dztxcl(kt) + ts(kb,2)*dztxcl(kb)
            tmx(2) = tsm(2)/zsm

            do while (chk_lb .or. chk_la)

!             check for an unstable level (lb) below kb
              if (kb .ge. kbo) chk_lb = .false.
              do while (chk_lb)
                chk_lb = .false.
                lb = kb + 1
                ru = dens (tmx(1)-to(lb),tmx(2)-so(lb),lb)
                rl = dens (ts(lb,1)-to(lb),ts(lb,2)-so(lb),lb)
                if (ru .gt. rl) then
!                 add new level to sums
                  kb = lb
                  zsm = zsm + dztxcl(kb)
                  tsm(1) = tsm(1) + ts(kb,1)*dztxcl(kb)
                  tmx(1) = tsm(1)/zsm
                  tsm(2) = tsm(2) + ts(kb,2)*dztxcl(kb)
                  tmx(2) = tsm(2)/zsm
                  chk_la = .true.
                  if (kb .lt. kbo) chk_lb = .true.
                endif
              enddo

!             check for an unstable level (la) above kt
!             chk_la = .true.
              if (kt .le. 1) chk_la = .false.
              do while (chk_la)
                chk_la = .false.
                la = kt - 1
                ru = dens (ts(la,1)-to(kt),ts(la,2)-so(kt),kt)
                rl = dens (tmx(1)-to(kt),tmx(2)-so(kt),kt)
                if (ru .gt. rl) then
!                 add new level to sums
                  kt = la
                  zsm = zsm + dztxcl(kt)
                  tsm(1) = tsm(1) + ts(kt,1)*dztxcl(kt)
                  tmx(1) = tsm(1)/zsm
                  tsm(2) = tsm(2) + ts(kt,2)*dztxcl(kt)
                  tmx(2) = tsm(2)/zsm
                  chk_lb = .true.
                endif
              enddo
            enddo

!           mix all tracers from kt to kb
            do k=kt,kb
              ts(k,1) = tmx(1)
              ts(k,2) = tmx(2)
            enddo
            do n=3,nt
              tsm(3) = c0
              do k=kt,kb
                tsm(3) = tsm(3) + ts(k,n)*dztxcl(k)
              enddo
              tmx(3) = tsm(3)/zsm
              do k=kt,kb
                ts(k,n) = tmx(3)
              enddo
            enddo

            kt = kb + 1
          else
            kt = kb
          endif

!         continue the search for other unstable regions
          kb = kt + 1
        enddo
        
      enddo

      return
      end

      subroutine OTECplace

      implicit none

      include "size.h"
      include "param.h"
      include "scalar.h"
      include "coord.h"
      include "grdvar.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "mw.h"
      include "OTEC.h"
      include "tmngr.h"

      integer i, j, n, m
      integer numadd               ! number to add in a grid per check
      integer numsub               ! number to subtract in a grid per check
      integer numplc               ! number of plants placed
      integer numpek               ! number of plants at peak power (at yerpek)
      integer locpot               ! potential plant locations in a grid cell
      integer nummov               ! number of plants in a grid cell that need to move
      integer nr(imt*jmt), kr(10)  ! randomized index arrays

      real yr                      ! current year
      real grdchk                  ! minimum gradient to check for placement
      real grdinc                  ! increnet for checking temperature gradient
      real grdadj                  ! adjusted grid average temperature gradient
      real grdpro                  ! portion of unused temperature gradient 
      real rempro                  ! portion of plant locations remaining

# if defined O_OTEC_old_place
      return
# endif
      nmtplt = 0
      numpek = pwrpek/pwrplt
      yr = relyr + year0

# if defined O_OTEC_linear
!     calculate the number of available plants over time using a linear increase    
!     assume build capacity is ramped up to meet projected demand pwrpek in yerpek
      if (yr .lt. yerstr) then
        nmtblt = int(pwrstr/pwrplt)
      elseif (yr .lt. yerpek .and. yerpek .gt. yerstr) then
        nmtblt = numpek*(yr-yerstr)/(yerpek-yerstr)
      elseif (yr .lt. yerend .and. yerend .gt. yerstr) then
        nmtblt = numpek*(yr-yerend)/(yerpek-yerend)
      else
        nmtblt = int(pwrend/pwrplt)
      endif
# else
!     calculate the number of available plants over time using raised cosine     
!     assume build capacity is ramped up to meet projected demand pwrpek in yerpek
      if (yr .lt. yerstr) then
        nmtblt = int(pwrstr/pwrplt)
      elseif (yr .lt. yerpek .and. yerpek .gt. yerstr) then
        nmtblt = numpek*(0.5 - 0.5*cos(pi*(yr-yerstr)/(yerpek-yerstr)))
      elseif (yr .lt. yerend .and. yerend .gt. yerstr) then
        nmtblt = numpek*(0.5 + 0.5*cos(pi*(yr-yerpek)/(yerend-yerpek)))
      else
        nmtblt = int(pwrend/pwrplt)
      endif
# endif

      if (nmtblt .gt. 0) then
       
!       estimate plant demand from emissions (done in OTECdata)
        nmtdem = int(elcdem/pwrplt)

!       limit build to demand
        nmtblt = min(nmtblt, nmtdem)

!       find all possible locations of viable plants
        nmtplt = 0
        locrem(:,:) = 0.
        lctrem = 0
        lctpot = 0
        nmtmov = 0
        grdchk = 0.
        grdinc = 1.
        numplc = 0
        numadd = 1
        numsub = 1

# if defined O_OTEC_increment_place        
        nmtloc = 0
        do j=2,jmtm1
          do i=2,imtm1
            if (kmt(i,j) .ge. levupt) then
              if (grdavg(i,j).ge.grdmin .and. locmap(i,j).gt.1) then
!               find the maximum temperature gradient
                if (grdavg(i,j) .ge. grdchk) grdchk = grdavg(i,j)
!               save total area, number and locations of viable plant grid cells
                nmtloc = nmtloc + 1
                indplt(nmtloc,1) = i
                indplt(nmtloc,2) = j
                locpot = int(tgarea(i,j)/aremin)
                lctpot = lctpot + locpot
                locrem(i,j) = locpot - numplt(i,j)
#  if defined O_OTEC_adjust_gradient
!               move plants if the gradient is becoming depleted
                if (grdstr(i,j) .gt. 0. .and. float(locpot) .gt. 0.) then
                  grdpro = max(0.,(grdstr(i,j)-grdavg(i,j))/grdstr(i,j))
                  rempro = float(locrem(i,j))/float(locpot)
                  if (rempro .lt. grdpro) then
!                   move enough plants to make the portion of remaining area equal 
!                   to the portion of the remaining gradient
                    nummov = max(0,int(grdpro*float(locpot))-numplt(i,j))
                    nmtmov = nmtmov + nummov
                    numplt(i,j) = numplt(i,j) - nummov
                    locrem(i,j) = locrem(i,j) + nummov
                  endif
                endif
#  endif
                nmtplt = nmtplt + numplt(i,j)
                lctrem = lctrem + locrem(i,j)
              else
!               need to move due to inadequate gradient or location no longer allowed
                nmtmov = nmtmov + numplt(i,j)
                numplt(i,j) = 0
              endif
            endif
          enddo
        enddo
# endif


        if (nmtblt .gt. nmtplt) then 

          if (lctrem .gt. 0) then 
# if defined O_OTEC_increment_place
            grdchk = float(int(grdchk)) - grdinc ! starting temperature gradient for search
# endif
!           place plants in areas with the largest temperature gradient 
            m = 0
            n = 1
            do
# if defined O_OTEC_increment_place
#  if defined O_OTEC_random_place
!             scramble location indicies to randomize deployment
              if (n .eq. 1) call scramble(nmtloc, nr)
              i = indplt(nr(n),1)
              j = indplt(nr(n),2)
              call scramble(11, kr)
#  else
              i = indplt(n,1)
              j = indplt(n,2)
              kr(1) = 6
#  endif
!             check the average gradient is larger than the current grdchk
!             add noise to temperature gradient to eliminate hard boundries
#  if defined O_OTEC_scale_gradient
!             scale the useable gradient with the portion of area remaining
!             in order to reduce overloading a grid cell with plants
              grdadj = (grdavg(i,j) - grdmin)*locrem(i,j)
     &                 /float(int(tgarea(i,j)/aremin)) + grdmin
              if (grdadj-float(kr(1)-1)/10+0.5 .ge. grdchk) then
#  else
              if (grdavg(i,j)-float(kr(1)-1)/10+0.5 .ge. grdchk) then
#  endif
# endif
!               check if this grid also has enough area for additional plants
                if (locrem(i,j) .ge. numadd) then
!                 place a place and add to the total number and discharge area
                  nmtplt = nmtplt + numadd
                  numplt(i,j) = numplt(i,j) + numadd
                  locrem(i,j) = locrem(i,j) - numadd
                  lctrem = lctrem - numadd
                  numplc = numplc + 1
!                 exit if number of plants meets number built or demand or area is zero
                  if (nmtplt .ge. nmtblt) exit
                  if (nmtplt .ge. nmtdem) exit
                  if (lctrem .le. 0) exit
                endif
# if defined O_OTEC_increment_place
              endif
# endif
!             increment the grid index to check the next viable grid cell
              n = n + 1
              if (n .gt. nmtloc) then
                m = m + 1
                if (m .gt. 10000) then
                  print*, "=> Warning: problem placing plants in OTEC.f"
                  print*, "number of plants placed: ", numplc
                  print*, "number of plants to meet demand: ", nmtdem
                  print*, "number of plants built: ", nmtblt
                  print*, "number of plants in total: ", nmtplt
                  print*, "locations remaining: ", lctrem
                  exit
                endif
!               if at the last viable cell start seach at lower gradient
                n = 1
# if defined O_OTEC_increment_place
                grdchk = max(grdchk - grdinc, grdmin)            
# endif
              endif
            enddo
          endif

        elseif (nmtblt .lt. nmtplt) then 
         
!         decommission plants from suboptimal locations (smallest dt)
          grdchk = grdmin + grdinc ! starting temperature gradient for search
!         remove plants in areas with the smallest temperature gradient 
          m = 0
          n = nmtloc
          do
# if defined O_OTEC_increment_place
#  if defied O_OTEC_random_place
!           scramble location indicies to randomize decommissioning
            if (n .eq. 1) call scramble(nmtloc, nr)
            i = indplt(nr(n),1)
            j = indplt(nr(n),2)
            call scramble(11, kr)
#  else
            i = indplt(n,1)
            j = indplt(n,2)
            kr(1) = 6
#  endif
!           check the average gradient is smaller than the current grdchk
!           add noise to temperature gradient to eliminate hard boundries
            if (grdavg(i,j)-float(kr(1)-1)/9+0.5 .le. grdchk) then
# endif
!             check if this grid has enough plants to remove
              if (numplt(i,j) .ge. numsub) then
!               remove a plant and subtract from the total number and discharge area
                nmtplt = nmtplt - numsub
                numplt(i,j) = numplt(i,j) - numsub
                locrem(i,j) = locrem(i,j) + numsub
                lctrem = lctrem + numsub
                numplc = numplc - numsub
                nmtmov = nmtmov + numsub
!               exit if number of plants meets number built
                if (nmtplt .le. nmtblt) exit
              endif
# if defined O_OTEC_increment_place
            endif
# endif
!           increment the grid index to check the next viable grid cell
            n = n - 1
            if (n .le. 1) then
              m = m + 1
              if (m .gt. 10000) then
                print*, "=> Warning: problem removing plants in OTEC.f"
                print*, "number of plants removed: ", -numplc
                print*, "number of plants to meet demand: ", nmtdem
                print*, "number of plants built: ", nmtblt
                print*, "number of plants in total: ", nmtplt
                print*, "locations remaining: ", lctrem
                exit
              endif
!             if at the last viable cell start seach at higher gradient
              n = nmtloc
# if defined O_OTEC_increment_place
              grdchk = max(grdchk + grdinc, grdmin)            
# endif
            endif
          enddo
        endif
      endif

      return
      end

      subroutine scramble(num, array)

      integer i, j, k, m, n, num, temp
      integer array(num)

      real u

      array = [(i,i=1,num)]
      n = 1
      m = num
      do k=1,2
         do i=1,m
            call random_number(u)
            j = n + floor((m+1-n)*u)
            temp=array(j)
            array(j) = array(i)
            array(i) = temp
         enddo
      enddo
#endif

      return
      end
