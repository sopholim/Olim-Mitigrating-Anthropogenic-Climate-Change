      subroutine tracer (joff, js, je, is, ie)

#if defined O_mom
!=======================================================================
!     compute tracers at "tau+1" for rows js through je in the MW.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name

      integer istrt, iend, i, k, j, ip, kr, jq, n, jp, jrow, iou, js
      integer je, limit, joff, is, ie, kmx, m, kb, idiag, index
      integer it(10), iu(10), ib(10), ic(10), nfnpzd, mfnpzd, mxfnpzd
      integer id_time, id_xt, id_yt, id_zt

      logical inqvardef, exists

      real rctheta, declin, gl, impo, expo, npp, time
      real remi, excr, graz, morp, morpt, morz, temp, swr, dayfrac
      real phin, dz, prca, dprca, nud, bct, tap, fo2, so2, ai, hi, hs
      real npp_D, graz_D, morp_D, no3flag, deni, nfix
      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, zmax, cont, drho
      real drhom1, wt, ahbi_cstr, ahbi_csu_dyur, gamma, rrstd, fy, fyz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "cregin.h"
      include "csbc.h"
      include "emode.h"
      include "grdvar.h"
      include "hmixc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "timeavgs.h"
      include "tmngr.h"
      include "vmixc.h"
# if defined O_save_convection || defined O_carbon_14
      include "diaga.h"
# endif
# if defined O_ice
#  if defined O_ice_cpts
      include "cpts.h"
#  endif
      include "ice.h"
# endif
# if defined O_npzd
#  if defined O_embm
      include "atm.h"
#  endif
      include "npzd.h"
# endif
# if defined O_carbon_fnpzd
      include "calendar.h"
# endif

      parameter (istrt=2, iend=imt-1)

      real twodt(km)
# if defined O_OTEC
      include "OTEC.h"

      real wti, wtm
      real inipro(km,nt), mixpro(km,nt)
      
      integer locpot

#  if defined O_OTEC_old_place
      integer NotDone, iplace, numloc

      real yr, tda, nmtpek, arepot
#  else
      integer nmtpek
#  endif
# endif
# if defined O_plume
      real  work(imt,km,jsmw:jemw)
# endif
# if defined O_npzd
      real snpzd(ntnpzd), tnpzd(ntnpzd)
# endif
# if defined O_npzd || defined O_carbon_14
      real src(imt,km,jsmw:jemw,nsrc)
# endif
# if defined O_carbon_fnpzd
      real tmpijk(imtm2,jmtm2,km), tmpi(imtm2), tmpj(jmtm2)
      real, allocatable :: fnpzd(:,:,:,:)
      save fnpzd, nfnpzd, mfnpzd, mxfnpzd
# endif
# if defined O_isopycmix
      include "isopyc.h"
# endif
      
      include "fdift.h"

# if defined O_carbon_fnpzd
!-----------------------------------------------------------------------
!     create file for writing npzd 3d fluxes of carbon and alkalinity
!-----------------------------------------------------------------------
      if (.not. allocated (fnpzd)) then
        allocate ( fnpzd(imt,jmt,km,2) )
        fnpzd(:,:,:,:) = 0.
!       initialize time record counter
        mfnpzd = 0
!       initialize records written counter
        nfnpzd = 0
!       set the repeating time interval to one year
        mxfnpzd = nint(daylen*yrlen/dtts)
        fname = new_file_name ("O_fnpzd.nc")
        inquire (file=trim(fname), exist=exists)
        if (exists) then
!         limit time interval to what exists
          call openfile (fname, iou)
          call getdimlen ('time', iou, nfnpzd)
          mxfnpzd = min(mxfnpzd, nfnpzd)
#  if defined O_carbon
          if (inqvardef('O_fnpzd1', iou)) nfnpzd = mxfnpzd
#  endif
#  if defined O_npzd_alk
          if (inqvardef('O_fnpzd2', iou)) nfnpzd = mxfnpzd
#  endif
        else
!         define file if not there
          call openfile (fname, iou)
          call redef (iou)
          call defdim ('time', iou, 0, id_time)
          call defdim ('longitude', iou, imtm2, id_xt)
          call defdim ('latitude', iou, jmtm2, id_yt)
          call defdim ('depth', iou, km, id_zt)
          it(1) = id_time
          call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &,     'time', 'time', 'years since 0-1-1')
          call putatttext (iou, 'time', 'calendar', calendar)
          it(1) = id_xt
          call defvar ('longitude', iou, 1, it, c0, c0, 'X', 'D'
     &,     'longitude', 'longitude', 'degrees')
          it(1) = id_yt
          call defvar ('latitude', iou, 1, it, c0, c0, 'Y', 'D'
     &,     'latitude', 'latitude', 'degrees')
          it(1) = id_zt
          call defvar ('depth', iou, 1, it, c0, c0, 'Z', 'D'
     &,     'depth', 'depth', 'm')
          it(1) = id_xt
          it(2) = id_yt
          it(3) = id_zt
          it(4) = id_time
#  if defined O_carbon
          call defvar ('O_fnpzd1', iou, 4, it, c0, c0, ' ', 'D'
     &,       ' ', ' ', ' ')
#  endif
#  if defined O_npzd_alk
          call defvar ('O_fnpzd2', iou, 4, it, c0, c0, ' ', 'D'
     &,       ' ', ' ', ' ')
#  endif
          call enddef (iou)
          ib(:) = 1
          ic(:) = imtm2
          tmpi(1:imtm2) = xt(2:imtm1)
          call putvara ('longitude', iou, imtm2, ib, ic, xt, c1, c0)
          ib(:) = 1
          ic(:) = jmtm2
          tmpj(1:jmtm2) = yt(2:jmtm1)
          call putvara ('latitude', iou, jmtm2, ib, ic, yt, c1, c0)
          ib(:) = 1
          ic(:) = km
          call putvara ('depth', iou, km, ib, ic, zt, c1, c0)
        endif
      endif
# endif
!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     build coefficients to minimize advection and diffusion computation
!-----------------------------------------------------------------------

#if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic || defined O_isopycmix
      limit = min(je+1+joff,jmt) - joff
      do j=js,limit
# else
      do j=js,je
# endif
        jrow = j + joff
        do i=istrt-1,iend
          cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
          cstdxt2r(i,j)   = cstr(jrow)*dxtr(i)*p5
          cstdxur(i,j)    = cstr(jrow)*dxur(i)
# if defined O_consthmix && !defined O_bryan_lewis_horizontal && !defined O_biharmonic
          ah_cstdxur(i,j) = diff_cet*cstr(jrow)*dxur(i)
# endif
        enddo
      enddo
# if defined O_plume

!-----------------------------------------------------------------------
!     find depth of penetration for the llnl simple plume model
!-----------------------------------------------------------------------

!     set parameters
      zmax = 160.0e2 ! use constant depth penetration (if cont=0)
!      cont = 0.0    ! use rho contrast for penetration (if cont>0)
      cont = 0.4e-3  ! 0.4 kg/m3 = 0.4e-3g/cm3

      do j=js,je
        jrow = j + joff
        do i=is,ie
          subz(i,jrow) = zmax
        enddo
        if (cont .gt. c0) then
          call state_ref (t(1,1,1,1,tau), t(1,1,1,2,tau)
     &,                   work(1,1,jsmw), j, j, 1, imt, 1)
          do i=is,ie
            kmx = kmt(i,jrow)
            subz(i,jrow) = zw(1)
            do k=2,kmx
              drho = work(i,k,j) - work(i,1,j)
              drhom1 = work(i,k-1,j) - work(i,1,j)
              if (drho .ge. cont .and. drhom1 .lt. cont)
     &          subz(i,jrow) = zt(k-1) + (zt(k) - zt(k-1))
     &                      *(cont - drhom1)/(drho - drhom1)
            enddo
            if (drho .le. cont) subz(i,jrow) = zw(kmx)
          enddo
        endif

!-----------------------------------------------------------------------
!     calculate "3d stf" multiplier for the llnl simple plume model
!-----------------------------------------------------------------------

        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            subz(i,jrow) = max(min(subz(i,jrow),zw(kmt(i,jrow))),zw(1))
            work(i,1,j) = c1/subz(i,jrow)
            do k=2,km
              wt =  min(max(c0,(subz(i,jrow) - zw(k-1))),dzt(k))*dztr(k)
              work(i,k,j) = wt*work(i,1,j)
            enddo
          else
            do k=1,km
              work(i,k,j) = c0
            enddo
          endif
        enddo
      enddo
# endif
# if defined O_npzd

!-----------------------------------------------------------------------
!     calculation of biological interactions
!-----------------------------------------------------------------------

      declin = sin((mod(relyr,1.) - 0.22)*2.*pi)*0.4   ! declination

      do k=1,km
        twodt(k) = c2dtts*dtxcel(k)
        nbio(k) = twodt(k)/dtnpzd
        dtbio(k) = twodt(k)/nbio(k)
        rdtts(k) = 1./twodt(k)
        rnbio(k) = 1./nbio(k)
      enddo
      tap = 2.*alpha*par

      do j=js,je
        jrow = j + joff

        do i=is,ie

          if (kmt(i,jrow) .gt. 0) then

#  if defined O_ice
#   if defined O_ice_cpts
            ai = 0.
            hi = 0.
            hs = 0.
            do n=1,ncat
              ai =  ai + A(i,jrow,2,n)
              hi =  hi + heff(i,jrow,2,n)
              hs =  hs + hseff(i,jrow,2,n)
            enddo
#   else
            ai = aice(i,jrow,2)
            hi = hice(i,jrow,2)
            hs = hsno(i,jrow,2)
#   endif
#  else
            ai = 0.
            hi = 0.
            hs = 0.
#  endif
!           calculate day fraction and incoming solar
!           angle of incidence = lat - declin, refraction index = 1.33
            rctheta = max(-1.5, min(1.5, tlat(i,jrow)/radian - declin))
            rctheta = kw/sqrt(1. - (1. - cos(rctheta)**2.)/1.33**2.)
            dayfrac = min( 1., -tan(tlat(i,jrow)/radian)*tan(declin))
            dayfrac = max(1e-12, acos(max(-1., dayfrac))/pi)
#  if defined O_embm
            swr = dnswr(i,jrow)*1e-3*(1. + ai*(exp(-ki*(hi + hs)) - 1.))
#  else
            swr = 200.
#  endif
            expo = 0.0
            impo = 0.0
            phin = 0.0 ! integrated phytoplankton
            prca = 0.0 ! integrated production of calcite

            kmx = min(kmt(i,jrow), kpzd)
            do k=1,kmx

!-----------------------------------------------------------------------
!             limit tracers to positive values
!-----------------------------------------------------------------------

              tnpzd(1) = max(t(i,k,j,ipo4,taum1), trcmin)
              tnpzd(2) = max(t(i,k,j,iphyt,taum1), trcmin)
              tnpzd(3) = max(t(i,k,j,izoop,taum1), trcmin)
              tnpzd(4) = max(t(i,k,j,idetr,taum1), trcmin)
#  if defined O_npzd_nitrogen
              tnpzd(5) = max(t(i,k,j,ino3,taum1), trcmin)
              tnpzd(6) = max(t(i,k,j,idiaz,taum1), trcmin)
#  endif

              swr = swr*exp(-kc*phin)
              phin = phin + tnpzd(2)*dzt(k)
              gl = tap*swr*exp(ztt(k)*rctheta)
              impo = expo*dztr(k)
              bct = bbio**(cbio*t(i,k,j,itemp,taum1))
#  if defined O_npzd_o2
!             decrease remineralisation rate in oxygen minimum zone
              nud = nud0*(0.65+0.35*tanh(t(i,k,j,io2,taum1)*1000.-6.))
#  else
              nud = nud0
#  endif

!-----------------------------------------------------------------------
!             call the npzd model
!-----------------------------------------------------------------------

              call npzd_src (tnpzd, nbio(k), dtbio(k), gl, bct, impo
     &,                      dzt(k), dayfrac, wd(k), rkwz(k), nud
     &,                      snpzd, expo, graz, morp, morz
#  if defined O_save_npzd
     &,                      npp, morpt, remi, excr
#   if defined O_npzd_nitrogen
     &,                      npp_D, graz_D, morp_D, nfix
#   endif
#  endif
     &                      )

              snpzd(1:4) = snpzd(1:4)*rdtts(k)
#  if defined O_npzd_nitrogen
              snpzd(5:6) = snpzd(5:6)*rdtts(k)
#  endif
              expo = expo*rnbio(k)
#  if defined O_save_npzd
              rexpo(i,k,j) = expo
              rgraz(i,k,j) = graz*rnbio(k)
              rmorp(i,k,j) = morp*rnbio(k)
              rmorz(i,k,j) = morz*rnbio(k)
              rnpp(i,k,j) = npp*rnbio(k)
              rmorpt(i,k,j) = morpt*rnbio(k)
              rremi(i,k,j) = remi*rnbio(k)
              rexcr(i,k,j) = excr*rnbio(k)
#   if defined O_npzd_nitrogen
              rnpp_D(i,k,j) = npp_D*rnbio(k)
              rgraz_D(i,k,j) = graz_D*rnbio(k)
              rmorp_D(i,k,j) = morp_D*rnbio(k)
              rnfix(i,k,j) = nfix*rnbio(k)
#   endif
#  endif
!-----------------------------------------------------------------------
!             calculate detritus at the bottom and remineralize
!-----------------------------------------------------------------------

              if (k .eq. kmt(i,jrow)) then
#  if defined O_sed
                if (addflxo .and. eots)
     &            sbc(i,jrow,irorg) = sbc(i,jrow,irorg)
     &                              + expo*redctn*twodt(1)*dzt(k)
#  endif
#  if defined O_save_npzd
                rremi(i,k,j) = rremi(i,k,j) + expo
#  endif
                snpzd(1) = snpzd(1) + redptn*expo
#  if defined O_npzd_nitrogen
                snpzd(5) = snpzd(5) + expo
#  endif
              endif

!-----------------------------------------------------------------------
!             set source terms
!-----------------------------------------------------------------------

              src(i,k,j,ispo4) = snpzd(1)
              src(i,k,j,isphyt) = snpzd(2)
              src(i,k,j,iszoop) = snpzd(3)
              src(i,k,j,isdetr) = snpzd(4)
#  if defined O_npzd_nitrogen
              src(i,k,j,isno3) = snpzd(5)
              src(i,k,j,isdiaz) = snpzd(6)
#  endif
!             production of calcite
              dprca = (morp+morz+graz*(1.-gamma1))*capr*redctn*rnbio(k)
              prca = prca + dprca*dzt(k)
#  if defined O_carbon
              src(i,k,j,isdic) = (snpzd(1)*redctp - dprca)
#  endif
#  if defined O_npzd_alk
              src(i,k,j,isalk) = (-snpzd(1)*redntp*1.e-3 - 2.*dprca)
#  endif
#  if defined O_time_averages && defined O_save_npzd

!-----------------------------------------------------------------------
!             accumulate time averages
!-----------------------------------------------------------------------

              if (timavgperts .and. .not. euler2) then
                ta_rnpp(i,k,jrow) = ta_rnpp(i,k,jrow) + rnpp(i,k,j)
                ta_rgraz(i,k,jrow) = ta_rgraz(i,k,jrow) + rgraz(i,k,j)
                ta_rmorp(i,k,jrow) = ta_rmorp(i,k,jrow) + rmorp(i,k,j)
                ta_rmorpt(i,k,jrow)= ta_rmorpt(i,k,jrow) + rmorpt(i,k,j)
                ta_rmorz(i,k,jrow) = ta_rmorz(i,k,jrow) + rmorz(i,k,j)
                ta_rexcr(i,k,jrow) = ta_rexcr(i,k,jrow) + rexcr(i,k,j)
#   if defined O_npzd_nitrogen
                ta_rnpp_D(i,k,jrow) = ta_rnpp_D(i,k,jrow)
     &                              + rnpp_D(i,k,j)
                ta_rgraz_D(i,k,jrow) = ta_rgraz_D(i,k,jrow)
     &                               + rgraz_D(i,k,j)
                ta_rmorp_D(i,k,jrow) = ta_rmorp_D(i,k,jrow)
     &                               + rmorp_D(i,k,j)
                ta_rnfix(i,k,jrow) = ta_rnfix(i,k,jrow) + rnfix(i,k,j)
#   endif
              endif
#  endif
!             calculate total export to get total import for next layer
              expo = expo*dzt(k)

            enddo

#  if defined O_npzd_o2
            kmx = kmt(i,jrow)
            do k=1,kmx
!             limit oxygen consumption below concentrations of
!             5umol/kg as recommended in OCMIP
              fo2 = 0.5*tanh(t(i,k,j,io2,taum1)*1000. - 5.)
!             sink of oxygen
              so2 = src(i,k,j,ispo4)*redotp
              src(i,k,j,iso2) = -so2*(0.5 + fo2)
#   if defined O_npzd_nitrogen
!             add denitrification as source term for NO3
              no3flag = 0.5+sign(0.5,t(i,k,j,ino3,taum1)-trcmin)
!             800 = 0.8*1000 = (elec/mol O2)/(elec/mol NO3)*(mmol/mol)
              deni = 800.*no3flag*so2*(0.5 - fo2)
              src(i,k,j,isno3) = src(i,k,j,isno3) - deni
#    if defined O_save_npzd
              rdeni(i,k,jrow) = deni
#    endif
#   endif
            enddo
#  endif

!-----------------------------------------------------------------------
!           remineralize calcite
!-----------------------------------------------------------------------
            kmx = kmt(i,jrow)
            do k=1,kmx-1
#  if defined O_carbon
              src(i,k,j,isdic) = src(i,k,j,isdic) + prca*rcak(k)
#  endif
#  if defined O_npzd_alk
              src(i,k,j,isalk) = src(i,k,j,isalk) + 2.*prca*rcak(k)
#  endif
            enddo
#  if defined O_sed
            if (addflxo .and. eots)
     &        sbc(i,jrow,ircal) = sbc(i,jrow,ircal) + (prca*rcab(kmx)
     &                          - prca*rcak(kmx))*twodt(1)*dzt(kmx)
#  endif
#  if defined O_carbon
            src(i,kmx,j,isdic) = src(i,kmx,j,isdic) + prca*rcab(kmx)
#  endif
#  if defined O_npzd_alk
            src(i,kmx,j,isalk) = src(i,kmx,j,isalk) + 2.*prca*rcab(kmx)
#  endif
#  if defined O_time_averages && defined O_save_npzd

!-----------------------------------------------------------------------
!           accumulate time averages for full depth variables
!-----------------------------------------------------------------------
            if (timavgperts .and. .not. euler2) then
              kmx = kmt(i,jrow)
              expo = prca
              ta_rprocal(i,jrow) =  ta_rprocal(i,jrow) + prca
              do k=1,kmx
                expo = expo*dztr(k)
                ta_rremi(i,k,jrow) = ta_rremi(i,k,jrow) + rremi(i,k,j)
                ta_rexpo(i,k,jrow) = ta_rexpo(i,k,jrow) + rexpo(i,k,j)
                expo = expo - prca*rcak(k)
                ta_rexpocal(i,k,jrow) = ta_rexpocal(i,k,jrow) + expo
                expo = expo*dzt(k)
#   if defined O_npzd_nitrogen && defined O_npzd_o2
                ta_rdeni(i,k,jrow) = ta_rdeni(i,k,jrow) + rdeni(i,k,j)
#   endif
              enddo
            endif
#  endif

          endif

        enddo

      enddo
# endif
# if defined O_carbon_fnpzd

!-----------------------------------------------------------------------
!    option to write and read fixed biological fluxes for carbon and
!    alkalinity. code assumes that restarts line up with the first
!    record. the repeating time period and restarts should probably
!    line up with the start of the year.
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff

!       adjust counters and read fluxes, if first row of slab
        if (jrow .eq. 2) then
          mfnpzd = mfnpzd + 1
          if (mfnpzd .gt. mxfnpzd) mfnpzd = 1
          if (nfnpzd .le. mxfnpzd) nfnpzd = nfnpzd + 1
!         read fluxes, if enough records written
          if (nfnpzd .gt. mxfnpzd) then
            fname = new_file_name ("O_fnpzd.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              ib(:) = 1
              ic(:) = 1
              ic(1) = imtm2
              ic(2) = jmtm2
              ic(3) = km
              ib(4) = mfnpzd
              print*, "Reading fnpzd from: ", trim(fname),
     &          " record:", mfnpzd
              call openfile (fname, iou)
#  if defined O_carbon
              call getvara ('O_fnpzd1', iou, imtm2*jmtm2*km, ib, ic
     &,         tmpijk, c1, c0)
              fnpzd(2:imtm1,2:jmtm1,1:km,1)=tmpijk(1:imtm2,1:jmtm2,1:km)
#  endif
#  if defined O_npzd_alk
              call getvara ('O_fnpzd2', iou, imtm2*jmtm2*km, ib, ic
     &,         tmpijk, c1, c0)
              fnpzd(2:imtm1,2:jmtm1,1:km,2)=tmpijk(1:imtm2,1:jmtm2,1:km)
#  endif
            endif
          endif
        endif

        if (nfnpzd .le. mxfnpzd) then

!         store fluxes, if enough records written
          do i=is,ie
            do k=1,km
#  if defined O_carbon
              fnpzd(i,jrow,k,1) = src(i,k,j,isdic)
#  endif
#  if defined O_npzd_alk
              fnpzd(i,jrow,k,2) = src(i,k,j,isalk)
#  endif
            enddo
          enddo
!         write fluxes, if last row of slab
          if (jrow .eq. jmt-1) then
            fname = new_file_name ("O_fnpzd.nc")
            call openfile (fname, iou)
            print*, "Writing fnpzd to: ", trim(fname),
     &          " record:", mfnpzd
            ib(:) = 1
            ic(:) = 1
            ic(1) = imtm2
            ic(2) = jmtm2
            ic(3) = km
            ib(4) = mfnpzd
            time = relyr
            call putvars ('time', iou, mfnpzd, time, c1, c0)
#  if defined O_carbon
            tmpijk(1:imtm2,1:jmtm2,1:km) = fnpzd(2:imtm1,2:jmtm1,1:km,1)
            call putvara ('O_fnpzd1', iou, imtm2*jmtm2*km, ib, ic
     &,       tmpijk, c1, c0)
#  endif
#  if defined O_npzd_alk
            tmpijk(1:imtm2,1:jmtm2,1:km) = fnpzd(2:imtm1,2:jmtm1,1:km,2)
            call putvara ('O_fnpzd2', iou, imtm2*jmtm2*km, ib, ic
     &,       tmpijk, c1, c0)
#  endif
          endif

        else

!         set fluxes to source terms
          do i=is,ie
            do k=1,km
#  if defined O_carbon
              src(i,k,j,isdic) = fnpzd(i,jrow,k,1)
#  endif
#  if defined O_npzd_alk
              src(i,k,j,isalk) = fnpzd(i,jrow,k,2)
#  endif
            enddo
          enddo

        endif

      enddo
# endif
# if defined O_carbon && defined O_carbon_14

!-----------------------------------------------------------------------
!     set source for c14
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
#  if defined O_npzd
              src(i,k,j,isc14) = src(i,k,j,isdic)*rstd
     &                         - 3.836e-12*t(i,k,j,ic14,taum1)
#  else
              src(i,k,j,isc14) = - 3.836e-12*t(i,k,j,ic14,taum1)
#  endif
            enddo
          endif
        enddo
      enddo
# endif
# if defined O_sed && !defined O_sed_uncoupled

!-----------------------------------------------------------------------
!     set source terms from sediment model
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            k = kmt(i,jrow)
#  if defined O_carbon && defined O_npzd
            src(i,k,j,isdic) = src(i,k,j,isdic) + sbc(i,j,ibdicfx)
#   if defined O_global_sums
            if (addflxo .and. eots) dtoic = dtoic - sbc(i,j,ibdicfx)
     &        *c2dtts*dtxcel(k)*dxt(i)*dyt(jrow)*cst(jrow)*dzt(k)
#   endif
#  endif
#  if defined O_npzd_alk
            src(i,k,j,isalk) = src(i,k,j,isalk) + sbc(i,j,ibalkfx)
#  endif
          endif
        enddo
      enddo
# endif

!-----------------------------------------------------------------------
!     solve for one tracer at a time
!     n = 1 => temperature
!     n = 2 => salinity
!     n > 2 => other tracers (if applicable)
!-----------------------------------------------------------------------

      do n=1,nt

!-----------------------------------------------------------------------
!       calculate advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!       calculate diffusive flux across eastern and northern faces
!       of "T" cells due to various parameterizations for diffusion.
!-----------------------------------------------------------------------

# if defined O_consthmix
#  if !defined O_biharmonic || defined O_bryan_lewis_horizontal

!       diffusive flux on eastern face of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) =
#   if defined O_bryan_lewis_horizontal
     &                         diff_cet(k)*cstdxur(i,j)*
#   else
     &                         ah_cstdxur(i,j)*
#   endif
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#   if defined O_isopycmix

!       diffusive flux on northern face of "T" cells
!       (background for isopycnal mixing)

        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) =
#    if defined O_bryan_lewis_horizontal
     &                         diff_cnt(k)*
#    else
     &                         diff_cnt*
#    endif
     &           csu_dyur(jrow)*(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#   else

!       diffusive flux on northern face of "T" cells is not calculated.
!       instead, it is incorporated into DIFF_Ty for performance issues.

#   endif
#  else

!-----------------------------------------------------------------------
!       diffusive flux across eastern face of "T" cell
!       diffusive flux across northern face of "T" cell
!-----------------------------------------------------------------------

        m = 2 + n
        do j=js,je
          jrow = j + joff
          ahbi_cstr = diff_cet*cstr(jrow)
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = ahbi_cstr*dxur(i)*
     &                         (del2(i+1,k,j,m)-del2(i,k,j,m))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          ahbi_csu_dyur = diff_cnt*csu(jrow)*dyur(jrow)
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = ahbi_csu_dyur*
     &                        (del2(i,k,j+1,m) - del2(i,k,j,m))
            enddo
          enddo
        enddo
#  endif
# else
#  if defined O_smagnlmix

!       diffusive flux on eastern and northern faces of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = diff_cet(i,k,j)*cstdxur(i,j)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = diff_cnt(i,k,j)*csu_dyur(jrow)*
     &                     (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#  endif
# endif
!-----------------------------------------------------------------------
!       calculate diffusive flux across bottom face of "T" cells
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km-1
            do i=istrt,iend
              diff_fb(i,k,j) = diff_cbt(i,k,j)*dzwr(k)*
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo

# if defined O_isopycmix

!-----------------------------------------------------------------------
!       compute isopycnal diffusive flux through east, north,
!       and bottom faces of T cells.
!-----------------------------------------------------------------------

        call isoflux (joff, js, je, is, ie, n)
# endif

!-----------------------------------------------------------------------
!       set surface and bottom vert b.c. on "T" cells for diffusion
!       and advection. for isopycnal diffusion, set adiabatic boundary
!       conditions.
!       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
!             However, it is not set when k=km so it is set below.
!             adv_fb(i,km,j) is always zero (to within roundoff).
!-----------------------------------------------------------------------

        do j=js,je
          jrow   = j + joff
          do i=istrt,iend
            kb              = kmt(i,jrow)
# if defined O_replacst
            diff_fb(i,0,j)  = c0
# else
            diff_fb(i,0,j)  = stf(i,j,n)
# endif
            diff_fb(i,kb,j) = btf(i,j,n)
            adv_fb(i,0,j)   = adv_vbt(i,0,j)*(t(i,1,j,n,tau) +
     &                                        t(i,1,j,n,tau))
            adv_fb(i,km,j)  = adv_vbt(i,km,j)*t(i,km,j,n,tau)
          enddo
        enddo

# if defined O_source_term || defined O_npzd || defined O_carbon_14

!-----------------------------------------------------------------------
!       set source term for "T" cells
!-----------------------------------------------------------------------

        source(:,:,:) = c0

#  if defined O_npzd || defined O_carbon_14
        if (itrc(n) .ne. 0) then
          do j=js,je
            do k=1,km
              do i=istrt,iend
                source(i,k,j) = src(i,k,j,itrc(n))
              enddo
            enddo
          enddo
        endif
#  endif
#  if defined O_shortwave

!-----------------------------------------------------------------------
!       incorporate short wave penetration into source
!-----------------------------------------------------------------------

        if (n .eq. 1) then
          call swflux0 (joff, js, je, istrt, iend, source)
        endif
#  endif
# endif

!-----------------------------------------------------------------------
!       solve for "tau+1" tracer using statement functions to represent
!       each component of the calculation
!-----------------------------------------------------------------------

!       1st: solve using all components which are treated explicitly

        do j=js,je
          jrow   = j + joff
          do k=1,km
            twodt(k) = c2dtts*dtxcel(k)
            do i=istrt,iend
              t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &          DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j)
     &          - ADV_Tx(i,k,j) -  ADV_Ty(i,k,j,jrow,n) -  ADV_Tz(i,k,j)
# if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
     &          - ADV_Txiso(i,k,j,n) - ADV_Tyiso(i,k,j,jrow,n)
     &          - ADV_Tziso(i,k,j)
# endif
# if defined O_source_term || defined O_npzd || defined O_carbon_14
     &          + source(i,k,j)
# endif
# if defined O_plume
     &          + work(i,k,j)*subflux(i,jrow,n)
# endif
     &          )*tmask(i,k,j)
            enddo
          enddo
        enddo
# if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion

!       2nd: add in portion of vertical diffusion handled implicitly

        call ivdift (joff, js, je, istrt, iend, n, twodt)
# endif

        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo

!-----------------------------------------------------------------------
!       construct diagnostics associated with tracer "n"
!-----------------------------------------------------------------------

        call diagt1 (joff, js, je, istrt, iend, n, twodt)

!-----------------------------------------------------------------------
!       end of tracer component "n" loop
!-----------------------------------------------------------------------

      enddo

# if defined O_replacst

      gamma = zw(1)/twodt(1)
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          stf(i,j,1) = gamma*(sbc(i,jrow,isst) - t(i,1,j,1,taup1))
          stf(i,j,2) = gamma*(sbc(i,jrow,isss) - t(i,1,j,2,taup1))
          t(i,1,j,1,taup1) = sbc(i,jrow,isst)
          t(i,1,j,2,taup1) = sbc(i,jrow,isss)
        enddo
      enddo
# endif
# if defined O_convect_brine

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!     convect brine rejected under all ice categories
!-----------------------------------------------------------------------

      call convect_brine (joff, js, je, is, ie)
# else
#  if !defined O_implicitvmix || defined O_isopycmix
#   if defined O_OTEC

!-----------------------------------------------------------------------
!     Ocean Thermal Energy Conversion code
!-----------------------------------------------------------------------

#    if defined O_OTEC_old_place
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          if (kmt(i,jrow) .ge. levupt) then
            grdavg(i,j) = t(i,levdis,j,1,taup1) - t(i,levupt,j,1,taup1)
          endif
        enddo
      enddo

      nmtplt = 0. 
      numplt(:,:) = 0.
      nmtpek = pwrpek/pwrplt  ! total at yerpek (number of plants)
      yr = relyr + year0      ! current year

!     calculate the number of available plants over time      
!     assume build capacity is ramped up to meet demand in yerpek (nmtpek)
      if (yr .lt. yerstr) then
        nmtblt = int(pwrstr/pwrplt)
      elseif (yr .lt. yerpek .and. yerpek .gt. yerstr) then
        nmtblt = nmtpek*(0.5 - 0.5*cos(pi*(yr-yerstr)/(yerpek-yerstr)))
      elseif (yr .lt. yerend .and. yerend .gt. yerstr) then
        nmtblt = nmtpek*(0.5 + 0.5*cos(pi*(yr-yerpek)/(yerend-yerpek)))
      else
        nmtblt = int(pwrend/pwrplt)
      endif

!     estimate plant demand from emissions (done in co2data)
      nmtdem = elcdem/pwrplt

!     limit build to demand
      nmtblt = min(nmtblt, nmtdem)

!     find all possible locations of viable plants
      arepot = 0.
      numloc = 0
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          if (kmt(i,jrow) .ge. levupt) then
            if (grdavg(i,jrow) .ge. grdmin) then
              arepot = arepot + tgarea(i,jrow)
              numloc = numloc + 1
              indplt(numloc,1) = i
              indplt(numloc,2) = j
            endif
          endif
        enddo
      enddo

!     place or relocate plants to optimal locations (largest dtdu)
      nmtplt = 0
      artplt = 0.
      if (arepot .gt. 0) then 
        tda = pi*disrad**2
        n = 1
        m = 26
        iplace = 0
        NotDone = 1
        n = 1
        do while (NotDone .ne. 0)
          i = indplt(n,1)
          j = indplt(n,2)
          jrow = j + joff
          if (grdavg(i,jrow) .ge. m) then
            if (tgarea(i,jrow) - tda .gt. numplt(i,jrow)*tda) then
              iplace = iplace + 1
              numplt(i,jrow) = numplt(i,jrow) + 1
              nmtplt = nmtplt + 1
              artplt = artplt + tda
            endif
          endif
!         exit if number of plants meet number built or demand
          if (nmtplt .ge. nmtblt .or. nmtplt .ge. nmtdem) exit
          n = n + 1
          if (n .eq. numloc) then
            n = 1
            m = m - 1
            if (m .lt. grdmin) m = grdmin
            if (iplace .eq. 0) then
              m = m - 1
              if (m .ge. grdmin) then
                NotDone = 1
              else
                NotDone = 0
              endif
            endif
          endif
        enddo
      endif
      
      if (nmtblt .gt. 0) then
        nmtplt = 0.
        pwrtot = 0.
#    else
!     accumulate temperature gradient for averaging
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          if (kmt(i,jrow) .ge. levupt) then
            grdacc(i,jrow) = grdacc(i,jrow) + t(i,levdis,j,1,taup1)
     &                  - t(i,levupt,j,1,taup1)
          endif
        enddo
      enddo
      nmtacc = nmtacc + 1

      if (eoyear .and. eots .and. je+joff .ge. jmtm1) then
!       update annual average temperature gradient
        do j=2,jmtm1
          do i=2,imtm1
            grdavg(i,j) = grdacc(i,j)/nmtacc
            if (int(relyr + year0) .eq. int(yerstr)) then
              grdstr(i,jrow) = grdavg(i,jrow)
            endif
            grdacc(i,j) = 0.
          enddo
        enddo
        nmtacc = 0.

#     if defined O_OTEC_sort_place
        nmtloc = 0
        do j=2,jmtm1
          do i=2,imtm1
            if (kmt(i,j) .ge. levupt) then
              if (grdavg(i,j).ge.grdmin .and. locmap(i,j).gt.1) then
                nmtloc = nmtloc + 1
!               save temperature gradient for sorting
                locpot = int(tgarea(i,j)/aremin)
                lctpot = lctpot + locpot
                locrem(i,j) = locpot - numplt(i,j)
#      if defined O_OTEC_adjust_gradient
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
#      endif
#      if defined O_OTEC_scale_gradient
!               scale the useable gradient with the portion of area remaining
!               in order to reduce overloading a grid cell with plants
                grdavg(i,j) = (grdavg(i,j) - grdmin)*locrem(i,j)
     &                   /float(int(tgarea(i,j)/aremin)) + grdmin
#      endif
                grdsor(nmtloc) = grdavg(i,j)
                indplt(nmtloc,1) = i
                indplt(nmtloc,2) = j
              else
!               need to move due to inadequate gradient or location no longer allowed
                nmtmov = nmtmov + numplt(i,j)
                numplt(i,j) = 0
              endif
            endif
          enddo
        enddo
!       sort temperature gradient and indices for placement
        call quicksort(grdsor(1:nmtloc), indplt(1:nmtloc,1)
     &,                indplt(1:nmtloc,2), 1, nmtloc)
#     endif
      endif

      if (nmtblt .gt. 0) then
        pwrtot = nmtplt*pwrplt
        artplt = nmtplt*aremin
#    endif
!       count number of averaging times
        if (timavgperts .and. eots .and. joff .eq. 0) 
     &    nta_OTEC = nta_OTEC + 1      

        do j=js,je
          jrow = j + joff
          do i=istrt,iend

            if (numplt(i,jrow) .ne. 0) then
            
!             save data for time averaging 
              if (timavgperts .and. eots) then
                ta_numplt(i,jrow) = ta_numplt(i,jrow) + numplt(i,j)
                ta_locrem(i,jrow) = ta_locrem(i,jrow) + locrem(i,j)
                ta_grdavg(i,jrow) = ta_grdavg(i,jrow) + grdavg(i,j)
              endif

!             initialise the profiles
              do k=1,kmt(i,jrow)
                do n=1,nt
                  inipro(k,n) = t(i,k,j,n,taup1)
#    if !defined O_OTEC_old_mix
                mixpro(k,n) = 0.
#    endif
                enddo
              enddo

              do m=1,nmtdis-1

                call OTEC (inipro(1,1), i, j, joff, m, kmt)

                if (timavgperts .and. eots) then
                  wtm = disare(m)/dissum(nmtdis)
                  ta_totlev(i,jrow) = ta_totalk(i,jrow)+totalk(i,j)*wtm
                  ta_vendth(i,jrow) = ta_vendth(i,jrow)+vdepth(i,j)*wtm
                  ta_poterg(i,jrow) = ta_poterg(i,jrow)+pe(i,j)*wtm
                endif

#    if defined O_OTEC_old_mix
                wti = disare(m)*numplt(i,j)/tgarea(i,jrow)
#    else
                wtm = disare(m)/(dissum(m))
#     if defined O_OTEC_mixrings
                wti = disare(m+1)/(dissum(m+1))
#     endif
#    endif
                do k=1,kmt(i,jrow)
                  do n=1,nt
!                   mix the previous profiles
#    if defined O_OTEC_old_mix
                    t(i,k,j,n,taup1) = t(i,k,j,n,taup1)*(1.-wti)
     &                               + inipro(k,n)*wti
                    inipro(k,n) = t(i,k,j,n,taup1)
#    else
                    mixpro(k,n) = mixpro(k,n)*(1-wtm)
     &                          + inipro(k,n)*wtm
#     if defined O_OTEC_mixrings
!                   set the profile for the next ring of discharge
                    inipro(k,n) = mixpro(k,n)*(1.-wti)
     &                          + t(i,k,j,n,taup1)*wti
#     else
                    inipro(k,n) = t(i,k,j,n,taup1)
#     endif
#    endif
                  enddo
                enddo

              enddo

!             calculate the final profile
              call OTEC (inipro(1,1), i, j, joff, m, kmt)

              if (timavgperts .and. eots) then
                wtm = disare(nmtdis)/dissum(nmtdis)
                ta_totlev(i,jrow) = ta_totalk(i,jrow)+totalk(i,j)*wtm
                ta_vendth(i,jrow) = ta_vendth(i,jrow)+vdepth(i,j)*wtm
                ta_poterg(i,jrow) = ta_poterg(i,jrow)+pe(i,j)*wtm
              endif

#    if defined O_OTEC_old_mix
              wti = disare(nmtdis)*numplt(i,j)/tgarea(i,jrow)
#    else
              wtm = disare(nmtdis)/(dissum(nmtdis))
              wti = dissum(nmtdis)*numplt(i,j)/tgarea(i,jrow)
#    endif
              do k=1,kmt(i,jrow)
                do n=1,nt
!                 mix all discharge over the entire grid cell 
#    if defined O_OTEC_old_mix
                  t(i,k,j,n,taup1) = t(i,k,j,n,taup1)*(1.-wti)
     &                             + inipro(k,n)*wti
#    else
                  mixpro(k,n) = mixpro(k,n)*(1-wtm)
     &                        + inipro(k,n)*wtm
                  t(i,k,j,n,taup1) = t(i,k,j,n,taup1)*(1.-wti)
     &                             + mixpro(k,n)*wti
#    endif
                enddo
              enddo

#    if defined O_OTEC_old_place
              pwrtot = pwrtot + numplt(i,jrow)*pwrplt
              nmtplt = nmtplt + numplt(i,j)
#    endif
            endif

          enddo
        enddo

#    if defined O_OTEC_old_eff
        pwrtot = pwrtot*0.363
#    endif
#    if defined O_OTEC_old_out
        if (tsits .and. eots) then
#     if defined O_OTEC_old_place
          print*, "OTEC: number of viable locations: ",numloc
#     endif
          print*, "OTEC: estimate of power demand (TW): ",elcdem/1e12
          print*, "OTEC: number required for power demand: ",nmtdem
          print*, "OTEC: number of plants: ",nmtplt
#     if defined O_OTEC_old_place
          print*, "OTEC: percent of area used: ",100.*artplt/arepot
#     endif
          print*, "OTEC: emissions saved (PgC): ", emtred*3.156e-5
          print*, "OTEC: power produced:(TW)", pwrtot/1e12
        endif
#    endif

      endif
#   endif
!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!-----------------------------------------------------------------------

#   if defined O_fullconvect
      call convct2 (t(1,1,1,1,taup1), joff, js, je, is, ie, kmt)
      do j=js,je
        do n=1,nt
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
#   else
      call convct (t(1,1,1,1,taup1), ncon, joff, js, je, istrt, iend
     &,            kmt)
#   endif
#  endif
# endif
# if defined O_save_convection
      if (timavgperts .and. eots) then
        if (joff .eq. 0) nta_conv = nta_conv + 1
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            ta_totalk(i,jrow) = ta_totalk(i,jrow) + totalk(i,j)
            ta_vdepth(i,jrow) = ta_vdepth(i,jrow) + vdepth(i,j)
            ta_pe(i,jrow) = ta_pe(i,jrow) + pe(i,j)
          enddo
        enddo
      endif
# endif

!-----------------------------------------------------------------------
!     construct diagnostics after convection
!-----------------------------------------------------------------------

      idiag = 10
      call diagt2 (joff, js, je, istrt, iend, idiag)

# if defined O_fourfil || defined O_firfil

!-----------------------------------------------------------------------
!     filter tracers at high latitudes
!-----------------------------------------------------------------------

      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filt (joff, js, je)
      else
        write (stdout,'(a)')
     &  'Error: filtering requires is=2 and ie=imt-1 in tracer'
        stop '=>tracer'
      endif
# endif
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct diagnostics after filtering (for total dT/dt)
!-----------------------------------------------------------------------

      idiag = 1
      call diagt2 (joff, js, je, istrt, iend, idiag)

# if !defined O_replacst

!-----------------------------------------------------------------------
!     if needed, construct the Atmos S.B.C.(surface boundary conditions)
!     averaged over this segment
!     eg: SST and possibly SSS
!-----------------------------------------------------------------------

      call asbct (joff, js, je, istrt, iend, isst, itemp)
      call asbct (joff, js, je, istrt, iend, isss, isalt)
# endif
# if defined O_carbon
      call asbct (joff, js, je, istrt, iend, issdic, idic)
#  if defined O_carbon_14
      call asbct (joff, js, je, istrt, iend, issc14, ic14)
#  endif
# endif
# if defined O_npzd_alk
      call asbct (joff, js, je, istrt, iend, issalk, ialk)
# endif
# if defined O_npzd_o2
      call asbct (joff, js, je, istrt, iend, isso2, io2)
# endif
# if defined O_npzd
      call asbct (joff, js, je, istrt, iend, isspo4, ipo4)
#  if !defined O_npzd_no_vflux
      call asbct (joff, js, je, istrt, iend, issphyt, iphyt)
      call asbct (joff, js, je, istrt, iend, isszoop, izoop)
      call asbct (joff, js, je, istrt, iend, issdetr, idetr)
#  endif
#  if defined O_npzd_nitrogen
      call asbct (joff, js, je, istrt, iend, issno3, ino3)
#   if !defined O_npzd_no_vflux
      call asbct (joff, js, je, istrt, iend, issdiaz, idiaz)
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      call asbct (joff, js, je, istrt, iend, isscfc11, icfc11)
      call asbct (joff, js, je, istrt, iend, isscfc12, icfc12)
# endif
# if defined O_sed
!-----------------------------------------------------------------------
!     accumulate bottom tracer values for the sediment model
!-----------------------------------------------------------------------

      if (addflxo .and. eots) then
        if (joff .eq. 0) atsed = atsed + twodt(1)
        do j=js,je
          jrow = j + joff
          do i=is,ie
            if (kmt(i,jrow) .gt. 0) then
              k = kmt(i,jrow)
              sbc(i,jrow,ibtemp) = sbc(i,jrow,ibtemp)
     &                           + t(i,k,j,itemp,taup1)*twodt(1)
              sbc(i,jrow,ibsalt) = sbc(i,jrow,ibsalt)
     &                           + t(i,k,j,isalt,taup1)*twodt(1)
#  if defined O_carbon
              sbc(i,jrow,ibdic) = sbc(i,jrow,ibdic)
     &                          + t(i,k,j,idic,taup1)*twodt(1)
#  endif
#  if defined O_npzd_alk
              sbc(i,jrow,ibalk) = sbc(i,jrow,ibalk)
     &                          + t(i,k,j,ialk,taup1)*twodt(1)
#  endif
#  if defined O_npzd_o2
              sbc(i,jrow,ibo2) = sbc(i,jrow,ibo2)
     &                         + t(i,k,j,io2,taup1)*twodt(1)
#  endif
            endif
          enddo
        enddo
      endif
# endif
# if defined O_carbon && defined O_carbon_14

!-----------------------------------------------------------------------
!     calculate diagnostic delta carbon 14
!-----------------------------------------------------------------------

      if (tsiperts .or. timavgperts) then
        rrstd = 1000./rstd
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              dc14(i,k,j) = (rrstd*t(i,k,j,ic14,taup1)
     &                      /(t(i,k,j,idic,taup1) + epsln) - 1000.)
     &                      *tmask(i,k,j)
            enddo
          enddo
        enddo
      endif
      if (tsiperts .and. eots) then
        if (js+joff .eq. 2) dc14bar = 0.
        do j=js,je
          jrow = j + joff
          fy = cst(jrow)*dyt(jrow)
          do k=1,km
            fyz = fy*dzt(k)
            do i=istrt,iend
              dc14bar = dc14bar + dc14(i,k,j)*dxt(i)*fyz*tmask(i,k,j)
            enddo
          enddo
        enddo
      endif
      if (timavgperts .and. .not. euler2) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              ta_dc14(i,k,jrow) = ta_dc14(i,k,jrow) + dc14(i,k,j)
            enddo
          enddo
        enddo
      endif
# endif

      return
      end

      subroutine diagt1 (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     construct diagnostics associated with tracer component "n"

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = (1,2) = (u,v) velocity component
!       twodt = (2*dtts,dtts) on (leapfrog,mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, js, je, joff, is, ie
      integer mask, m

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, dtdx, dtdy, dtdz
      real r2dt, cosdyt, fx, darea, boxar, rtwodt, sumdx, delx
      real sumdxr, dxdy, dxdydz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "cregin.h"
      include "csbc.h"
# if defined O_tracer_averages
      include "ctavg.h"
# endif
      include "diag.h"
      include "diaga.h"
      include "emode.h"
      include "grdvar.h"
      include "hmixc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"

# if defined O_meridional_tracer_budget
      include "ctmb.h"
      real stor(imt,km), div(imt,km), sorc(imt,km), dif(imt,km)
      real t1(imt), t2(imt), t3(imt), t4(imt)
# endif
# if defined O_time_step_monitor
      real temp1(imt,km), temp2(imt,km), temp3(imt,km)
# endif
      real twodt(km)

# if defined O_isopycmix
      include "isopyc.h"
# endif
      include "fdift.h"

# if defined O_save_mixing_coeff

!-----------------------------------------------------------------------
!     diagnostic: estimate mixing coefficients on east, north, and
!                 bottom face of T cells from the flux
!-----------------------------------------------------------------------

      if (cmixts .and. n .eq. 1 .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=2,imt-1
              dtdx = (t(i+1,k,j,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i+1,k,j)*tmask(i,k,j)
     &                *cstr(jrow)*dxur(i) + epsln
              ce(i,k,j,2) = diff_fe(i,k,j)/dtdx
#  if !defined O_consthmix || defined O_biharmonic || defined O_isopycmix
              dtdy = (t(i,k,j+1,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i,k,j+1)*tmask(i,k,j)
     &                *dyur(jrow) + epsln
              cn(i,k,j,2) = csur(jrow)*diff_fn(i,k,j)/dtdy
#  else
              cn(i,k,j,2) = csur(jrow)*ah
#  endif
            enddo
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              dtdz = (t(i,k,j,1,taum1)-t(i,k+1,j,1,taum1))
     &                *tmask(i,k,j)*tmask(i,k+1,j)
     &                *dzwr(k) + epsln
              cb(i,k,j,2) = diff_fb(i,k,j)/dtdz
#  if defined O_isopycmix
     &                      + diff_fbiso(i,k,j)/dtdz
#  endif
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,2) = 0.0
          enddo
        enddo

        do j=js,je
          call setbcx (ce(1,1,j,2), imt, km)
          call setbcx (cn(1,1,j,2), imt, km)
          call setbcx (cb(1,1,j,2), imt, km)
        enddo
      endif
# endif

# if defined O_save_convection_full

!-----------------------------------------------------------------------
!     diagnostic: initialize temperature before convection
!-----------------------------------------------------------------------

      if (exconvts .and. n .eq. 1 .and. eots) then
        do j=js,je
          do k=1,km
            do i=1,imt
              excnv0(i,k,j) = t(i,k,j,1,taup1)
            enddo
          enddo
        enddo
      endif

# endif

# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     diagnostic: integrate |d(tracer)/dt|  and tracer variance on "tau"
!                 globally
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) then
        do j=js,je
          jrow = j + joff
          r2dt    = c1/c2dtts
          cosdyt  = cst(jrow)*dyt(jrow)
          do k=1,km
            fx = r2dt/dtxcel(k)
            do i=is,ie
              darea      = dzt(k)*dxt(i)*cosdyt*tmask(i,k,j)
              temp3(i,k) = t(i,k,j,n,tau)*darea
              temp1(i,k) = t(i,k,j,n,tau)**2*darea
              temp2(i,k) = abs(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))*
     &                     darea*fx
            enddo
            do i=is,ie
              tbar(k,n,jrow)   = tbar(k,n,jrow) + temp3(i,k)
              travar(k,n,jrow) = travar(k,n,jrow) + temp1(i,k)
              dtabs(k,n,jrow)  = dtabs(k,n,jrow) + temp2(i,k)
            enddo
          enddo
        enddo
      endif
# endif

# if defined O_tracer_averages

!-----------------------------------------------------------------------
!     diagnostic: accumulate tracers for averages under horizontal
!                 regions (use units of meters, rather than cm)
!-----------------------------------------------------------------------

      if (tavgts .and. eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            mask = mskhr(i,jrow)
            if (mask .ne. 0) then
              boxar = cst(jrow)*dxt(i)*dyt(jrow)*tmask(i,1,j)*0.0001
              sumbf(mask,n) = sumbf(mask,n) + stf(i,j,n)*boxar
              do k=1,km
                sumbk(mask,k,n) = sumbk(mask,k,n) + t(i,k,j,n,tau)
     &                             *boxar*dzt(k)*tmask(i,k,j)*0.01
             enddo
            endif
          enddo
        enddo
      endif
# endif

# if defined O_tracer_yz

!----------------------------------------------------------------------
!     diagnostic: integrate tracer equations in longitude
!----------------------------------------------------------------------

      if (tyzts .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            rtwodt = c1/twodt(k)
            sumdx = c0
            do m=1,5
              tyz(jrow,k,n,m) = c0
            enddo
            do i=is,ie
              delx  = dxt(i)*tmask(i,k,j)
              sumdx = sumdx + delx
              tyz(jrow,k,n,1) = tyz(jrow,k,n,1)+ delx*t(i,k,j,n,tau)
              tyz(jrow,k,n,2) = tyz(jrow,k,n,2) + rtwodt*delx*
     &          (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              tyz(jrow,k,n,3) =  tyz(jrow,k,n,3)
     &          - delx*(ADV_Ty(i,k,j,jrow,n) + ADV_Tz(i,k,j))
              tyz(jrow,k,n,4) =  tyz(jrow,k,n,4)
     &          + delx*(DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j))
#  if defined O_source_term || defined O_npzd || defined O_carbon_14
              tyz(jrow,k,n,5) = tyz(jrow,k,n,5) + delx*source(i,k,j)
#  endif
            enddo
            sumdxr = c1/(sumdx + epsln)
            do m=1,5
              tyz(jrow,k,n,m) = tyz(jrow,k,n,m)*sumdxr
            enddo
          enddo
        enddo
        if (js+joff .eq. 2) then
          do m=1,5
            do k=1,km
              tyz(1,k,n,m)   = c0
              tyz(jmt,k,n,m) = c0
            enddo
          enddo
        endif
      endif
# endif

# if defined O_meridional_tracer_budget

!----------------------------------------------------------------------
!     diagnostic: integrate equations in depth, lon, and time for
!                 each basin and jrow. basin #  0 is used to catch
!                 values in land areas
!----------------------------------------------------------------------

      if (tmbperts .and. eots) then
        do j=js,je
          jrow = j + joff
          if (jrow .eq. 2 .and. n .eq. 1) numtmb = numtmb + 1
          cosdyt = cst(jrow)*dyt(jrow)
          do k=1,km
            do i=is,ie
              stor(i,k) = 0.0
              div(i,k)  = 0.0
              sorc(i,k) = 0.0
              dif(i,k)  = 0.0
            enddo
          enddo
          do k=1,km
            rtwodt = c1/twodt(k)
            do i=is,ie
              dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
              dxdydz          = dxdy*dzt(k)
              stor(i,k) = stor(i,k) + rtwodt*dxdydz*
     &                            (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              div(i,k)  = div(i,k)  - dxdydz*ADV_Ty(i,k,j,jrow,n)
#  if defined O_source_term || defined O_npzd || defined O_carbon_14
              sorc(i,k) = sorc(i,k) + dxdydz*source(i,k,j)
#  endif
              dif(i,k)  = dif(i,k)  + dxdydz*DIFF_Ty(i,k,j,jrow,n)
            enddo
          enddo

!         the accumulation is done this way for issues of speed

          do i=is,ie
            t1(i) = stor(i,1)
            t2(i) = div(i,1)
            t3(i) = sorc(i,1)
            t4(i) = dif(i,1)
          enddo

          do k=2,km
            do i=is,ie
              t1(i) = t1(i) + stor(i,k)
              t2(i) = t2(i)  + div(i,k)
              t3(i) = t3(i) + sorc(i,k)
              t4(i) = t4(i)  + dif(i,k)
            enddo
          enddo

          do i=is,ie
            m               = msktmb(i,jrow)
            tstor(jrow,n,m) = tstor(jrow,n,m) + t1(i)
            tdiv(jrow,n,m)  = tdiv(jrow,n,m)  + t2(i)
            tsorc(jrow,n,m) = tsorc(jrow,n,m) + t3(i)
            tdif(jrow,n,m)  = tdif(jrow,n,m)  + t4(i)
          enddo

          k = 1
          do i=is,ie
            m               = msktmb(i,jrow)
            dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
            tflux(jrow,n,m) = tflux(jrow,n,m) + dxdy*stf(i,j,n)
          enddo
          if (n .eq. 1) then
            do k=1,km
              do i=is,ie
                m               = msktmb(i,jrow)
                dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
                dxdydz          = dxdy*dzt(k)
                smdvol(jrow,m)  = smdvol(jrow,m)  + dxdydz
              enddo
            enddo
          endif
        enddo
      endif
# endif

# if defined O_gyre_components

!-----------------------------------------------------------------------
!     diagnostic: compute the northward transport components of
!                 each tracer
!-----------------------------------------------------------------------

      if (gyrets .and. eots)  call gyre (joff, js, je, is, ie, n)
# endif

# if defined O_term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate r.h.s. terms in the tracer equations
!                 over specified regional volumes.
!-----------------------------------------------------------------------

      if (trmbts .and. eots)  call ttb1 (joff, js, je, is, ie, n)
# endif

# if defined O_xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s. terms in the tracer equations
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call txbt1 (joff, js, je, n)
# endif
# if defined O_mom_tbt

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s. terms in the tracer equations
!-----------------------------------------------------------------------

      if (tbtperts .and. eots) call tbt1 (joff, js, je, n)
# endif

      return
      end

      subroutine diagt2 (joff, js, je, is, ie, idiag)

!-----------------------------------------------------------------------
!     construct d(tracer)/dt diagnostics

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       idiag = 1  => total tracer change
!       idiag = 10 => change of tracer due to filtering(also convection)
!-----------------------------------------------------------------------

      implicit none

      integer idiag, j, js, je, k, i, joff, iocv, jrow, is, ie

      real rdt, reltim, period

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "diaga.h"
      include "iounit.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "timeavgs.h"

# if defined O_save_convection_full

!-----------------------------------------------------------------------
!     diagnostic: save tracer time change due to explicit convection
!     idiag = 10 signifies diagnostics immediately after convection
!-----------------------------------------------------------------------

      if (exconvts .and. idiag .eq. 10 .and. eots) then
        rdt = c1/c2dtts

!       excnv1 = epsln over land cells and d(convect)/dt over ocean

        do j=js,je
          do k=1,km
            do i=1,imt
              excnv1(i,k,j) = tmask(i,k,j)*
     &         (t(i,k,j,1,taup1)-excnv0(i,k,j))*rdt
     &                       + (c1-tmask(i,k,j))*epsln
            enddo
          enddo
        enddo

        reltim = relyr
        if (joff + js .eq. 2) then
          write (stdout,*) ' =>Writing explicit convection at ts=',itt
     &  , ' ',stamp
          call getunit (iocv, 'cvct.dta'
     &,                'unformatted sequential append ieee')

          period = 0.0
          iotext = 'read(iocv) reltim, period, imt, jmt, km, flag'
          write (iocv) stamp, iotext, expnam
          write (iocv) reltim, period, imt, jmt, km, epsln

          iotext = 'read(iocv) (xt(i),i=1,imt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, xt, imt)

          iotext = 'read(iocv) (yt(j),j=1,jmt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, yt, jmt)

          iotext = 'read(iocv) (zt(k),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, zt, km)

          call relunit (iocv)
        endif

        call getunit (iocv, 'cvct.dta'
     &,                'unformatted sequential append ieee')

        do j=js,je
          jrow = j+joff

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocv) (convect(i,k),i=1,imt),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, excnv1(1,1,j), imt*km)

        enddo
        call relunit(iocv)
      endif
# endif

# if defined O_term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate d/dt(tracer) over specified regional volumes
!                  after convection and filtering
!-----------------------------------------------------------------------

      if (trmbts .and. eots) call ttb2 (joff, js, je, is, ie, idiag)
# endif

# if defined O_xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate d/dt(tracer) after convection
!                 and filtering
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call txbt2 (joff, js, je, idiag)
# endif
# if defined O_mom_tbt

!-----------------------------------------------------------------------
!     diagnostic: accumulate d/dt(tracer) after convection
!                 and filtering

!-----------------------------------------------------------------------

      if (tbtperts .and. eots) call tbt2 (joff, js, je, idiag)
# endif

      return
      end

      subroutine asbct (joff, js, je, is, ie, isbc, itr)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C. (surface boundary conditions)

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       isbc  = index for sbc
!       itr   = index for tracer
!-----------------------------------------------------------------------

      implicit none

      integer isbc, itr, j, js, je, jrow, joff, i, is, ie

      real rts

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

!     initialize the Atmos S.B.C. at the start of each ocean segment
!     (do not alter values in land)

      if (isbc .le. 0 .or. itr .le. 0) return

      if (eots .and. osegs) then
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) sbc(i,jrow,isbc) = c0
          enddo
        enddo
      endif

!     accumulate surface tracers for the Atmos S.B.C. every time step

      if (eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            sbc(i,jrow,isbc) = sbc(i,jrow,isbc)+t(i,1,j,itr,taup1)
          enddo
        enddo
      endif

!     average the surface tracers for the Atmos S.B.C. at the end of
!     each ocean segment. (do not alter values in land)

      if (eots .and. osege) then
        rts = c1/ntspos
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0)
     &        sbc(i,jrow,isbc) = rts*sbc(i,jrow,isbc)
          enddo
        enddo
      endif

      return
      end

      subroutine ivdift (joff, js, je, is, ie, n, twodt)

# if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
!-----------------------------------------------------------------------
!     solve vertical diffusion of tracers implicitly

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = tracer component
!       twodt = (2*dtts, dtts) on (leapfrog, mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer j, js, je, k, i, is, ie, n, joff

      real rc2dt

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "mw.h"
      include "switch.h"
      include "vmixc.h"

      real twodt(km)

!     store terms to compute implicit vertical mixing on
!     diagnostic time steps

#  if defined O_xbts || defined O_mom_tbt
      if (eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#  else
#   if defined O_term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#   endif
#  endif

      call invtri (t(1,1,1,n,taup1), stf(1,1,n), btf(1,1,n)
     &, diff_cbt(1,1,jsmw), twodt, kmt, tmask(1,1,1), is, ie
     &, joff, js, je)

!     compute residual implicit vertical mixing

#  if defined O_xbts || defined O_mom_tbt
      if (eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
#  else
#   if defined O_term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
#   endif
#  endif
# endif

      return
      end

      subroutine swflux0 (joff, js, je, is, ie, source)

# if defined O_mom && defined O_shortwave
!=======================================================================
!       incorporate short wave penetration via the "source"
!       term. note that the divergence of shortwave for the first
!       level "divpen(1)" is compensating for the effect of having
!       the shortwave component already included in the total
!       surface tracer flux "stf(i,j,1)"

!       incorporating the penetrative shortwave radiative flux into
!       the vertical diffusive flux term directly is not correct when
!       vertical diffusion is solved implicitly.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

!     output:
!       source = source term penetrative short wave
!=======================================================================

      implicit none

      integer j, js, je, jrow, joff, k, i, is, ie

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "cshort.h"

      real source(imt,km,jsmw:jemw)

      if (iswr .ne. 0) then
        do j=js,je
          jrow   = j + joff
          do k=1,km
            do i=is,ie
              source(i,k,j) = source(i,k,j)
     &                       + sbc(i,jrow,ipsw)*divpen(k)
            enddo
          enddo
        enddo
      endif
# endif
#endif

      return
      end
      
      recursive subroutine quicksort(a, b, c, first, last)
        implicit none
        real  a(*), x, t
        integer b(*), c(*), first, last
        integer i, j

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
           do while (a(i) < x)
              i=i+1
           end do
           do while (x < a(j))
              j=j-1
           end do
           if (i >= j) exit
           t = a(i);  a(i) = a(j);  a(j) = t
           t = b(i);  b(i) = b(j);  b(j) = t
           t = c(i);  c(i) = c(j);  c(j) = t
           i=i+1
           j=j-1
        end do
        if (first < i-1) call quicksort(a, b, c, first, i-1)
        if (j+1 < last)  call quicksort(a, b, c, j+1, last)
      end subroutine quicksort
