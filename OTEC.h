!======================= include file "OTEC.h" =========================

!     variables for use with OTEC power generation

!     mxndis = maximum number of discrete discharge volumes 
!     nmtdis = total number of discrete discharge volumes
!     levdis = level of discharge
!     levupt = level of uptake
!     numplt = number of plants in a grid cell
!     nmtplt = total number of plants
!     nmtacc = total number of gradient accumulations (for averaging)
!     nmtpip = total number of pipes per plant
!     nmtdem = total number of plants to reach current demand
!     nmtblt = total number of plants built
!     nmtmov = total number of plants needed to be moved
!     nmtloc = total number of grid locations
!     locrem = number of locations remaining in a grid cell
!     lctrem = total number of locations remaining in total
!     lctpot = total number of potential locations
!     locmap = map of potential locations
!     indplt = grid indices of possible plant locations

!     dptdis = depth of discharge (cm)
!     dptupt = depth of uptake (cm)
!     piprad = radius of pipe (cm)
!     pipare = area of pipe (cm2)
!     pipvol = volume of discharge from a pipe (cm2)
!     pipvel = velocity of water the in pipe (cm2 s-1)
!     dismax = maximum volume of dischage (cm3 s-1)
!     disrad = radius of discharge(cm)
!     disare = area of discharge (cm2)
!     disvol = volume of discharge (cm3)
!     disvel = velocity of discharge (cm2 s-1)
!     pwrstr = starting power (zero before yerstr)(W)
!     pwrpek = peak power at maximum deployment (W)
!     pwrend = ending power (continues after yerend) (W)
!     pwrplt = power produced by a single plant (W)
!     pwrtot = total power produced by all plants (W)
!     yerstr = start year (year)
!     yerpek = peak year (year)
!     yerend = final year (year)
!     grdmin = minimum tempeature gradient (K or C)
!     grdacc = temperature gradient accumulator (K or C)
!     grdavg = annual average temperature gradient (K or C)
!     grdstr = starting average temperature gradient at yerstr (K or C)
!     grdsor = sorted temperature gradient (K or C)
!     artplt = total area of plants (cm2)
!     aremin = minimum area for a plant (cm2)
!     OTEC_yr = OTEC forcing year (year)

      integer mxndis
      parameter (mxndis=50)

# if defined O_OTEC_old

      integer nmtdis, levdis, levupt, numplt, nmtplt
      integer nmtacc, nmtpip, locrem
      integer nmtmov, nmtloc, lctrem, lctpot, locmap
      integer indplt
      
      common /OTEC_i/ nmtdis, levdis, levupt, numplt(imt,jmt), nmtplt
      common /OTEC_i/ nmtacc, nmtpip, locrem(imt,jmt)
      common /OTEC_i/ nmtmov, nmtloc, lctrem, lctpot, locmap(imt,jmt)
      common indplt(imtm2*jmtm2,2)

      real nmtdem, nmtblt
      common /OTEC_r/ nmtdem, nmtblt
      
# else
      integer nmtdis, levdis, levupt, numplt, nmtplt
      integer nmtacc, nmtpip, locrem, nmtdem, nmtblt
      integer nmtmov, nmtloc, lctrem, lctpot, locmap
      integer indplt
      
      common /OTEC_i/ nmtdis, levdis, levupt, numplt(imt,jmt), nmtplt
      common /OTEC_i/ nmtacc, nmtpip, locrem(imt,jmt), nmtdem, nmtblt
      common /OTEC_i/ nmtmov, nmtloc, lctrem, lctpot, locmap(imt,jmt)
      common indplt(imtm2*jmtm2,2)
#endif

      real dptdis, dptupt, piprad, pipare, pipvol, pipvel      
      real disare, disvol, disvel
      real dissum, dismax, disrad, pwrstr, pwrpek
      real pwrend, pwrplt, pwrtot, yerstr, yerpek, yerend
      real emtred, elcdem, grdmin, grdsor
      real grdacc, grdavg, grdstr
      real artplt, aremin, OTEC_yr

      common /OTEC_r/ dptdis, dptupt, piprad, pipare, pipvol, pipvel      
      common /OTEC_r/ disare(mxndis), disvol(mxndis), disvel(mxndis)
      common /OTEC_r/ dissum(mxndis), dismax, disrad, pwrstr, pwrpek
      common /OTEC_r/ pwrend, pwrplt, pwrtot, yerstr, yerpek, yerend
      common /OTEC_r/ emtred, elcdem, grdmin, grdsor(imtm2*jmtm2)
      common /OTEC_r/ grdacc(imt,jmt), grdavg(imt,jmt), grdstr(imt,jmt)
      common /OTEC_r/ artplt, aremin, OTEC_yr

!     time step integral variables

      real tai_nmtplt, tai_lctrem, tai_lctpot, tai_nmtblt
      real tai_nmtdem, tai_nmtmov, tai_pwrtot, tai_elcdem
      real tai_emtred
      
      common /OTEC_r/ tai_nmtplt, tai_lctrem, tai_lctpot, tai_nmtblt
      common /OTEC_r/ tai_nmtdem, tai_nmtmov, tai_pwrtot, tai_elcdem
      common /OTEC_r/ tai_emtred