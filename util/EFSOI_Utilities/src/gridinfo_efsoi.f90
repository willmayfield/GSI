module gridinfo_efsoi
!$$$  module documentation block
!
! module: gridinfo_efsoi               read horizontal (lons, lats) and
!                                      vertical (pressure) information from
!                                      ensemble mean first guess file.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
! prgmmr: groff            org: emc                    date: 2018-05-24
!
! abstract: This module reads gfg_YYYYMMDDHH_fhr06_ensmean, and
! extracts information about the analysis grid, including the
! longitudes and latitudes of the analysis grid points and
! the pressure on each grid point/vertical level.
!
! Public Subroutines:
!   getgridinfo_efsoi: read latitudes, longitudes, pressures and orography for analysis grid,
!    broadcast to each task. Compute spherical cartesian coordinate values
!    for each analysis horizontal grid point.
!   gridinfo_cleanup_efsoi: deallocate allocated module variables.
!
! Public Variables:
!   npts: number of analysis grid points in the horizontal (from module params).
!   nlevs: number of analysis vertical levels (from module params).
!   ntrac: number of 'tracer' model state variables (3 for GFS,
!    specific humidity, ozone and cloud condensate).
!   ptop: (real scalar) pressure (hPa) at top model layer interface.
!   lonsgrd(npts): real array of analysis grid longitudes (radians).
!   latsgrd(npts): real array of analysis grid latitudes (radians).
!   logp(npts,ndim):  -log(press) for all 2d analysis grids. Assumed invariant
!   in assimilation window, computed fro ensemble mean at middle of window.
!   gridloc(3,npts): spherical cartesian coordinates (x,y,z) for analysis grid.
!
! Modules Used: mpisetup, params, kinds
!
! program history log:
!   2009-02-23  Initial version.
!   2016-05-02: shlyaeva: Modification for reading state vector from table
!   2016-04-20  Modify to handle the updated nemsio sig file (P, DP & DPDT removed)
!   2018-05-24  Groff: Adapting EnKF gridinfo_efsoi module to EFSOI needs.
!
! attributes:
!   language: f95
!
!$$$

use mpisetup, only: nproc, mpi_integer, mpi_real4, mpi_comm_world
use params, only: datapath,nlevs,nlons,nlats,fgfileprefixes,nvars
use kinds, only: r_kind, i_kind, r_double, r_single
use constants, only: one,zero,pi,cp,rd,grav,rearth,max_varname_length
use specmod, only: sptezv_s, sptez_s, init_spec_vars, isinitialized, asin_gaulats, &
    ndimspec => nc
use reducedgrid_mod, only: reducedgrid_init, regtoreduced, reducedtoreg,&
                           nptsred, lonsred, latsred
use mpeu_util, only: getindex
!use statevec_efsoi, only: clevels !, cvars2d, cvars3d, nc3d
implicit none
private
public :: getgridinfo_efsoi, gridinfo_cleanup_efsoi
integer(i_kind),public :: nlevs_pres, idvc
real(r_single),public :: ptop
real(r_single),public, allocatable, dimension(:) :: lonsgrd, latsgrd
! arrays passed to kdtree2 routines must be single
real(r_single),public, allocatable, dimension(:,:) :: gridloc
real(r_single),public, allocatable, dimension(:,:) :: logp
integer,public :: npts
integer,public :: ntrunc
integer(i_kind),public, allocatable, dimension(:) :: id_u, id_v, id_t, id_q
integer(i_kind),public :: id_ps, ncdim
! supported variable names in anavinfo
character(len=max_varname_length),public, dimension(10) :: vars3d_supported = (/'u   ', 'v   ', 'tv  ', 'q   ', 'oz  ', 'cw  ', 'tsen', 'prse', 'ql  ', 'qi  '/)
character(len=max_varname_length),public, dimension(3)  :: vars2d_supported = (/'ps ', 'pst', 'sst' /)
! supported variable names in anavinfo
contains

subroutine getgridinfo_efsoi(fileprefix, reducedgrid, clevels, nc3d)
! read latitudes, longitudes and pressures for analysis grid,
! broadcast to each task.
use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,&
                         nemsio_getfilehead,nemsio_getheadvar,&
                         nemsio_readrecv,nemsio_init, nemsio_realkind
implicit none

character(len=120), intent(in) :: fileprefix
logical, intent(in)            :: reducedgrid
integer, intent(in) :: nc3d
integer, dimension(0:nc3d), intent(in) :: clevels
integer(i_kind) nlevsin, ierr, iunit, k, nn, idvc 
character(len=500) filename
integer(i_kind) iret,i,j,nlonsin,nlatsin,u_ind,v_ind,tv_ind,q_ind,ps_ind
real(r_kind), allocatable, dimension(:) :: ak,bk,spressmn,tmpspec
real(r_kind), allocatable, dimension(:,:) :: pressimn,presslmn
real(r_single),allocatable,dimension(:,:,:) :: nems_vcoord
real(r_kind) kap,kapr,kap1
real(nemsio_realkind), dimension(nlons*nlats) :: nems_wrk
type(nemsio_gfile) :: gfile

iunit = 77
kap = rd/cp
kapr = cp/rd
kap1 = kap + one
nlevs_pres=nlevs+1
if (nproc .eq. 0) then
   filename = trim(adjustl(datapath))//trim(adjustl(fileprefix))//"ensmean"
   call nemsio_init(iret=iret)
   if(iret/=0) then
      write(6,*)'grdinfo: gfs model: problem with nemsio_init, iret=',iret, ', file: ', trim(filename)
      call stop2(23)
   end if
   call nemsio_open(gfile,filename,'READ',iret=iret)
   if (iret/=0) then
      write(6,*)'grdinfo: gfs model: problem with nemsio_open, iret=',iret, ', file: ', trim(filename)
      call stop2(23)
   endif
   call nemsio_getfilehead(gfile,iret=iret, dimx=nlonsin, dimy=nlatsin,&
                           dimz=nlevsin,jcap=ntrunc,idvc=idvc)
   ! set ntrunc to nlats if missing
   ! (only used for inflation smoothing and mass balance adjustment if use_gfsnemsio = T)
   ! FV3GFS write component does not include JCAP, infer from nlatsin
   if (ntrunc < 0) ntrunc = nlatsin-2
   if (iret/=0) then
      write(6,*)'grdinfo: gfs model: problem with nemsio_getfilehead, iret=',iret, ', file: ', trim(filename)
      call stop2(23)
   endif
   print *,'ntrunc = ',ntrunc
   if (nlons /= nlonsin .or. nlats /= nlatsin .or. nlevs /= nlevsin) then
      print *,'incorrect dims in nemsio file'
      print *,'expected',nlons,nlats,nlevs
      print *,'got',nlonsin,nlatsin,nlevsin
      call stop2(23)
   end if
endif
call mpi_bcast(ntrunc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! initialize spectral module on all tasks.
if (.not. isinitialized) call init_spec_vars(nlons,nlats,ntrunc,4)

if (nproc .eq. 0) then
   ! get pressure, lat/lon information from ensemble mean file.
   allocate(presslmn(nlons*nlats,nlevs))
   allocate(pressimn(nlons*nlats,nlevs+1))
   allocate(spressmn(nlons*nlats))
   call nemsio_readrecv(gfile,'pres','sfc',1,nems_wrk,iret=iret)
   if (iret/=0) then
      write(6,*)'grdinfo: gfs model: problem with nemsio_readrecv(ps), iret=',iret
      call stop2(23)
   endif

!     Extract vertical coordinate descriptions nems_vcoord.
!     nems_vcoord(gfshead%levs+1,3,2) dimension is hardwired here.
!     Present NEMSIO modules do not allow flexibility of 2nd and 3rd
!     array dimension for nems_vcoord, for now, it is hardwired as
!     (levs,3,2) If NEMS changes the setting of vcoord dimension,
!     GSI needs to update its setting of nems_vcoord accordingly.

   if (allocated(nems_vcoord))     deallocate(nems_vcoord)
   allocate(nems_vcoord(nlevs_pres,3,2))
   call nemsio_getfilehead(gfile,iret=iret,vcoord=nems_vcoord)
   if ( iret /= 0 ) then
      write(6,*)' gridinfo_efsoi:  ***ERROR*** problem reading header ', &
         'vcoord, Status = ',iret
      call stop2(99)
   endif

   spressmn = 0.01_r_kind*nems_wrk ! convert ps to millibars.
   !print *,'min/max spressmn = ',minval(spressmn),maxval(spressmn)

   allocate(ak(nlevs+1),bk(nlevs+1))

   if ( idvc == 0 ) then                         ! sigma coordinate, old file format.
      ak = zero
      bk = nems_vcoord(1:nlevs+1,1,1)
   elseif ( idvc == 1 ) then                     ! sigma coordinate
      ak = zero
      bk = nems_vcoord(1:nlevs+1,2,1)
   elseif ( idvc == 2 .or. idvc == 3 ) then      ! hybrid coordinate
      ak = 0.01_r_kind*nems_vcoord(1:nlevs+1,1,1) ! convert to mb
      bk = nems_vcoord(1:nlevs+1,2,1)
   else
      write(6,*)'gridinfo_efsoi:  ***ERROR*** INVALID value for idvc=',idvc
      call stop2(85)
   endif

   ! pressure at interfaces
   do k=1,nlevs+1
      pressimn(:,k) = ak(k)+bk(k)*spressmn(:)
   enddo
   call nemsio_close(gfile, iret=iret)
   ptop = ak(nlevs+1)
   deallocate(ak,bk)

   if (reducedgrid) then
      call reducedgrid_init(nlons,nlats,asin_gaulats)
      npts = nptsred 
   else
      npts = nlons*nlats
   end if
   allocate(latsgrd(npts),lonsgrd(npts))
   allocate(logp(npts,nlevs_pres)) ! log(ens mean first guess press) on mid-layers
   allocate(gridloc(3,npts))
   !==> pressure at interfaces.
   if (reducedgrid) then
      lonsgrd(:) = lonsred(:)
      latsgrd(:) = latsred(:)
   else
      nn = 0
      do j=1,nlats
         do i=1,nlons
            nn = nn + 1
            lonsgrd(nn) = 2._r_single*pi*float(i-1)/nlons
            latsgrd(nn) = asin_gaulats(j)
         enddo
      enddo
   endif
   do k=1,nlevs
     ! layer pressure from Phillips vertical interpolation.
     presslmn(:,k) = ((pressimn(:,k)**kap1-pressimn(:,k+1)**kap1)/&
                      (kap1*(pressimn(:,k)-pressimn(:,k+1))))**kapr
   end do
   print *,'ensemble mean first guess surface pressure:'
   print *,minval(spressmn),maxval(spressmn)
   !do k=1,nlevs
   !   print *,'min/max ens mean press level',&
   !   k,'=',minval(presslmn(:,k)),maxval(presslmn(:,k))
   !   print *,'min/max ens mean press interface',&
   !   k,'=',minval(pressimn(:,k)),maxval(pressimn(:,k))
   !enddo
   ! logp holds log(pressure) or pseudo-height on grid, for each level/variable.
   do k=1,nlevs
      ! all variables to be updated are on mid-layers, not layer interfaces.
      if (reducedgrid) then
         call regtoreduced(presslmn(:,k),logp(:,k))
         logp(:,k) = -log(logp(:,k))
      else
         logp(:,k) = -log(presslmn(:,k))
      endif
      !print *,'min/max presslmn',k,minval(presslmn(:,k)),maxval(presslmn(:,k)),minval(logp(:,k)),maxval(logp(:,k))
   end do
   if (reducedgrid) then
      call regtoreduced(spressmn,logp(:,nlevs_pres))
      logp(:,nlevs_pres) = -log(logp(:,nlevs_pres))
   else
      logp(:,nlevs_pres) = -log(spressmn(:))
   endif
   deallocate(spressmn,presslmn,pressimn)
end if
call mpi_bcast(npts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (nproc .ne. 0) then
   ! allocate arrays on other (non-root) tasks
   allocate(latsgrd(npts),lonsgrd(npts))
   allocate(logp(npts,nlevs_pres)) ! log(ens mean first guess press) on mid-layers
   allocate(gridloc(3,npts))
   ! initialize reducedgrid_mod on other tasks.
   if (reducedgrid) then
      call reducedgrid_init(nlons,nlats,asin_gaulats)
   end if
endif
!call mpi_bcast(logp,npts*nlevs_pres,mpi_real4,0,MPI_COMM_WORLD,ierr)
do k=1,nlevs_pres
  call mpi_bcast(logp(1,k),npts,mpi_real4,0,MPI_COMM_WORLD,ierr)
enddo
call mpi_bcast(lonsgrd,npts,mpi_real4,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(latsgrd,npts,mpi_real4,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(ptop,1,mpi_real4,0,MPI_COMM_WORLD,ierr)
  
!==> precompute cartesian coords of analysis grid points.
do nn=1,npts
   gridloc(1,nn) = cos(latsgrd(nn))*cos(lonsgrd(nn))
   gridloc(2,nn) = cos(latsgrd(nn))*sin(lonsgrd(nn))
   gridloc(3,nn) = sin(latsgrd(nn))
end do

! Identify EFSOI relevant state variable indices
u_ind   = getindex(vars3d_supported, 'u')   !< indices in the state var arrays
v_ind   = getindex(vars3d_supported, 'v')   ! U and V (3D)
tv_ind  = getindex(vars3d_supported, 'tv')  ! Tv (3D)
q_ind   = getindex(vars3d_supported, 'q')   ! Q (3D)
ps_ind  = getindex(vars2d_supported, 'ps')  ! Ps (2D)

! Index of each elements
allocate(id_u(nlevs))
allocate(id_v(nlevs))
allocate(id_t(nlevs))
allocate(id_q(nlevs))
do k=1,nlevs
   id_u(k) = clevels(u_ind-1) + k
   id_v(k) = clevels(v_ind-1) + k
   id_t(k) = clevels(tv_ind-1) + k
   id_q(k) = clevels(q_ind-1) + k
end do
id_ps = clevels(nc3d) + ps_ind

end subroutine getgridinfo_efsoi

subroutine gridinfo_cleanup_efsoi()
if (allocated(lonsgrd)) deallocate(lonsgrd)
if (allocated(latsgrd)) deallocate(latsgrd)
if (allocated(logp)) deallocate(logp)
if (allocated(gridloc)) deallocate(gridloc)
if (allocated(id_u)) deallocate(id_u) 
if (allocated(id_v)) deallocate(id_v) 
if (allocated(id_t)) deallocate(id_t) 
if (allocated(id_q)) deallocate(id_q) 
end subroutine gridinfo_cleanup_efsoi

end module gridinfo_efsoi
