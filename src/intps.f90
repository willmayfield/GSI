module intpsmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   intpsmod    module for intps and its tangent linear intps_tl
!   prgmmr:
!
! abstract: module for intps and its tangent linear intps_tl
!
! program history log:
!   2005-05-12  Yanqiu zhu - wrap intps and its tangent linear intps_tl into one module
!   2005-11-16  Derber - remove interfaces
!   2008-11-26  Todling - remove intps_tl; add interface back
!   2009-08-13  lueken - update documentation
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - implemented obs adjoint test  
!   2013-10-28  todling - rename p3d to prse
!   2014-04-12       su - add non linear qc from Purser's scheme
!
! subroutines included:
!   sub intps_
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC intps

interface intps; module procedure &
          intps_
end interface

contains

subroutine intps_(pshead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intps       apply nonlin qc obs operator for ps
!   prgmmr: derber           org: np23                date: 1991-02-26
!
! abstract: apply observation operator and adjoint for ps observations
!           with nonlinear qc operator
!
! program history log:
!   1991-02-26  derber
!   1997-12-12  weiyu yang
!   1999-08-24  derber, j., yang, w., first frozen mpp version
!   2004-08-02  treadon - add only to module use, add intent in/out
!   2004-10-08  parrish - add nonlinear qc option
!   2005-03-01  parrish - nonlinear qc change to account for inflated obs error
!   2005-04-11  treadon - merge intps and intps_qc into single routine
!   2005-08-02  derber  - modify for variational qc parameters for each ob
!   2005-09-28  derber  - consolidate location and weight arrays
!   2005-10-21  su      - modify for variational qc
!   2006-07-28  derber  - modify to use new inner loop obs data structure
!   2007-02-15  rancic - add foto
!   2007-03-19  tremolet - binning of observations
!   2007-06-05  tremolet - use observation diagnostics structure
!   2007-07-09  tremolet - observation sensitivity
!   2008-06-02  safford - rm unused vars
!   2008-01-04  tremolet - Don't apply H^T if l_do_adjoint is false
!   2008-11-28  todling  - turn FOTO optional; changed ptr%time handle
!   2010-05-13  todling  - update to use gsi_bundlemod; update interface
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - introduced ladtest_obs         

!
!   input argument list:
!     pshead  - obs type pointer to obs structure
!     sp      - ps increment in grid space
!     rp
!
!   output argument list:
!     rp      - ps results from observation operator 
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: half,one,tiny_r_kind,cg_term,r3600,two
  use obsmod, only: ps_ob_type,lsaveobsens,l_do_adjoint
  use qcmod, only: nlnqc_iter,varqc_iter
  use gridmod, only: latlon1n1
  use jfunc, only: jiter,l_foto,xhat_dt,dhat_dt
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_4dvar, only: ladtest_obs
  implicit none

! Declare passed variables
  type(ps_ob_type),pointer,intent(in   ) :: pshead
  type(gsi_bundle),        intent(in   ) :: sval
  type(gsi_bundle),        intent(inout) :: rval

! Declare local variables
  integer(i_kind) ier,istatus
  integer(i_kind) j1,j2,j3,j4
! real(r_kind) penalty
  real(r_kind) cg_ps,val,p0,grad,wnotgross,wgross,ps_pg
  real(r_kind) w1,w2,w3,w4,time_ps
  real(r_kind),pointer,dimension(:) :: xhat_dt_prse
  real(r_kind),pointer,dimension(:) :: dhat_dt_prse
  real(r_kind),pointer,dimension(:) :: sp
  real(r_kind),pointer,dimension(:) :: rp
  type(ps_ob_type), pointer :: psptr

!  If no ps data return
  if(.not. associated(pshead))return
! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'prse',sp,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'prse',rp,istatus);ier=istatus+ier
  if(l_foto) then
     call gsi_bundlegetpointer(xhat_dt,'prse',xhat_dt_prse,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(dhat_dt,'prse',dhat_dt_prse,istatus);ier=istatus+ier
  endif
  if(ier/=0)return

  psptr => pshead
  do while (associated(psptr))
     j1=psptr%ij(1)
     j2=psptr%ij(2)
     j3=psptr%ij(3)
     j4=psptr%ij(4)
     w1=psptr%wij(1)
     w2=psptr%wij(2)
     w3=psptr%wij(3)
     w4=psptr%wij(4)
     

!    Forward model
     val=w1* sp(j1)+w2* sp(j2)+w3* sp(j3)+w4* sp(j4)
     if (l_foto) then
        time_ps=psptr%time*r3600
        val=val+&
          (w1*xhat_dt_prse(j1)+w2*xhat_dt_prse(j2)+ &
           w3*xhat_dt_prse(j3)+w4*xhat_dt_prse(j4))*time_ps
     endif

     if (lsaveobsens) then
        psptr%diags%obssen(jiter) = val*psptr%raterr2*psptr%err2
     else
        if (psptr%luse) psptr%diags%tldepart(jiter)=val
     endif
  
     if (l_do_adjoint) then
        if (lsaveobsens) then
           grad = psptr%diags%obssen(jiter)
  
        else
           if( .not. ladtest_obs)   val=val-psptr%res
!          gradient of nonlinear operator
           if (nlnqc_iter .and. psptr%pg > tiny_r_kind .and.  &
                                psptr%b  > tiny_r_kind) then
              ps_pg=psptr%pg*varqc_iter
              cg_ps=cg_term/psptr%b                           ! b is d in Enderson
              wnotgross= one-ps_pg                            ! pg is A in Enderson
              wgross =ps_pg*cg_ps/wnotgross                   ! wgross is gama in Enderson
              p0=wgross/(wgross+exp(-half*psptr%err2*val**2)) ! p0 is P in Enderson
              val=val*(one-p0)                                ! term is Wqc in Enderson
           endif
           if ( psptr%jb  > tiny_r_kind) then
              val=sqrt(two*psptr%jb)*tanh(sqrt(psptr%err2*psptr%raterr2)*val/sqrt(two*psptr%jb))
           endif
           if ( psptr%jb  > tiny_r_kind) then
              grad = val*sqrt(psptr%raterr2*psptr%err2)
              if( ladtest_obs) then
                 grad = val
              else
                 grad = val*psptr%raterr2*psptr%err2
              end if
           endif
        endif

!       Adjoint
        rp(j1)=rp(j1)+w1*grad
        rp(j2)=rp(j2)+w2*grad
        rp(j3)=rp(j3)+w3*grad
        rp(j4)=rp(j4)+w4*grad


        if (l_foto) then
           grad=grad*time_ps
           dhat_dt_prse(j1)=dhat_dt_prse(j1)+w1*grad
           dhat_dt_prse(j2)=dhat_dt_prse(j2)+w2*grad
           dhat_dt_prse(j3)=dhat_dt_prse(j3)+w3*grad
           dhat_dt_prse(j4)=dhat_dt_prse(j4)+w4*grad
        endif

     endif

     psptr => psptr%llpoint
  end do
  return
end subroutine intps_

end module intpsmod
