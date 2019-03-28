!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE pwinterfacemod
  !----------------------------------------------------------------------------
  !
  ! ... This ``interface" holds temporary data like gvector maps and data
  ! ... partitioning among mpi processes. 
  ! ... performs
  USE klist,            ONLY : nelec,nelup,neldw,nkstot,nks, xk
  USE wvfct,            ONLY : nbnd, igk, g2kin, npw, ecutwfc, npw, npwx
  USE gvect,            ONLY : ngm, g
  USE cell_base,        ONLY : tpiba2, at
  USE io_files,         ONLY : nwordwfc, iunwfc
  USE wavefunctions_module, ONLY : evc
  USE fft_base,             ONLY : dffts
  USE mp_pools, ONLY: inter_pool_comm, intra_pool_comm, nproc_pool, me_pool
  !
  IMPLICIT NONE
  SAVE

  INTEGER :: mesh(3)
  INTEGER, ALLOCATABLE :: indx(:), igtog(:), igtomin(:) 
  INTEGER :: current_kpt=-1 !

  INTEGER :: ngtot   !the total number of k+g vectors < ecut on a process.  
  INTEGER :: ngtot_g !the total number of k+g vectors < ecut on ALL processes.
  INTEGER :: npw_g   !
  INTEGER :: npwx_g 
  LOGICAL :: gvec_set=.false. !this flag determines if gvector maps have been
                              !initialized.

  INTEGER, ALLOCATABLE :: ng_proc(:), ng_disp(:)
  !ng_proc:  1:nproc array with ng_proc(i)=ngtot, where i is the rank of a
  !process.
  
  !ng_disp:  1:nproc array.  Cumulative sum of ng_proc.  ng_disp(1)=0.
  !          this is for MPI_Gatherv, when number of gvectors owned by process
  !          are not the same.

  PUBLIC  :: getwfn_info, get_gvecs, get_band
  PRIVATE :: create_g_maps, get_bands_at_kpt 
  
 CONTAINS

  SUBROUTINE getwfn_info(nbands_, nktot_, nkloc_, mesh_, ngtot_, npw_, npwx_)
      USE klist,            ONLY : nelec,nelup,neldw,nkstot,nks, xk
      USE wvfct,            ONLY : nbnd, igk, g2kin, npw, ecutwfc, npw, npwx
      USE gvect,            ONLY : ngm, g
      USE cell_base,        ONLY : tpiba2
      USE io_files,         ONLY : nwordwfc, iunwfc
      USE fft_base,         ONLY : dffts
      USE mp_world,         ONLY : world_comm
      USE mp,               ONLY : mp_sum, mp_rank, mp_barrier
      USE io_global, ONLY: stdout
     INTEGER, INTENT(OUT) :: nbands_, nktot_, nkloc_, ngtot_, npw_, npwx_
     INTEGER, INTENT(OUT) :: mesh_(3)

     nbands_=nbnd
     nktot_=nkstot
     nkloc_=nks
!     WRITE(stdout,*) 'rank = ',mp_rank(world_comm)
!     WRITE(stdout,*) 'nbands = ',nbands_
!     WRITE(stdout,*) 'nktot = ',nktot_
!     WRITE(stdout,*) 'nkloc = ',nkloc_
!     WRITE(stdout,*) 'gvec_set = ',gvec_set
!     WRITE(stdout,*) 'mpi_comm_world = ',world_comm
!     WRITE(stdout,*) 'inter_pool_comm = ',inter_pool_comm
!     WRITE(stdout,*) 'intra_pool_comm = ',intra_pool_comm 
     call mp_barrier(world_comm)
     call flush_unit(stdout)
     IF( gvec_set .eqv. .false. ) CALL create_g_maps()
     mesh_(1) = mesh(1)
     mesh_(2) = mesh(2)
     mesh_(3) = mesh(3)

     ngtot_g=ngtot
     npw_g=npw
     npwx_g=npwx
    
     !sum the local ngtot, npw's, and npwx's over all processes.  results stored
     !in X_g.
     CALL mp_sum(ngtot_g,intra_pool_comm) 
     CALL mp_sum(npw_g, intra_pool_comm)
     CALL mp_sum(npwx_g, intra_pool_comm)

     ngtot_=ngtot_g  
     npw_=npw_g
     npwx_=npwx_g

  END SUBROUTINE getwfn_info

 
 SUBROUTINE create_g_maps()
   USE mp_world,         ONLY : world_comm
   USE io_global, ONLY: stdout, ionode, ionode_id
   USE mp,        ONLY: mp_gather, mp_barrier, mp_bcast
   INTEGER :: ik, ig
   ALLOCATE (indx (4*npwx) )
   ALLOCATE (igtog (4*npwx) )
   ALLOCATE (igtomin(4*npwx) )
    
   ALLOCATE(ng_proc(nproc_pool))
   ALLOCATE(ng_disp(nproc_pool))
   
   indx=0
   igtog=0
   igtomin=0

   ng_proc=0
   ng_disp=0

   mesh(1) = dffts%nr1
   mesh(2) = dffts%nr2
   mesh(3) = dffts%nr3 

   
   !Taken from pw2qmcpack.  This creates a map from the gvector list to energy sorted k+g vectors in
   !the array igk(:).  
   !indx notes which gvectors are used among all kvectors.
   !
   !NOTE:  MAPS & gvecs are local to this mpiprocess.  gvectors are distributed before we
   !get here.  
   DO ik = 1, nkstot
      CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      DO ig =1, npw
         IF( igk(ig) > 4*npwx ) & 
              CALL errore ('pwinterface','increase allocation of index', ig)
         indx( igk(ig) ) = 1 
      ENDDO
   ENDDO
    
   ngtot = 0
   
   DO ig = 1, 4*npwx
     IF( indx(ig) == 1 ) THEN
       ngtot = ngtot + 1          
       igtog(ngtot) = ig   !in k-collected gvector list, find where ig lies in it.  
       igtomin(ig) = ngtot
     ENDIF
   ENDDO 
   
   gvec_set=.true.
!   WRITE(*,*) "calling mp_gather in create_g_maps"
!   WRITE(*,*) "ngtot=",ngtot," ng_proc=",ng_proc," ionode_id=",ionode_id," intra_pool_comm=",intra_pool_comm 
   CALL mp_gather(ngtot, ng_proc, ionode_id, intra_pool_comm) !collect all local
   CALL mp_bcast(ng_proc, ionode_id, intra_pool_comm)         !broadcast ng_proc to all processes
   !Set up the displacement vector for gatherv.  
   ng_disp(0)=0         
   DO ik=2,nproc_pool
     ng_disp(ik)=ng_disp(ik-1)+ng_proc(ik-1)
   ENDDO

 END SUBROUTINE create_g_maps 

 SUBROUTINE get_gvecs(gvecs_)
   USE io_global, ONLY: stdout, ionode, ionode_id
   USE mp,        ONLY: mp_gather, mp_bcast, mp_barrier
   INTEGER :: disp(nproc_pool), proc_dat(nproc_pool)
   INTEGER :: ig, ig_c
   INTEGER :: gvecs_(3,ngtot_g)
   INTEGER :: gvecs_local(3,ngtot)

   !We have local versions of ng_disp and ng_proc, except modified
   ! to account for 3xn vectors being sent.  Multiply ng_proc and disp by 3
   ! to handle this.     
   disp(1:nproc_pool) = 3*ng_disp(1:nproc_pool)
   proc_dat(1:nproc_pool)=3*ng_proc(1:nproc_pool)
   IF (gvec_set .eqv. .false.) CALL create_g_maps()
 
   DO ig=1, ngtot
     ig_c =igtog(ig)
     gvecs_local(1,ig)=NINT(at(1,1)*g(1,ig_c)+at(2,1)*g(2,ig_c)+at(3,1)*g(3,ig_c))
     gvecs_local(2,ig)=NINT(at(1,2)*g(1,ig_c)+at(2,2)*g(2,ig_c)+at(3,2)*g(3,ig_c))
     gvecs_local(3,ig)=NINT(at(1,3)*g(1,ig_c)+at(2,3)*g(2,ig_c)+at(3,3)*g(3,ig_c))
    ENDDO

  CALL mp_gather(gvecs_local(1:3,:), gvecs_(1:3,:), ng_proc, ng_disp, ionode_id, intra_pool_comm)
  CALL mp_bcast(gvecs_, ionode_id, intra_pool_comm)  !again, collect gvecs to ionode, and then broadcast to all 

 END SUBROUTINE get_gvecs

 SUBROUTINE get_bands_at_kpt(ik)
   USE wavefunctions_module, ONLY : evc
   USE cell_base, ONLY: omega, alat, tpiba2, at, bg
   USE gvect, ONLY: ngm, g
   USE gvecs, ONLY : nls, nlsm  
   USE klist , ONLY: nks, nelec, nelup, neldw, wk, xk, nkstot
   USE wvfct, ONLY: npw, npwx, nbnd, igk, g2kin, wg, et, ecutwfc
   USE buffers,              ONLY : get_buffer
   USE io_global, ONLY: stdout, ionode, ionode_id

   INTEGER :: ik

   CALL gk_sort (xk (1:3, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
   CALL get_buffer (evc, nwordwfc, iunwfc, ik)
   
 END SUBROUTINE get_bands_at_kpt

 SUBROUTINE get_band(cg_, ibnd_, ik_)
   USE wavefunctions_module, ONLY : evc
   USE io_global, ONLY: stdout, ionode, ionode_id
   USE mp,        ONLY: mp_gather, mp_bcast

   COMPLEX*16, INTENT(OUT) :: cg_(ngtot_g) ! plane wave coefficients that are returned.
   COMPLEX*16 :: cg_local(ngtot)           ! the local pw coefficients corresponding to local gvecs. 
   INTEGER, INTENT(IN) :: ik_, ibnd_       ! kpoint id and band number respectively.  
!   WRITE(*,*) "pwinterfacemod::get_band(cg, ibnd, ik)"
!   WRITE(*,*) "   ibnd=",ibnd_, " ik=", ik_
   
   IF (gvec_set .eqv. .false.) CALL create_g_maps()
   
   IF (ik_ /= current_kpt) THEN            !only read and grab PW data if it hasn't be done for this kpoint already.  
     CALL get_bands_at_kpt(ik_)            !  e.g. you have to read all the bands at a kpoint
     current_kpt = ik_
   ENDIF
 
   cg_local(:)=(0.d0,0.d0)
   cg_local(igtomin(igk(1:npw))) = evc(1:npw,ibnd_)
    
   CALL mp_gather(cg_local, cg_,ng_proc, ng_disp, ionode_id, intra_pool_comm) 
   CALL mp_bcast(cg_, ionode_id, intra_pool_comm)

 END SUBROUTINE get_band

 SUBROUTINE get_atom_pos(R)
     USE ions_base,       ONLY : nat, tau
     USE cell_base,       ONLY : alat
     USE kinds,             ONLY:  DP     
      !
     REAL(DP), INTENT(OUT) :: R(3,nat)
      !
      R     = tau*alat

 END SUBROUTINE get_atom_pos


 SUBROUTINE set_atom_pos(R)
     USE ions_base,       ONLY : nat, tau
     USE cell_base,       ONLY : alat
     USE kinds,             ONLY:  DP     
      !
     REAL(DP), INTENT(IN) :: R(3,nat)
      !
      tau     = R/alat

 END SUBROUTINE set_atom_pos


 SUBROUTINE get_atom_types(ityp_)
     USE ions_base,       ONLY : nat, ityp
      !
     INTEGER, INTENT(OUT) :: ityp_(nat)
     !
     ityp_=ityp
 END SUBROUTINE get_atom_types

 SUBROUTINE set_atom_types(ityp_)
     USE ions_base,       ONLY : nat, ityp
      !
     INTEGER, INTENT(IN) :: ityp_(nat)
     !
     ityp=ityp_
   
 END SUBROUTINE set_atom_types

END MODULE pwinterfacemod
