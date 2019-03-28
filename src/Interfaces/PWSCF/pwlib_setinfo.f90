    !
    !------------------------------------------------------------------------
    !SUBROUTINE pwlib_getinfo(nat_,charge_,nup_,ndown_,nktot_,nkloc_,box_,  &
     ! kpts_,wkp_,ngrid_,nbds_,ngridx_)
    SUBROUTINE pwlib_getinfo(nat_,charge_,nup_,ndown_,nkloc_,box_,nktot_,kpts_)
      !------------------------------------------------------------------------
      !
      USE ions_base,        ONLY : nat, nsp
      USE klist,            ONLY : xk,wk,nelec,nelup,neldw,nkstot,nks,tot_charge 
      USE wvfct,            ONLY : nbnd  
      USE fft_base,         ONLY : dffts
      USE cell_base,        ONLY : at,alat,tpiba
     USE kinds,             ONLY:  DP     
      !
      IMPLICIT NONE
      !
      !
      INTEGER, INTENT(OUT) :: nat_,charge_,nup_,ndown_,nktot_,nkloc_
      INTEGER :: jk
      REAL(DP), INTENT(OUT) :: box_(3,3), kpts_(3,nks) 
      !
      nat_ = nat 
      charge_ = tot_charge 
      nup_ = NINT(nelup)      
      ndown_ = NINT(neldw)
      nktot_ = nkstot
      nkloc_ = nks
      box_ = at*alat

      kpts_(1:3,1:nks) = xk(1:3,1:nks)*tpiba
      
      ! 
    END SUBROUTINE pwlib_getinfo
    !
   SUBROUTINE pwlib_setbox_data(box_)
     USE cell_base,       ONLY : at, alat
     USE kinds,             ONLY:  DP     
 
     IMPLICIT NONE
     REAL(DP), INTENT(OUT) :: box_(3,3)
     
     box_=at*alat 

   END SUBROUTINE pwlib_getbox_data

   

   SUBROUTINE pwlib_setatom_info(nat_, nsp_)
      !------------------------------------------------------------------------
      !
      USE ions_base,       ONLY : nat, nsp
      !
      INTEGER, INTENT(OUT) :: nat_, nsp_
      !
      nat_ = nat
      nsp_ = nsp

   END SUBROUTINE pwlib_getatom_info

!! Get/Set Routines for ion data.  
   !Returns R in bohr.  

   SUBROUTINE pwlib_setatom_data(R, ityp_)
      !------------------------------------------------------------------------
      !
      USE ions_base,       ONLY : nat, nsp, ityp, tau
      USE cell_base,       ONLY : alat
     USE kinds,             ONLY:  DP     
      !
      REAL(DP), INTENT(OUT) :: R(3,nat)
      INTEGER, INTENT(OUT) :: ityp_(nat)
      !
      R     = tau*alat
      ityp_ = ityp

   END SUBROUTINE pwlib_setatom_data

   SUBROUTINE pwlib_setatom_pos(R)
     USE pwinterfacemod, ONLY : get_atom_pos
     USE ions_base,  ONLY : nat

     REAL(DP), INTENT(OUT) :: R(3,nat)

     CALL set_atom_pos(R)
  
    END SUBROUTINE pwlib_setatom_pos


   SUBROUTINE pwlib_setatom_types(ityp_)
     USE pwinterfacemod, ONLY : get_atom_types
     USE ions_base,  ONLY : nat

     INTEGER, INTENT(IN) :: ityp_(nat)

     CALL set_atom_types(ityp_)

   END SUBROUTINE pwlib_setatom_types
  
!   SUBROUTINE pwlib_setatom_data(R,ityp, nat_, nsp_)
!      !------------------------------------------------------------------------
!      USE ions_base,      ONLY nat, 
!   END SUBROUTINE pwlib_setatom_data
    ! anums=atomic numbers, vcharge=valence charges, masses=ion masses in amu,
    ! names are element names.  
   SUBROUTINE pwlib_setspecies_data(anums, vcharge, masses, names)
      !------------------------------------------------------------------------
      !
      USE ions_base,       ONLY : nsp, atm, zv, amass
     USE kinds,             ONLY:  DP     
      !
      INTEGER, EXTERNAL :: atomic_number
      INTEGER, INTENT(OUT) :: anums(nsp), vcharge(nsp)
      REAL(DP), INTENT(OUT) :: masses(nsp)
      CHARACTER(LEN=3), INTENT(OUT) :: names(nsp)
      INTEGER :: ik
      !
      do ik=1,nsp
        anums(ik)=atomic_number(TRIM(atm(ik)))
      end do 

      vcharge = zv
      masses=amass
      names=atm
      
   END SUBROUTINE pwlib_setspecies_data

  SUBROUTINE pwlib_setelectron_info(nelec_, nup_, ndown_, nktot_, nkloc_)
      USE klist,            ONLY : nelec,nelup,neldw 
   !   USE fft_base,         ONLY : dffts
        
      INTEGER, INTENT(OUT) :: nelec_, nup_, ndown_
      
      nelec_  = nelec 
      nup_    = nelup
      ndown_  = neldw
      
  END SUBROUTINE pwlib_setelectron_info
 
  SUBROUTINE pwlib_setwfn_info(nbands_, nktot_, nkloc_, mesh_, ngtot_, npw_, npwx_)
     USE pwinterfacemod,    ONLY : getwfn_info 
     USE mp,                ONLY : mp_rank 
           USE mp_world,         ONLY : world_comm
      USE io_global, ONLY: stdout

     INTEGER, INTENT(OUT) :: nbands_, nktot_, nkloc_, ngtot_, npw_, npwx_
     INTEGER, INTENT(OUT) :: mesh_(3)
     

     WRITE(*,*) 'b4 getwfn_info...'
          WRITE(stdout,*) 'rank = ',mp_rank(world_comm)

     CALL setwfn_info(nbands_, nktot_, nkloc_, mesh_, ngtot_, npw_, npwx_)
    
     
  END SUBROUTINE pwlib_getwfn_info

  SUBROUTINE pwlib_setwfn_gvecs(gvecs_)
   USE pwinterfacemod,    ONLY: ngtot, igtog, get_gvecs

   INTEGER :: gvecs_(3,ngtot)
 
   CALL get_gvecs(gvecs_)
 
  END SUBROUTINE pwlib_setwfn_gvecs

  SUBROUTINE pwlib_setwfn_band(cg_, iband_, ik_)
     USE pwinterfacemod,    ONLY:  ngtot_g, igtomin, igtog, get_band
     USE kinds,             ONLY:  DP     
     INTEGER, INTENT(IN) :: iband_, ik_
     COMPLEX(KIND=DP), INTENT(OUT) :: cg_(ngtot_g)

     CALL get_band(cg_, iband_, ik_)
  
  END SUBROUTINE pwlib_setwfn_band

  SUBROUTINE pwlib_setwfn_kpoints(klist_, weights_)
      USE klist,           ONLY  : nkstot, xk,wk
      USE cell_base,       ONLY  : tpiba, at
     USE kinds,             ONLY:  DP     

      REAL(DP), INTENT(OUT) :: klist_(3,nkstot), weights_(nkstot)
      klist_(1,:)=at(1,1)*xk(1,:)+at(2,1)*xk(2,:)+at(3,1)*xk(3,:)
      klist_(2,:)=at(1,2)*xk(1,:)+at(2,2)*xk(2,:)+at(3,2)*xk(3,:)
      klist_(3,:)=at(1,3)*xk(1,:)+at(2,3)*xk(2,:)+at(3,3)*xk(3,:)
      weights_=wk

  END SUBROUTINE pwlib_setwfn_kpoints
  
  SUBROUTINE pwlib_setwfn_eigenvals(eigval_, ik)
      USE wvfct,          ONLY  : nbnd, et
     USE kinds,             ONLY:  DP     

      REAL(DP), INTENT(OUT) :: eigval_(nbnd)
      INTEGER, INTENT(IN) :: ik
      eigval_(1:nbnd)=0.5*et(1:nbnd, ik)
      
  END SUBROUTINE pwlib_setwfn_eigenvals

  SUBROUTINE pwlib_setmpifft(nnr_,mesh_)
    USE fft_base,         ONLY : dffts, dfftp
    INTEGER, INTENT(OUT) :: nnr_, mesh_(3)
   
    nnr_ = dffts%nnr

    mesh_(1)=dffts%nr1x
    mesh_(2)=dffts%nr2x
    mesh_(3)=dffts%nr3x
    

  END SUBROUTINE pwlib_setmpifft

  
