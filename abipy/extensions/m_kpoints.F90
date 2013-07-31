module m_kpoints

 implicit none

 public

 real(8),private,parameter :: TOL_KDIFF=0.0001

contains

!----------------------------------------------------------------------

!!****f* m_bz_mesh/map_bz2ibz
!! NAME
!! map_bz2ibz
!!
!! FUNCTION
!!
!! INPUTS
!!  nkibz=Number of k points in IBZ.
!!  kibz(3,nkibz)=Coordinates of k-points in the IBZ.
!!  nkbz=Number of k points in the BZ.
!!  kbz(3,nkbz)= k-points in the whole BZ
!!  nsym=Number of symmetry operations.
!!  timrev=2 if time reversal symmetry can be used; 1 otherwise.
!!  symrec(3,3,nsym)=Symmetry operation matrices in reciprocal space.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations.
!!
!! OUTPUT
!!  ktab(nkbz)= table giving for each k-point in the BZ (array kbz), the corresponding irreducible point in the array (kibz)
!!    k_BZ= (IS) kIBZ where S is one of the symrec operations and I is the inversion or the identity
!!    where k_BZ = (IS) k_IBZ and S = \transpose R^{-1} 
!!  ktabi(nkbz)= for each k-point in the BZ defines whether inversion has to be 
!!   considered in the relation k_BZ=(IS) k_IBZ (1 => only S; -1 => -S)  
!!  ktabo(nkbz)= the symmetry operation S that takes k_IBZ to each k_BZ
!!
!! SOURCE

subroutine map_bz2ibz(nkibz,kibz,nkbz,kbz,nsym,timrev,symrec,symafm,bz2ibz,ktabi,ktabo,ierr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkibz,nkbz,nsym,timrev
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym)
 integer,intent(out) :: bz2ibz(nkbz),ktabi(nkbz),ktabo(nkbz)
 real(8),intent(in) :: kibz(3,nkibz),kbz(3,nkbz)

!Local variables ------------------------------
!scalars
 integer :: ik_bz,ik_ibz,isym,itim,ans
!arrays
 integer :: g0(3)
 real(8) :: knew(3) 

! *************************************************************************

 ierr=0

 bz2ibz=-1; ktabi=-1; ktabo=-1
 
 do ik_bz=1,nkbz
   ibz_loop: do ik_ibz=1,nkibz
     !
     ! === Loop over time-reversal I and symmetry operations S  ===
     ! * Use spatial inversion instead of time reversal whenever possible.
     do itim=1,timrev
       do isym=1,nsym
         if (symafm(isym)==-1) CYCLE
         !
         ! * Form IS k
         knew = (3-2*itim) * MATMUL(symrec(:,:,isym),kibz(:,ik_ibz))
         !
         ! * Check whether it has already been found (to within a RL vector).
         call isamek(knew,kbz(:,ik_bz),ans,g0) 
         if (ans==1) then 
           bz2ibz (ik_bz) = ik_ibz - 1
           ktabo(ik_bz) = isym - 1
           ktabi(ik_bz) = 3-2*itim
           EXIT ibz_loop
         end if
         !
        end do 
     end do 
     !
   end do ibz_loop
 end do

 if (ANY(bz2ibz==-1)) then
   ierr = 1
 end if

end subroutine map_bz2ibz
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/isamek
!! NAME
!! isamek
!!
!! FUNCTION
!! Test two k-points for equality. 
!! Return .TRUE. is they are equal within a reciprocal lattice vector G0.
!!
!! INPUTS
!!  k1(3),k2(3)=The two k points to be compared.
!!
!! OUTPUT
!! Return .TRUE. if they are the same within a RL vector,
!!        .FALSE. if they are different.
!! G0(3)=if .TRUE. G0(3) is the reciprocal lattice vector such as k1=k2+G0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine isamek(k1,k2,ans,g0)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ans
!arrays
 integer,intent(out) :: g0(3)
 real(8),intent(in) :: k1(3),k2(3)

!Local variables-------------------------------
!scalars
 real(8) :: kdiff(3)

! *************************************************************************

 kdiff = k1-k2

 if (is_integer(kdiff,TOL_KDIFF)) then
   ans = 1
   g0=NINT(kdiff)
 else 
   ans = 0
   g0=HUGE(1)
 end if

end subroutine isamek
!!***

!!****f* m_numeric_tools/is_integer
!! NAME
!!  is_integer
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE
function is_integer(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars
 real(8),intent(in) :: tol
 logical :: ans
!arrays
 real(8),intent(in) :: rr(:)

!Local variables-------------------------------
!scalars
! *************************************************************************

 ans=ALL((ABS(rr-NINT(rr))<tol))

end function is_integer
!!***

end module m_kpoints

!----------------------------------------------------------------------
