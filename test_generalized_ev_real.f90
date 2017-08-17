! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This program tests ELSI eigensolver.
!!
!!
! alvarovm.github.io : added Matrix Market format reader 
!

program test_generalized_ev_real

   use ELSI_PRECISION, only: r8,i4
   use ELSI
   use elpa_utilities, only : error_unit

   implicit none

   include "mpif.h"

   character(5) :: m_storage
   character(3) :: m_operation
   character(128) :: arg1
   character(128) :: arg2

   integer(kind=i4) :: n_proc,nprow,npcol,myid
   integer(kind=i4) :: myprow,mypcol
   integer(kind=i4) :: mpi_comm_global,mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: BLACS_CTXT
   integer(kind=i4) :: sc_desc(9)
   integer(kind=i4) :: n_basis,n_states
   integer(kind=i4) :: matrix_size,supercell(3)
   integer(kind=i4) :: solver
   integer(kind=i4) :: local_row,local_col,ldm
   integer(kind=i4) :: n
   integer(kind=i4) :: i,j
   integer(kind=i4) :: info

   real(kind=r8) :: n_electrons,frac_occ,sparsity,orb_r_cut
   real(kind=r8) :: k_point(3)
   real(kind=r8) :: e_test,e_ref,e_tol
   real(kind=r8) :: t1,t2
   integer(kind=i4), external :: numroc

   ! VY: Reference value from calculations on Apr 5, 2017.
   real(kind=r8), parameter :: e_elpa  = -126.817462901838_r8

   real(kind=r8), allocatable :: mat_a(:,:),mat_b(:,:),mat_tmp(:,:),e_vec(:,:)
   real(kind=r8), allocatable :: e_val(:)

#ifdef WITH_OPENMP
   integer :: omp_get_max_threads,  required_mpi_thread_level, provided_mpi_thread_level
#endif

   type(elsi_handle) :: elsi_h

   integer na, nev ,nb
   integer*8 natoms, natoms2, natoms3
   character*256 dum
   integer*8 na8
   character*256 filenameA
   character*256 filenameB
   integer lenarg


   ! Initialize MPI
#ifndef WITH_OPENMP
   call MPI_Init(mpierr)
#else
   required_mpi_thread_level = MPI_THREAD_MULTIPLE

   call mpi_init_thread(required_mpi_thread_level,     &
                        provided_mpi_thread_level, mpierr)

   if (required_mpi_thread_level .ne. provided_mpi_thread_level) then
     write(error_unit,*) "MPI ERROR: MPI_THREAD_MULTIPLE is not provided on this system"
     write(error_unit,*) "           only ", mpi_thread_level_name(provided_mpi_thread_level), " is available"
     call exit(77)
   endif

#endif
   mpi_comm_global = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm_global,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm_global,myid,mpierr)

   ! Read command line arguments

   if(myid == 0) then
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      write(*,'("  Now start testing  elsi_ev_real + ELPA")')
      write(*,*)
   endif
   solver = 1

   if(myid==0) then
      call get_command_argument(1,filenameA,lenarg,info)
      call get_command_argument(2,filenameB,lenarg,info)
      if(info/=0) then
         write(error_unit,*) 'Usage: test_real matrix_file'
         call mpi_abort(mpi_comm_global,1,mpierr)
      endif
      OPEN (10, FILE = filenameA, STATUS = 'OLD')
      OPEN (20, FILE = filenameB, STATUS = 'OLD')

      if(info/=0) then
         write(error_unit,*) 'Error: Unable to open ',trim(filenameA)
         call mpi_abort(mpi_comm_global,1,mpierr)
      endif
   endif

   ! Read matrix size
   if(myid==0) read(10,*) dum
   if(myid==0) read(10,*) natoms, natoms2, natoms3
   if(myid==0) na= natoms


   ! Read matrix size
   if(myid==0) read(20,*) dum
   if(myid==0) read(20,*) natoms, natoms2, natoms3
   if(myid==0) nb= natoms

   call mpi_bcast(nb, 1, mpi_integer, 0, mpi_comm_global, mpierr)
   call mpi_bcast(na, 1, mpi_integer, 0, mpi_comm_global, mpierr)
   call mpi_barrier(mpi_comm_global, mpierr) ! Just for safety

   matrix_size = na
   !
   ! Solve half spectrum
   !
   n_states = matrix_size/2
   !n_states = 10
   if(myid==0) print *,'Matrix size: ',na
   if(myid==0) print *,'N_states size: ',n_states

   if(myid==0.and.na.ne.nb) then
         write(6,*) 'Size ofA .ne. B'
         call mpi_abort(mpi_comm_global,1,mpierr)
   endif

   e_ref = e_elpa
   e_tol = 1e-10_r8

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(n_proc))),2,-1
      if(mod(n_proc,npcol) == 0) exit
   enddo
   nprow = n_proc/npcol

   ! Set block size
   blk = 32

   ! Set up BLACS
   BLACS_CTXT = mpi_comm_global
   call BLACS_Gridinit(BLACS_CTXT,'r',nprow,npcol)
   call BLACS_Gridinfo(BLACS_CTXT,nprow,npcol,myprow,mypcol)

   local_row = numroc(matrix_size,blk,myprow,0,nprow)
   local_col = numroc(matrix_size,blk,mypcol,0,npcol)

   ldm = max(local_row,1)

   call descinit(sc_desc,matrix_size,matrix_size,blk,blk,&
                 0,0,BLACS_CTXT,ldm,info)

   allocate(mat_a(local_row,local_col))
   allocate(mat_b(local_row,local_col))
   allocate(e_vec(local_row,local_col))
   allocate(e_val(matrix_size))

   t1 = MPI_Wtime()
   call read_matrix(10, matrix_size, mat_a, ubound(mat_a,1), blk, myprow, mypcol, nprow, npcol)
   if(myid==0) close(10)
   t2 = MPI_Wtime()
   if(myid==0) print *,' A IO time = ', t2-t1
   t1 = MPI_Wtime()
   call read_matrix(20, matrix_size, mat_b, ubound(mat_b,1), blk, myprow, mypcol, nprow, npcol)
   if(myid==0) close(20)
   t2 = MPI_Wtime()
   if(myid==0) print *,' B IO time = ', t2-t1

   ! Initialize ELSI
  !n_electrons = 2.0_r8*n_states
    n_electrons = 0.0_r8

   call elsi_init(elsi_h,solver,1,0,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(elsi_h,mpi_comm_global)
   call elsi_set_blacs(elsi_h,BLACS_CTXT,blk)


   ! Customize ELSI
   call elsi_customize(elsi_h,print_detail=.true.)
   !!call elsi_customize(elsi_h,uplo=2)
   !call elsi_customize(elsi_h,singularity_tolerance=1.0_r8)
   !call elsi_customize(elsi_h,no_singularity_check=.true.)
   !!call elsi_customize(elsi_h,overlap_is_unit=.true.)
   
   t1 = MPI_Wtime()

   call elsi_ev_real(elsi_h,mat_a,mat_b,e_val,e_vec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished test program")')
      write(*,'("  | Total computation time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(elsi_h)

  if(myid == 0) then
      if( n_states .lt. 1025) then
        do i=1, n_states
         print  * , "e=" ,  i  , e_val(i)
        enddo
      endif
  endif
   deallocate(e_val)
   deallocate(e_vec)
   deallocate(mat_b)
   deallocate(mat_a)

   call MPI_Finalize(mpierr)

end program
!-------------------------------------------------------------------------------
subroutine read_matrix(iunit, na, a, lda, nblk, my_prow, my_pcol, np_rows, np_cols)

   implicit none
   include 'mpif.h'

   integer, intent(in) :: iunit, na, lda, nblk, my_prow, my_pcol, np_rows, np_cols
   real*8, intent(out) :: a(lda, *)

   integer i, j, lr, lc, myid, mpierr
   integer, allocatable :: l_row(:), l_col(:)

   real*8, allocatable :: col(:)
   integer*8 :: ii,jj
   real*8 ::dtmp
   integer :: i_code

   ! allocate and set index arrays

   allocate(l_row(na))
   allocate(l_col(na))

   ! Mapping of global rows/cols to local

   l_row(:) = 0
   l_col(:) = 0

   lr = 0 ! local row counter
   lc = 0 ! local column counter

   do i = 1, na

     if( MOD((i-1)/nblk,np_rows) == my_prow) then
       ! row i is on local processor
       lr = lr+1
       l_row(i) = lr
     endif

     if( MOD((i-1)/nblk,np_cols) == my_pcol) then
       ! column i is on local processor
       lc = lc+1
       l_col(i) = lc
     endif

   enddo

   call mpi_comm_rank(MPI_COMM_WORLD,myid,mpierr)
   allocate(col(na))

   read(iunit,*, iostat=i_code) ii, jj, dtmp
   do i=1,na
      col(1:i)=0.0
      if(myid==0)  then
          do while (i.eq.ii)
            col(jj)=dtmp
            read(iunit,*, iostat=i_code) ii, jj, dtmp
            if(i_code<0) goto 100
            if(i_code>0) stop ' error' 
          enddo
 100  continue
      endif
      call mpi_bcast(col,i,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      if(l_col(i) > 0) then
         do j=1,i
            if(l_row(j)>0) a(l_row(j),l_col(i)) = col(j)
         enddo
      endif
      if(l_row(i) > 0) then
         do j=1,i-1
            if(l_col(j)>0) a(l_row(i),l_col(j)) = col(j)
         enddo
      endif
   enddo

   deallocate(l_row, l_col, col)

end subroutine read_matrix
