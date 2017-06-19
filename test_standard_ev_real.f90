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
! alvarovm.github.io : added Matrix Market format reader 
!
program test_standard_ev_real

   use ELSI_PRECISION, only: r8,i4
   use ELSI
   use elpa_utilities, only : error_unit

   implicit none

   include "mpif.h"

   character(128) :: arg1
   character(128) :: arg2
   character(128) :: arg3

   integer(kind=i4) :: n_proc,nprow,npcol,myid,myprow,mypcol
   integer(kind=i4) :: mpi_comm_global,mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: BLACS_CTXT
   integer(kind=i4) :: sc_desc(9)
   integer(kind=i4) :: info
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: n_states
   integer(kind=i4) :: solver
   integer(kind=i4) :: local_row,local_col,ldm
   integer(kind=i4) :: n
   integer(kind=i4) :: i,j
   integer(kind=i4), external :: numroc

   real(kind=r8) :: t1,t2
   real(kind=r8), allocatable :: mat_a(:,:),mat_b(:,:),mat_tmp(:,:),e_vec(:,:)
   real(kind=r8), allocatable :: e_val(:)

#ifdef WITH_OPENMP
   integer :: omp_get_max_threads,  required_mpi_thread_level, provided_mpi_thread_level
#endif

   type(elsi_handle) :: elsi_h

   integer na, nev
   integer*8 natoms, natoms2, natoms3
   character*256 dum
   integer*8 na8
   character*256 filename
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

   solver = 1 
   if(myid == 0) then
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      write(*,'("  This test program performs the following computational steps:")')
      write(*,*)
      write(*,*)
      write(*,'("  Now start testing  elsi_ev_real + ELPA")')
      write(*,*)
   endif

   if(myid==0) then
      call get_command_argument(1,filename,lenarg,info)
      if(info/=0) then
         write(error_unit,*) 'Usage: test_real matrix_file'
         call mpi_abort(mpi_comm_global,1,mpierr)
      endif
      OPEN (10, FILE = filename, STATUS = 'OLD')

  !   open(10,file=filename, form="unformatted", action='READ',status='OLD',iostat=info)
      if(info/=0) then
         write(error_unit,*) 'Error: Unable to open ',trim(filename)
         call mpi_abort(mpi_comm_global,1,mpierr)
      endif
   endif
   call mpi_barrier(mpi_comm_global, mpierr) ! Just for safety



   ! Read matrix size
   if(myid==0) read(10,*) dum
   if(myid==0) read(10,*) natoms, natoms2, natoms3
   if(myid==0) na= natoms

   call mpi_bcast(na, 1, mpi_integer, 0, mpi_comm_global, mpierr)
   call mpi_barrier(mpi_comm_global, mpierr) ! Just for safety

   matrix_size = na
   n_states = matrix_size
   if(myid==0) print *,'Matrix size: ',na
   if(myid==0) print *,'N_states size: ',n_states




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

   ! Read matrix and distribute
   allocate(mat_a(local_row,local_col))

   t1 = MPI_Wtime()
   call read_matrix(10, matrix_size, mat_a, ubound(mat_a,1), blk, myprow, mypcol, nprow, npcol)
   if(myid==0) close(10)
   t2 = MPI_Wtime()
   if(myid==0) print *,'IO time = ', t2-t1


   ! Initialize ELSI
   call elsi_init(elsi_h,solver,1,0,matrix_size,0.0_r8,n_states)
   call elsi_set_mpi(elsi_h,mpi_comm_global)
   call elsi_set_blacs(elsi_h,BLACS_CTXT,blk)

   allocate(mat_b(1,1)) ! Dummy allocation
   allocate(e_vec(local_row,local_col))
   allocate(e_val(matrix_size))

   ! Customize ELSI
   call elsi_customize(elsi_h,print_detail=.true.)
   call elsi_customize(elsi_h,overlap_is_unit=.true.)
   
   t1 = MPI_Wtime()

   ! Solve problem
   call elsi_ev_real(elsi_h,mat_a,mat_b,e_val,e_vec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished test program")')
      write(*,'("  | Total computation time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(elsi_h)

  if(myid == 0 .and. n_states .le. 1024) then
      do i=1, n_states
      print  * , "e=" ,  i  , e_val(i)
      enddo
  endif
   deallocate(mat_a)
   deallocate(mat_b)
   deallocate(e_val)
   deallocate(e_vec)

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
      !if(myid==0)    read(iunit) col(1:i)
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
