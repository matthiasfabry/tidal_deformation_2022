! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_binary_extras

   use star_lib
   use star_def
   use const_def
   use chem_def
   use num_lib
   use binary_def
   use binary_lib
   use run_star_extras

   implicit none

   integer, parameter :: ilx_pre_ms = 1

contains

   ! extras_ functions
   subroutine extras_binary_controls(binary_id, ierr)
      integer :: binary_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      ! Set these function pointers to point to the functions you wish to use in
      ! your run_binary_extras. Any which are not set, default to a null_ version
      ! which does nothing.
      b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
      b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
      b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
      b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

      b% extras_binary_startup => extras_binary_startup
      b% extras_binary_start_step => extras_binary_start_step
      b% extras_binary_check_model => extras_binary_check_model
      b% extras_binary_finish_step => extras_binary_finish_step
      b% extras_binary_after_evolve => extras_binary_after_evolve

      ! Once you have set the function pointers you want, then uncomment this (or set it in your
      ! star_job inlist)
      ! to disable the printed warning message,
      b% warn_binary_extra = .false.

   end subroutine extras_binary_controls

   integer function extras_binary_startup(binary_id, restart, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      logical, intent(in) :: restart
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      if (.not. restart) then
         b% lxtra(ilx_pre_ms) = .true.
      end if
      ! b% s1% job% warn_run_star_extras = .false.
      extras_binary_startup = keep_going
   end function  extras_binary_startup

   integer function extras_binary_start_step(binary_id, ierr)
      type (binary_info), pointer :: b
      type (star_info), pointer :: s
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      call binary_ptr(binary_id, b, ierr)
      extras_binary_start_step = keep_going

      if (ierr /= 0) return  ! failure in  binary_ptr

      extras_binary_start_step = keep_going
   end function  extras_binary_start_step

   integer function extras_binary_check_model(binary_id)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      extras_binary_check_model = keep_going

   end function extras_binary_check_model

   integer function extras_binary_finish_step(binary_id)
      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      type (binary_info), pointer :: b
      type (star_info), pointer :: s
      integer, intent(in) :: binary_id
      integer :: ierr
      extras_binary_finish_step = keep_going
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return ! failure in  binary_ptr
      call star_ptr(b% star_ids(1), s, ierr)
      if(ierr/=0) return

      if (s% model_number > 1000 .and. b% lxtra(ilx_pre_ms)) then
         extras_binary_finish_step = terminate
         write(*, *) "Terminate due to pre-MS evolution taking too long"
         return
      end if

      if (b% lxtra(ilx_pre_ms) .and. &
            abs(log10(abs(s% L_nuc_burn_total * Lsun / s% L(1)))) < 0.005 .and. &
            s% star_age > 1d2) then
         ! if here, primary reached thermal equilibrium (reached ZAMS), so activate RLOF
         ! this is the amount of overflow of a q=1 system at L2, anything more than this
         ! is too much
         b% lxtra(ilx_pre_ms) = .false.
         b% ignore_rlof_flag = .false.
         !b% s1% max_mdot_redo_cnt = 100
         s% scale_max_correction = 0.02d0
         s% ignore_species_in_max_correction = .true.
         if (b% point_mass_i /= 2) then
            !b% s2% max_mdot_redo_cnt = 100
            b% s2% scale_max_correction = 0.02d0
            b% s2% ignore_species_in_max_correction = .true.
         end if
         write(*, *) "ZAMS reached, allow RLOF!"
      else if (b% lxtra(ilx_pre_ms) .and. &
            (abs(log10(abs(s% L_nuc_burn_total * Lsun / s% L(1)))) > 0.005 .or. &
                  s% star_age < 1d2)) then
         ! if here, still not in ZAMS, keep period fixed

         b% period = b% initial_period_in_days * 86400
         b% separation = pow(s% cgrav(1)*(b% m(1)+b% m(2)) * b% period * b% period/(4*pi2), one_third)
         b% angular_momentum_j = b% m(1) * b% m(2) * sqrt(s% cgrav(1) * &
               b% separation / (b% m(1) + b% m(2)))
         b% rl(1) = binary_eval_rlobe(b% m(1), b% m(2), b% separation)
         b% rl(2) = binary_eval_rlobe(b% m(2), b% m(1), b% separation)
         b% rl_relative_gap(1) = (b% r(1) - b% rl(1)) / b% rl(1) ! gap < 0 means out of contact
         b% rl_relative_gap(2) = (b% r(2) - b% rl(2)) / b% rl(2) ! gap < 0 means out of contact

         !         keep stars synchronized
         if (b% point_mass_i /= 1 .and. s% rotation_flag) &
               call star_set_uniform_omega(b% s1% id, 2 * pi / b% period, ierr)
!               call star_set_uniform_omega(b% s2% id, 2 * pi / b% period, ierr)
      end if

      if(b% r(1) >= 0.99 * b% rl(1))then
         extras_binary_finish_step=terminate
         write(*, *) 'nearing L1 overflow, stopping'
      end if

   end function extras_binary_finish_step

   subroutine extras_binary_after_evolve(binary_id, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if

   end subroutine extras_binary_after_evolve

   ! functions for extra data
   integer function how_many_extra_binary_history_header_items(binary_id)
      use binary_def, only : binary_info
      integer, intent(in) :: binary_id
      how_many_extra_binary_history_header_items = 0
   end function how_many_extra_binary_history_header_items

   subroutine data_for_extra_binary_history_header_items(binary_id, n, names, vals, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id, n
      character (len = maxlen_binary_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
   end subroutine data_for_extra_binary_history_header_items

   integer function how_many_extra_binary_history_columns(binary_id)
      use binary_def, only : binary_info
      integer, intent(in) :: binary_id
      how_many_extra_binary_history_columns = 1
   end function how_many_extra_binary_history_columns

   subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id, n
      character (len = maxlen_binary_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      call binary_ptr(binary_id, b, ierr)
      if(ierr /= 0) return
      names(1) = 'Pre_zams'
      if(b% lxtra(ilx_pre_ms)) then
         vals(1) = 1
      else
         vals(1) = 0
      end if
   end subroutine data_for_extra_binary_history_columns

   ! custom physics functions


end module run_binary_extras
