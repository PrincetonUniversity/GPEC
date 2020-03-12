module cvode
   use, intrinsic :: iso_c_binding
   use dll_module
   implicit none

   contains
   
   function cxintgrnd(neq,x,y,ydot,energy_imaxis,ximag,xnutype,energy_we,energy_nuk,energy_leff,xnufac,energy_wb,energy_n, &
      energy_wd,xf0type,energy_wn,energy_wt,qt) result(output)
 
   ! interface to our shared lib
   abstract interface
      function xintgrnd(neq,x,y,ydot,energy_imaxis,ximag,xnutype,energy_we,energy_nuk,energy_leff,xnufac, energy_wb, energy_n, &
            energy_wd,xf0type,energy_wn,energy_wt,qt) result(output)
         use, intrinsic :: iso_c_binding
         implicit none
         integer, intent(in) ::  neq
         real*8, intent(in) :: x, y(neq), ydot(neq)
         real*8 ximag, energy_we, energy_nuk, energy_leff, xnufac, energy_wb, energy_wd, energy_wn, energy_wt
         logical :: energy_imaxis, qt
         real*8 :: output
         integer :: xnutype, energy_n, xf0type
      end function xintgrnd
   end interface

   abstract interface
      function test(x) result(output)
         use, intrinsic :: iso_c_binding
         implicit none
         real*8 :: x
         real*8 ximag, energy_we, energy_nuk, energy_leff, xnufac, energy_wb, energy_wd, energy_wn, energy_wt
         logical :: energy_imaxis, qt
         real*8 :: output
         integer :: xnutype, energy_n, xf0type
      end function test
   end interface
 
   type(os_type) :: os
   type(dll_type) :: dll
   integer :: errstat
   character(1024) :: errmsg
   type(c_funptr) :: cfun
   procedure(xintgrnd), pointer :: fproc
   logical :: exists, energy_imaxis, qt
   !something
        integer ::  neq
        real*8 x, y(neq), ydot(neq)
        real*8 ximag, energy_we, energy_nuk, energy_leff, xnufac, energy_wb, energy_wd, energy_wn, energy_wt
   integer :: xnutype, energy_n, xf0type
   real :: output
 
   call init_os_type(1,os)
   call init_dll(dll)
 
   dll%filename="/p/gpec/users/rgates/GPEC/cvode/libxintgrnd.so"
   dll%procname="xintgrnd"

        call load_dll(os, dll, errstat, errmsg )
   if (errstat == 0) then
        call c_f_procpointer(dll%procaddr,fproc)
        output = fproc(neq,x,y,ydot,energy_imaxis,ximag,xnutype,energy_we,energy_nuk,energy_leff,xnufac,energy_wb,energy_n, &
            energy_wd, xf0type, energy_wn, energy_wt, qt)
 
        call free_dll (os, dll, errstat, errmsg )
        if (errstat /= 0) then
            write(*,*)"load_dll: errstat=", errstat
            print *, errmsg
        end if
   else
        write(*,*)"load_dll: errstat=", errstat
        print *, errmsg
   end if
   end function
 
end module cvode
