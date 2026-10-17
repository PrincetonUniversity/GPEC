program gen_supplied_rates
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, exp, operator(*)
    use fortsym_diff, only: diff
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: logarithm, logarithm_derivative, roots(2)
    type(kernel_ir_t) :: ir
    type(kernel_emit_spec_t) :: spec
    type(str_t) :: source
    character(2048) :: output, revision
    character(:), allocatable :: message
    integer :: unit
    logical :: ok

    call get_command_argument(1,output)
    call get_command_argument(2,revision)
    if(len_trim(output)==0.or.len_trim(revision)==0) &
        error stop 'usage: gen_supplied_rates OUTPUT FORTSYM_REVISION'
    call arena%init()
    logarithm=sym(arena,'logarithm')
    logarithm_derivative=sym(arena,'logarithm_derivative')
    roots(1)=exp(logarithm)
    roots(2)=diff(roots(1),logarithm)*logarithm_derivative
    call lower_kernel_ir(roots,ir,ok,message)
    if(.not.ok) error stop message
    spec%name=str('supplied_rate_leaf')
    spec%args=[str('logarithm'),str('logarithm_derivative')]
    spec%outputs=[str('rate'),str('rate_derivative')]
    spec%generator=str('gen/gen_supplied_rates.f90')
    spec%generator_revision=str(trim(revision))
    spec%regenerate_command=str( &
        'gen_supplied_rates pentrc/supplied_rates.inc FORTSYM_REVISION')
    source=emit_fortran_kernel_ir(ir,spec,ok,message)
    if(.not.ok) error stop message
    open(newunit=unit,file=trim(output),status='replace')
    write(unit,'(a)',advance='no') chars(source)
    close(unit)
end program gen_supplied_rates
