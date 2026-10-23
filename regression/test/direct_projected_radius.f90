program direct_projected_radius
  use read_eq_mod
  implicit none
  integer :: j, idx
  real(r8) :: radius, angle, r, z, expected, residual

  call get_command_argument(1,eq_filename)
  mpsi=32
  mtheta=64
  psilow=.05_r8
  psihigh=.95_r8
  etol=1e-5_r8
  grid_type='ldp'
  jac_type='hamada'
  power_bp=0
  power_b=0
  power_r=0
  verbose=.false.
  open(out_unit,file='gpec.out',status='new')
  call read_eq_efit()

  idx=mpsi/2
  expected=psio*(1-sq%xs(idx))
  residual=0
  do j=0,mtheta-1
    call bicube_eval(rzphi,sq%xs(idx),real(j,r8)/mtheta,0)
    radius=sqrt(rzphi%f(1))
    angle=twopi*(real(j,r8)/mtheta+rzphi%f(2))
    r=ro+radius*cos(angle)
    z=zo+radius*sin(angle)
    call bicube_eval(psi_in,r,z,0)
    residual=max(residual,abs(psi_in%f(1)-expected)/psio)
  end do
  write(*,'(es25.16)') residual
end program
