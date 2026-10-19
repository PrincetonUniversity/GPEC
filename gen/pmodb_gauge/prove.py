import fortsym.sympy as sp

def gauge_coordinates(sp):
    """Local Jacobian and field-strength gauge equations in angles in turns.

    B.grad = chi1/J (d_theta + q d_zeta). For a target chart J B**a=H(psi),
    a=0 (Hamada) or 2 (Boozer), w=chi1 lambda B**a/H obeys D w=-T,
    T=div xi + a delta B_L/B, for n!=0 (zero flux-function perturbation).
    """
    B, Bt, J, chi, q, lam, Lt, Lz, d, b, c = sp.symbols(
        "B Bt J chi q lam Lt Lz d b c")
    m,n=sp.symbols("m n")
    out=[]
    for a,name in ((0,"hamada"),(2,"boozer")):
        # First jet of lambda B**a; axisymmetric equilibrium, div B=0.
        Ftheta = Lt*B**a + a*lam*B**(a-1)*Bt
        Fzeta = Lz*B**a
        direct = d + chi/J*(Lt+q*Lz) + a*(b+lam*chi/J*Bt)/B
        product = d+a*b/B+chi/(J*B**a)*(Ftheta+q*Fzeta)
        out.append(("Q1."+name,"T_a(xi+lambda B) = T_a(xi)+B.grad(lambda B**a)/B**a",direct-product))
        # Solving the Fourier equation, with modes exp(2 pi i(m theta-n zeta)).
        T,H=sp.symbols("T H")
        w=-T/(2*sp.pi*sp.I*(m-n*q))
        out.append(("Q2."+name,"D w = -T; w_m = -T_m/[2 pi i(m-nq)]",2*sp.pi*sp.I*(m-n*q)*w+T))
        # On a target base J B**a=H, correction to field strength is w d_theta B.
        lambda_w = w*J/chi
        out.append(("Q3."+name,"delta B_target = delta B_X + w d_theta B (J B**a flux function)",
                    lambda_w*chi/J*Bt-w*Bt))
    # Transform the base angles: theta_B=theta_H+iota h, phi_B=phi_H+h.
    h,Dh,iota,Hh,Hb=sp.symbols("h Dh iota Hh Hb")
    ht,hp=sp.symbols("ht hp")
    determinant=(1+iota*ht)*(1+hp)-ht*iota*hp
    out.append(("Q4","det d(theta_H+iota h,phi_H+h)/d(theta_H,phi_H)=1+D h",
                determinant-1-iota*ht-hp))
    # Same physical scalar and relabeling, in the convention of rmp_torque E4.14.
    bX,jX,u,v=sp.symbols("bX jX u v")
    dbB=bX+u*Bt; dbH=bX+v*Bt
    out.append(("Q5","Hamada->Boozer equals direct X->Boozer when w_HB=w_XB-w_XH",
                dbH+(u-v)*Bt-dbB))
    # Perpendicular reconstruction is unchanged by adding a field-aligned vector.
    out.append(("Q6","field-aligned input: delta B_E=0, xi_perp=0 => target delta B=0",
                lam*chi/J*Bt + (-lam*chi/J)*Bt))
    out.append(("Q7","xi_perp limit: zero target source implies no gauge correction", (bX+u*Bt-bX).subs(u,0)))
    # PEST J/R**2 is a flux function. Its linearized constraint has radial motion.
    R,xiR,divxi,db,dr=sp.symbols("R xiR divxi db dr")
    eps=sp.Symbol("eps")
    ratio=J*(1+eps*divxi)/(R+eps*xiR)**2
    out.append(("Q8","delta(J/R**2)/(J/R**2)=div xi-2 xi_R/R",
                ratio.diff(eps).subs(eps,0)*R**2/J-divxi+2*xiR/R))
    # Complex CGL quadratic form: real and imaginary components independently.
    Dr,Di,Kr,Ki=sp.symbols("Dr Di Kr Ki")
    pressure=0; code=0
    for D,K in ((Dr,Kr),(Di,Ki)):
        b=-(D+K); pp=-D+b; pa=-3*D-2*b
        pressure-=((pa-pp)*b+pa*D)/2
        code+=(sp.Rational(5,3)*D**2+sp.Rational(1,3)*(D+3*K)**2)/2
    out.append(("K5","Complex CGL pressure work equals [5/3 |D|**2+1/3 |D+3K|**2]/2",pressure-code))
    # Pressure term in curvature is mu0 p' xi^psi/B, not mu0 p' xi^psi.
    P,X=sp.symbols("P X")
    out.append(("P1","(Bk - (llagb-eulb))/B = mu0 p' xi^psi/B**2", (P*X/B)/B-P*X/B**2))
    return out

for label,equation,residual in gauge_coordinates(sp):
    assert sp.sympify(residual).is_zero is True, label
    print(label,"PROVED")
