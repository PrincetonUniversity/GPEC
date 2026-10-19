"""Symbolic identities behind the GPEC pmodb |B| columns.

Every function returns ``(label, equation, residual)`` triples; a residual must
vanish identically.  The module is written against a SymPy-shaped namespace
``sp`` so that the same derivation runs on FortSym (``fortsym.sympy``, the
authority used by ``derive.py``).  Labels refer to ``docs/pmodb_gauge.rst``.

GPEC/DCON coordinates (psi, theta, zeta): angles in turns, general
straight-field-line Jacobian ``J(psi, theta)`` (any ``jac_type``),

    B = chi1 grad(zeta - q theta) x grad(psi),   J B^theta = chi1,
    J B^zeta = chi1 q,                           B^psi = 0,

and the surface displacement (Glasser 2016 Eqs. 9 and 18)

    xi_s = (xi x B)_psi = chi1 (q xi^theta - xi^zeta).

Pointwise identities are proved on jets: every field is an affine (1-jet)
polynomial in (psi, theta, zeta) about the origin with independent symbolic
coefficients, and the residual is evaluated at the origin after
differentiation.  Each identity involves at most first derivatives of the
jetted fields, so it holds for all smooth fields iff it holds for all 1-jets.
The metric is general except for ``det g = J**2``, imposed by solving for
``g_33``.
"""

from __future__ import annotations

COORDS = ("psi", "theta", "zeta")


class Dual:
    """First-order jet at the origin: value ``v`` and derivatives ``d[i]``.

    ``d is None`` marks a value whose own derivatives were not carried (a
    derivative of a jet); differentiating it again raises.
    """

    def __init__(self, v, d=(0, 0, 0)):
        self.v, self.d = v, (None if d is None else tuple(d))

    @staticmethod
    def lift(x):
        return x if isinstance(x, Dual) else Dual(x)

    def _need(self):
        if self.d is None:
            raise ValueError("second derivative requested from a 1-jet")
        return self.d

    def __add__(self, o):
        o = Dual.lift(o)
        if self.d is None or o.d is None:
            return Dual(self.v + o.v, None)
        return Dual(self.v + o.v, [a + b for a, b in zip(self.d, o.d)])

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, None if self.d is None else [-a for a in self.d])

    def __sub__(self, o):
        return self + (-Dual.lift(o))

    def __rsub__(self, o):
        return Dual.lift(o) - self

    def __mul__(self, o):
        o = Dual.lift(o)
        if self.d is None or o.d is None:
            return Dual(self.v * o.v, None)
        return Dual(self.v * o.v, [self.v * b + a * o.v for a, b in zip(self.d, o.d)])

    __rmul__ = __mul__

    def inv(self):
        if self.d is None:
            return Dual(1 / self.v, None)
        return Dual(1 / self.v, [-a / self.v ** 2 for a in self.d])

    def __truediv__(self, o):
        return self * Dual.lift(o).inv()

    def __rtruediv__(self, o):
        return Dual.lift(o) * self.inv()

    def __pow__(self, k):
        assert isinstance(k, int) and k >= 1
        out = self
        for _ in range(k - 1):
            out = out * self
        return out


class Setup:
    """Jetted GPEC fields (1-jets at the origin) in a SymPy-shaped namespace."""

    def __init__(self, sp, metric=False):
        self.sp = sp
        self.u = (0, 1, 2)                         # psi, theta, zeta
        self.chi1 = Dual(sp.Symbol("chi1", positive=True))
        self.q = self.jet("q", (0,))
        self.J = self.jet("J", (0, 1))
        # Contravariant displacement components (general, incl. xi_parallel).
        self.xi = [self.jet(f"x{c[0]}", self.u) for c in COORDS]
        self.Bup = [Dual(0), self.chi1 / self.J, self.chi1 * self.q / self.J]
        if metric:
            g = {}
            for (i, j) in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2)):
                g[(i, j)] = self.jet(f"g_{i + 1}{j + 1}", (0, 1))
            minor = g[(0, 0)] * g[(1, 1)] - g[(0, 1)] ** 2
            rest = (-g[(0, 0)] * g[(1, 2)] ** 2 + 2 * g[(0, 1)] * g[(1, 2)] * g[(0, 2)]
                    - g[(1, 1)] * g[(0, 2)] ** 2)
            g[(2, 2)] = (self.J ** 2 - rest) / minor
            self.g = [[g[(min(i, j), max(i, j))] for j in range(3)] for i in range(3)]

    def jet(self, name, deps):
        sp = self.sp
        return Dual(sp.Symbol(f"{name}_0"),
                    [sp.Symbol(f"{name}_{COORDS[i]}") if i in deps else 0 for i in range(3)])

    def d(self, f, i):
        f = Dual.lift(f)
        return Dual(f._need()[i], None)

    def at0(self, e):
        return Dual.lift(e).v

    # --- tensor algebra in the (psi, theta, zeta) chart -----------------------
    @staticmethod
    def eps(i, j, k):
        return (i - j) * (j - k) * (k - i) // 2

    def cross_cov(self, a, b):
        """(a x b)_i = J eps_ijk a^j b^k for contravariant a, b."""
        return [self.J * sum(self.eps(i, j, k) * a[j] * b[k]
                             for j in range(3) for k in range(3)) for i in range(3)]

    def curl_up(self, A):
        """(curl A)^i = eps^ijk d_j A_k / J for covariant A."""
        return [sum(self.eps(i, j, k) * self.d(A[k], j)
                    for j in range(3) for k in range(3)) / self.J for i in range(3)]

    def div(self, v):
        return sum(self.d(self.J * v[i], i) for i in range(3)) / self.J

    def lower(self, v):
        return [sum(self.g[i][j] * v[j] for j in range(3)) for i in range(3)]

    def dot(self, a, b):
        return sum(a[i] * b[i] for i in range(3))

    def xi_s(self, v):
        return self.chi1 * (self.q * v[1] - v[2])


def induction(sp):
    """G1: GPEC's J dB^i (gpeq.f:107-109) is J curl(xi x B)^i for any Jacobian."""
    S = Setup(sp)
    xi, B = S.xi, S.Bup
    JdB = [S.J * c for c in S.curl_up(S.cross_cov(xi, B))]
    xs = S.xi_s(xi)
    d = S.d
    gpec = [S.chi1 * (d(xi[0], 1) + S.q * d(xi[0], 2)),
            d(xs, 2) - S.chi1 * d(xi[0], 0),
            -d(xs, 1) - S.chi1 * d(S.q * xi[0], 0)]
    out = [(f"G1.{i + 1}", f"J dB^{COORDS[i]} = GPEC bw{'ptz'[i]} (pointwise)",
            S.at0(JdB[i] - gpec[i])) for i in range(3)]
    # (xi x B)_psi = xi_s: the surface displacement of Glasser (2016) Eq. 9.
    out.append(("G1.4", "(xi x B)_psi = chi1 (q xi^theta - xi^zeta) = xi_s",
                S.at0(S.cross_cov(xi, B)[0] - xs)))
    # A field-aligned displacement changes neither xi^psi nor xi_s.
    lam = S.jet("lam", S.u)
    shifted = [xi[i] + lam * B[i] for i in range(3)]
    out.append(("G1.5", "xi_s(xi + lam B) = xi_s(xi)",
                S.at0(S.xi_s(shifted) - xs)))
    return out


def fourier_symbols(sp):
    """G1F: Fourier symbols of G1 for xi ~ exp(2 pi i (m theta - n zeta))."""
    I, pi = sp.I, sp.pi
    m, n, chi1 = sp.Symbol("m"), sp.Symbol("n"), sp.Symbol("chi1", positive=True)
    q, q1 = sp.Symbol("q"), sp.Symbol("q1")
    X, X1, Xs = sp.Symbol("xsp"), sp.Symbol("xsp1"), sp.Symbol("xss")
    dth, dze = 2 * pi * I * m, -2 * pi * I * n
    out = []
    out.append(("G1F.1", "bwp = chi1 (m - n q) 2 pi i xsp",
                chi1 * (dth + q * dze) * X - chi1 * (m - n * q) * 2 * pi * I * X))
    out.append(("G1F.2", "bwt = -(chi1 xsp1 + 2 pi i n xss)",
                (dze * Xs - chi1 * X1) + (chi1 * X1 + 2 * pi * I * n * Xs)))
    out.append(("G1F.3", "bwz = -(chi1 (q1 xsp + q xsp1) + 2 pi i m xss)",
                (-dth * Xs - chi1 * (q1 * X + q * X1))
                + (chi1 * (q1 * X + q * X1) + 2 * pi * I * m * Xs)))
    return out


def incompressible_closure(sp):
    """G2: GPEC's xi^theta, xi^zeta (gpeq.f gpeq_contra) follow from div xi = 0."""
    S = Setup(sp)
    xi, J, q, chi1, d = S.xi, S.J, S.q, S.chi1, S.d
    D = sum(d(J * xi[i], i) for i in range(3))          # J div xi
    xs = S.xi_s(xi)
    Jt, Jz = J * xi[1], J * xi[2]
    r1 = (d(Jt, 1) + q * d(Jt, 2)) - d(J * xs, 2) / chi1 + d(J * xi[0], 0) - D
    r2 = (d(Jz, 1) + q * d(Jz, 2)) + q * d(J * xi[0], 0) + d(J * xs, 1) / chi1 - q * D
    r3 = q * Jt - Jz - J * xs / chi1
    return [("G2.1", "(d_t + q d_z)(J xi^t) = d_z(J xi_s)/chi1 - d_psi(J xi^psi) + J div xi",
             S.at0(r1)),
            ("G2.2", "(d_t + q d_z)(J xi^z) = -q d_psi(J xi^psi) - d_t(J xi_s)/chi1"
             " + q J div xi", S.at0(r2)),
            ("G2.3", "q xwt - xwz = J xi_s / chi1", S.at0(r3))]


def columns(sp):
    """G3-G5: eulb, xi_perp components, the gauge of lagb, and llagb."""
    S = Setup(sp, metric=True)
    xi, B, J, q, chi1 = S.xi, S.Bup, S.J, S.q, S.chi1
    Bm = S.jet("Bm", S.u[:2])                     # |B|(psi, theta), axisymmetric
    dB = S.curl_up(S.cross_cov(xi, B))
    dBcov = S.lower(dB)
    B2 = S.dot(S.lower(B), B)
    xs = S.xi_s(xi)
    xv = S.lower(xi)                               # GPEC xvt, xvz (covariant)
    xmt, xmz = J * xi[1], J * xi[2]                # GPEC xmt, xmz (J xi^i)
    sigma = S.dot(xv, B) / B2
    xperp = [xi[i] - sigma * B[i] for i in range(3)]
    out = []
    out.append(("G3", "B.dB = chi1 (bvt + q bvz)/J  (eulb times |B|)",
                S.at0(S.dot(dBcov, B) - chi1 * (dBcov[1] + q * dBcov[2]) / J)))
    out.append(("G4.1", "J xi_perp^theta = xmt - (chi1/B)^2 (xvt + q xvz)/J",
                S.at0(J * xperp[1] - (xmt - chi1 ** 2 / B2 * (xv[1] + q * xv[2]) / J))))
    out.append(("G4.2", "J xi_perp^zeta = xmz - q (chi1/B)^2 (xvt + q xvz)/J",
                S.at0(J * xperp[2] - (xmz - q * chi1 ** 2 / B2 * (xv[1] + q * xv[2]) / J))))
    # B . xi_perp = 0 as a statement about arbitrary vectors: the metric-jet
    # form exceeds FortSym's polynomial-cancellation limit (issue #84), so the
    # projection is proved on independent covariant B_i, contravariant B^i, xi^i.
    Bl = [sp.Symbol(f"Bl_{i}") for i in range(3)]
    Bu = [sp.Symbol(f"Bu_{i}") for i in range(3)]
    xu = [sp.Symbol(f"xu_{i}") for i in range(3)]
    sig = sum(Bl[i] * xu[i] for i in range(3)) / sum(Bl[i] * Bu[i] for i in range(3))
    out.append(("G4.3", "B_i xi_perp^i = 0 for xi_perp^i = xi^i - (B_j xi^j / B_k B^k) B^i",
                sum(Bl[i] * (xu[i] - sig * Bu[i]) for i in range(3))))
    coef = xs / (chi1 * q) + xmz / (J * q) - chi1 ** 2 / (J ** 2 * B2) * (xv[1] + q * xv[2])
    out.append(("G4.4", "xi_perp^theta = xms/(chi1 q) + xmz/(J q)"
                " - (chi1/(J B))^2 (xvt + q xvz)  (curvfun, gpout.f)",
                S.at0(xperp[1] - coef)))
    # Gauge of lagb: xi_g = xi_perp + lam B with xi_g^zeta = 0.
    lam = -xperp[2] / B[2]
    xg = [xperp[i] + lam * B[i] for i in range(3)]
    out.append(("G5.1", "xi_g^zeta = 0 for lam = -xi_perp^zeta / B^zeta", S.at0(xg[2])))
    out.append(("G5.2", "xi_g^theta = xi_s/(chi1 q)", S.at0(xg[1] - xs / (chi1 * q))))
    out.append(("G5.2b", "xi_g^psi = xi^psi", S.at0(xg[0] - xi[0])))
    dBg = S.curl_up(S.cross_cov(xg, B))
    dBp = S.curl_up(S.cross_cov(xperp, B))
    out.append(("G5.3", "curl(xi_g x B) = curl(xi_perp x B)  (theta component)",
                S.at0(dBg[1] - dBp[1])))
    eul = S.dot(dBcov, B)                          # |B| * eulb
    grad = lambda v: sum(v[i] * S.d(Bm, i) for i in range(3))
    lagb = eul + Bm * (xi[0] * S.d(Bm, 0) + xs / (chi1 * q) * S.d(Bm, 1))
    lagb_g = eul + Bm * grad(xg)
    llagb = eul + Bm * grad(xperp)
    llagb_code = lagb + Bm * S.d(Bm, 1) * (xmz / (J * q) - chi1 ** 2 / (J ** 2 * B2)
                                          * (xv[1] + q * xv[2]))
    out.append(("G5.4", "lagb (gpout.f) = dB_L(xi_g), the xi^zeta = 0 gauge",
                S.at0(lagb - lagb_g)))
    out.append(("G5.5", "lagb - dB_L(xi_perp) = -(xi_perp^zeta/q) d_theta B",
                S.at0(lagb - llagb + Bm * xperp[2] / q * S.d(Bm, 1))))
    out.append(("G5.6", "llagb (gpdiag.f) = dB_L(xi_perp)", S.at0(llagb_code - llagb)))
    # Park (2009) Eq. 7 with the div-free xi (GPEC's own xmt, xmz).
    lagb_div_free = eul + Bm * (xi[0] * S.d(Bm, 0) + xi[1] * S.d(Bm, 1))
    out.append(("G5.7", "dB_L(xi, div-free) = lagb + (xmz/(J q)) d_theta B",
                S.at0(lagb_div_free - lagb - Bm * xmz / (J * q) * S.d(Bm, 1))))
    # Boozer-type closed form: B_i xi_perp^i = 0 and the xi_s definition.
    Bc = S.lower(B)
    out.append(("G5.8", "xi_perp^zeta = -(B_theta xi_s/chi1 + q B_psi xi^psi)"
                "/(B_theta + q B_zeta)",
                S.at0(xperp[2] + (Bc[1] * xs / chi1 + q * Bc[0] * xi[0])
                      / (Bc[1] + q * Bc[2]))))
    return out


def force_balance_identity(sp):
    """G6: |B| times (llagb + Bdivxprp + Bkxprp) = -xi_perp . F, F = (curl B) x B - mu0 grad p.

    With pmodb's Bkxprp = xi_perp.grad B + mu0 p' xi^psi / B (force-balance form
    of B kappa.xi_perp, gpout.f curvfun) and Bdivxprp = B div xi_perp, this is
    B.dB + xi_perp.grad(B^2) + B^2 div xi_perp + xi_perp.((curl B) x B) = 0.
    """
    S = Setup(sp, metric=True)
    xi, B = S.xi, S.Bup
    Bc = S.lower(B)
    B2 = S.dot(Bc, B)
    sigma = S.dot(S.lower(xi), B) / B2
    xperp = [xi[i] - sigma * B[i] for i in range(3)]
    dB = S.curl_up(S.cross_cov(xi, B))
    curlB = S.curl_up(Bc)
    jxb = S.cross_cov(curlB, B)                    # ((curl B) x B)_i
    res = (S.dot(S.lower(dB), B) + sum(xperp[i] * S.d(B2, i) for i in range(3))
           + B2 * S.div(xperp) + S.dot(jxb, xperp))
    return [("G6", "B.dB + xi_perp.grad B^2 + B^2 div xi_perp + xi_perp.(curl B x B) = 0",
             S.at0(res))]


def pentrc_kappa_matrices(sp):
    """G7: idcon.f fmodb sband/tband give J B^2 (Bkxprp/B) from xi^psi and xi^alpha.

    PENTRC set_peq (inputs.f90) and GPEC.jl set_perturbation_data! form
    jbb_kapx = smat xsp + tmat xi_alpha with xi_alpha = xms/chi1 (gpout_xclebsch).
    The statement is algebraic at a point (only the values of d_psi B^2 and
    d_theta B^2 enter), so the metric, Jacobian and derivatives are independent
    symbols here:
    sband = J (mu0 p' + d_psi B^2/2) - chi1^2 (d_theta B^2/(2 B^2)) (g12 + q g13)/J,
    tband = chi1^2 (d_theta B^2/(2 B^2)) (g23 + q g33)/J.
    """
    chi1, q, J, P, Dp, Dt = sp.symbols("chi1 q J P1 Dp Dt")
    # B^2 = B^i B_i = (chi1/J)(B_theta + q B_zeta) needs no det g = J^2, so the
    # metric components are independent here.
    g11, g12, g13, g22, g23, g33 = sp.symbols("g11 g12 g13 g22 g23 g33")
    g = [[g11, g12, g13], [g12, g22, g23], [g13, g23, g33]]
    xi = list(sp.symbols("xp xt xz"))
    B = [0, chi1 / J, chi1 * q / J]
    Bl = [sum(g[i][j] * B[j] for j in range(3)) for i in range(3)]
    B2 = sum(Bl[i] * B[i] for i in range(3))
    sigma = sum(Bl[i] * xi[i] for i in range(3)) / B2
    xperp = [xi[i] - sigma * B[i] for i in range(3)]
    xs = chi1 * (q * xi[1] - xi[2])
    lhs = J * (xperp[0] * (P + Dp) + xperp[1] * Dt)
    sband = J * (P + Dp) - chi1 ** 2 * (Dt / B2) * (g12 + q * g13) / J
    tband = chi1 ** 2 * (Dt / B2) * (g23 + q * g33) / J
    res = lhs - sband * xi[0] - tband * xs / chi1
    # The residual is linear in xi; proving it for each basis displacement keeps
    # every zero test within FortSym's 12-variable polynomial limit (issue #84).
    out = []
    for k, name in enumerate(("psi", "theta", "zeta")):
        r = res
        for j in range(3):
            r = r.subs(xi[j], 1 if j == k else 0)
        out.append((f"G7.{k + 1}", "J B^2 (xi_perp.grad B + mu0 p' xi^psi/B)/B"
                    f" = sband xsp + tband xms/chi1, xi = e_{name}", r))
    return out


def cartesian_kappa(sp):
    """C1, C2: Cartesian 1-jet forms of the kinematic and curvature identities."""
    x, y, z = sp.symbols("x y z")
    X = (x, y, z)
    beta = sp.Symbol("beta", positive=True)
    g = [[sp.Symbol(f"h_{i}{j}") for j in range(3)] for i in range(3)]
    B = [[0, 0, beta][i] + sum(g[i][j] * X[j] for j in range(3)) for i in range(3)]
    d = lambda f, c: f.diff(c) if hasattr(f, "diff") else 0
    dot = lambda a, b: sum(a[i] * b[i] for i in range(3))
    cross = lambda a, b: [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                          a[0] * b[1] - a[1] * b[0]]
    curl = lambda v: [d(v[2], y) - d(v[1], z), d(v[0], z) - d(v[2], x),
                      d(v[1], x) - d(v[0], y)]
    div = lambda v: sum(d(v[i], X[i]) for i in range(3))
    at0 = lambda e: e.subs(x, 0).subs(y, 0).subs(z, 0)
    Bm = sp.sqrt(dot(B, B))
    b = [c / Bm for c in B]
    kappa = [sum(b[j] * d(b[i], X[j]) for j in range(3)) for i in range(3)]
    a = [sp.Symbol(f"a_{i}") + sum(sp.Symbol(f"a_{i}{j}") * X[j] for j in range(3))
         for i in range(3)]
    xp = cross(B, a)
    gradB = lambda v: sum(v[j] * d(Bm, X[j]) for j in range(3))
    dBpar = dot(b, curl(cross(xp, B)))
    return [("C1", "b.curl(xi_perp x B) + xi_perp.grad B = -B(div xi_perp + kappa.xi_perp)",
             at0(dBpar + gradB(xp) + Bm * (div(xp) + dot(kappa, xp)))),
            ("C2", "B kappa.xi_perp = xi_perp.grad B + xi_perp.((curl B) x B)/B",
             at0(Bm * dot(kappa, xp) - gradB(xp) - dot(xp, cross(curl(B), B)) / Bm))]


def regularisation(sp):
    """R1-R3: Park-regularised components in GPEC (reg_flag) on one harmonic.

    Diagonal model: one harmonic, J independent of theta.  R = s^2/(s^2 + reg^2),
    s = m - n q.  gpeq_sol uses xmp1 = R xsp1 in bmt/bmz and xms; gpeq_contra
    builds xwt/xwz from xmp1 and multiplies by R again; gpout_pmodb uses
    R xsp1 explicitly in divxprp.
    """
    I, pi = sp.I, sp.pi
    m, n, q, chi1 = sp.symbols("m n q chi1")
    J0, J1 = sp.symbols("J_0 J_1")                 # J and dJ/dpsi (theta-independent)
    X, X1, Xs, R = sp.symbols("xsp xsp1 xms R")
    s = m - n * q
    xwt = -(J0 * R * X1 + J1 * X + 2 * pi * I * n / chi1 * J0 * Xs) / (2 * pi * I * s)
    xwz = -(q * J0 * R * X1 + q * J1 * X + 2 * pi * I * m / chi1 * J0 * Xs) / (2 * pi * I * s)
    xmt, xmz = R * xwt, R * xwz
    jdiv = J0 * R * X1 + J1 * X + 2 * pi * I * m * xmt - 2 * pi * I * n * xmz
    out = [("R1", "J div xi (reg) = (1 - R)(J R xsp1 + J' xsp)",
            jdiv - (1 - R) * (J0 * R * X1 + J1 * X))]
    jdiv_unreg = (J0 * X1 + J1 * X + 2 * pi * I * m * xwt.subs(R, 1)
                  - 2 * pi * I * n * xwz.subs(R, 1))
    out.append(("R2", "J div xi (R = 1) = 0", jdiv_unreg))
    # Divergence of the regularised Eulerian field: J div(dB) = d_psi bwp + ... .
    bwp1 = 2 * pi * I * chi1 * (s * X1 - n * sp.Symbol("q1") * X)
    bmt = -(chi1 * R * X1 + 2 * pi * I * n * Xs)
    bmz = -(chi1 * (sp.Symbol("q1") * X + q * R * X1) + 2 * pi * I * m * Xs)
    out.append(("R3", "J div dB (reg) = 2 pi i chi1 (m - n q)(1 - R) xsp1",
                bwp1 + 2 * pi * I * m * bmt - 2 * pi * I * n * bmz
                - 2 * pi * I * chi1 * s * (1 - R) * X1))
    return out


def cgl(sp):
    """K1-K2: the PENTRC/GPEC.jl fcgl integrand against the CGL closure.

    D = div xi_perp, K = kappa.xi_perp, dB/B = -(D + K) (C1).  CGL:
    dp_perp/p = -D + dB/B, dp_par/p = -3D - 2 dB/B (no xi_par).  Logan (2013)
    Eq. 4 energy density, -(1/2)[(dp_par - dp_perp) dB/B + dp_par D]/p.
    """
    D, K = sp.symbols("D K")
    b = -(D + K)
    dpe, dpa = -D + b, -3 * D - 2 * b
    w = -((dpa - dpe) * b + dpa * D) / 2
    code = lambda k: (sp.Rational(5, 3) * D ** 2 + sp.Rational(1, 3) * (D + 3 * k) ** 2) / 2
    return [("K1", "CGL energy density = (1/2)[(5/3) D^2 + (1/3)(D + 3K)^2]", w - code(K)),
            ("K2neg", "fcgl with kapx = -(dbob + divx)/2 equals the CGL density",
             w - code(-(b + D) / 2))]


def drive_gauge(sp):
    """K3, K4: gauge (in)variance of the Logan (2013) Eq. 5 bounce drive.

    Along a field line (arc length l, M = 1): v = sqrt(2(E - mu B)), integrand
    (1/v)[(2E - 3 mu B) dB/B + (2E - 2 mu B) div xi].  Under xi -> xi + lam B,
    dB -> dB + lam B dB/dl and div xi -> div xi + B dlam/dl (C1, Park 2011 Eq. 4).
    """
    l, E, mu = sp.Symbol("l"), sp.Symbol("E", positive=True), sp.Symbol("mu", positive=True)
    b0, b1, l0, l1 = sp.symbols("b_0 b_1 lam_0 lam_1")
    B = b0 + b1 * l
    lam = l0 + l1 * l
    v = sp.sqrt(2 * (E - mu * B))
    dBdl, dlam = B.diff(l), lam.diff(l)
    change = ((2 * E - 3 * mu * B) * lam * dBdl + (2 * E - 2 * mu * B) * B * dlam) / v
    total = (v * B * lam).diff(l)
    # pmodb_file path: dbob = lagb/B (gauge xi_g), divx = -(lagb + Bkxprp)/B.
    change_pmodb = ((2 * E - 3 * mu * B) * lam * dBdl - (2 * E - 2 * mu * B) * lam * dBdl) / v
    # Multiplied by v > 0 (away from turning points): FortSym leaves the
    # unmultiplied sqrt quotient UNKNOWN.
    return [("K3", "v [drive change under xi -> xi + lam B - d/dl (v B lam)] = 0",
             ((change - total) * v).subs(l, 0)),
            ("K4", "pmodb_file-path drive change = -mu B lam dB/dl / v",
             (change_pmodb + mu * B * lam * dBdl / v).subs(l, 0))]


def all_identities(sp):
    return (induction(sp) + fourier_symbols(sp) + incompressible_closure(sp)
            + columns(sp) + pentrc_kappa_matrices(sp)
            + cartesian_kappa(sp) + regularisation(sp)
            + [c for c in cgl(sp) if not c[0].endswith("neg")] + drive_gauge(sp))


def open_gaps(sp):
    """Identities FortSym does not decide; their content is covered elsewhere.

    G6 is the coordinate form of C1 + C2 (proved on Cartesian 1-jets) with
    det g = J^2 imposed through g_33.  The zero test exceeds the native
    polynomial-term limit (lazy-fortran/fortsym#84 class); it is recorded, not
    counted as proved.
    """
    return force_balance_identity(sp)


def negative_controls(sp):
    """Residuals that must NOT vanish (the claims they encode are false)."""
    S = Setup(sp, metric=True)
    xi, B, J, q, chi1 = S.xi, S.Bup, S.J, S.q, S.chi1
    Bm = S.jet("Bm", S.u[:2])
    B2 = S.dot(S.lower(B), B)
    sigma = S.dot(S.lower(xi), B) / B2
    xperp = [xi[i] - sigma * B[i] for i in range(3)]
    xs = S.xi_s(xi)
    lagb_adv = xi[0] * S.d(Bm, 0) + xs / (chi1 * q) * S.d(Bm, 1)
    perp_adv = sum(xperp[i] * S.d(Bm, i) for i in range(3))
    out = [("N1", "lagb is the xi_perp-only Lagrangian |B| (false)",
            S.at0(lagb_adv - perp_adv))]
    out.append(("N2", "lagb is Park 2009 Eq. 7 with the div-free xi (false)",
                S.at0(lagb_adv - (xi[0] * S.d(Bm, 0) + xi[1] * S.d(Bm, 1)))))
    out += [c for c in cgl(sp) if c[0].endswith("neg")]
    return out
