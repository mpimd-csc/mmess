Function handles for solving state-space transformed matrix equations
via Krylov subspace methods.
Assume the linear first-order system of the form

    Ex' = Ax + Bu,
     y  = Cx,

with E = LU invertible.
Then in the state-space transformed case, we are considering the system

    z' = (L\A/U)z + (L\B)u,
    y  = (C/U)z.
 