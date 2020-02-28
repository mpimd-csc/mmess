Dependencies of M-M.E.S.S. 2.0
==============================

Matlab R2014a and above, or Octave 4.0 and above.

Some functions can benefit from the Control Systems Toolbox in Matlab or the
Control package in Octave, but fallbacks exist in case those are not available.

Matlab R2017a has improved handling of negative definite matrices with in the
"backslash" operator. We recommend using this version or later ones for optimal
performance.

Note that Matlab R2017a and R2017b also contain a bug in "backslash" that can
cause extraordinarily slow computations with certain block structured matrices.
Use
```
spparms('usema57', 0);
```
to fix this, or upgrade to at least R2017b update 5.
