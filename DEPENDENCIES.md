Dependencies of M-M.E.S.S. 2.0
============================================================================

Matlab R2014a and above, or Octave 4.0 and above. 

Some functions can benefit from the Control Systems Toolbox in either
of the above, but fallbacks exist in case it is not available.

Matlab R2017a has improved handling negative definite matrices with "\". 
We recommend using this version or later ones for optimal performance.

Note that Matlab R2017a R2017b contain a bug in "\" that can cause
extraordinarily slow computations with certain block structured
matrices. Use spparms('usema57', 0); to fix this, or upgrade to at
least R2017b update 5.
