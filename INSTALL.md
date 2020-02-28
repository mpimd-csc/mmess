Install Instructions  						{#install}
====================

## Dependencies

Generally, the Software should run fine in plain MATLAB and
OCTAVE. However, for optimal performance of the projection based
approaches it can be advantageous to have the corresponding control
systems toolbox providing optimized solvers for dense Lyapunov and
Riccati equations.


## Installing M.E.S.S.

### ZIP-archive based downloads
Simply unpack the content of the archive to some folder, change to
that folder before you use it and run `mess_path`. This will add this
folder together with all required subfolders to your MATLAB/Octave search
path. If you do not want to do this every time you need to work with
M.E.S.S., simply run `mypath = mess_path;` and add the content of
`mypath` to the default search path permanently.

### mltbx downloads (MATLAB only, recommended for R2014a and above)
This is the Matlab toolbox file based download. Simply open the
downloaded file in the Matlab directory browser and it should
automatically install the toolbox and set the paths correctly.

### tar.gz downloads (Octave only, recommended for version 4.2.2 and above) 
This is the Octave package file intended for use with the Octave
package manager pkg. On the Octave prompt, change to the folder where
the download is located and run `pkg install mess-2.0.tar.gz`. 

Since this package depends on the `control` package for the solution
of small dense Lyapunov and Riccati equations, make sure to have that
package installed already. This can be done, e.g. by running `pkg
install -forge control` with an active internet connection, to
fetch and install the package from Octave-Forge directly.
