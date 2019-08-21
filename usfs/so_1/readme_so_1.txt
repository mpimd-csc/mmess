Function Handles for Second Order Equations, e.g.
Linear Beam 1D Example from 
http://morwiki.mpi-magdeburg.mpg.de/morwiki/index.php/Linear_1D_Beam

Second Order System
Mx" + D x' + K x 	= B u
y 	= C x

Transformed to First Order System
|-K 0||x' |= |0  -K||x | + |0|
|0  M||x''|  |-K -D||x'|   |B|u
   E   x'  =  A      x   +  B

Attention the Matrix M D K are symmetric and quadratic.
K is a fullrank Matrix.
The fieldnames have to end with _  to indicate that the Data 
are inputdata for the Algorithm.
eqn.M_
eqn.K_
eqn.E_
eqn.B
