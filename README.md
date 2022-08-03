[![paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=KYWUWS86GSFGL)

### Hydrodynamics

This is a port of my JavaScript hydrodynamics solver to C++, 
combined with my Tensor library so I can write every dimension case at once without using runtime-specified loop sizes, for the hopes of unrolling and template inlining for faster results 

It depends on:
- GLApp
- Tensor
- Profiler
- Parallel

TO DO:
- external forces of any kind crash - due to singular matrix - most likely due to negative sqrts, most likely due to total energy inconsistency / negative potential energies computed from the external force
- analytic flux jacobian and eignevector inverses
- cut down on all those lookups? pointer offsetting, as ugly as it is
- add arbitrary boundaries
- mouse interaction?
- configuration scripts, especially as I add different EOS's with new variables 

Sources:
	the "Hydrodynamics II" book found online here: http://www.mpia-hd.mpg.de/homes/dullemon/lectures/hydrodynamicsII/ 
	http://www.cfdbooks.com/cfdcodes.html 
	"Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro
	http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf
