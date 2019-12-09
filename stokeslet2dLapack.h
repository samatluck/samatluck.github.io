/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 
 Use Lapack
*/

#include "stokeslet2d.h"

class Rstokeslet2DLapack: public Rstokeslet2D {
public:
	Rstokeslet2DLapack(int numberOfParticles);
	virtual ~Rstokeslet2DLapack();

	virtual bool solve();
};
