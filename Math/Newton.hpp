//
//  Newton.hpp
//  Math
//
//  Created by Edward Janne on 5/22/21.
//

#ifndef Newton_hpp
#define Newton_hpp

#include <math.h>
#include <stdio.h>

#include "Functor.hpp"

double newtonRaphsonSolve(double target, double hint, const Functor &f, const Functor &d, int maxSteps = 100, double tol = pow(10.0, -4.0), double epsilon = pow(10.0, -4.0));

#endif /* Newton_hpp */
