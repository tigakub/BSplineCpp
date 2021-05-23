//
//  Legendre.hpp
//  Math
//
//  Created by Edward Janne on 5/22/21.
//

#ifndef Legendre_hpp
#define Legendre_hpp

#include <stdio.h>

#include "Functor.hpp"

double legendreIntegrate(int order, double from, double to, const Functor &f);

#endif /* Legendre_hpp */
