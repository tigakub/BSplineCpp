//
//  Functor.hpp
//  Math
//
//  Created by Edward Janne on 5/22/21.
//

#ifndef Functor_hpp
#define Functor_hpp

#include <stdio.h>

class Functor
{
    public:
        Functor() { }
        virtual ~Functor() { }
    
        virtual double operator()(double x) const = 0;
};

#endif /* Functor_hpp */
