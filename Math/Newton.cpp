//
//  Newton.cpp
//  Math
//
//  Created by Edward Janne on 5/22/21.
//

#include "Newton.hpp"

double newtonRaphsonSolve(double target, double hint, const Functor &f, const Functor &d, int maxSteps, double tol, double epsilon)
{
    double x0 = hint;
    double x1 = hint;
    int deltaCount = 0;
    double lastDelta = 0.0;
    int n = 0;
    while(maxSteps--) {
        double yp = d(x0);
        if(abs(yp) < epsilon) return x1;
        
        double newTry = f(x0);
        x1 = x0 - (newTry - target) / yp;
        
        if((abs(x1 - x0) / abs(x1)) < tol) return x1;
        
        float delta = abs(x0 - x1);
        if(delta == lastDelta) {
            if(deltaCount > 10) {
                x1 = (x0 + x1) * 0.5;
                return x1;
            }
            deltaCount++;
        }
        lastDelta = delta;
        
        x0 = x1;
        n++;
    }
    
    return x1;
}
