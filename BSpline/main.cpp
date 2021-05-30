//
//  main.cpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#include <iostream>
#include "BSpline.hpp"
#include "Parametizer.hpp"

float cps[] =
{
    -2.0, -1.0,
    -1.0, 1.0,
    -0.25, 1.0,
    0.25, -1.0,
    1.0, -1.0,
    2.0, 1.0
};

float knots[10];

int main(int argc, const char * argv[])
{
    BSpline spline(cps, knots, 2, 6);
    spline.init(6);
    
    Parametizer param(spline);
    param.init();
    
    vector<float> parametization = param.parametizeSigmoidal(10);
    
    cout << "Spline length: " << param.length << endl;
    for(int i = 1; i < parametization.size(); i++) {
        float t0 = parametization[i-1];
        float t1 = parametization[i];
        float a0 = param.arcLength(t0);
        float a1 = param.arcLength(t1);
        cout << "Arc " << i << " t: " << t1 << " l: " << a1 << " d: " << a1 - a0 << endl;
    }
    
    return 0;
}
