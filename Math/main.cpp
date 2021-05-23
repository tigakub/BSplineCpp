//
//  main.cpp
//  Math
//
//  Created by Edward Janne on 5/22/21.
//

#include <iostream>

#include "WayPoint.hpp"
#include "BSpline.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    vector<WayPoint> wv;
    float p0[] = { -2.0, -1.0, 0.0 };
    float p1[] = { -1.0, 1.0, 0.0 };
    float p2[] = { -0.25, 1.0, 0.0 };
    float p3[] = { 0.25, -1.0, 0.0 };
    float p4[] = { 1.0, -1.0, 0.0 };
    float p5[] = { 2.0, 1.0, 0.0 };
    
    
    WayPoint w0(p0, 3), w1(p1, 3), w2(p2, 3), w3(p3, 3), w4(p4, 3), w5(p5,3);
    wv.push_back(w0);
    wv.push_back(w1);
    wv.push_back(w2);
    wv.push_back(w3);
    wv.push_back(w4);
    wv.push_back(w5);
    
    BSpline<WayPoint> spline(wv);
    vector<float> uniformt = spline.parameterizeLinear(10);
    
    cout << "Spline length: " << spline.length() << endl << endl;
    cout << "Linear interpolation" << endl;
    for(int i = 1; i < int(uniformt.size()); i++) {
        float t0 = uniformt[i - 1];
        float t1 = uniformt[i];
        float a0 = spline.arcLength(t0);
        float a1 = spline.arcLength(t1);
        cout  << "Arc " << i << " t: " << t1 << " dt: " << t1 - t0 << " l: " << a1 << " d: " << (a1 - a0) << endl;
    }
    cout << endl;
    
    vector<float> sigmoidt = spline.parameterizeSigmoid(10);
    
    cout << "Sigmoid interpolation" << endl;
    for(int i = 1; i < int(sigmoidt.size()); i++) {
        float t0 = sigmoidt[i - 1];
        float t1 = sigmoidt[i];
        float a0 = spline.arcLength(t0);
        float a1 = spline.arcLength(t1);
        cout  << "Arc " << i << " t: " << t1 << " dt: " << t1 - t0 << " l: " << a1 << " d: " << (a1 - a0) << endl;
    }
    cout << endl;
    
    for(int i = 0; i <= 10; i++) {
        float t = 0.1 * float(i);
        WayPoint p(spline.evaluate(t));
        cout << "t: " << t << ", p: " << p << endl;
    }
    
    return 0;
}
