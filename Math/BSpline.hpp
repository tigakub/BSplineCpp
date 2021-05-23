//
//  BSpline.hpp
//  Math
//
//  Created by Edward Janne on 5/22/21.
//

#ifndef BSpline_hpp
#define BSpline_hpp

#include <stdio.h>
#include <vector>

#include "Functor.hpp"
#include "Newton.hpp"
#include "Legendre.hpp"

using namespace std;

template <class T>
class BSpline
{
    protected:
        vector<T> cp;
        vector<float> knots;
        int order;
        int degree;
        
    public:
        BSpline();
        BSpline(vector<T> &iCP, int iOrder = 4);
        BSpline(const BSpline<T> &iSpline);
        virtual ~BSpline();
        
        T p(float t);
        
    protected:
        float basis(int i, int k, float t);
        void calcKnots();
};

template <class T>
BSpline<T>::BSpline()
: cp(), knots(), order(0), degree(0)
{ }

template <class T>
BSpline<T>::BSpline(vector<T> &iCP, int iOrder)
: cp(iCP), knots(), order(iOrder), degree(iOrder - 1)
{
    calcKnots();
}

template <class T>
BSpline<T>::BSpline(const BSpline<T> &iSpline)
: cp(iSpline.cp), knots(iSpline.knots), order(iSpline.order), degree(iSpline.degree)
{ }

template <class T>
BSpline<T>::~BSpline()
{ }

template <class T>
T BSpline<T>::p(float t)
{
    if(t < 0.0) t = 0.0;
    float max = knots.back();
    if(t > max) t = max;
    
    T r;
    int cpSize = int(cp.size());
    int n = order - 1;
    
    for(int i = 0; i < cpSize; i++) {
        r += cp[i] * basis(i, n, t);
    }
    
    return r;
}

template <class T>
float BSpline<T>::basis(int i, int k, float t)
{
    float r = 0.0;
    switch(k) {
        case 0:
            if((knots[i] <= t) && (t <= knots[i+1]))
                r = 1.0;
            else
                r = 0.0;
            break;
        default:
            float n0 = t - knots[i];
            float d0 = knots[i + k] - knots[i];
            float b0 = basis(i, k - 1, t);
            
            float n1 = knots[i + k + 1] - t;
            float d1 = knots[i + k + 1] - knots[i + 1];
            float b1 = basis(i + 1, k - 1, t);
            
            float left = 0.0;
            float right = 0.0;
            
            if((b0 != 0.0) && (d0 != 0.0))
                left = n0 * b0 / d0;
            if((b1 != 0.0) && (d1 != 0.0))
                right = n1 * b1 / d1;
                
            r = left + right;
            break;
    }
    return r;
}

template <class T>
void BSpline<T>::calcKnots()
{
    int cpSize = int(cp.size());
    int knotCount = cpSize + order;
    
    knots.clear();
    
    int i = order;
    while(i--) {
        knots.push_back(0.0);
    }
    
    float k = 0.0;
    for(i = order; i < cpSize; i++) {
        k += 1.0;
        knots.push_back(float(k));
    }
    
    k += 1.0;
    for(i = cpSize; i < knotCount; i++) {
        knots.push_back(float(k));
    }
}

#endif /* BSpline_hpp */
