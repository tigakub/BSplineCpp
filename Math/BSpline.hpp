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
        typedef float (BSpline<T>::*callback)(float);
        
        class MagPFunctor : public Functor {
            protected:
                BSpline<T> &bSpline;
                
            public:
                MagPFunctor(BSpline<T> &iSpline)
                : Functor(), bSpline(iSpline) { }
                
                virtual double operator()(double x) const {
                    return double(bSpline.p(float(x)).mag());
                }
        };
        
        class MagDFunctor : public Functor {
            protected:
                BSpline<T> &bSpline;
                
            public:
                MagDFunctor(BSpline<T> &iSpline)
                : Functor(), bSpline(iSpline) { }
                
                virtual double operator()(double x) const {
                    return double(bSpline.d(float(x)).mag());
                }
        };
        
        class SegmentArcLengthFunctor : public Functor {
            protected:
                BSpline<T> &bSpline;
                int segment;
                
            public:
                SegmentArcLengthFunctor(BSpline<T> &iSpline, int iSegment)
                : Functor(), bSpline(iSpline), segment(iSegment) { }
                
                virtual double operator()(double x) const {
                    return double(bSpline.segmentArcLength(segment, float(x)));
                }
        };
        
        class SegmentArcLengthPrimeFunctor : public Functor {
            protected:
                BSpline<T> &bSpline;
                int segment;
                
            public:
                SegmentArcLengthPrimeFunctor(BSpline<T> &iSpline, int iSegment)
                : Functor(), bSpline(iSpline), segment(iSegment) { }
                
                virtual double operator()(double x) const {
                    return double(bSpline.segmentArcLengthPrime(segment, float(x)));
                }
        };
        
    protected:
        float mLength;
        vector<float> spanLengths;
        vector<T> cp;
        vector<float> knots;
        int order;
        int degree;
        MagPFunctor magP;
        MagDFunctor magD;
        
    public:
        BSpline();
        BSpline(vector<T> &iCP, int iOrder = 4);
        BSpline(const BSpline<T> &iSpline);
        virtual ~BSpline();
        
        vector<float> parameterizeLinear(int iDivisions);
        vector<float> parameterizeSigmoid(int iDivisions);
        T p(float t);
        T d(float t);
        T evaluate(float t);
        float arcLength(float t);
        float timeForArc(float iArc);
        float segmentArcLength(int iSegment, float t);
        float segmentArcLengthPrime(int iSegment, float t);
        float timeForSegmentArc(int iSegment, float iArg);
        float length() const { return mLength ;}
        
    protected:
        float basis(int i, int k, float t);
        void calcKnots();
        float cacheSpanLengths();
};

template <class T>
BSpline<T>::BSpline()
: mLength(0.0), spanLengths(), cp(), knots(), order(0), degree(0), magP(*this), magD(*this)
{ }

template <class T>
BSpline<T>::BSpline(vector<T> &iCP, int iOrder)
: mLength(0.0), spanLengths(), cp(iCP), knots(), order(iOrder), degree(iOrder - 1), magP(*this), magD(*this)
{
    calcKnots();
    cacheSpanLengths();
}

template <class T>
BSpline<T>::BSpline(const BSpline<T> &iSpline)
: mLength(iSpline.mLength), spanLengths(iSpline.spanLengths), cp(iSpline.cp), knots(iSpline.knots), order(iSpline.order), degree(iSpline.degree), magP(*this), magD(*this)
{ }

template <class T>
BSpline<T>::~BSpline()
{ }

template <class T>
vector<float> BSpline<T>::parameterizeLinear(int iDivisions)
{
    vector<float> times;
    
    float lastKnot = knots.back();
    times.push_back(0.0);
    float step = mLength / float(iDivisions);
    float arc = step;
    for(int i = 1; i < iDivisions; i++) {
        times.push_back(float(timeForArc(arc)));
        arc += step;
    }
    times.push_back(lastKnot);
    
    return times;
}

template <class T>
vector<float> BSpline<T>::parameterizeSigmoid(int iDivisions)
{
    vector<float> times;
    
    float lastKnot = knots.back();
    times.push_back(0.0);
    float step = 1.0 / float(iDivisions);
    for(float i = step; i < 1.0; i += step) {
        float arc = mLength / (1 + exp(-14.0 * (i - 0.5)));
        times.push_back(float(timeForArc(arc)));
    }
    times.push_back(lastKnot);
    
    return times;
}

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
T BSpline<T>::d(float t)
{
    if(t < 0.0) t = 0.0;
    float max = knots.back();
    if(t > max) t = max;
    
    T r;
    int cpSize = int(cp.size());
    int n = order - 1;
    
    for(int i = 0; i < (cpSize - 1); i++) {
        float u0 = knots[i + n + 1];
        float u1 = knots[i + 1];
        float fn = (float(n) / (u0 - u1));
        T p((cp[i + 1] - cp[i]) * fn);
        float b = basis(i + 1, n - 1, t);
        r += p * b;
    }
    
    return r;
}

template <class T>
T BSpline<T>::evaluate(float t)
{
    if(t < 0.0) t = 0.0;
    if(t > 1.0) t = 1.0;
    float n = timeForArc(t * mLength);
    return p(n);
}

template <class T>
float BSpline<T>::arcLength(float t)
{
    int start = order - 1;
    float arcLen = 0.0;
    int i = start;
    
    while(t < knots[i]) {
        arcLen += spanLengths[i - start];
        i++;
    }
    
    float spanLength = legendreIntegrate(64, knots[i], t, magD);
    
    return arcLen + spanLength;
}

template <class T>
float BSpline<T>::timeForArc(float iArc)
{
    if(iArc >= mLength) return knots[knots.size() - order];
    
    int i = 0;
    float len = 0.0;
    
    while((len + spanLengths[i]) < iArc) {
        len += spanLengths[i];
        i++;
    }
    
    return timeForSegmentArc(i, iArc - len);
}

template <class T>
float BSpline<T>::segmentArcLength(int iSegment, float t)
{
    float t0 = knots[iSegment + order - 1];
    
    return legendreIntegrate(64, t0, t0 + t, magD);
}

template <class T>
float BSpline<T>::segmentArcLengthPrime(int iSegment, float t)
{
    float t0 = knots[iSegment + order - 1];
    T deriv(d(t0 + t));
    
    return deriv.mag();
}

template <class T>
float BSpline<T>::timeForSegmentArc(int iSegment, float iArc)
{
    int start = order - 1;
    float t0 = knots[iSegment + start];
    float hint = 0.0;
    
    return t0 + newtonRaphsonSolve(iArc, hint, SegmentArcLengthFunctor(*this, iSegment), SegmentArcLengthPrimeFunctor(*this, iSegment));
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

template <class T>
float BSpline<T>::cacheSpanLengths()
{
    mLength = 0.0;
    int start = order - 1;
    int end = int(knots.size()) - order;
    
    spanLengths.clear();
    
    for(int i = start; i < end; i++) {
        float t0 = knots[i];
        float t1 = knots[i + 1];
        float spanLength = legendreIntegrate(64, t0, t1, magD);
        spanLengths.push_back(spanLength);
        mLength += spanLength;
    }
    
    return mLength;
}

#endif /* BSpline_hpp */
