//
//  ParametizingBSpline.hpp
//  Math
//
//  Created by Edward Janne on 5/23/21.
//

#ifndef ParametizingBSpline_hpp
#define ParametizingBSpline_hpp

#include <stdio.h>
#include "BSpline.hpp"

template <class T>
class ParametizingBSpline : public BSpline<T>
{
    protected:
        class MagPFunctor : public Functor {
            protected:
                ParametizingBSpline<T> &bSpline;
                
            public:
                MagPFunctor(ParametizingBSpline<T> &iSpline)
                : Functor(), bSpline(iSpline) { }
                
                virtual double operator()(double x) const {
                    return double(bSpline.p(float(x)).mag());
                }
        };
        
        class MagDFunctor : public Functor {
            protected:
                ParametizingBSpline<T> &bSpline;
                
            public:
                MagDFunctor(ParametizingBSpline<T> &iSpline)
                : Functor(), bSpline(iSpline) { }
                
                virtual double operator()(double x) const {
                    return double(bSpline.d(float(x)).mag());
                }
        };
        
        class SegmentArcLengthFunctor : public Functor {
            protected:
                ParametizingBSpline<T> &bSpline;
                int segment;
                
            public:
                SegmentArcLengthFunctor(ParametizingBSpline<T> &iSpline, int iSegment)
                : Functor(), bSpline(iSpline), segment(iSegment) { }
                
                virtual double operator()(double x) const {
                    return double(bSpline.segmentArcLength(segment, float(x)));
                }
        };
        
        class SegmentArcLengthPrimeFunctor : public Functor {
            protected:
                ParametizingBSpline<T> &bSpline;
                int segment;
                
            public:
                SegmentArcLengthPrimeFunctor(ParametizingBSpline<T> &iSpline, int iSegment)
                : Functor(), bSpline(iSpline), segment(iSegment) { }
                
                virtual double operator()(double x) const {
                    return double(bSpline.segmentArcLengthPrime(segment, float(x)));
                }
        };
        
    protected:
        float mLength;
        vector<float> spanLengths;
        MagPFunctor magP;
        MagDFunctor magD;
    
    public:
        ParametizingBSpline();
        ParametizingBSpline(vector<T> &iCP, int iOrder = 4);
        ParametizingBSpline(const BSpline<T> &iSpline);
        
        vector<float> parameterizeLinear(int iDivisions);
        vector<float> parameterizeSigmoid(int iDivisions);
        T d(float t);
        T evaluate(float t);
        float arcLength(float t);
        float timeForArc(float iArc);
        float segmentArcLength(int iSegment, float t);
        float segmentArcLengthPrime(int iSegment, float t);
        float timeForSegmentArc(int iSegment, float iArg);
        float length() const { return mLength ;}
        float cacheSpanLengths();
};

template <class T>
ParametizingBSpline<T>::ParametizingBSpline()
: BSpline<T>(), mLength(0.0), spanLengths(), magP(*this), magD(*this)
{ }

template <class T>
ParametizingBSpline<T>::ParametizingBSpline(vector<T> &iCP, int iOrder)
: BSpline<T>(iCP, iOrder), mLength(0.0), spanLengths(), magP(*this), magD(*this)
{
    cacheSpanLengths();
}

template <class T>
ParametizingBSpline<T>::ParametizingBSpline(const BSpline<T> &iSpline)
: BSpline<T>(iSpline), mLength(), spanLengths(), magP(*this), magD(*this)
{
    cacheSpanLengths();
}

template <class T>
vector<float> ParametizingBSpline<T>::parameterizeLinear(int iDivisions)
{
    vector<float> times;
    
    float lastKnot = BSpline<T>::knots.back();
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
vector<float> ParametizingBSpline<T>::parameterizeSigmoid(int iDivisions)
{
    vector<float> times;
    
    float lastKnot = BSpline<T>::knots.back();
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
T ParametizingBSpline<T>::d(float t)
{
    if(t < 0.0) t = 0.0;
    float max = BSpline<T>::knots.back();
    if(t > max) t = max;
    
    T r;
    int cpSize = int(BSpline<T>::cp.size());
    int n = BSpline<T>::order - 1;
    
    for(int i = 0; i < (cpSize - 1); i++) {
        float u0 = BSpline<T>::knots[i + n + 1];
        float u1 = BSpline<T>::knots[i + 1];
        float fn = (float(n) / (u0 - u1));
        T p((BSpline<T>::cp[i + 1] - BSpline<T>::cp[i]) * fn);
        float b = BSpline<T>::basis(i + 1, n - 1, t);
        r += p * b;
    }
    
    return r;
}

template <class T>
T ParametizingBSpline<T>::evaluate(float t)
{
    if(t < 0.0) t = 0.0;
    if(t > 1.0) t = 1.0;
    float n = timeForArc(t * mLength);
    return BSpline<T>::p(n);
}

template <class T>
float ParametizingBSpline<T>::arcLength(float t)
{
    int start = BSpline<T>::order - 1;
    float arcLen = 0.0;
    int i = start;
    
    while(t < BSpline<T>::knots[i]) {
        arcLen += spanLengths[i - start];
        i++;
    }
    
    float spanLength = legendreIntegrate(64, BSpline<T>::knots[i], t, magD);
    
    return arcLen + spanLength;
}

template <class T>
float ParametizingBSpline<T>::timeForArc(float iArc)
{
    if(iArc >= mLength) return BSpline<T>::knots[BSpline<T>::knots.size() - BSpline<T>::order];
    
    int i = 0;
    float len = 0.0;
    
    while((len + spanLengths[i]) < iArc) {
        len += spanLengths[i];
        i++;
    }
    
    return timeForSegmentArc(i, iArc - len);
}

template <class T>
float ParametizingBSpline<T>::segmentArcLength(int iSegment, float t)
{
    float t0 = BSpline<T>::knots[iSegment + BSpline<T>::order - 1];
    
    return legendreIntegrate(64, t0, t0 + t, magD);
}

template <class T>
float ParametizingBSpline<T>::segmentArcLengthPrime(int iSegment, float t)
{
    float t0 = BSpline<T>::knots[iSegment + BSpline<T>::order - 1];
    T deriv(d(t0 + t));
    
    return deriv.mag();
}

template <class T>
float ParametizingBSpline<T>::timeForSegmentArc(int iSegment, float iArc)
{
    int start = BSpline<T>::order - 1;
    float t0 = BSpline<T>::knots[iSegment + start];
    float hint = 0.0;
    
    return t0 + newtonRaphsonSolve(iArc, hint, SegmentArcLengthFunctor(*this, iSegment), SegmentArcLengthPrimeFunctor(*this, iSegment));
}

template <class T>
float ParametizingBSpline<T>::cacheSpanLengths()
{
    mLength = 0.0;
    int start = BSpline<T>::order - 1;
    int end = int(BSpline<T>::knots.size()) - BSpline<T>::order;
    
    spanLengths.clear();
    
    for(int i = start; i < end; i++) {
        float t0 = BSpline<T>::knots[i];
        float t1 = BSpline<T>::knots[i + 1];
        float spanLength = legendreIntegrate(64, t0, t1, magD);
        spanLengths.push_back(spanLength);
        mLength += spanLength;
    }
    
    return mLength;
}

#endif /* ParametizingBSpline_hpp */
