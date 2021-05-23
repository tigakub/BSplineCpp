//
//  WayPoint.cpp
//  Math
//
//  Created by Edward Janne on 5/22/21.
//

#include <math.h>
#include "WayPoint.hpp"

WayPoint::WayPoint()
: mSize(0), mArray(nullptr)
{ }

WayPoint::WayPoint(vector<float> &values)
: mSize(int(values.size())), mArray(mSize ? new float[mSize] : nullptr)
{
    int i = mSize;
    while(i--) mArray[i] = values[i];
}

WayPoint::WayPoint(float *values, int length)
: mSize(length), mArray(mSize ? new float[mSize] : nullptr)
{
    int i = mSize;
    while(i--) mArray[i] = values[i];
}

WayPoint::WayPoint(const WayPoint &p)
: mSize(p.mSize), mArray(mSize ? new float[mSize] : nullptr)
{
    int i = mSize;
    while(i--)  mArray[i] = p[i];
}

WayPoint::WayPoint(WayPoint &&p)
: mSize(p.mSize), mArray(p.mArray)
{
    p.mSize = 0;
    p.mArray = nullptr;
}

WayPoint::~WayPoint()
{
    if(mArray) delete [] mArray;
}

void WayPoint::addHocZero(int s)
{
    if(!s || mArray) return;
    mSize = s;
    mArray = new float[mSize];
    int i = mSize;
    while(i--) mArray[i] = 0.0;
}

float WayPoint::operator[](int i) const
{
    if(!mArray) return 0.0;
    if(i >= mSize) return 0.0;
    return mArray[i];
}

float &WayPoint::operator[](int i)
{
    return mArray[i];
}

int WayPoint::size() const
{
    return mSize;
}

WayPoint &WayPoint::operator+=(const WayPoint &p)
{
    addHocZero(p.mSize);
    int i = mSize;
    while(i--) mArray[i] += p[i];
    return *this;
}

WayPoint &WayPoint::operator-=(const WayPoint &p)
{
    addHocZero(p.mSize);
    int i = mSize;
    while(i--) mArray[i] -= p[i];
    return *this;
}

WayPoint &WayPoint::operator*=(const WayPoint &p)
{
    addHocZero(p.mSize);
    int i = mSize;
    while(i--) mArray[i] *= p[i];
    return *this;
}

WayPoint &WayPoint::operator/=(const WayPoint &p)
{
    addHocZero(p.mSize);
    int i = mSize;
    while(i--) mArray[i] /= p[i];
    return *this;
}

WayPoint WayPoint::operator+(const WayPoint &p) const
{
    WayPoint t(*this);
    t += p;
    return t;
}

WayPoint WayPoint::operator-(const WayPoint &p) const
{
    WayPoint t(*this);
    t -= p;
    return t;
}

WayPoint WayPoint::operator*(const WayPoint &p) const
{
    WayPoint t(*this);
    t *= p;
    return t;
}

WayPoint WayPoint::operator/(const WayPoint &p) const
{
    WayPoint t(*this);
    t /= p;
    return t;
}


WayPoint &WayPoint::operator+=(float s)
{
    if(!mArray) return *this;
    int i = mSize;
    while(i--) mArray[i] += s;
    return *this;
}

WayPoint &WayPoint::operator-=(float s)
{
    if(!mArray) return *this;
    int i = mSize;
    while(i--) mArray[i] -= s;
    return *this;
}

WayPoint &WayPoint::operator*=(float s)
{
    if(!mArray) return *this;
    int i = mSize;
    while(i--) mArray[i] *= s;
    return *this;
}

WayPoint &WayPoint::operator/=(float s)
{
    return (*this) *= (1.0 / s);
}

WayPoint WayPoint::operator+(float s) const
{
    WayPoint t(*this);
    t += s;
    return t;
}

WayPoint WayPoint::operator-(float s) const
{
    WayPoint t(*this);
    t -= s;
    return t;
}

WayPoint WayPoint::operator*(float s) const
{
    WayPoint t(*this);
    t *= s;
    return t;
}

WayPoint WayPoint::operator/(float s) const
{
    WayPoint t(*this);
    t *= (1.0 / s);
    return t;
}

float WayPoint::mag()
{
    if(!mArray) return 0.0;
    int i = mSize;
    float r = 0;
    while(i--) {
        r += mArray[i] * mArray[i];
    }
    return sqrt(r);
}

ostream &operator<<(ostream &os, const WayPoint &p)
{
    os << "(";
    for(int i = 0; i < p.size(); i++) {
        os << (i ? ", " : "") << p[i];
    }
    os << ")";
    return os;
}
