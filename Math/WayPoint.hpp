//
//  WayPoint.hpp
//  Math
//
//  Created by Edward Janne on 5/22/21.
//

#ifndef WayPoints_hpp
#define WayPoints_hpp

#include <iostream>
#include <stdio.h>
#include <vector>

using namespace std;

class WayPoint
{
    protected:
        int mSize;
        float *mArray;
        
        void addHocZero(int s);
        
    public:
        WayPoint();
        WayPoint(vector<float> &values);
        WayPoint(float *values, int length);
        WayPoint(const WayPoint &p);
        WayPoint(WayPoint &&p);
        virtual ~WayPoint();
        
        virtual float operator[](int i) const;
        virtual float &operator[](int i);
        virtual int size() const;
        
        virtual WayPoint &operator+=(const WayPoint &p);
        virtual WayPoint &operator-=(const WayPoint &p);
        virtual WayPoint &operator*=(const WayPoint &p);
        virtual WayPoint &operator/=(const WayPoint &p);
        
        virtual WayPoint operator+(const WayPoint &p) const;
        virtual WayPoint operator-(const WayPoint &p) const;
        virtual WayPoint operator*(const WayPoint &p) const;
        virtual WayPoint operator/(const WayPoint &p) const;
        
        virtual WayPoint &operator+=(float s);
        virtual WayPoint &operator-=(float s);
        virtual WayPoint &operator*=(float s);
        virtual WayPoint &operator/=(float s);
        
        virtual WayPoint operator+(float s) const;
        virtual WayPoint operator-(float s) const;
        virtual WayPoint operator*(float s) const;
        virtual WayPoint operator/(float s) const;
        
        virtual float mag();
};

ostream &operator<<(ostream &os, const WayPoint &p);

#endif /* WayPoints_hpp */
