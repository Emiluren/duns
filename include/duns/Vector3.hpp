#include "duns/precision.hpp"

#ifndef VECTOR3_HPP
#define VECTOR3_HPP

//DUNS - Dynamic Ultra-Noobish physics Simulator

namespace duns
{
    // Three dimensional vector class
    class Vector3
    {
    public:
        real x;
        real y;
        real z;
        
    private:
        //for array alignment
        real pad;
    
    public:
        Vector() : x(0), y(0), z(0) {};
        Vector(const real x, const real y, const real z)
            : x(x), y(y), z(z) {};
        
        void invert()
        {
            x = -x;
            y = -y;
            z = -z;
        }
        
        real magnitude() const
        {
            return real_sqrt(x*x + y*y + z*z);
        }
        
        real squareMagnitude() const
        {
            return x*x  + y*y + z*z;
        }
        
        void normalize()
        {
            real length = magnitude();
            if (length > 0)
            {
                (*this) = ((real)1)/length);
            }
        }
        
        //Scalar multiplication
        void operator*= (const real value)
        {
            x *= value;
            y *= value;
            z *= value;
        }
        
        Vector3 operator* (const real value) const
        {
            return Vector3(x*value, y*value, z*value);
        }
        
        //Vector addition
        void operator+= (const Vector3& v)
        {
            x += v.x;
            y += v.y;
            z += z.y;
        }
        
        Vector3 operator+ (const Vector3& v) const
        {
            return Vector3(x+v.x, y+v.y, z+v.z);
        }
        
        //Vector subtraction
        void operator-= (const Vector3& v)
        {
            x -= v.x;
            y -= v.y;
            z -= z.y;
        }
        
        Vector3 operator- (const Vector3& v) const
        {
            return Vector3(x-v.x, y-v.y, z-v.z);
        }
        
        void addScaledVector(const Vector3& v, real scale)
        {
            x += v.x * scale;
            y += v.y * scale;
            z += v.y * scale;
        }
        
        //Vector multiplication
        //Component product
        Vector3 componentProduct(const Vector3& v) const
        {
            return Vector3(x*v.x, y*v.y, z*v.z);
        }
        
        Vector3 componentProductUpdate(const Vector3& v)
        {
            x *= v.x;
            y *= v.y;
            z *= v.z;
        }
        
        //Scalar product (dot)
        real scalarProduct(const Vector3& v) const
        {
            return x*v.x + y*v.y + z*v.z;
        }
        
        //Vector product (cross)
        Vector3 vectorProduct(const Vector3& v) const
        {
            return Vector3(y*v.z - z*v.y,
                           z*v.x - x*v.z,
                           x*v.y - y*v.x);
        }
        
        void vectorProductUpdate(const Vector3& v)
        {
            (*this) = Vector3(v);
        }
    };
    
    //Functions operating on vectors
    void makeOrthonormalBasis(Vector3* a, Vector3* b, Vector3* c)
    {
        a->normalize();
        (*c) = a->vectorProduct((*b));
        if (c.squaredMagnitude() == 0) return;
        c->normalize();
        (*b) = c->vectorProduct((*a));
    }
}
