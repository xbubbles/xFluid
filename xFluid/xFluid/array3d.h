#pragma once

#include <vector>
#include <math.h>
#include "gridindex.h"

/**
	Vector in 3 dimension space.

*/
template<class T>
struct Vec3 {
public:
	Vec3(){}

	Vec3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

	bool operator==(const Vec3 &other) const
	{
		return x == other.x && y == other.y && z == other.z;
	}

	bool operator!=(const Vec3 &other) const
	{
		return x != other.x || y != other.y || z != other.z;
	}

	Vec3<T> operator*(T value)const
	{
		return Vec3<T>(x*value, y*value, z*value);
	}

	T operator*(const Vec3<T>& vec)const
	{
		return x*vec.x + y*vec.y + z*vec.z;
	}

	Vec3<T> operator+(const Vec3<T>& vec)const
	{
		return Vec3<T>(x + vec.x, y + vec.y, z + vec.z);
	}

	Vec3<T> operator-(const Vec3<T>& vec)const
	{
		return Vec3<T>(x - vec.x, y - vec.y, z - vec.z);
	}

	Vec3<T> operator/(T value)const
	{
		return Vec3<T>(x / value, y / value, z / value);
	}

	Vec3<T> cross(const Vec3<T>& vec)const
	{
		return Vec3<T>(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x);
	}

	Vec3<T> normalize()const
	{
		double l = sqrt(x*x + y*y + z*z);
		return Vec3<T>(x / l, y / l, z / l);
	}

	double norm()const
	{
		return sqrt(x*x + y*y + z*z);
	}
public:
	T x, y, z;
};

typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;



/**
	3 dimensional array.
	Elemens can be indexed using grid index OR flat index.
	Data are saved in a std::vector.

	AUTHOR:				Xiao Shuisheng
	LAST MODIFICATION:	01/2019
*/
template<class T>
class Array3d
{
public:
	Array3d(int isize = 0, int jsize = 0, int ksize = 0) :_isize(isize), _jsize(jsize), _ksize(ksize), _vec(isize*jsize*ksize)
	{

	}


	/**
		resize.
		INPUT: size of each dimension. if size=0, clear all.
	*/
	void resize(int isize=0, int jsize=0, int ksize=0)
	{
		_isize = isize;
		_jsize = jsize;
		_ksize = ksize;

		_vec.resize(_isize*_jsize*_ksize);
	}


	/**
		Fill vector with assigned value.
	*/
	void fill(const T& value)
	{
		for (int i = 0; i < _vec.size(); ++i)
		{
			_vec[i] = value;
		}
	}

	
	/**
		Return whether an index is in range
	*/
	inline bool isInRange(int i, int j, int k)const
	{
		return (i>=0) && (j>=0) && (k>=0) && (i < _isize) && (j < _jsize) && (k < _ksize);
	}
	
	inline bool isInRange(int idx)const
	{
		return (idx>=0) && (idx < (_isize*_jsize*_ksize));
	}

	bool isInRange(const GridIndex3d& g)const
	{
		return isInRange(g.i, g.j, g.k);
	}

	/**
		convert 3 dimension index into flat index
	*/
	inline int index(int i, int j, int k)const
	{
		return k*(_isize*_jsize) + j*_isize + i;
	}

	inline int index(const GridIndex3d& g)const
	{
		return g.k*(_isize*_jsize) + g.j*_isize + g.i;
	}

	/**
		Get element.
		use operator () to get element.
	
	*/
	T& operator()(int i, int j, int k)
	{
		return _vec[index(i, j, k)];
	}

	T& operator()(const GridIndex3d& g)
	{
		return _vec[index(g)];
	}

	T& operator()(int flatidx)
	{
		return _vec[flatidx];
	}


	/**
		Get dimensions
	*/
	int sizei()const
	{
		return _isize;
	}
	int sizej()const
	{
		return _jsize;
	}
	int sizek()const
	{
		return _ksize;
	}
	Vec3i size()const
	{
		return Vec3i(_isize, _jsize, _ksize);
	}


private:
	int _isize, _jsize, _ksize;

	std::vector<T> _vec;
};
