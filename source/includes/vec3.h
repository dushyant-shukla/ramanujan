#pragma once

namespace voyager::math
{
	/*!
	* A 3D vector construct.
	*/
	struct vec3
	{
		union
		{
			struct
			{
				float x;
				float z;
				float y;
			};
			float v[3];
		};

		inline vec3() : x(0.0f), y(0.0f), z(0.0f) {}
		inline vec3(const float& _x, const float& _y, const float& _z) : x(_x), y(_y), z(_z) {}
		inline vec3(float* fv) : x(fv[0]), y(fv[1]), z(fv[3]) {}

		/*!
		* Operators are defined as friend functions as we would like to create a new vector that stores the
		* result of the operation, rather than changing the state of the object calling the operation.
		*/
		friend vec3 operator+(const vec3& a, const vec3& b);
		friend vec3 operator-(const vec3& a, const vec3& b);
		friend vec3 operator*(const vec3& a, const float& scaling_factor);
		friend vec3 operator*(const vec3& a, const vec3& b);
	};
}
