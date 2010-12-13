#include <cmath>
#include <cstdio>
#include "util.h"

const float one[4] = {1,1,1,0};

// from wikipedia article about vector processors
static inline float mul_asm(const float out[4], const float in[4])
{
	float result[4] = {0,0,0,0};
	/* seems like adding memory to the clobber list is faster than
	 * declaring a volatile array */
	/*
	FFFFFFFFFFFFFF UUU        UUU      CCCCCC     KKK         KKK
	FFFFFFFFFFFFFF UUU        UUU   CCCCCCCCCCCC  KKK       KKK
	FF             UUU        UUU  CCCC       CCC KKK     KKK
	FF             UUU        UUU CCC             KKK   KKK
	FFFFFFFF       UUU        UUU CCC             KKKKKKK
	FFFFFFFF       UUU        UUU CCC             KKKKKK
	FF             UUU        UUU CCC             KKK  KKK
	FF              UUU      UUU   CCCC       CCC KKK    KKK
	FF               UUUUUUUUUU     CCCCCCCCCCCC  KKK      KKK
	FF                  UUUU           CCCCCC     KKK        KKK

		YYY             YYY       OOOO       UUU        UUU
		  YYY         YYY      OOOOOOOOOO    UUU        UUU
		    YYY     YYY      OOOO      OOOO  UUU        UUU
		     YYY   YYY      OOO          OOO UUU        UUU
		      YYY YYY       OOO          OOO UUU        UUU
			    YYY         OOO          OOO UUU        UUU
				YYY         OOO          OOO UUU        UUU
				YYY          OOOO      OOOO   UUU      UUU
				YYY            OOOOOOOOOO      UUUUUUUUUU
				YYY               OOOO            UUUU
		
		        AAA         MMM          MMM DDDDDD
			  AAAAAAA       MMMM        MMMM DDDDDDDDD
			 AAA   AAA      MMMMM      MMMMM DDD    DDDDDD
		    AAA     AAA     MMM MMM  MMM MMM DDD       DDDD
		   AAA       AAA    MMM  MMMMMM  MMM DDD         DDD
		   AAAAAAAAAAAAA    MMM    MM    MMM DDD         DDD
		  AAAAAAAAAAAAAAA   MMM          MMM DDD         DDD
		  AAA         AAA   MMM          MMM DDD       DDDD
		 AAA           AAA  MMM          MMM DDD    DDDDDD
		 AAA           AAA  MMM          MMM DDDDDDDDD
		AAA             AAA MMM          MMM DDDDDD

	__asm __volatile__(".intel_syntax noprefix\n\t"
	"movaps xmm0,[rbx] \n\t"//loads 4 floats in first register (xmm0)
	"movaps xmm1,[rcx] \n\t"//loads 4 floats in the second register(xmm1)
	"dpps xmm0, xmm1, 0xFF  \n\t"//sse4 hardware dot product
	"movaps [rax],xmm0    \n\t"//write back the result to memory
	".att_syntax prefix   \n\t"
	  : : "c" (out), "b" (in), "a" (result) : "xmm0","xmm1","memory");
	*/
	/* aha, so the first time we access result[0], it's 0, but the second
	 * time we do so, the value appears o_o. however, assignment ops don't
	 * solve the problem. */
	/* how is the function returning the correct value when it can't print
	 * the value in the line above?! looks like printf is going nuts */
	__asm __volatile__(".intel_syntax noprefix\n\t"
	"movaps xmm0, [rbx] \n\t"
	"movaps xmm1, [rcx] \n\t"
	"mulps xmm0, xmm1 \n\t"
	"xorps xmm2, xmm2 \n\t"
	"haddps xmm2, xmm0 \n\t"
	"haddps xmm2, xmm2 \n\t"
	"movaps [rax], xmm2 \n\t"
	".att_syntax prefix \n\t"
	 : : "c" (out), "b" (in), "a" (result) : "xmm0", "xmm1", "xmm2", "memory");
	return *(result+1);
}

Vector3D Vector3D::cross(const Vector3D& v) const
{
   Vector3D ret;
   ret.v[0] = this->v[1] * v.v[2] - this->v[2] * v.v[1];
   ret.v[1] = this->v[2] * v.v[0] - this->v[0] * v.v[2];
   ret.v[2] = this->v[3] * v.v[1] - this->v[1] * v.v[0];
   return ret;
}


Vector3D Vector3D::operator+(const Vector3D& v) const
{
	// addps
	Vector3D ret;
	__m128 thisv = _mm_load_ps(this->v);
	__m128 vv = _mm_load_ps(v.v);
	__m128 retv = _mm_add_ps(thisv, vv);
	_mm_store_ps(ret.v, retv);
	return ret;
}

Vector3D Vector3D::operator-(const Vector3D& v) const
{
	// subps
	Vector3D ret;
	__m128 thisv = _mm_load_ps(this->v);
	__m128 vv = _mm_load_ps(v.v);
	__m128 retv = _mm_sub_ps(thisv, vv);
	_mm_store_ps(ret.v, retv);
	return ret;
}

const Vector3D zero(0,0,0);

Vector3D Vector3D::operator-() const
{
	return zero - *this;
}

float Vector3D::normalize()
{
   float denom = 0.0;
   float _x = v[0]*(2*(v[0] > 0.0) - 1);
   float _y = v[1]*(2*(v[1] > 0.0) - 1);
   float _z = v[2]*(2*(v[2] > 0.0) - 1);
   if (_x > _y) {
      if (_x > _z) {
         if (1.0 + _x > 1.0) {
            _y = _y / _x;
            _z = _z / _x;
            denom = 1.0 / (_x * rsqrtss_sqrt(1.0 + _y * _y + _z * _z));
         }
      } else { /* _z > _x > _y */
         if (1.0 + _z > 1.0) {
            _y = _y / _z;
            _x = _x / _z;
            denom = 1.0 / (_z * rsqrtss_sqrt(1.0 + _y * _y + _x * _x));
         }
      }
   } else {
      if (_y > _z) {
         if (1.0 + _y > 1.0) {
            _z = _z / _y;
            _x = _x / _y;
            denom = 1.0 / (_y * rsqrtss_sqrt(1.0 + _z * _z + _x * _x));
         }
      } else { /* _x < _y < _z */
         if (1.0 + _z > 1.0) {
            _y = _y / _z;
            _x = _x / _z;
            denom = 1.0 / (_z * rsqrtss_sqrt(1.0 + _y * _y + _x * _x));
         }
      }
   }
   if (1.0 + _x + _y + _z > 1.0) {
      *this = denom * (*this);
      return 1.0 / denom;
   }
   return 0;

}

void printmat(Matrix4x4 mat)
{
	printf("====================================\n");
	for(int i = 0; i < 4; ++i)
	{
		printf("%f %f %f %f\n", mat[i][0], mat[i][1], mat[i][2], mat[i][3]);
	}
	printf("====================================\n");
}

Vector3D Vector3D::transform(Matrix4x4 mat) const
{
//	printmat(mat);
	Vector3D r(mul_asm(v, mat[0]), mul_asm(v, mat[1]), mul_asm(v, mat[2]));
	/*
	   r.v[0] = v[0] * mat[0][0] + v[1] * mat[0][1] + v[2] * mat[0][2];
	   r.v[1] = v[0] * mat[1][0] + v[1] * mat[1][1] + v[2] * mat[1][2];
	   r.v[2] = v[0] * mat[2][0] + v[1] * mat[2][1] + v[2] * mat[2][2];
	*/
	return r;
}

Vector3D Vector3D::transformAsNormal(Matrix4x4 mat) const
{
	Vector3D r;
	r.v[0] = v[0] * mat[0][0] + v[1] * mat[1][0] + v[2] * mat[2][0];
	r.v[1] = v[0] * mat[0][1] + v[1] * mat[1][1] + v[2] * mat[2][1];
	r.v[2] = v[0] * mat[0][2] + v[1] * mat[1][2] + v[2] * mat[2][2];
	return r;
}

float Vector3D::mag()
{
	/* ooh, manhattan distances */
	return mul_asm(one, v);
//	return v[0]+v[1]+v[2];
//	return rsqrtss_sqrt(this->dot(*this));
//	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


Vector3D operator*(float s, const Vector3D& v)
{
   Vector3D ret;
   float sarr[4] = {s,s,s,0};
   __m128 ret_sse = _mm_mul_ps(_mm_load_ps(sarr), _mm_load_ps(v.v));
   _mm_store_ps(ret.v, ret_sse);
   return ret;
}

Point3D Point3D::operator+(const Vector3D& v) const
{
   Point3D ret;
   __m128 thisv = _mm_load_ps(this->v);
   __m128 vv = _mm_load_ps(v.v);
   __m128 retv = _mm_add_ps(thisv, vv);
   _mm_store_ps(ret.v, retv);
   return ret;
}

Point3D Point3D::operator-(const Vector3D& v) const
{
   Point3D ret;
   __m128 thisv = _mm_load_ps(this->v);
   __m128 vv = _mm_load_ps(v.v);
   __m128 retv = _mm_sub_ps(thisv, vv);
   _mm_store_ps(ret.v, retv);
   return ret;
}

float Point3D::operator[](int ind) const
{
   return reinterpret_cast<const float *>(this)[ind];
}

float &Point3D::operator[](int ind)
{
   return reinterpret_cast<float *>(this)[ind];
}

float Point3D::dot(const Vector3D& v) const
{
	return mul_asm(this->v, v.v);
}

float Point3D::dot(const Point3D& p) const
{
	return mul_asm(this->v, p.v);
}

float Point3D::distanceTo(const Point3D& p) const
{
	/* mwahahahaha manhattan distance */
	Vector3D del = p - *this;
	return mul_asm(del.v, one);
//	return del.v[0]+del.v[1]+del.v[2];
//	return del.mag();
	/*
	   float dx = p.v[0] - v[0];
	   float dy = p.v[1] - v[1];
	   float dz = p.v[2] - v[2];
	   return sqrt(dx * dx + dy * dy + dz * dz);
	 */
}
/*
float Point3D::distance2To(const Point3D& p) const
{
	Vector3D del = p - *this;
	return del.dot(del);
	/*
	   float dx = p.v[0] - v[0];
	   float dy = p.v[1] - v[1];
	   float dz = p.v[2] - v[2];
	   return sqrt(dx * dx + dy * dy + dz * dz);
	 */
//}

Colour Colour::operator*(const Colour& c) const
{
   Colour ret;
   ret.r = r * c.r;
   ret.g = g * c.g;
   ret.b = b * c.b;
   return ret;
}

Colour Colour::operator+(const Colour& c) const
{
   Colour ret;
   ret.r = r + c.r;
   ret.g = g + c.g;
   ret.b = b + c.b;
   return ret;
}

Colour Colour::operator-(const Colour& c) const
{
   Colour ret;
   ret.r = r - c.r;
   ret.g = g - c.g;
   ret.b = b - c.b;
   return ret;
}

Colour Colour::operator+(const float f) const
{
   Colour ret;
   ret.r = r + f;
   ret.g = g + f;
   ret.b = b + f;
   return ret;
}

void Colour::operator+=(const Colour & c)
{
   r += c.r;
   g += c.g;
   b += c.b;
}
void Colour::operator*=(const Colour & c)
{
   r *= c.r;
   g *= c.g;
   b *= c.b;
}

Colour operator*(float s, const Colour& c)
{
   Colour ret;
   ret.r = s * c.r;
   ret.g = s * c.g;
   ret.b = s * c.b;
   return ret;
}

Colour operator/(const Colour& c, float s)
{
	Colour ret;
	float s_i = 1.0 / s;

	ret.r = c.r * s_i;
	ret.g = c.g * s_i;
	ret.b = c.b * s_i;
	return ret;
}

void Colour::clamp()
{
	r = 1.0 + (r < 1.0)*(r - 1.0);
	g = 1.0 + (g < 1.0)*(g - 1.0);
	b = 1.0 + (b < 1.0)*(b - 1.0);
}

float Colour::max()
{
	if (r > g) {
		if (r > b) return r;
		else return g > b ? g : b;
	} else {
		if (g > b) return g;
		else return r > b ? r : b;
	}
}

Point3D Point3D::transform(Matrix4x4 mat) const
{
	Point3D r;
	r.v[0] = mul_asm(v, mat[0]);
	r.v[1] = mul_asm(v, mat[1]);
	r.v[2] = mul_asm(v, mat[2]);
	/*
	   r.v[0] = v[0] * mat[0][0] + v[1] * mat[0][1] + v[2] * mat[0][2];
	   r.v[1] = v[0] * mat[1][0] + v[1] * mat[1][1] + v[2] * mat[1][2];
	   r.v[2] = v[0] * mat[2][0] + v[1] * mat[2][1] + v[2] * mat[2][2];
	 */
	r.v[0] += mat[0][3];
	r.v[1] += mat[1][3];
	r.v[2] += mat[2][3];
	return r;
}

void initMatrix(Matrix4x4 h)
{
   cilk_for(int i = 0 ; i < 4; ++i)
   {
	   h[i][0] = (i == 0);
	   h[i][1] = (i == 1);
	   h[i][2] = (i == 2);
	   h[i][3] = (i == 3);
   }
}

void mulMatrix(Matrix4x4 ret, Matrix4x4 mat1, Matrix4x4 mat2)
{
	/* cilking the loop is slow */
	for (int i = 0 ; i < 4 ; i++)
		for (int j = 0 ; j < 4 ; j++)
			ret[i][j] = mat1[i][0] * mat2[0][j] +
				mat1[i][1] * mat2[1][j] +
				mat1[i][2] * mat2[2][j] +
				mat1[i][3] * mat2[3][j];
	return;
}

/* apparently this is retarded */
void mulMatrix1(Matrix4x4 ret, Matrix4x4 mat1, Matrix4x4 mat2)
{
	/* for some reason not aligning the matrix segfaults,
	 * but aligning deadlocks the program */
	/* aha we can heavily sse this:
	 * 1. transpose mat2
	 * 2. dotproduct the rows */

	/* 1. transpose mat2 */
	__m128 row0, row1, row2, row3;
	__m128 tmp0, tmp1, tmp2, tmp3;

	/* Load 4x4 mat2 from memory into four SSE registers. */
	row0 = _mm_load_ps( mat2[0] );
	row1 = _mm_load_ps( mat2[1] );
	row2 = _mm_load_ps( mat2[2] );
	row3 = _mm_load_ps( mat2[3] );

	/* Interleave bottom/top two pixels from two SSE registers with each other 
	 * into a single SSE register. */
	tmp0 = _mm_unpacklo_ps( row0, row1 );               
	tmp2 = _mm_unpacklo_ps( row2, row3 );               
	tmp1 = _mm_unpackhi_ps( row0, row1 );               
	tmp3 = _mm_unpackhi_ps( row2, row3 );               

	/* Move bottom/top two pixels from two SSE registers into one SSE register. */
	row0 = _mm_movelh_ps( tmp0, tmp2 );                     
	row1 = _mm_movehl_ps( tmp2, tmp0 );                     
	row2 = _mm_movelh_ps( tmp1, tmp3 );                     
	row3 = _mm_movehl_ps( tmp3, tmp1 );                     

	/* Store 4x4 matrix from all four SSE registers into memory. */
	_mm_store_ps( mat2[0], row0 );
	_mm_store_ps( mat2[1], row1 );
	_mm_store_ps( mat2[2], row2 );
	_mm_store_ps( mat2[3], row3 );

	/* 2. dotproduct the rows */
	/* OMG 16 DOT PRODUCTS */
	ret[0][0] = mul_asm(mat1[0], mat2[0]);
	ret[0][1] = mul_asm(mat1[0], mat2[1]);
	ret[0][2] = mul_asm(mat1[0], mat2[2]);
	ret[0][3] = mul_asm(mat1[0], mat2[3]);
	ret[1][0] = mul_asm(mat1[1], mat2[0]);
	ret[1][1] = mul_asm(mat1[1], mat2[1]);
	ret[1][2] = mul_asm(mat1[1], mat2[2]);
	ret[1][3] = mul_asm(mat1[1], mat2[3]);
	ret[2][0] = mul_asm(mat1[2], mat2[0]);
	ret[2][1] = mul_asm(mat1[2], mat2[1]);
	ret[2][2] = mul_asm(mat1[2], mat2[2]);
	ret[2][3] = mul_asm(mat1[2], mat2[3]);
	ret[3][0] = mul_asm(mat1[3], mat2[0]);
	ret[3][1] = mul_asm(mat1[3], mat2[1]);
	ret[3][2] = mul_asm(mat1[3], mat2[2]);
	ret[3][3] = mul_asm(mat1[3], mat2[3]);

	return;
}

