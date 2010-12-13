#include "util.h"

#include "gtest/gtest.h"

namespace {

// gtest is smart enough to compare doubles that are "almost" equal, but it
// doesn't work for arrays of doubles, so we provide this helper.  You'd think
// this should be a function, but gtest macros use "return", so that doesn't
// work.
#define EXPECT_MAT_EQ(expected, actual) \
   do { \
      for (int i = 0; i < 4; i++) { \
         for (int j = 0; j < 4; j++) { \
            EXPECT_EQ(expected[i][j], actual[i][j]); \
         } \
      } \
   } while (0)

/* print out debug info */
TEST(UtilTest, dumbTest)
{
	printf("Size of Vector3D = %lu\n", sizeof(Vector3D));
	printf("Size of Point3D = %lu\n", sizeof(Point3D));
	EXPECT_EQ(1,1);
}

/* sqrt test */
TEST(UtilTest, rsqrtss_sqrtTest)
{
	printf("sqrt(1000000) = %f\n", rsqrtss_sqrt(1000000));
	printf("sqrt(100) = %f\n", rsqrtss_sqrt(100));
	printf("sqrt(0) = %f\n", rsqrtss_sqrt(0));
	EXPECT_EQ(1, 1);
}

/* colour test */
TEST(UtilTest, colorClampTest)
{
	Colour c(0.5, 0.5, 0.5);
	Colour d(0.75, 0.75, 0.75);
	c += d;
	c.clamp();
	EXPECT_EQ(c.r, 1.);
	EXPECT_EQ(c.g, 1.);
	EXPECT_EQ(c.b, 1.);
}

/* dot product */
TEST(UtilTest, vecVecDotProduct)
{
	Vector3D a(0,1,2);
	Vector3D b(1,2,3);
	float result = 8.;
	float dot = a.dot(b);
//	printf("dot = %.6f\n", dot);
	EXPECT_EQ(dot, result);
}

TEST(UtilTest, vecPtDotProduct)
{
	Vector3D a(0,1,2);
	Point3D b(1,2,3);
	float result = 8.;
	float dot = b.dot(a);
//	printf("dot = %.6f\n", dot);
	EXPECT_EQ(dot, result);
}

TEST(UtilTest, ptPtDotProduct)
{
	Point3D a(0,1,2);
	Point3D b(1,2,3);
	float result = 8.;
	float dot = a.dot(b);
//	printf("dot = %.6f\n", dot);
	EXPECT_EQ(dot, result);
}

/* vec arithmetic */
TEST(UtilTest, vecAdd)
{
	Vector3D a(0,1,2);
	Vector3D b(1,2,3);
	Vector3D c = a+b;
	EXPECT_EQ(c.v[0], 1);
	EXPECT_EQ(c.v[1], 3);
	EXPECT_EQ(c.v[2], 5);
}

TEST(UtilTest, vecSub)
{
	Vector3D a(0,1,2);
	Vector3D b(1,2,3);
	Vector3D c = a-b;
	EXPECT_EQ(c.v[0], -1);
	EXPECT_EQ(c.v[1], -1);
	EXPECT_EQ(c.v[2], -1);
}

/* pt arithmetic */
TEST(UtilTest, ptAdd)
{
	Point3D a(0,1,2);
	Vector3D b(1,2,3);
	Point3D c = a+b;
	EXPECT_EQ(c.v[0], 1);
	EXPECT_EQ(c.v[1], 3);
	EXPECT_EQ(c.v[2], 5);
}

TEST(UtilTest, ptVecSub)
{
	Point3D a(0,1,2);
	Vector3D b(1,2,3);
	Point3D c = a-b;
	EXPECT_EQ(c.v[0], -1);
	EXPECT_EQ(c.v[1], -1);
	EXPECT_EQ(c.v[2], -1);
}

TEST(UtilTest, ptPtSub)
{
	Point3D a(0,1,2);
	Point3D b(1,2,3);
	Vector3D c = a-b;
	EXPECT_EQ(c.v[0], -1);
	EXPECT_EQ(c.v[1], -1);
	EXPECT_EQ(c.v[2], -1);
}

TEST(UtilTest, matrixMulIdentity) {
   Matrix4x4 A;
   Matrix4x4 B;
   Matrix4x4 C;
   initMatrix(A);
   initMatrix(B);
   // C = A * B
   mulMatrix(C, A, B);
   Matrix4x4 expected = {
      {1, 0, 0, 0},
      {0, 1, 0, 0},
      {0, 0, 1, 0},
      {0, 0, 0, 1},
   };
   EXPECT_MAT_EQ(expected, C);
}

// Tests that Foo does Xyz.
TEST(UtilTest, matrixMulSimple) {
   Matrix4x4 A = {
      {1, 2, 3, 4},
      {1, 2, 3, 4},
      {1, 2, 3, 4},
      {1, 2, 3, 4},
   };
   Matrix4x4 B;
   memcpy(B, A, sizeof(A));
   Matrix4x4 C;
   // C = A * B
   mulMatrix(C, A, B);
   Matrix4x4 expected = {
      {10, 20, 30, 40},
      {10, 20, 30, 40},
      {10, 20, 30, 40},
      {10, 20, 30, 40},
   };
   EXPECT_MAT_EQ(expected, C);
}

}  // namespace
