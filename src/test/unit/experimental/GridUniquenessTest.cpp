/* ---------------------------------------------------------------------
 * Numenta Platform for Intelligent Computing (NuPIC)
 * Copyright (C) 2018, Numenta, Inc.  Unless you have an agreement
 * with Numenta, Inc., for a separate license for this software code, the
 * following terms and conditions apply:
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero Public License version 3 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero Public License for more details.
 *
 * You should have received a copy of the GNU Affero Public License
 * along with this program.  If not, see http://www.gnu.org/licenses.
 *
 * http://numenta.org/licenses/
 * ----------------------------------------------------------------------
 */

/** @file
 * Unit tests for Grid Uniqueness
 */

#include <nupic/utils/Log.hpp>

#include <nupic/experimental/GridUniqueness.hpp>
#include "gtest/gtest.h"

#include <vector>
#include <cmath>

using namespace nupic;
using namespace nupic::experimental::grid_uniqueness;
using std::vector;
using std::pair;

namespace {
  TEST(GridUniquenessTest, ZeroInsideSquare)
  {
    const vector<double> scales = {2, 3, 6, 7, 21}; // Factors of 42
    vector<vector<vector<double>>> domainToPlaneByModule;
    vector<vector<vector<double>>> latticeBasisByModule;
    for (double scale : scales)
    {
      domainToPlaneByModule.push_back({
          {1/scale, 0},
          {0, 1/scale},
        });
      latticeBasisByModule.push_back({
          {1, 0},
          {0, 1},
        });
    }

    // Zero in center of square.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule,
                       {41.0, 41.0}, {2.0, 2.0}, 0.01));

    // Zero at bottom-left of square.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule,
                       {42.0, 42.0}, {2.0, 2.0}, 0.01));

    // Zero at top-right of square.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule,
                       {40.0, 40.0}, {2.0, 2.0}, 0.01));

    // Zero at top-left of square.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule,
                       {42.0, 40.0}, {2.0, 2.0}, 0.01));

    // Zero at bottom-right of square.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule,
                       {40.0, 42.0}, {2.0, 2.0}, 0.01));

    // Zero in bottom-left quadrant.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule,
                       {40.5, 40.5}, {2.0, 2.0}, 0.01));
  }


  TEST(GridUniquenessTest, ZeroOutsideSquare)
  {
    const vector<double> scales = {2, 3, 6, 7, 21}; // Factors of 42
    vector<vector<vector<double>>> domainToPlaneByModule;
    vector<vector<vector<double>>> latticeBasisByModule;
    for (double scale : scales)
    {
      domainToPlaneByModule.push_back({
          {1/scale, 0},
          {0, 1/scale},
        });
      latticeBasisByModule.push_back({
          {1, 0},
          {0, 1},
        });
    }

    // Zero out-of-bounds of square (bottom-left).
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {41.0, 41.0}, {0.5, 0.5}, 0.1));

    // Zero out-of-bounds of square (top-right).
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {42.5, 42.5}, {0.5, 0.5}, 0.1));

    // Zero out-of-bounds of square (bottom-right).
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {42.5, 41.0}, {0.5, 0.5}, 0.1));

    // Zero out-of-bounds of square (top-left).
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {41.0, 42.5}, {0.5, 0.5}, 0.1));
  }

  /**
   * Test rectangles just outside and just inside the radius of location zero.
   * Specifically focus on the area that would be covered by a square around
   * zero but not by a circle.
   */
  TEST(GridUniquenessTest, ZeroHasACircularRadius)
  {
    const vector<vector<vector<double>>> domainToPlaneByModule = {
      {{1, 0},
       {0, 1}}
    };
    const vector<vector<vector<double>>> latticeBasisByModule = {
      {{1, 0},
       {0, 1}}
    };

    const double phaseResolution = 0.2;

    double d = phaseResolution / (2 * sqrt(2));

    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {d + phaseResolution/100, d + phaseResolution/100},
                               {0.2, 0.2}, phaseResolution));

    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {d - phaseResolution/100, d - phaseResolution/100},
                               {0.2, 0.2}, phaseResolution));
  }

  TEST(GridUniquenessTest, 1DSanityCheck)
  {
    const vector<double> scales = {2, 3, 6, 7, 21}; // Factors of 42
    vector<vector<vector<double>>> domainToPlaneByModule;
    vector<vector<vector<double>>> latticeBasisByModule;
    for (double scale : scales)
    {
      // Add two modules at each scale.
      domainToPlaneByModule.push_back({
          {1/scale},
          {0},
        });
      latticeBasisByModule.push_back({
          {1, 0},
          {0, 1},
        });
      domainToPlaneByModule.push_back({
          {0},
          {1/scale},
        });
      latticeBasisByModule.push_back({
          {1, 0},
          {0, 1},
        });
    }

    // Zero in center of range.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {41.0}, {2.0}, 0.01));

    // Zero at left of range.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {42.0}, {2.0}, 0.01));

    // Zero at right of range.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {40.0}, {2.0}, 0.01));

    // Zero in left half.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {40.5}, {2.0}, 0.01));

    // Zero out-of-bounds of square (bottom-left).
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {41.0}, {0.5}, 0.1));

    // Zero out-of-bounds of square (top-right).
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {42.5}, {0.5}, 0.1));
  }

  TEST(GridUniquenessTest, NegativeNumbersTest)
  {
    const vector<double> scales = {2, 3, 6, 7, 21}; // Factors of 42
    vector<vector<vector<double>>> domainToPlaneByModule;
    vector<vector<vector<double>>> latticeBasisByModule;
    for (double scale : scales)
    {
      domainToPlaneByModule.push_back({
          {1/scale, 0},
          {0, 1/scale},
        });
      latticeBasisByModule.push_back({
          {1, 0},
          {0, 1},
        });
    }

    // Zero at a slightly negative number
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {-1.0, 41.0}, {0.999, 2.0}, 0.1));

    // Zero at a very negative number
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {-43.0, -43.0}, {2.0, 2.0}, 0.01));

    // Zero out-of-bounds of negative square
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {-1.0, 41.0}, {0.5, 2.0}, 0.1));

    // Zero out-of-bounds of negative square, both axes, bottom left
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {-43.0, -43.0}, {0.5, 0.5}, 0.1));

    // Zero out-of-bounds of negative square, both axes, top right
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {-41.5, -41.5}, {0.5, 0.5}, 0.1));
  }

  TEST(GridUniquenessTest, Project3DOnto2D)
  {
    // Choose scales that will have phase 0.5 at 42.
    const vector<double> scales = {4, 12, 28, 84};
    vector<vector<vector<double>>> domainToPlaneByModule;
    vector<vector<vector<double>>> latticeBasisByModule;
    for (double scale : scales)
    {
      // Shift the phase by 0.1 for every unit z.
      domainToPlaneByModule.push_back({
          {1/scale, 0, 0.1},
          {0, 1/scale, 0.1},
        });
      latticeBasisByModule.push_back({
          {1, 0},
          {0, 1},
        });
    }

    // Verify the z coordinate shifts the phase at 84.0 away from 0
    ASSERT_FALSE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {83.5, 83.5, 4.5}, {1.0, 1.0, 1.0}, 0.1));

    // Verify the z coordinate shifts the phase at 42.0 from 0.5 to 0.
    ASSERT_TRUE(
      findGridCodeZero(domainToPlaneByModule, latticeBasisByModule, {41.5, 41.5, 4.5}, {1.0, 1.0, 1.0}, 0.01));
  }

  vector<vector<double>> invert2DMatrix(const vector<vector<double>>& M)
  {
    const double detInv = 1 / (M[0][0]*M[1][1] - M[0][1]*M[1][0]);
    return {{detInv*M[1][1], -detInv*M[0][1]},
            {-detInv*M[1][0], detInv*M[0][0]}};
  }

  /**
   * Combine a square grid and hexagonal grid to place a common zero at a
   * particular point.
   *
   * Unrotated square and hexagonal lattices have the same lattice points along
   * the x axis. Rotate and scale this axis to place the first collision at any
   * arbitrary point.
   */
  vector<vector<vector<double>>> getBasisWithNearestZeroAt(double x, double y)
  {
    const double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan2(y, x);

    return {
      // Square
      invert2DMatrix(
        {{r*cos(theta), r*cos(theta + M_PI/4)},
         {r*sin(theta), r*sin(theta + M_PI/4)}}),
      // Hexagon
      invert2DMatrix(
        {{r*cos(theta), r*cos(theta + M_PI/3)},
         {r*sin(theta), r*sin(theta + M_PI/3)}}),
    };
  }

  TEST(GridUniquenessTest, ComputeGridUniquenessHypercubeTestPositive)
  {
    const vector<vector<vector<double>>> latticeBasisByModule = {
      {{1, 0},
       {0, 1}},
      {{1, 0},
       {0, 1}}
    };

    // Zero to the right of the ignored area.
    EXPECT_EQ(12,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(12.5, 0.25),
                      latticeBasisByModule,
                      0.01, 0.5).first));

    // Zero upper-right of the ignored area.
    EXPECT_EQ(6,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(6.5, 6.5),
                      latticeBasisByModule,
                      0.01, 0.5).first));

    // Zero above the ignored area.
    EXPECT_EQ(12,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(0.25, 12.5),
                      latticeBasisByModule,
                      0.01, 0.5).first));
  }

  TEST(GridUniquenessTest, ComputeGridUniquenessHypercubeTestNegative)
  {
    const vector<vector<vector<double>>> latticeBasisByModule = {
      {{1, 0},
       {0, 1}},
      {{1, 0},
       {0, 1}}
    };

    // Zero to the left of the ignored area.
    EXPECT_EQ(12,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(-12.5, 0.25),
                      latticeBasisByModule,
                      0.01, 0.5).first));

    // Zero upper-left of the ignored area.
    EXPECT_EQ(6,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(-6.5, 6.5),
                      latticeBasisByModule,
                      0.01, 0.5).first));

    // Zero above the ignored area.
    EXPECT_EQ(12,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(-0.25, 12.5),
                      latticeBasisByModule,
                      0.01, 0.5).first));
  }


  /**
   * Test 1: Upper right region
   * Test 2: Upper left region
   *
   * Test A: Test an area that's outside both the circle around zero
   * on the plane and the circle around zero on the torus.
   *
   * Test B: Test an area that's outside the circle around zero on the
   * plane but that intersects the circle around zero on the torus.
   *
   * Test C: Test an area that intersects the circle around zero on
   * the plane but is outside the circle around zero on the torus.
   *
   * Test D: Test an area that intersects the circle around zero on
   * the plane and the circle around zero on the torus.
   *
   * (Tests 1C and 2B are impossible.)
   */
  TEST(GridUniquenessTest, DistanceOnPlaneNotTorus)
  {
    const vector<vector<vector<double>>> domainToPlaneByModule = {
      {{1, 0},
       {0, 1}}
    };

    const vector<vector<vector<double>>> torusToPlaneByModule = {
      invert2DMatrix({{cos(0), cos(M_PI/3)},
                      {sin(0), sin(M_PI/3)}})
    };

    const double sideLength = 0.1;
    const double readoutResolution = 0.1;
    const double r = readoutResolution/2;

    //
    // Test 1
    //
    const double planeTheta1 = M_PI/6;
    const vector<double> planePoint1 = {r*cos(planeTheta1), r*sin(planeTheta1)};

    const double torusTheta1 = M_PI/4;
    const vector<double> torusPoint1 = {r*cos(torusTheta1), r*sin(torusTheta1)};
    const vector<double> torusPoint1OnPlane = {
      torusPoint1[0]*cos(0) + torusPoint1[1]*cos(M_PI/3),
      torusPoint1[0]*sin(0) + torusPoint1[1]*sin(M_PI/3),
    };

    // Test 1A
    const vector<double> corner1A = {
      torusPoint1OnPlane[0] + 0.01,
      torusPoint1OnPlane[1] + 0.01,
    };
    ASSERT_FALSE(findGridCodeZero(domainToPlaneByModule, torusToPlaneByModule, corner1A,
                                  {sideLength, sideLength}, readoutResolution));

    // Test 1B
    const vector<double> corner1B = {
      planePoint1[0] + (torusPoint1OnPlane[0] - planePoint1[0])/2,
      planePoint1[1] + (torusPoint1OnPlane[1] - planePoint1[1])/2
    };
    ASSERT_FALSE(findGridCodeZero(domainToPlaneByModule, torusToPlaneByModule, corner1B,
                                  {sideLength, sideLength}, readoutResolution));

    // Test 1D
    const vector<double> corner1D = {
      planePoint1[0] - 0.01,
      planePoint1[1] - 0.01,
    };
    ASSERT_TRUE(findGridCodeZero(domainToPlaneByModule, torusToPlaneByModule, corner1D,
                                 {sideLength, sideLength}, readoutResolution));

    //
    // Test 2
    //
    const double planeTheta2 = 2*M_PI/3;

    const vector<double> planePoint2 = {r*cos(planeTheta2), r*sin(planeTheta2)};

    const double torusTheta2 = 3*M_PI/4;
    const vector<double> torusPoint2 = {r*cos(torusTheta2), r*sin(torusTheta2)};
    const vector<double> torusPoint2OnPlane = {
      torusPoint2[0]*cos(0) + torusPoint2[1]*cos(M_PI/3),
      torusPoint2[0]*sin(0) + torusPoint2[1]*sin(M_PI/3),
    };

    // Test 2A
    const vector<double> corner2A = {
      planePoint2[0] - 0.01 - sideLength,
      planePoint2[1] + 0.01,
    };
    ASSERT_FALSE(findGridCodeZero(domainToPlaneByModule, torusToPlaneByModule, corner2A,
                                  {sideLength, sideLength}, readoutResolution));

    // Test 2C
    const vector<double> corner2C = {
      planePoint2[0] + (torusPoint2OnPlane[0] - planePoint2[0])/2 - sideLength,
      planePoint2[1] + (torusPoint2OnPlane[1] - planePoint2[1])/2
    };
    ASSERT_TRUE(findGridCodeZero(domainToPlaneByModule, torusToPlaneByModule, corner2C,
                                 {sideLength, sideLength}, readoutResolution));

    // Test 2D
    const vector<double> corner2D = {
      torusPoint2OnPlane[0] + 0.01 - sideLength,
      torusPoint2OnPlane[1] - 0.01,
    };
    ASSERT_TRUE(findGridCodeZero(domainToPlaneByModule, torusToPlaneByModule, corner2D,
                                 {sideLength, sideLength}, readoutResolution));
  }

  TEST(GridUniquenessTest, DeterministicDespiteMultithreading)
  {
    // This is a previously randomly-generated matrix that triggers multiple
    // threads to find zeros, some nearer, some further. Further away zeros
    // should not override nearby zeros, and when a thread finds a zero, it
    // shouldn't cancel threads that are searching for nearer zeros.
    const vector<vector<vector<double>>> domainToPlaneByModule = {
      {{0.4088715361390395, -0.9999112506968285, 0.8109731922797785, 0.25590203028822855},
       {-0.9125919498523434, -0.013322564689428938, 0.5850833115066139, 0.9667027210545974}},
      {{-0.704978485994098, -0.016909658985638815, 0.41508560377373277, 0.19893559514770887},
       {-0.05482092926492031, -0.7069045645863304, 0.5724543139323832, 0.6785459667430253}}};
    const vector<vector<vector<double>>> latticeBasisByModule = {
      {{1, 0},
       {0, 1}},
      {{1, 0},
       {0, 1}}
    };

    for (int i = 0; i < 100; i++)
    {
      ASSERT_EQ(0.5,
                computeGridUniquenessHypercube(domainToPlaneByModule,
                                               latticeBasisByModule,
                                               0.2, 0.5).first);
    }
  }
}
