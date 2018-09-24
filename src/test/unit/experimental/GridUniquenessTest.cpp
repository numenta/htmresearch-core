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
    vector<vector<vector<double>>> A;
    for (double scale : scales)
    {
      A.push_back({
          {1/scale, 0},
          {0, 1/scale},
        });
    }

    vector<double> displacement(2);

    // Zero in center of square.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {41.0, 41.0}, {2.0, 2.0}, 0.01, displacement));

    // Zero at bottom-left of square.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {42.0, 42.0}, {2.0, 2.0}, 0.01, displacement));

    // Zero at top-right of square.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {40.0, 40.0}, {2.0, 2.0}, 0.01, displacement));

    // Zero at top-left of square.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {42.0, 40.0}, {2.0, 2.0}, 0.01, displacement));

    // Zero at bottom-right of square.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {40.0, 42.0}, {2.0, 2.0}, 0.01, displacement));

    // Zero in bottom-left quadrant.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {40.5, 40.5}, {2.0, 2.0}, 0.01, displacement));
  }


  TEST(GridUniquenessTest, ZeroOutsideSquare)
  {
    const vector<double> scales = {2, 3, 6, 7, 21}; // Factors of 42
    vector<vector<vector<double>>> A;
    for (double scale : scales)
    {
      A.push_back({
          {1/scale, 0},
          {0, 1/scale},
        });
    }

    vector<double> displacement(2);

    // Zero out-of-bounds of square (bottom-left).
    ASSERT_FALSE(
      findGridDisplacementZero(A, {41.0, 41.0}, {0.5, 0.5}, 0.1, displacement));

    // Zero out-of-bounds of square (top-right).
    ASSERT_FALSE(
      findGridDisplacementZero(A, {42.5, 42.5}, {0.5, 0.5}, 0.1, displacement));

    // Zero out-of-bounds of square (bottom-right).
    ASSERT_FALSE(
      findGridDisplacementZero(A, {42.5, 41.0}, {0.5, 0.5}, 0.1, displacement));

    // Zero out-of-bounds of square (top-left).
    ASSERT_FALSE(
      findGridDisplacementZero(A, {41.0, 42.5}, {0.5, 0.5}, 0.1, displacement));
  }

  /**
   * Test rectangles just outside and just inside the radius of location zero.
   * Specifically focus on the area that would be covered by a square around
   * zero but not by a circle.
   */
  TEST(GridUniquenessTest, ZeroHasACircularRadius)
  {
    const vector<vector<vector<double>>> A = {
      {{1, 0},
       {0, 1}}
    };

    const double phaseResolution = 0.2;

    double d = phaseResolution / (2 * sqrt(2));

    vector<double> displacement(2);
    ASSERT_FALSE(
      findGridDisplacementZero(A, {d + phaseResolution/100, d + phaseResolution/100},
                               {0.2, 0.2}, phaseResolution, displacement));

    ASSERT_TRUE(
      findGridDisplacementZero(A, {d - phaseResolution/100, d - phaseResolution/100},
                               {0.2, 0.2}, phaseResolution, displacement));
  }

  TEST(GridUniquenessTest, 1DSanityCheck)
  {
    const vector<double> scales = {2, 3, 6, 7, 21}; // Factors of 42
    vector<vector<vector<double>>> A;
    for (double scale : scales)
    {
      // Add two modules at each scale.
      A.push_back({
          {1/scale},
          {0},
        });
      A.push_back({
          {0},
          {1/scale},
        });
    }

    vector<double> displacement(1);

    // Zero in center of range.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {41.0}, {2.0}, 0.01, displacement));

    // Zero at left of range.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {42.0}, {2.0}, 0.01, displacement));

    // Zero at right of range.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {40.0}, {2.0}, 0.01, displacement));

    // Zero in left half.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {40.5}, {2.0}, 0.01, displacement));

    // Zero out-of-bounds of square (bottom-left).
    ASSERT_FALSE(
      findGridDisplacementZero(A, {41.0}, {0.5}, 0.1, displacement));

    // Zero out-of-bounds of square (top-right).
    ASSERT_FALSE(
      findGridDisplacementZero(A, {42.5}, {0.5}, 0.1, displacement));
  }

  TEST(GridUniquenessTest, NegativeNumbersTest)
  {
    const vector<double> scales = {2, 3, 6, 7, 21}; // Factors of 42
    vector<vector<vector<double>>> A;
    for (double scale : scales)
    {
      A.push_back({
          {1/scale, 0},
          {0, 1/scale},
        });
    }

    vector<double> displacement(2);
    // Zero at a slightly negative number
    ASSERT_TRUE(
      findGridDisplacementZero(A, {-1.0, 41.0}, {0.999, 2.0}, 0.1, displacement));

    // Zero at a very negative number
    ASSERT_TRUE(
      findGridDisplacementZero(A, {-43.0, -43.0}, {2.0, 2.0}, 0.01, displacement));

    // Zero out-of-bounds of negative square
    ASSERT_FALSE(
      findGridDisplacementZero(A, {-1.0, 41.0}, {0.5, 2.0}, 0.1, displacement));

    // Zero out-of-bounds of negative square, both axes, bottom left
    ASSERT_FALSE(
      findGridDisplacementZero(A, {-43.0, -43.0}, {0.5, 0.5}, 0.1, displacement));

    // Zero out-of-bounds of negative square, both axes, top right
    ASSERT_FALSE(
      findGridDisplacementZero(A, {-41.5, -41.5}, {0.5, 0.5}, 0.1, displacement));
  }

  TEST(GridUniquenessTest, Project3DOnto2D)
  {
    // Choose scales that will have phase 0.5 at 42.
    const vector<double> scales = {4, 12, 28, 84};
    vector<vector<vector<double>>> A;
    for (double scale : scales)
    {
      // Shift the phase by 0.1 for every unit z.
      A.push_back({
          {1/scale, 0, 0.1},
          {0, 1/scale, 0.1},
        });
    }

    vector<double> displacement(3);

    // Verify the z coordinate shifts the phase at 84.0 away from 0
    ASSERT_FALSE(
      findGridDisplacementZero(A, {83.5, 83.5, 4.5}, {1.0, 1.0, 1.0}, 0.1, displacement));

    // Verify the z coordinate shifts the phase at 42.0 from 0.5 to 0.
    ASSERT_TRUE(
      findGridDisplacementZero(A, {41.5, 41.5, 4.5}, {1.0, 1.0, 1.0}, 0.01, displacement));
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
    // Zero to the right of the ignored area.
    EXPECT_EQ(12,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(12.5, 0.25),
                      0.01, 0.5).first));

    // Zero upper-right of the ignored area.
    EXPECT_EQ(6,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(6.5, 6.5),
                      0.01, 0.5).first));

    // Zero above the ignored area.
    EXPECT_EQ(12,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(0.25, 12.5),
                      0.01, 0.5).first));
  }

  TEST(GridUniquenessTest, ComputeGridUniquenessHypercubeTestNegative)
  {
    // Zero to the left of the ignored area.
    EXPECT_EQ(12,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(-12.5, 0.25),
                      0.01, 0.5).first));

    // Zero upper-left of the ignored area.
    EXPECT_EQ(6,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(-6.5, 6.5),
                      0.01, 0.5).first));

    // Zero above the ignored area.
    EXPECT_EQ(12,
              floor(computeGridUniquenessHypercube(
                      getBasisWithNearestZeroAt(-0.25, 12.5),
                      0.01, 0.5).first));
  }

  TEST(GridUniquenessTest, DeterministicDespiteMultithreading)
  {
    // This is a previously randomly-generated matrix that triggers multiple
    // threads to find zeros, some nearer, some further. Further away zeros
    // should not override nearby zeros, and when a thread finds a zero, it
    // shouldn't cancel threads that are searching for nearer zeros.
    const vector<vector<vector<double>>> A = {
      {{0.4088715361390395, -0.9999112506968285, 0.8109731922797785, 0.25590203028822855},
       {-0.9125919498523434, -0.013322564689428938, 0.5850833115066139, 0.9667027210545974}},
      {{-0.704978485994098, -0.016909658985638815, 0.41508560377373277, 0.19893559514770887},
       {-0.05482092926492031, -0.7069045645863304, 0.5724543139323832, 0.6785459667430253}}};

    for (int i = 0; i < 100; i++)
    {
      ASSERT_EQ(0.5,
                computeGridUniquenessHypercube(A, 0.2, 0.5).first);
    }
  }
}
