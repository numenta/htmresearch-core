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
 * Implementations for GridUniqueness.hpp
 */

#include <nupic/experimental/GridUniqueness.hpp>
#include <nupic/utils/Log.hpp>

#include <math.h>
#include <signal.h>
#include <time.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using std::vector;
using std::pair;

static std::atomic<bool> quitting(false);

pair<double,double> transform2D(const vector<vector<double>>& M,
                                pair<double,double> p)
{
  return {M[0][0]*p.first + M[0][1]*p.second,
          M[1][0]*p.first + M[1][1]*p.second};
}

vector<vector<double>> invert2DMatrix(const vector<vector<double>>& M)
{
  const double detInv = 1 / (M[0][0]*M[1][1] - M[0][1]*M[1][0]);
  return {{detInv*M[1][1], -detInv*M[0][1]},
          {-detInv*M[1][0], detInv*M[0][0]}};
}

pair<double,double> transformND(const vector<vector<double>>& M,
                                const double p[])
{
  pair<double,double> result = {0, 0};

  for (size_t col = 0; col < M[0].size(); col++)
  {
    result.first += M[0][col]*p[col];
  }

  for (size_t col = 0; col < M[1].size(); col++)
  {
    result.second += M[1][col]*p[col];
  }

  return result;
}

/**
 * Enumerate the points of a lattice near or within a specified rectangle. This
 * is equivalent to checking whether any circles centered on the points of a
 * lattice overlap the rectangle.
 */
class LatticePointEnumerator
{
public:
  LatticePointEnumerator(const vector<vector<double>>& latticeBasis,
                         const vector<vector<double>>& inverseLatticeBasis,
                         double x0, double y0, double width, double height, double r)
    :latticeBasis_(latticeBasis), x0_(x0), y0_(y0), width_(width),
     height_(height), rSquared_(pow(r, 2))
  {
    // Find the bounding box of the rectangle in the lattice's basis.
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    pair<double, double> ij;
    ij = transform2D(inverseLatticeBasis, {x0 - r, y0 - r});
    xmin = std::min(xmin, ij.first);
    xmax = std::max(xmax, ij.first);
    ij = transform2D(inverseLatticeBasis, {x0 + width + r, y0 - r});
    xmin = std::min(xmin, ij.first);
    xmax = std::max(xmax, ij.first);
    ij = transform2D(inverseLatticeBasis, {x0 - r, y0 + height + r});
    xmin = std::min(xmin, ij.first);
    xmax = std::max(xmax, ij.first);
    ij = transform2D(inverseLatticeBasis, {x0 + width + r, y0 + height + r});
    xmin = std::min(xmin, ij.first);
    xmax = std::max(xmax, ij.first);

    iMin_ = ceil(xmin);
    iMax_ = floor(xmax);

    ij = transform2D(inverseLatticeBasis, {x0, y0});
    iStart_ = floor(ij.first);
    jStart_ = floor(ij.second);

    this->restart();
  }

  bool getNext(pair<double,double> *out)
  {
    bool foundContainedPoint = false;

    while (!foundContainedPoint && i_ <= iMax_)
    {
      const pair<double, double> p = transform2D(latticeBasis_, {i_, j_});

      const double nearestX = std::max(x0_,
                                       std::min(p.first,
                                                x0_ + width_));
      const double nearestY = std::max(y0_,
                                       std::min(p.second,
                                                y0_ + height_));

      const double dSquared = (pow(p.first - nearestX, 2) +
                               pow(p.second - nearestY, 2));

      if (j_ == j0_)
      {
        dSquared_j0_ = dSquared;
      }

      if (dSquared < innerSweepMin_)
      {
        innerSweepMin_ = dSquared;
        jForInnerSweepMin_ = j_;
      }

      foundContainedPoint = (dSquared <= rSquared_);
      if (foundContainedPoint)
      {
        *out = p;
      }

      // If we're moving away from the box, end this part of the inner sweep.
      const bool endInnerSweep = (!foundContainedPoint &&
                                  dSquared >= dSquaredPrev_);

      if (endInnerSweep)
      {
        if (sweepingDown_)
        {
          sweepingDown_ = false;

          j_ = j0_ + 1;
          dSquaredPrev_ = dSquared_j0_;
        }
        else
        {
          if (i_ == iStart_)
          {
            jForInnerSweepMin_i0_ = jForInnerSweepMin_;
          }

          if (sweepingLeft_)
          {
            if (--i_ < iMin_)
            {
              sweepingLeft_ = false;
              i_ = iStart_ + 1;

              j0_ = jForInnerSweepMin_i0_;
              j_ = j0_;
            }
            else
            {
              j0_ = jForInnerSweepMin_;
              j_ = j0_;
            }
          }
          else
          {
            ++i_;
            j0_ = jForInnerSweepMin_;
            j_ = j0_;
          }

          sweepingDown_ = true;
          innerSweepMin_ = std::numeric_limits<double>::max();
          dSquaredPrev_ = std::numeric_limits<double>::max();
        }
      }
      else
      {
        if (sweepingDown_)
        {
          j_ -= 1;
        }
        else
        {
          j_ += 1;
        }

        dSquaredPrev_ = dSquared;
      }
    }

    return foundContainedPoint;
  }

  void restart()
  {
    dSquaredPrev_ = std::numeric_limits<double>::max();
    innerSweepMin_ = std::numeric_limits<double>::max();
    i_ = iStart_;
    j_ = j0_ = jStart_;
    finished_ = false;
    sweepingLeft_ = true;
    sweepingDown_ = true;
  }

private:

  const vector<vector<double>> &latticeBasis_;
  const double x0_;
  const double y0_;
  const double width_;
  const double height_;
  const double rSquared_;

  double dSquaredPrev_;
  double innerSweepMin_;
  long long jForInnerSweepMin_;
  long long jForInnerSweepMin_i0_;
  double dSquared_j0_;

  bool finished_;
  bool sweepingLeft_;
  bool sweepingDown_;

  long long iStart_;
  long long jStart_;
  long long i_;
  long long j_;
  long long j0_;
  long long iMin_;
  long long iMax_;
};

/**
 * Enumerate the vertices of a hyperrectangle by incrementing an integer and
 * converting each binary representation into a vertex.
 */
class HyperrectangleVertexEnumerator
{
public:
  HyperrectangleVertexEnumerator(const double x0[],
                                 const double dims[],
                                 size_t numDims)
    :x0_(x0), dims_(dims), numDims_(numDims), upper_(pow(2, numDims)),
     bitvector_(0x0)
  {
  }

  bool getNext(double *out)
  {
    if (bitvector_ >= upper_)
    {
      return false;
    }

    for (size_t bit = 0; bit < numDims_; bit++)
    {
      out[bit] = x0_[bit];
      if (bitvector_ & (0x1 << bit))
      {
        out[bit] += dims_[bit];
      }
    }

    bitvector_++;

    return true;
  }

  void restart()
  {
    bitvector_ = 0x0;
  }

private:
  const double *x0_;
  const double *dims_;
  const size_t numDims_;
  const double upper_;
  unsigned int bitvector_;
};

/**
 * Quickly check a few points in this hyperrectangle to see if they have grid
 * code zero.
 */
bool tryFindGridCodeZero(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<vector<vector<double>>>& latticeBasisByModule,
  const vector<vector<vector<double>>>& inverseLatticeBasisByModule,
  size_t numDims,
  const double x0[],
  const double dims[],
  double readoutResolution,
  double vertexBuffer[])
{
  // Add a small epsilon to handle situations where floating point math causes a
  // vertex to be non-zero-overlapping here and zero-overlapping in
  // tryProveGridCodeZeroImpossible. With this addition, anything
  // zero-overlapping in tryProveGridCodeZeroImpossible is guaranteed to be
  // zero-overlapping here, so the program won't get caught in infinite
  // recursion.
  const double r = readoutResolution/2 + 0.000000001;
  const double rSquared = pow(r, 2);

  HyperrectangleVertexEnumerator vertices(x0, dims, numDims);
  while (vertices.getNext(vertexBuffer))
  {
    bool vertexDisqualified = false;

    for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
    {
      const pair<double, double> pointOnPlane =
        transformND(domainToPlaneByModule[iModule], vertexBuffer);

      LatticePointEnumerator latticePoints(latticeBasisByModule[iModule],
                                           inverseLatticeBasisByModule[iModule],
                                           pointOnPlane.first,
                                           pointOnPlane.second, 0, 0, r);

      bool isZero = false;

      pair<double, double> latticePoint;
      while (!isZero && latticePoints.getNext(&latticePoint))
      {
        isZero = (pow(latticePoint.first - pointOnPlane.first, 2) +
                  pow(latticePoint.second - pointOnPlane.second, 2) <= rSquared);
      }

      if (!isZero)
      {
        vertexDisqualified = true;
        break;
      }
    }

    if (!vertexDisqualified)
    {
      return true;
    }
  }

  return false;
}

/**
 * Quickly check whether this hyperrectangle excludes grid code zero
 * in any individual module.
 */
bool tryProveGridCodeZeroImpossible(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<vector<vector<double>>>& latticeBasisByModule,
  const vector<vector<vector<double>>>& inverseLatticeBasisByModule,
  size_t numDims,
  const double x0[],
  const double dims[],
  double readoutResolution,
  double vertexBuffer[])
{
  const double r = readoutResolution/2;

  HyperrectangleVertexEnumerator vertices(x0, dims, numDims);
  for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
  {
    // Figure out which lattice points we need to check.
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();
    vertices.restart();
    while (vertices.getNext(vertexBuffer))
    {
      const pair<double,double> p = transformND(domainToPlaneByModule[iModule],
                                                vertexBuffer);
      xmin = std::min(xmin, p.first);
      xmax = std::max(xmax, p.first);
      ymin = std::min(ymin, p.second);
      ymax = std::max(ymax, p.second);
    }
    LatticePointEnumerator latticePoints(latticeBasisByModule[iModule],
                                         inverseLatticeBasisByModule[iModule],
                                         xmin, ymin, (xmax - xmin), (ymax - ymin),
                                         r);

    pair<double, double> latticePoint;
    bool foundLatticeCollision = false;
    while (!foundLatticeCollision && latticePoints.getNext(&latticePoint))
    {
      // At this point, there might not actually be a lattice collision. The
      // bounding box collides with a lattice point, but it might not collide
      // with the actual polygon. However, we can just lazily assume that it
      // collided. If it didn't actually, this function will get called again
      // with a smaller box, and eventually the space will be broken into
      // sufficiently small boxes that each of them have no overlap with the
      // lattice of at least one module.
      foundLatticeCollision = true;
    }

    if (!foundLatticeCollision)
    {
      // This module never gets near grid code zero for the provided range of
      // locations. So this range can't possibly contain grid code zero.
      return true;
    }
  }

  return false;
}

/**
 * Temporarily set a value.
 */
class SwapValueRAII
{
public:
  SwapValueRAII(double* value, double newValue)
    :value_(value), initial_(*value)
  {
    *value_ = newValue;
  }

  ~SwapValueRAII()
  {
    *value_ = initial_;
  }

private:
  double *value_;
  double initial_;
};

/**
 * Helper function that doesn't allocate any memory, so it's much better for
 * recursion.
 */
bool findGridCodeZeroHelper(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<vector<vector<double>>>& latticeBasisByModule,
  const vector<vector<vector<double>>>& inverseLatticeBasisByModule,
  size_t numDims,
  double x0[],
  double dims[],
  double readoutResolution,
  double vertexBuffer[],
  std::atomic<bool>& shouldContinue)
{
  if (!shouldContinue)
  {
    return false;
  }

  if (tryFindGridCodeZero(domainToPlaneByModule, latticeBasisByModule,
                          inverseLatticeBasisByModule, numDims, x0, dims,
                          readoutResolution, vertexBuffer))
  {
    return true;
  }

  if (tryProveGridCodeZeroImpossible(domainToPlaneByModule,
                                     latticeBasisByModule,
                                     inverseLatticeBasisByModule, numDims, x0,
                                     dims, readoutResolution, vertexBuffer))
  {
    return false;
  }

  size_t iWidestDim = std::distance(dims,
                                    std::max_element(dims, dims + numDims));
  {
    SwapValueRAII swap1(&dims[iWidestDim], dims[iWidestDim] / 2);

    if (findGridCodeZeroHelper(domainToPlaneByModule, latticeBasisByModule,
                               inverseLatticeBasisByModule, numDims, x0, dims,
                               readoutResolution, vertexBuffer, shouldContinue))
    {
      return true;
    }

    {
      SwapValueRAII swap2(&x0[iWidestDim], x0[iWidestDim] + dims[iWidestDim]);
      return findGridCodeZeroHelper(domainToPlaneByModule, latticeBasisByModule,
                                    inverseLatticeBasisByModule, numDims, x0,
                                    dims, readoutResolution, vertexBuffer,
                                    shouldContinue);
    }
  }
}

bool nupic::experimental::grid_uniqueness::findGridCodeZero(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<vector<vector<double>>>& latticeBasisByModule,
  const vector<double>& x0,
  const vector<double>& dims,
  double readoutResolution,
  vector<double>* pointWithGridCodeZero)
{
  // Avoid doing any allocations in each recursion.
  vector<double> x0Copy(x0);
  vector<double> dimsCopy(dims);
  std::atomic<bool> shouldContinue(true);

  vector<double> defaultPointBuffer;

  if (pointWithGridCodeZero != nullptr)
  {
    NTA_ASSERT(pointWithGridCodeZero->size() == dims.size());
  }
  else
  {
    defaultPointBuffer.resize(dims.size());
    pointWithGridCodeZero = &defaultPointBuffer;
  }

  NTA_ASSERT(domainToPlaneByModule[0].size() == 2);

  vector<vector<vector<double>>> inverseLatticeBasisByModule;
  for (const vector<vector<double>>& latticeBasis : latticeBasisByModule)
  {
    inverseLatticeBasisByModule.push_back(invert2DMatrix(latticeBasis));
  }

  return findGridCodeZeroHelper(
    domainToPlaneByModule, latticeBasisByModule, inverseLatticeBasisByModule,
    dimsCopy.size(), x0Copy.data(), dimsCopy.data(), readoutResolution,
    pointWithGridCodeZero->data(), shouldContinue);
}


struct GridUniquenessState {
  // Constants (thread-safe)
  const vector<vector<vector<double>>>& domainToPlaneByModule;
  const vector<vector<vector<double>>>& latticeBasisByModule;
  const vector<vector<vector<double>>>& inverseLatticeBasisByModule;
  const double readoutResolution;
  const size_t numDims;

  // Task management
  double baselineRadius;
  double expansionRadiusGoal;
  vector<double> expansionProgress;
  size_t expandingDim;
  bool positiveExpand;
  bool continueExpansion;

  // Results
  vector<double> pointWithGridCodeZero;
  double foundPointBaselineRadius;

  // Thread management
  std::mutex& mutex;
  std::condition_variable& finished;
  size_t numActiveThreads;
  vector<double> threadBaselineRadius;
  vector<vector<double>> threadQueryX0;
  vector<vector<double>> threadQueryDims;
  vector<std::atomic<bool>> threadShouldContinue;
  vector<bool> threadRunning;
};

std::string vecs(const vector<double>& v)
{
  std::ostringstream oss;
  oss << "[";

  for (size_t i = 0; i < v.size(); i++)
  {
    oss << v[i];
    if (i < v.size() - 1)
    {
      oss << ", ";
    }
  }

  oss << "]";

  return oss.str();
}

void recordResult(size_t iThread, GridUniquenessState& state,
                  const vector<double>& pointWithGridCodeZero)
{
  state.continueExpansion = false;
  if (state.threadBaselineRadius[iThread] < state.foundPointBaselineRadius)
  {
    state.foundPointBaselineRadius = state.threadBaselineRadius[iThread];
    state.pointWithGridCodeZero = pointWithGridCodeZero;

    // Notify all others that they should stop unless they're checking a lower
    // base width.
    for (size_t iOtherThread = 0;
         iOtherThread < state.threadBaselineRadius.size();
         iOtherThread++)
    {
      if (iOtherThread != iThread &&
          state.threadShouldContinue[iOtherThread] &&
          (state.threadBaselineRadius[iOtherThread] >=
           state.foundPointBaselineRadius))
      {
        state.threadShouldContinue[iOtherThread] = false;
      }
    }
  }
}

void claimNextTask(size_t iThread, GridUniquenessState& state)
{
  state.threadBaselineRadius[iThread] = state.baselineRadius;

  vector<double>& x0 = state.threadQueryX0[iThread];
  vector<double>& dims = state.threadQueryDims[iThread];

  // Determine all but the final dimension.
  for (size_t iDim = 0; iDim < state.numDims - 1; iDim++)
  {
    dims[iDim] = 2*state.expansionProgress[iDim];
    x0[iDim] = -state.expansionProgress[iDim];
  }

  // Optimization: for the final dimension, don't go negative. Half of the
  // hypercube will be equal-and-opposite phases of the other half, so we ignore
  // the lower half of the final dimension.
  dims[state.numDims - 1] = state.expansionProgress[state.numDims - 1];
  x0[state.numDims - 1] = 0;

  // Make the changes specific to this query.
  dims[state.expandingDim] = state.expansionRadiusGoal - state.baselineRadius;
  x0[state.expandingDim] = state.positiveExpand
    ? state.baselineRadius
    : -state.expansionRadiusGoal;

  // Update.
  if (state.positiveExpand &&
      // Optimization: Don't check the negative for the final
      // dimension. (described above)
      state.expandingDim < state.numDims - 1)
  {
    state.positiveExpand = false;
  }
  else
  {
    state.positiveExpand = true;
    state.expansionProgress[state.expandingDim] = state.expansionRadiusGoal;
    if (++state.expandingDim >= state.numDims)
    {
      state.baselineRadius = state.expansionRadiusGoal;
      state.expansionRadiusGoal *= 1.01;
      state.expandingDim = 0;
    }
  }
}

void findGridCodeZeroThread(size_t iThread, GridUniquenessState& state)
{
  bool foundGridCodeZero = false;
  vector<double> x0(state.numDims);
  vector<double> dims(state.numDims);
  vector<double> pointWithGridCodeZero(state.numDims);

  while (!quitting)
  {
    // Modify the shared state. Record the results, decide the next task,
    // volunteer to do it.
    {
      std::lock_guard<std::mutex> lock(state.mutex);

      if (foundGridCodeZero)
      {
        recordResult(iThread, state, pointWithGridCodeZero);
      }

      if (!state.continueExpansion)
      {
        break;
      }

      // Select task params.
      claimNextTask(iThread, state);

      // Make an unshared copy that findGridCodeZeroHelper can modify.
      x0 = state.threadQueryX0[iThread];
      dims = state.threadQueryDims[iThread];
    }

    // Perform the task.
    foundGridCodeZero = findGridCodeZeroHelper(
      state.domainToPlaneByModule, state.latticeBasisByModule,
      state.inverseLatticeBasisByModule, state.numDims, x0.data(), dims.data(),
      state.readoutResolution, pointWithGridCodeZero.data(),
      state.threadShouldContinue[iThread]);
  }

  // This thread is exiting.
  {
    std::lock_guard<std::mutex> lock(state.mutex);
    if (--state.numActiveThreads == 0)
    {
      state.finished.notify_all();
    }
    state.threadRunning[iThread] = false;
  }
}

pair<double,double> rotateClockwise(double theta, double x, double y)
{
  return {cos(theta)*x + sin(theta)*y,
          -sin(theta)*x + cos(theta)*y};
}

/**
 * Change the matrices so that movement tends to be along the x axis on the
 * plane. The end result is the same, but this improves performance because we
 * draw bounding boxes around the projected hyperrectangles. This is especially
 * beneficial in 1D, totally eliminating diagonal motion so that the bounding
 * box perfectly encloses the projected line.
 */
void optimizeMatrices(vector<vector<vector<double>>> *domainToPlaneByModule,
                      vector<vector<vector<double>>> *latticeBasisByModule)
{
  for (size_t iModule = 0; iModule < domainToPlaneByModule->size(); iModule++)
  {
    vector<vector<double>> &domainToPlane = (*domainToPlaneByModule)[iModule];
    vector<vector<double>> &latticeBasis = (*latticeBasisByModule)[iModule];

    size_t iLongest = (size_t) -1;
    double dLongest = std::numeric_limits<double>::max();
    for (size_t iColumn = 0; iColumn < domainToPlane[0].size(); iColumn++)
    {
      double length = sqrt(pow(domainToPlane[0][iColumn], 2) +
                           pow(domainToPlane[1][iColumn], 2));
      if (length < dLongest)
      {
        dLongest = length;
        iLongest = iColumn;
      }
    }

    const double theta = atan2(domainToPlane[1][iLongest],
                               domainToPlane[0][iLongest]);
    for (size_t iColumn = 0; iColumn < domainToPlane[0].size(); iColumn++)
    {
      const pair<double, double> newColumn = rotateClockwise(theta,
                                                             domainToPlane[0][iColumn],
                                                             domainToPlane[1][iColumn]);
      domainToPlane[0][iColumn] = newColumn.first;
      domainToPlane[1][iColumn] = newColumn.second;
    }
    for (size_t iColumn = 0; iColumn < latticeBasis[0].size(); iColumn++)
    {
      const pair<double, double> newColumn = rotateClockwise(theta,
                                                             latticeBasis[0][iColumn],
                                                             latticeBasis[1][iColumn]);
      latticeBasis[0][iColumn] = newColumn.first;
      latticeBasis[1][iColumn] = newColumn.second;
    }
  }
}

pair<double,vector<double>>
nupic::experimental::grid_uniqueness::computeGridUniquenessHypercube(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<vector<vector<double>>>& latticeBasisByModule,
  double readoutResolution,
  double ignoredCenterDiameter)
{
  typedef std::chrono::steady_clock Clock;

  // Manually handle interrupts so that they're handled when running in a
  // Jupyter notebook, and to make the threads return cleanly.
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler =
    [](int s) { quitting = true; };
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;
  sigaction(SIGINT, &sigIntHandler, NULL);

  NTA_CHECK(domainToPlaneByModule.size() == latticeBasisByModule.size())
    << "The two arrays of matrices must be the same length (one per module) "
    << "Actual: " << domainToPlaneByModule.size()
    << " " << latticeBasisByModule.size();

  NTA_CHECK(domainToPlaneByModule[0].size() == 2)
    << "Each matrix should have two rows -- the modules are two-dimensional. "
    << "Actual: " << domainToPlaneByModule[0].size();

  NTA_CHECK(latticeBasisByModule[0][0].size() == 2)
    << "There should be two lattice basis vectors. "
    << "Actual: " << latticeBasisByModule[0][0].size();

  const size_t numDims = domainToPlaneByModule[0][0].size();
  NTA_CHECK(numDims < sizeof(int)*8)
    << "Unsupported number of dimensions: " << numDims;

  vector<vector<vector<double>>> domainToPlaneByModule2(domainToPlaneByModule);
  vector<vector<vector<double>>> latticeBasisByModule2(latticeBasisByModule);
  optimizeMatrices(&domainToPlaneByModule2, &latticeBasisByModule2);

  vector<vector<vector<double>>> inverseLatticeBasisByModule;
  for (const vector<vector<double>>& latticeBasis : latticeBasisByModule2)
  {
    inverseLatticeBasisByModule.push_back(invert2DMatrix(latticeBasis));
  }

  // Use condition_variables to enable periodic logging while waiting for the
  // threads to finish.
  std::mutex stateMutex;
  std::condition_variable finished;

  size_t numThreads = std::thread::hardware_concurrency();

  GridUniquenessState state = {
    domainToPlaneByModule2,
    latticeBasisByModule2,
    inverseLatticeBasisByModule,
    readoutResolution,
    numDims,

    ignoredCenterDiameter,
    ignoredCenterDiameter * 2,
    vector<double>(numDims, ignoredCenterDiameter),
    0,
    true,
    true,

    vector<double>(numDims),
    std::numeric_limits<double>::max(),

    stateMutex,
    finished,
    0,
    vector<double>(numThreads, std::numeric_limits<double>::max()),
    vector<vector<double>>(numThreads, vector<double>(numDims)),
    vector<vector<double>>(numThreads, vector<double>(numDims)),
    vector<std::atomic<bool>>(numThreads),
    vector<bool>(numThreads, true),
  };

  for (size_t i = 0; i < numThreads; i++)
  {
    state.threadShouldContinue[i] = true;
  }

  {
    std::unique_lock<std::mutex> lock(stateMutex);
    for (size_t i = 0; i < numThreads; i++)
    {
      std::thread(findGridCodeZeroThread, i,  std::ref(state)).detach();
      state.numActiveThreads++;
    }

    const auto tStart = Clock::now();
    auto tNextPrint = tStart + std::chrono::seconds(10);

    bool processingQuit = false;

    while (true)
    {
      if (state.finished.wait_until(lock,
                                    tNextPrint) == std::cv_status::timeout)
      {
        NTA_INFO << "";
        NTA_INFO << domainToPlaneByModule.size() << " modules, " << numDims
                 << " dimensions, "
                 << std::chrono::duration_cast<std::chrono::seconds>(
                   Clock::now() - tStart).count() << " seconds elapsed";

        if (state.foundPointBaselineRadius < std::numeric_limits<double>::max())
        {
          NTA_INFO << "**Hypercube side length upper bound: "
                   << state.foundPointBaselineRadius << "**";
          NTA_INFO << "**Grid code zero found at: "
                   << vecs(state.pointWithGridCodeZero) << "**";
        }

        tNextPrint = Clock::now() + std::chrono::seconds(10);

        for (size_t iThread = 0; iThread < state.threadBaselineRadius.size();
             iThread++)
        {
          if (state.threadRunning[iThread])
          {
            if (state.threadShouldContinue[iThread])
            {
              NTA_INFO << "  Thread " << iThread
                       << " assuming hypercube side length lower bound "
                       << state.threadBaselineRadius[iThread] << ", querying x0 "
                       << vecs(state.threadQueryX0[iThread]) << " and dims "
                       << vecs(state.threadQueryDims[iThread]);
            }
            else
            {
              NTA_INFO << "  Thread " << iThread << " has been ordered to stop.";
            }
          }
          else
          {
            NTA_INFO << "  Thread " << iThread << " is finished.";
          }
        }
      }
      else if (quitting && !processingQuit)
      {
        // The condition_variable returned due to an interrupt. We still need to
        // wait for threads to exit.
        processingQuit = true;
        for (size_t iThread = 0; iThread < state.threadBaselineRadius.size();
             iThread++)
        {
          state.threadShouldContinue[iThread] = false;
        }
      }
      else
      {
        break;
      }
    }
  }

  if (quitting)
  {
    // The process might not be ending, the caller (e.g. the Python shell) is
    // likely to catch this exception and continue, so prepare to run again.
    quitting = false;
    NTA_THROW << "Caught interrupt signal";
  }

  return {state.foundPointBaselineRadius, state.pointWithGridCodeZero};
}
