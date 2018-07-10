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

#include <limits>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <chrono>
#include <atomic>
#include <algorithm>

using std::vector;
using std::pair;

static bool quitting = false;

pair<double,double> computePhaseDisplacement(
  const vector<vector<double>>& A,
  const double dx[])
{
  pair<double,double> result = {0, 0};

  for (size_t col = 0; col < A[0].size(); col++)
  {
    result.first += A[0][col]*dx[col];
  }

  for (size_t col = 0; col < A[1].size(); col++)
  {
    result.second += A[1][col]*dx[col];
  }

  return result;
}


/**
 * Checks whether a 2D range overlaps any squares of width w centered on points
 * in the square lattice.
 */
bool overlapsAnySquareLatticeSquare(double left, double right, double bottom,
                                    double top, double w)
{
  const double r = w/2;

  // Find the square lattice "hull" of the 2D range. Expand the range outward to
  // capture any squares that are centered near the range, then collapse it
  // inward (via ceil and floor) to find the first and last square along each
  // axis.
  const double hullLeft = ceil(left - r);
  const double hullRight = floor(right + r);
  const double hullBottom = ceil(bottom - r);
  const double hullTop = floor(top + r);

  for (double x = hullLeft; x <= hullRight; x += 1.0)
  {
    for (double y = hullBottom; y <= hullTop; y += 1.0)
    {
      if (right < x - r ||
          x + r < left ||
          top < y - r ||
          y + r < bottom)
      {
        // They don't overlap
      }
      else
      {
        return true;
      }
    }
  }

  return false;
}

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
 * displacement zero.
 */
bool tryFindGridDisplacementZero(
  const vector<vector<vector<double>>>& A,
  size_t numDims,
  const double x0[],
  const double dims[],
  double phaseResolution,
  double vertexBuffer[])
{
  // Add a small epsilon to handle situations where floating point math causes a
  // vertex to be non-zero-overlapping here and zero-overlapping in
  // tryProveNegative. With this addition, anything zero-overlapping in
  // tryProveNegative is guaranteed to be zero-overlapping here, so the program
  // won't get caught in infinite recursion.
  const double r = phaseResolution/2 + 0.000000001;

  HyperrectangleVertexEnumerator vertices(x0, dims, numDims);
  while (vertices.getNext(vertexBuffer))
  {
    bool vertexDisqualified = false;
    for (const vector<vector<double>>& module : A)
    {
      pair<double,double> dPhase = computePhaseDisplacement(module,
                                                            vertexBuffer);
      dPhase.first -= floor(dPhase.first);
      dPhase.second -= floor(dPhase.second);
      if (dPhase.first > 0.5) { dPhase.first = 1.0 - dPhase.first; }
      if (dPhase.second > 0.5) { dPhase.second = 1.0 - dPhase.second; }
      if (dPhase.first > r || dPhase.second > r)
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
 * Quickly check whether this hyperrectangle excludes grid displacement zero in
 * any individual module.
 */
bool tryProveGridDisplacementZeroImpossible(
  const vector<vector<vector<double>>>& A,
  size_t numDims,
  const double x0[],
  const double dims[],
  double phaseResolution,
  double vertexBuffer[])
{
  HyperrectangleVertexEnumerator vertices(x0, dims, numDims);
  for (const vector<vector<double>>& module : A)
  {
    double projection_left = std::numeric_limits<double>::max();
    double projection_right = std::numeric_limits<double>::lowest();
    double projection_bottom = std::numeric_limits<double>::max();
    double projection_top = std::numeric_limits<double>::lowest();
    vertices.restart();
    while (vertices.getNext(vertexBuffer))
    {
      pair<double,double> phase = computePhaseDisplacement(module,
                                                           vertexBuffer);
      projection_left = std::min(projection_left, phase.first);
      projection_right = std::max(projection_right, phase.first);
      projection_bottom = std::min(projection_bottom, phase.second);
      projection_top = std::max(projection_top, phase.second);
    }

    if (!overlapsAnySquareLatticeSquare(projection_left, projection_right,
                                        projection_bottom, projection_top,
                                        phaseResolution))
    {
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
bool findDisplacementZeroHelper(
  const vector<vector<vector<double>>>& A,
  size_t numDims,
  double x0[],
  double dims[],
  double phaseResolution,
  double vertexBuffer[],
  std::atomic<bool>& shouldContinue)
{
  if (!shouldContinue)
  {
    return false;
  }

  if (tryFindGridDisplacementZero(A, numDims, x0, dims, phaseResolution,
                                  vertexBuffer))
  {
    return true;
  }

  if (tryProveGridDisplacementZeroImpossible(A, numDims, x0, dims,
                                             phaseResolution, vertexBuffer))
  {
    return false;
  }

  size_t iWidestDim = std::distance(dims,
                                    std::max_element(dims, dims + numDims));
  {
    SwapValueRAII swap1(&dims[iWidestDim], dims[iWidestDim] / 2);

    if (findDisplacementZeroHelper(A, numDims, x0, dims, phaseResolution,
                                   vertexBuffer, shouldContinue))
    {
      return true;
    }

    {
      SwapValueRAII swap2(&x0[iWidestDim], x0[iWidestDim] + dims[iWidestDim]);
      return findDisplacementZeroHelper(A, numDims, x0, dims, phaseResolution,
                                        vertexBuffer, shouldContinue);
    }
  }
}

bool nupic::experimental::grid_uniqueness::findGridDisplacementZero(
  const vector<vector<vector<double>>>& A,
  const vector<double>& x0,
  const vector<double>& dims,
  double phaseResolution,
  vector<double>& displacementWithGridCodeZero)
{
  // Avoid doing any allocations in each recursion.
  vector<double> x0Copy(x0);
  vector<double> dimsCopy(dims);
  std::atomic<bool> shouldContinue(true);

  NTA_ASSERT(A[0].size() == 2);
  NTA_ASSERT(displacementWithGridCodeZero.size() == dims.size());

  return findDisplacementZeroHelper(
    A, dimsCopy.size(), x0Copy.data(), dimsCopy.data(), phaseResolution,
    displacementWithGridCodeZero.data(), shouldContinue);
}


struct GridUniquenessState {
  // Constants (thread-safe)
  const vector<vector<vector<double>>>& A;
  const double phaseResolution;
  const size_t numDims;

  // Task management
  double baselineRadius;
  double expansionRadiusGoal;
  vector<double> expansionProgress;
  size_t expandingDim;
  bool positiveExpand;
  bool continueExpansion;

  // Results
  vector<double> displacementWithGridCodeZero;
  double displacementBaselineRadius;

  // Thread management
  std::mutex& mutex;
  std::condition_variable& finished;
  size_t numActiveThreads;
  vector<double> threadBaselineRadius;
  vector<vector<double>> threadQueryX0;
  vector<vector<double>> threadQueryDims;
  vector<std::atomic<bool>> threadShouldContinue;
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
                  const vector<double>& displacementWithGridCodeZero)
{
  state.continueExpansion = false;
  if (state.threadBaselineRadius[iThread] < state.displacementBaselineRadius)
  {
    state.displacementBaselineRadius = state.threadBaselineRadius[iThread];
    state.displacementWithGridCodeZero = displacementWithGridCodeZero;

    // Notify all others that they should stop unless they're checking a lower
    // base width.
    for (size_t iOtherThread = 0;
         iOtherThread < state.threadBaselineRadius.size();
         iOtherThread++)
    {
      if (iOtherThread != iThread &&
          state.threadShouldContinue[iOtherThread] &&
          (state.threadBaselineRadius[iOtherThread] >=
           state.displacementBaselineRadius))
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

void findDisplacementZeroThread(size_t iThread, GridUniquenessState& state)
{
  bool foundGridDisplacementZero = false;
  vector<double> x0(state.numDims);
  vector<double> dims(state.numDims);
  vector<double> displacementWithGridCodeZero(state.numDims);

  while (!quitting)
  {
    // Modify the shared state. Record the results, decide the next task,
    // volunteer to do it.
    {
      std::lock_guard<std::mutex> lock(state.mutex);

      if (foundGridDisplacementZero)
      {
        recordResult(iThread, state, displacementWithGridCodeZero);
      }

      if (!state.continueExpansion)
      {
        break;
      }

      // Select task params.
      claimNextTask(iThread, state);

      // Make an unshared copy that findDisplacementZeroHelper can modify.
      x0 = state.threadQueryX0[iThread];
      dims = state.threadQueryDims[iThread];
    }

    // Perform the task.
    foundGridDisplacementZero = findDisplacementZeroHelper(
      state.A, state.numDims, x0.data(), dims.data(),
      state.phaseResolution,
      displacementWithGridCodeZero.data(), state.threadShouldContinue[iThread]);
  }

  // This thread is exiting.
  {
    std::lock_guard<std::mutex> lock(state.mutex);
    if (--state.numActiveThreads == 0)
    {
      state.finished.notify_all();
    }
  }
}


pair<double,vector<double>>
nupic::experimental::grid_uniqueness::computeGridUniquenessHypercube(
  const vector<vector<vector<double>>>& A,
  double phaseResolution,
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

  NTA_CHECK(A[0].size() == 2)
    << "Each matrix should have two rows -- the modules are two-dimensional. "
    << "Actual: " << A[0].size();

  const size_t numDims = A[0][0].size();
  NTA_CHECK(numDims < sizeof(int)*8)
    << "Unsupported number of dimensions: " << numDims;

  // Use condition_variables to enable periodic logging while waiting for the
  // threads to finish.
  std::mutex stateMutex;
  std::condition_variable finished;

  size_t numThreads = std::thread::hardware_concurrency();

  GridUniquenessState state = {
    A,
    phaseResolution,
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
    vector<std::atomic<bool>>(numThreads)
  };

  for (size_t i = 0; i < numThreads; i++)
  {
    state.threadShouldContinue[i] = true;
  }

  {
    std::unique_lock<std::mutex> lock(stateMutex);
    for (size_t i = 0; i < numThreads; i++)
    {
      std::thread(findDisplacementZeroThread, i,  std::ref(state)).detach();
      state.numActiveThreads++;
    }

    const auto tStart = Clock::now();
    auto tNextPrint = tStart + std::chrono::seconds(10);

    while (true)
    {
      if (state.finished.wait_until(lock,
                                    tNextPrint) == std::cv_status::timeout)
      {
        NTA_INFO << "";
        NTA_INFO << A.size() << " modules, " << numDims << " dimensions, "
                 << std::chrono::duration_cast<std::chrono::seconds>(
                   Clock::now() - tStart).count() << " seconds elapsed, "
                 << "currently querying radius " << state.expansionRadiusGoal;
        tNextPrint = Clock::now() + std::chrono::seconds(10);

        for (size_t iThread = 0; iThread < state.threadBaselineRadius.size();
             iThread++)
        {
          if (state.threadShouldContinue[iThread])
          {
            NTA_INFO << "  Thread " << iThread << " is currently querying x0 "
                     << vecs(state.threadQueryX0[iThread]) << " and dims "
                     << vecs(state.threadQueryDims[iThread]);
          }
          else
          {
            NTA_INFO << " Thread " << iThread << " has been ordered to stop.";
          }
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

  return {state.displacementBaselineRadius, state.displacementWithGridCodeZero};
}
