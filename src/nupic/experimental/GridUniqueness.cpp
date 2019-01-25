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

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/geometries/register/point.hpp>

using std::vector;
using std::pair;
namespace bg = boost::geometry;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)


template<typename T>
class ThreadSafeQueue {
public:
  void put(T v)
  {
    std::unique_lock<std::mutex> lock(mutex_);
    queue_.push(v);
    cv_.notify_one();
  }

  T take()
  {
    std::unique_lock<std::mutex> lock(mutex_);
    while (queue_.empty())
    {
      cv_.wait(lock);
    }

    T ret = queue_.front();
    queue_.pop();
    return ret;
  }

private:
  std::mutex mutex_;
  std::condition_variable cv_;
  std::queue<T> queue_;
};


class ScheduledTask {
public:
  template <typename T, typename F>
  ScheduledTask(T timeout, F f)
  {
    timerThread_ = std::thread(
      // Don't use default capture [&]. It causes the addresses of the captured
      // variables of f to get mangled, though this effect goes away when this
      // method/constructor is inlined so it only occurs on debug builds or with
      // __attribute__((noinline)). This occurs on clang and gcc, maybe others.
      [this, timeout, f]() { this->waiterThread_(timeout, f); });
  }

  ~ScheduledTask()
  {
    {
      std::unique_lock<std::mutex> lock(timerMutex_);
      timerCancelled_ = true;
      timerCondition_.notify_one();
    }
    timerThread_.join();
  }

private:

  template <typename T, typename F>
  void waiterThread_(T tTimeout, F f)
  {
    std::unique_lock<std::mutex> lock(this->timerMutex_);

    while (!this->timerCancelled_)
    {
      if (this->timerCondition_.wait_until(lock, tTimeout) ==
          std::cv_status::timeout)
      {
        f();
        this->timerCancelled_ = true;
      }
    }
  }

  std::mutex timerMutex_;
  std::condition_variable timerCondition_;
  bool timerCancelled_ = false;
  std::thread timerThread_;
};


enum Message {
  Interrupt,
  Timeout,
  Exiting
};


static size_t g_captureInterruptsCounter = 0;
static std::mutex g_messageQueuesMutex;
static vector<ThreadSafeQueue<Message>*> g_messageQueues;
static void (*g_prevHandler)(int) = nullptr;


// Custom interrupt processing is particularly necessary in Jupyter notebooks.
void processInterrupt(int sig)
{
  {
    std::unique_lock<std::mutex> lock(g_messageQueuesMutex);
    for (ThreadSafeQueue<Message>* messageQueue : g_messageQueues)
    {
      messageQueue->put(Message::Interrupt);
    }
  }
}

class CaptureInterruptsRAII
{
public:
  CaptureInterruptsRAII(ThreadSafeQueue<Message>* messages)
    : messages_(messages)
  {
    std::unique_lock<std::mutex> lock(g_messageQueuesMutex);

    if (g_captureInterruptsCounter++ == 0)
    {
      g_prevHandler = signal(SIGINT, processInterrupt);
    }

    g_messageQueues.push_back(messages);
  }

  ~CaptureInterruptsRAII()
  {
    std::unique_lock<std::mutex> lock(g_messageQueuesMutex);

    if (--g_captureInterruptsCounter == 0)
    {
      signal(SIGINT, g_prevHandler);
      g_prevHandler = nullptr;
    }

    for (auto it = g_messageQueues.begin(); it != g_messageQueues.end(); it++)
    {
      if (*it == messages_)
      {
        g_messageQueues.erase(it);
        break;
      }
    }
  }

private:
  ThreadSafeQueue<Message>* messages_;
};

template<typename T>
struct SquareMatrix2D {
  T v00;
  T v01;
  T v10;
  T v11;
};

pair<double,double> transform2D(const SquareMatrix2D<double>& M,
                                pair<double,double> p)
{
  return {M.v00*p.first + M.v01*p.second,
          M.v10*p.first + M.v11*p.second};
}

SquareMatrix2D<double> invert2DMatrix(const vector<vector<double>>& M)
{
  const double detInv = 1 / (M[0][0]*M[1][1] - M[0][1]*M[1][0]);
  return {detInv*M[1][1], -detInv*M[0][1],
          -detInv*M[1][0], detInv*M[0][0]};
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


struct LatticeBox {
  double xmin;
  double xmax;
  pair<double,double> middle;
};

/**
 * Enumerate the points of a lattice near or within a specified rectangle. This
 * is equivalent to checking whether any circles centered on the points of a
 * lattice overlap the rectangle.
 *
 * This enumeration works by first guessing at an initial lattice point (i, j),
 * then performing a series of nested sweeps. It sweeps downward, decrementing j
 * until decrementing j causes the distance of the lattice point from the
 * rectangle to increase. Then it sweeps upward using the same rules. It sweeps
 * leftward, decrementing i and repeating these downward and upward sweeps at
 * each i until it passes the lowest i possible for this rectangle. It repeats
 * the process sweeping right. As it sweeps left and right it attempts to choose
 * a good starting j value by considering results from the previous downward and
 * upward sweeps.
 */
class LatticePointEnumerator
{
public:
  LatticePointEnumerator(const SquareMatrix2D<double>& latticeBasis,
                         const SquareMatrix2D<double>& inverseLatticeBasis,
                         const LatticeBox& cachedLatticeBox,
                         const pair<double,double>& shift,
                         double left, double right, double bottom, double top,
                         double rSquared)
    :latticeBasis_(latticeBasis), left_(left), right_(right), bottom_(bottom),
     top_(top), rSquared_(rSquared)
  {
    const pair<double,double> ijShift = transform2D(inverseLatticeBasis, shift);
    iMin_ = ceil(cachedLatticeBox.xmin + ijShift.first);
    iMax_ = floor(cachedLatticeBox.xmax + ijShift.first);
    i_ = iStart_ = floor(cachedLatticeBox.middle.first + ijShift.first);
    j_ = j0_ = jStart_ = floor(cachedLatticeBox.middle.second + ijShift.second);
  }

  LatticePointEnumerator(const SquareMatrix2D<double>& latticeBasis,
                         const SquareMatrix2D<double>& inverseLatticeBasis,
                         double left, double right, double bottom, double top,
                         double r, double rSquared)
    :latticeBasis_(latticeBasis), left_(left), right_(right), bottom_(bottom),
     top_(top), rSquared_(rSquared)
  {
    // Find the bounding box of the rectangle in the lattice's basis.
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();

    const double paddedLeft = left_ - r;
    const double paddedRight = right_ + r;
    const double paddedBottom = bottom_ - r;
    const double paddedTop = top_ + r;

    pair<double, double> ij;
    ij = transform2D(inverseLatticeBasis, {paddedLeft, paddedBottom});
    xmin = std::min(xmin, ij.first);
    xmax = std::max(xmax, ij.first);
    ij = transform2D(inverseLatticeBasis, {paddedRight, paddedBottom});
    xmin = std::min(xmin, ij.first);
    xmax = std::max(xmax, ij.first);
    ij = transform2D(inverseLatticeBasis, {paddedLeft, paddedTop});
    xmin = std::min(xmin, ij.first);
    xmax = std::max(xmax, ij.first);
    ij = transform2D(inverseLatticeBasis, {paddedRight, paddedTop});
    xmin = std::min(xmin, ij.first);
    xmax = std::max(xmax, ij.first);

    iMin_ = ceil(xmin);
    iMax_ = floor(xmax);

    ij = transform2D(inverseLatticeBasis, {left, bottom});
    iStart_ = floor(ij.first);
    jStart_ = floor(ij.second);

    i_ = iStart_;
    j_ = j0_ = jStart_;
  }

  bool getNext(pair<double,double> *out)
  {
    bool foundContainedPoint = false;

    while (!foundContainedPoint && (sweepingLeft_ || i_ <= iMax_))
    {
      const pair<double, double> p = transform2D(latticeBasis_, {i_, j_});

      const double nearestX = std::max(left_,
                                       std::min(p.first,
                                                right_));
      const double nearestY = std::max(bottom_,
                                       std::min(p.second,
                                                top_));

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
          --j_;
        }
        else
        {
          ++j_;
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

  const SquareMatrix2D<double>& latticeBasis_;
  const double left_;
  const double right_;
  const double bottom_;
  const double top_;
  const double rSquared_;

  double dSquaredPrev_ = std::numeric_limits<double>::max();
  double innerSweepMin_ = std::numeric_limits<double>::max();
  long long jForInnerSweepMin_ = std::numeric_limits<long long>::lowest();
  long long jForInnerSweepMin_i0_ = std::numeric_limits<long long>::lowest();
  double dSquared_j0_;

  bool finished_ = false;
  bool sweepingLeft_ = true;
  bool sweepingDown_ = true;

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
  HyperrectangleVertexEnumerator(const double dims[], size_t numDims)
    :dims_(dims), numDims_(numDims), upper_(pow(2, numDims)), bitvector_(0x0)
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
      out[bit] = (bitvector_ & (0x1 << bit))
        ? dims_[bit]
        : 0;
    }

    bitvector_++;

    return true;
  }

  void restart()
  {
    bitvector_ = 0x0;
  }

private:
  const double *dims_;
  const size_t numDims_;
  const double upper_;
  unsigned int bitvector_;
};

/**
 * Compute d % 1.0, returning a value within the range [-0.5, 0.5]
 */
double mod1_05(double d)
{
  double ret = d - floor(d);
  if (ret > 0.5)
  {
    ret -= 1.0;
  }

  return ret;
}

/**
 * Quickly check a few points in this hyperrectangle to see if they have grid
 * code zero.
 */
bool tryFindGridCodeZero(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<SquareMatrix2D<double>>& latticeBasisByModule,
  const vector<SquareMatrix2D<double>>& inverseLatticeBasisByModule,
  size_t numDims,
  const double x0[],
  const double dims[],
  double rSquared,
  double vertexBuffer[])
{
  for (size_t iDim = 0; iDim < numDims; iDim++)
  {
    vertexBuffer[iDim] = x0[iDim] + (dims[iDim]/2);
  }

  for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
  {
    const pair<double, double> pointOnPlane =
      transformND(domainToPlaneByModule[iModule], vertexBuffer);

    const pair<double, double> pointOnUnrolledTorus =
      transform2D(inverseLatticeBasisByModule[iModule], pointOnPlane);

    const pair<double, double> pointOnTorus = {
      mod1_05(pointOnUnrolledTorus.first),
      mod1_05(pointOnUnrolledTorus.second)
    };

    const pair<double, double> pointOnPlaneNearestZero =
      transform2D(latticeBasisByModule[iModule], pointOnTorus);

    if (pow(pointOnPlaneNearestZero.first, 2) +
        pow(pointOnPlaneNearestZero.second, 2) > rSquared)
    {
      return false;
    }
  }

  return true;
}

bool lineSegmentIntersectsCircle2(pair<double, double> start,
                                  pair<double, double> end,
                                  pair<double, double> center,
                                  pair<double, double> unitVector,
                                  double lineLength,
                                  double rSquared)
{
  const double centerDistanceAlongLine =
    (unitVector.first*(center.first - start.first) +
     unitVector.second*(center.second - start.second));

  pair<double,double> nearestPointOnLine;
  if (centerDistanceAlongLine <= 0)
  {
    nearestPointOnLine = start;
  }
  else if (centerDistanceAlongLine < lineLength)
  {
    nearestPointOnLine = {
      start.first + unitVector.first * centerDistanceAlongLine,
      start.second + unitVector.second * centerDistanceAlongLine,
    };
  }
  else
  {
    nearestPointOnLine = end;
  }

  return (pow(center.first - nearestPointOnLine.first, 2) +
          pow(center.second - nearestPointOnLine.second, 2) <= rSquared);
}

bool lineSegmentIntersectsCircle(pair<double, double> start,
                                 pair<double, double> end,
                                 pair<double, double> center,
                                 double rSquared)
{
  pair<double,double> unitVector = {end.first - start.first,
                                    end.second - start.second};
  const double lineLength = sqrt(pow(unitVector.first, 2) +
                                 pow(unitVector.second, 2));
  unitVector.first /= lineLength;
  unitVector.second /= lineLength;

  return lineSegmentIntersectsCircle2(start, end, center, unitVector, lineLength,
                                      rSquared);
}

struct LineInfo2D {
  pair<double,double> unitVector;
  double length;
};

bool latticePointOverlapsShadow(
  pair<double, double> latticePoint,
  const vector<pair<double, double>>& convexHullUnshifted,
  const vector<LineInfo2D>& lines,
  pair<double, double> shift,
  size_t numDims,
  double rSquared)
{
  bool foundLatticeCollision = false;

  // If the circle around a lattice point overlaps the 2D shadow of this box,
  // then one of the following is true:
  //   A. It intersects one of the edges.
  //   B. It is totally contained within the shadow. (not possible in 1D)
  //
  // To detect possibility A, we check the distance between the lattice point
  // and each edge. If it's within the radius of the circle, then the two
  // intersect.
  //
  // To detect possibility B, we extend two rays leftward and rightward from the
  // lattice point. If both rays intersect an edge, then the lattice point is
  // contained within the shadow.

  bool leftRayCollided = false;
  bool rightRayCollided = false;

  for (size_t iPoint = 0; iPoint < convexHullUnshifted.size() - 1; iPoint++)
  {
    const pair<double,double> p1 = {
      convexHullUnshifted[iPoint].first + shift.first,
      convexHullUnshifted[iPoint].second + shift.second,
    };
    const pair<double,double> p2 = {
      convexHullUnshifted[iPoint+1].first + shift.first,
      convexHullUnshifted[iPoint+1].second + shift.second,
    };

    if (lineSegmentIntersectsCircle2(p1, p2, latticePoint,
                                     lines[iPoint].unitVector,
                                     lines[iPoint].length, rSquared))
    {
      foundLatticeCollision = true;
      break;
    }

    if (p1.second == p2.second)
    {
      // Special logic to avoid dividing by zero.
      if (p1.second == latticePoint.second)
      {
        // The ray passes straight through the edge.
        if (p1.first < latticePoint.first)
        {
          leftRayCollided = true;
        }
        else
        {
          rightRayCollided = true;
        }
      }
    }
    else if ((p1.second < latticePoint.second &&
              latticePoint.second < p2.second) ||
             (p2.second < latticePoint.second &&
              latticePoint.second < p1.second))
    {
      // Determine the point where it crosses.
      const double xCross = p1.first +
        (p2.first - p1.first) * ((latticePoint.second - p1.second) /
                                 (p2.second - p1.second));

      if (xCross < latticePoint.first)
      {
        leftRayCollided = true;
      }
      else
      {
        rightRayCollided = true;
      }
    }

    if (leftRayCollided && rightRayCollided)
    {
      foundLatticeCollision = true;
      break;
    }
  }

  return foundLatticeCollision;
}

vector<pair<double,double>> getShadowConvexHull(
  const vector<vector<double>>& domainToPlane,
  size_t numDims,
  const double dims[],
  double vertexBuffer[])
{
  if (numDims == 2)
  {
    // Optimization: in 2D we already know the convex hull.
    const double point1[2] = {0, 0};
    const double point2[2] = {0, dims[1]};
    const double point3[2] = {dims[0], dims[1]};
    const double point4[2] = {dims[0], 0};

    const pair<double,double> p1 =  transformND(domainToPlane, point1);

    return {p1,
            transformND(domainToPlane, point2),
            transformND(domainToPlane, point3),
            transformND(domainToPlane, point4),
            p1};
  }

  typedef boost::tuple<float, float> point;
  typedef bg::model::polygon<point> polygon;

  polygon poly;

  HyperrectangleVertexEnumerator vertices(dims, numDims);
  while (vertices.getNext(vertexBuffer))
  {
    const pair<double,double> p = transformND(domainToPlane, vertexBuffer);
    bg::append(poly, point(p.first, p.second));
  }

  polygon hull;
  bg::convex_hull(poly, hull);

  vector<pair<double, double>> convexHull;
  convexHull.reserve(hull.outer().size());

  for (auto it = hull.outer().begin(); it != hull.outer().end(); ++it)
  {
    convexHull.push_back({bg::get<0>(*it),
                          bg::get<1>(*it)});
  }

  return convexHull;
}

struct BoundingBox2D {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

/**
 * Quickly check whether this hyperrectangle excludes grid code zero
 * in any individual module.
 */
bool tryProveGridCodeZeroImpossible_1D(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<SquareMatrix2D<double>>& latticeBasisByModule,
  const vector<SquareMatrix2D<double>>& inverseLatticeBasisByModule,
  size_t numDims,
  const double x0[],
  const double dims[],
  double r,
  double rSquared)
{
  const double point1 = x0[0];
  const double point2 = x0[0] + dims[0];

  for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
  {
    const pair<double,double> p1 = transformND(domainToPlaneByModule[iModule],
                                               &point1);
    const pair<double,double> p2 = transformND(domainToPlaneByModule[iModule],
                                               &point2);

    // Figure out which lattice points we need to check.
    const double xmin = std::min(p1.first, p2.first);
    const double xmax = std::max(p1.first, p2.first);
    const double ymin = std::min(p1.second, p2.second);
    const double ymax = std::max(p1.second, p2.second);
    LatticePointEnumerator latticePoints(latticeBasisByModule[iModule],
                                         inverseLatticeBasisByModule[iModule],
                                         xmin, xmax, ymin, ymax, r, rSquared);

    pair<double, double> latticePoint;
    bool foundLatticeCollision = false;
    while (!foundLatticeCollision && latticePoints.getNext(&latticePoint))
    {
      foundLatticeCollision =
        lineSegmentIntersectsCircle(p1, p2, latticePoint, rSquared);
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

const double shadowWidthThreshold = 0.5;

BoundingBox2D computeBoundingBox(const vector<pair<double,double>>& shadow)
{
  BoundingBox2D boundingBox = {
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::lowest(),
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::lowest(),
  };

  for (const pair<double,double>& p : shadow)
  {
    boundingBox.xmin = std::min(boundingBox.xmin, p.first);
    boundingBox.xmax = std::max(boundingBox.xmax, p.first);
    boundingBox.ymin = std::min(boundingBox.ymin, p.second);
    boundingBox.ymax = std::max(boundingBox.ymax, p.second);
  }

  return boundingBox;
}

vector<LineInfo2D> computeShadowLines(const vector<pair<double,double>>& shadow)
{
  vector<LineInfo2D> lines;
  lines.reserve(shadow.size() - 1);

  for (size_t iPoint = 0; iPoint < shadow.size() - 1; iPoint++)
  {
    const pair<double,double> start = shadow[iPoint];
    const pair<double,double> end = shadow[iPoint+1];

    pair<double,double> unitVector = {end.first - start.first,
                                      end.second - start.second};
    const double lineLength = sqrt(pow(unitVector.first, 2) +
                                   pow(unitVector.second, 2));
    unitVector.first /= lineLength;
    unitVector.second /= lineLength;

    lines.push_back({unitVector, lineLength});
  }

  return lines;
}

LatticeBox computeLatticeBox(
  const BoundingBox2D& boundingBox,
  const SquareMatrix2D<double>& inverseLatticeBasis,
  double r)
{
  const double paddedLeft = boundingBox.xmin - r;
  const double paddedRight = boundingBox.xmax + r;
  const double paddedBottom = boundingBox.ymin - r;
  const double paddedTop = boundingBox.ymax + r;

  double xmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::lowest();

  pair<double, double> ij;
  ij = transform2D(inverseLatticeBasis, {paddedLeft, paddedBottom});
  xmin = std::min(xmin, ij.first);
  xmax = std::max(xmax, ij.first);
  ij = transform2D(inverseLatticeBasis, {paddedRight, paddedBottom});
  xmin = std::min(xmin, ij.first);
  xmax = std::max(xmax, ij.first);
  ij = transform2D(inverseLatticeBasis, {paddedLeft, paddedTop});
  xmin = std::min(xmin, ij.first);
  xmax = std::max(xmax, ij.first);
  ij = transform2D(inverseLatticeBasis, {paddedRight, paddedTop});
  xmin = std::min(xmin, ij.first);
  xmax = std::max(xmax, ij.first);

  return {xmin, xmax,
          transform2D(inverseLatticeBasis, {(paddedRight + paddedLeft) / 2,
                                            (paddedTop + paddedBottom) / 2})};
}

/**
 * Quickly check whether this hyperrectangle excludes grid code zero
 * in any individual module.
 */
bool tryProveGridCodeZeroImpossible(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<SquareMatrix2D<double>>& latticeBasisByModule,
  const vector<SquareMatrix2D<double>>& inverseLatticeBasisByModule,
  size_t numDims,
  const double x0[],
  const double dims[],
  double r,
  double rSquared,
  double vertexBuffer[],
  vector<vector<vector<pair<double,double>>>>& cachedShadows,
  vector<vector<vector<LineInfo2D>>>& cachedShadowLines,
  vector<vector<BoundingBox2D>>& cachedShadowBoundingBoxes,
  vector<vector<LatticeBox>>& cachedLatticeBoxes,
  size_t frameNumber)
{
  if (numDims == 1)
  {
    return tryProveGridCodeZeroImpossible_1D(
      domainToPlaneByModule, latticeBasisByModule, inverseLatticeBasisByModule,
      numDims, x0, dims, r, rSquared);
  }

  NTA_ASSERT(frameNumber <= cachedShadowBoundingBoxes.size());

  if (frameNumber == cachedShadowBoundingBoxes.size())
  {
    vector<vector<pair<double,double>>> shadowByModule;
    shadowByModule.reserve(domainToPlaneByModule.size());

    vector<BoundingBox2D> boundingBoxByModule;
    boundingBoxByModule.reserve(domainToPlaneByModule.size());

    vector<LatticeBox> latticeBoxByModule;
    latticeBoxByModule.reserve(domainToPlaneByModule.size());

    vector<vector<LineInfo2D>> linesByModule;
    linesByModule.reserve(domainToPlaneByModule.size());

    for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
    {
      const vector<pair<double, double>> shadow =
        getShadowConvexHull(domainToPlaneByModule[iModule], numDims, dims,
                            vertexBuffer);

      const BoundingBox2D boundingBox = computeBoundingBox(shadow);;
      boundingBoxByModule.push_back(boundingBox);

      latticeBoxByModule.push_back(
        computeLatticeBox(boundingBox, inverseLatticeBasisByModule[iModule],
                          r));

      if (boundingBox.xmax - boundingBox.xmin > shadowWidthThreshold ||
          boundingBox.ymax - boundingBox.ymin > shadowWidthThreshold)
      {
        shadowByModule.push_back({});
        linesByModule.push_back({});
      }
      else
      {
        shadowByModule.push_back(shadow);
        linesByModule.push_back(computeShadowLines(shadow));
      }
    }

    cachedShadows.push_back(shadowByModule);
    cachedShadowBoundingBoxes.push_back(boundingBoxByModule);
    cachedLatticeBoxes.push_back(latticeBoxByModule);
    cachedShadowLines.push_back(linesByModule);
  }

  for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
  {
    // Figure out which lattice points we need to check.
    const pair<double,double> shift =
      transformND( domainToPlaneByModule[iModule], x0);
    const BoundingBox2D& boundingBox =
      cachedShadowBoundingBoxes[frameNumber][iModule];
    const double xmin = boundingBox.xmin + shift.first;
    const double xmax = boundingBox.xmax + shift.first;
    const double ymin = boundingBox.ymin + shift.second;
    const double ymax = boundingBox.ymax + shift.second;

    LatticePointEnumerator latticePoints(
      latticeBasisByModule[iModule], inverseLatticeBasisByModule[iModule],
      cachedLatticeBoxes[frameNumber][iModule], shift, xmin, xmax, ymin, ymax,
      rSquared);

    pair<double, double> latticePoint;
    bool foundLatticeCollision = false;
    while (!foundLatticeCollision && latticePoints.getNext(&latticePoint))
    {
      // At this point, the bounding box collides with a lattice point, but it
      // might not collide with the actual polygon formed by the 2D shadow of
      // the box. We don't actually need to check whether it collides with the
      // polygon; we can just say that it might. If it doesn't, this function
      // will get called again with a smaller box, and eventually the space will
      // be broken into sufficiently small boxes that each have no overlap with
      // the lattice of at least one module.
      //
      // With high dimensional boxes, this approach can be slow. It checks a
      // large number of tiny polygons that are just outside the range of a
      // lattice point. It has to divide each of these tiny polygons until the
      // bounding box doesn't touch the lattice point. With a thorough approach,
      // it has to perform much fewer checks, but each check is slower.
      //
      // To get the best of both worlds, we do non-thorough checks when the
      // shadow is large, and begin doing thorough checks when the shadow is
      // small.
      if (xmax - xmin > shadowWidthThreshold ||
          ymax - ymin > shadowWidthThreshold)
      {
        // Rely on the bounding box check.
        foundLatticeCollision = true;
      }
      else
      {
        foundLatticeCollision =
          latticePointOverlapsShadow(latticePoint,
                                     cachedShadows[frameNumber][iModule],
                                     cachedShadowLines[frameNumber][iModule],
                                     shift, numDims, rSquared);
      }
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
  const vector<SquareMatrix2D<double>>& latticeBasisByModule,
  const vector<SquareMatrix2D<double>>& inverseLatticeBasisByModule,
  size_t numDims,
  double x0[],
  double dims[],
  double r,
  double rSquaredPositive,
  double rSquaredNegative,
  double vertexBuffer[],
  vector<vector<vector<pair<double,double>>>>& cachedShadows,
  vector<vector<vector<LineInfo2D>>>& cachedShadowLines,
  vector<vector<BoundingBox2D>>& cachedShadowBoundingBoxes,
  vector<vector<LatticeBox>>& cachedLatticeBoxes,
  size_t frameNumber,
  std::atomic<bool>& shouldContinue)
{
  if (!shouldContinue)
  {
    return false;
  }

  if (tryProveGridCodeZeroImpossible(domainToPlaneByModule,
                                     latticeBasisByModule,
                                     inverseLatticeBasisByModule, numDims, x0,
                                     dims, r, rSquaredNegative, vertexBuffer,
                                     cachedShadows, cachedShadowLines,
                                     cachedShadowBoundingBoxes,
                                     cachedLatticeBoxes,
                                     frameNumber))
  {
    return false;
  }

  if (tryFindGridCodeZero(domainToPlaneByModule, latticeBasisByModule,
                          inverseLatticeBasisByModule, numDims, x0, dims,
                          rSquaredPositive, vertexBuffer))
  {
    return true;
  }

  size_t iWidestDim = std::distance(dims,
                                    std::max_element(dims, dims + numDims));
  {
    SwapValueRAII swap1(&dims[iWidestDim], dims[iWidestDim] / 2);
    if (findGridCodeZeroHelper(
          domainToPlaneByModule, latticeBasisByModule,
          inverseLatticeBasisByModule, numDims, x0, dims, r, rSquaredPositive,
          rSquaredNegative, vertexBuffer, cachedShadows, cachedShadowLines,
          cachedShadowBoundingBoxes, cachedLatticeBoxes,
          frameNumber + 1, shouldContinue))
    {
      return true;
    }

    {
      SwapValueRAII swap2(&x0[iWidestDim], x0[iWidestDim] + dims[iWidestDim]);
      return findGridCodeZeroHelper(
        domainToPlaneByModule, latticeBasisByModule,
        inverseLatticeBasisByModule, numDims, x0, dims, r, rSquaredPositive,
        rSquaredNegative, vertexBuffer, cachedShadows, cachedShadowLines,
        cachedShadowBoundingBoxes, cachedLatticeBoxes,
        frameNumber + 1, shouldContinue);
    }
  }
}

struct GridUniquenessState {
  // Constants (thread-safe)
  const vector<vector<vector<double>>>& domainToPlaneByModule;
  const vector<SquareMatrix2D<double>>& latticeBasisByModule;
  const vector<SquareMatrix2D<double>>& inverseLatticeBasisByModule;
  const double readoutResolution;
  const double meanScaleEstimate;
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
  std::condition_variable& finishedCondition;
  bool finished;
  size_t numActiveThreads;
  vector<double> threadBaselineRadius;
  vector<vector<double>> threadQueryX0;
  vector<vector<double>> threadQueryDims;
  vector<std::atomic<bool>> threadShouldContinue;
  std::atomic<bool>& quitting;
  vector<bool> threadRunning;
};

template<typename T>
std::string vecs(const vector<T>& v)
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

  vector<long long> numBinsByDim(state.numDims);

  // Add a small epsilon to handle situations where floating point math causes
  // a vertex to be non-zero-overlapping here and zero-overlapping in
  // tryProveGridCodeZeroImpossible. With this addition, anything
  // zero-overlapping in tryProveGridCodeZeroImpossible is guaranteed to be
  // zero-overlapping here, so the program won't get caught in infinite
  // recursion.
  const double rSquaredPositive = pow(state.readoutResolution/2 + 0.000000001, 2);
  const double rSquaredNegative = pow(state.readoutResolution/2, 2);

  while (!state.quitting)
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

    vector<vector<vector<pair<double,double>>>> cachedShadows;
    vector<vector<vector<LineInfo2D>>> cachedShadowLines;
    vector<vector<BoundingBox2D>> cachedShadowBoundingBoxes;
    vector<vector<LatticeBox>> cachedLatticeBoxes;

    // Optimization: if the box is large, break it into small chunks rather than
    // relying completely on the divide-and-conquer to break into
    // reasonable-sized chunks.

    // Use a longer bin size for 1D. A 1D slice of a 2D plane can be relatively
    // long before it has high probability of colliding with a lattice point in
    // every module.
    const double scalesPerBin = (state.numDims == 1)
      ? 2.5
      : 0.55;

    for (size_t iDim = 0; iDim < state.numDims; iDim++)
    {
      numBinsByDim[iDim] = ceil(dims[iDim] / (scalesPerBin *
                                              state.meanScaleEstimate));
      dims[iDim] /= numBinsByDim[iDim];
    }

    const vector<double>& x0_orig = state.threadQueryX0[iThread];
    vector<long long> currentBinByDim(state.numDims, 0);
    while (state.threadShouldContinue[iThread])
    {
      for (size_t iDim = 0; iDim < state.numDims; iDim++)
      {
        x0[iDim] = x0_orig[iDim] + currentBinByDim[iDim]*dims[iDim];
      }

      foundGridCodeZero = findGridCodeZeroHelper(
        state.domainToPlaneByModule, state.latticeBasisByModule,
        state.inverseLatticeBasisByModule, state.numDims, x0.data(),
        dims.data(), state.readoutResolution/2, rSquaredPositive,
        rSquaredNegative, pointWithGridCodeZero.data(), cachedShadows,
        cachedShadowLines, cachedShadowBoundingBoxes, cachedLatticeBoxes, 0,
        state.threadShouldContinue[iThread]);

      if (foundGridCodeZero) break;

      // Increment as little endian arithmetic with a varying base.
      bool overflow = true;
      for (size_t iDigit = 0; iDigit < state.numDims; iDigit++)
      {
        overflow = ++currentBinByDim[iDigit] == numBinsByDim[iDigit];
        if (!overflow) break;
        currentBinByDim[iDigit] = 0;
      }

      if (overflow) break;
    }
  }

  // This thread is exiting.
  {
    std::lock_guard<std::mutex> lock(state.mutex);
    if (--state.numActiveThreads == 0)
    {
      state.finished = true;
      state.finishedCondition.notify_all();
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
    double dLongest = std::numeric_limits<double>::lowest();
    for (size_t iColumn = 0; iColumn < domainToPlane[0].size(); iColumn++)
    {
      double length = sqrt(pow(domainToPlane[0][iColumn], 2) +
                           pow(domainToPlane[1][iColumn], 2));
      if (length > dLongest)
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

  vector<vector<vector<double>>> domainToPlaneByModule2(domainToPlaneByModule);
  vector<vector<vector<double>>> latticeBasisByModule2(latticeBasisByModule);
  optimizeMatrices(&domainToPlaneByModule2, &latticeBasisByModule2);

  vector<SquareMatrix2D<double>> latticeBasisByModule3;
  vector<SquareMatrix2D<double>> inverseLatticeBasisByModule;
  for (const vector<vector<double>>& latticeBasis : latticeBasisByModule2)
  {
    latticeBasisByModule3.push_back({
        latticeBasis[0][0], latticeBasis[0][1],
        latticeBasis[1][0], latticeBasis[1][1]});
    inverseLatticeBasisByModule.push_back(invert2DMatrix(latticeBasis));
  }

  vector<vector<vector<pair<double,double>>>> cachedShadows;
  vector<vector<vector<LineInfo2D>>> cachedShadowLines;
  vector<vector<BoundingBox2D>> cachedShadowBoundingBoxes;
  vector<vector<LatticeBox>> cachedLatticeBoxes;

  // Add a small epsilon to handle situations where floating point math causes a
  // vertex to be non-zero-overlapping here and zero-overlapping in
  // tryProveGridCodeZeroImpossible. With this addition, anything
  // zero-overlapping in tryProveGridCodeZeroImpossible is guaranteed to be
  // zero-overlapping here, so the program won't get caught in infinite
  // recursion.
  const double rSquaredPositive = pow(readoutResolution/2 + 0.000000001, 2);
  const double rSquaredNegative = pow(readoutResolution/2, 2);

  return findGridCodeZeroHelper(
    domainToPlaneByModule2, latticeBasisByModule3, inverseLatticeBasisByModule,
    dimsCopy.size(), x0Copy.data(), dimsCopy.data(), readoutResolution/2,
    rSquaredPositive, rSquaredNegative, pointWithGridCodeZero->data(),
    cachedShadows, cachedShadowLines, cachedShadowBoundingBoxes,
    cachedLatticeBoxes, 0, shouldContinue);
}

pair<double,vector<double>>
nupic::experimental::grid_uniqueness::computeGridUniquenessHypercube(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<vector<vector<double>>>& latticeBasisByModule,
  double readoutResolution,
  double ignoredCenterDiameter)
{
  typedef std::chrono::steady_clock Clock;

  enum ExitReason {
    Timeout,
    Interrupt,
    Completed
  };

  std::atomic<ExitReason> exitReason(ExitReason::Completed);

  ThreadSafeQueue<Message> messages;
  CaptureInterruptsRAII captureInterrupts(&messages);

  std::atomic<bool> quitting(false);
  std::thread messageThread(
    [&]() {
      while (true)
      {
        switch (messages.take())
        {
          case Message::Interrupt:
            quitting = true;
            exitReason = ExitReason::Interrupt;
            break;
          case Message::Timeout:
            quitting = true;
            exitReason = ExitReason::Timeout;
            break;
          case Message::Exiting:
            return;
        }
      }
    });

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

  vector<SquareMatrix2D<double>> latticeBasisByModule3;
  vector<SquareMatrix2D<double>> inverseLatticeBasisByModule;
  for (const vector<vector<double>>& latticeBasis : latticeBasisByModule2)
  {
    latticeBasisByModule3.push_back({
        latticeBasis[0][0], latticeBasis[0][1],
        latticeBasis[1][0], latticeBasis[1][1]});
    inverseLatticeBasisByModule.push_back(invert2DMatrix(latticeBasis));
  }

  double meanScaleEstimate = 0.0;

  for (const vector<vector<double>>& domainToPlane : domainToPlaneByModule2)
  {
    double longestDisplacementSquared = std::numeric_limits<double>::min();

    for (size_t iDim = 0; iDim < numDims; iDim++)
    {
      longestDisplacementSquared = std::max(longestDisplacementSquared,
                                            pow(domainToPlane[0][iDim], 2) +
                                            pow(domainToPlane[1][iDim], 2));
    }

    const double scaleEstimate = 1 / sqrt(longestDisplacementSquared);
    meanScaleEstimate += scaleEstimate;
  }
  meanScaleEstimate /= domainToPlaneByModule2.size();

  // Use condition_variables to enable periodic logging while waiting for the
  // threads to finish.
  std::mutex stateMutex;
  std::condition_variable finishedCondition;

  size_t numThreads = std::thread::hardware_concurrency();

  GridUniquenessState state = {
    domainToPlaneByModule2,
    latticeBasisByModule3,
    inverseLatticeBasisByModule,
    readoutResolution,

    meanScaleEstimate,
    numDims,

    ignoredCenterDiameter,
    ignoredCenterDiameter * 1.01,
    vector<double>(numDims, ignoredCenterDiameter),
    0,
    true,
    true,

    vector<double>(numDims),
    std::numeric_limits<double>::max(),

    stateMutex,
    finishedCondition,
    false,
    0,
    vector<double>(numThreads, std::numeric_limits<double>::max()),
    vector<vector<double>>(numThreads, vector<double>(numDims)),
    vector<vector<double>>(numThreads, vector<double>(numDims)),
    vector<std::atomic<bool>>(numThreads),
    quitting,
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

    bool printedInitialStatement = false;

    while (!state.finished)
    {
      if (state.finishedCondition.wait_until(
            lock, tNextPrint) == std::cv_status::timeout)
      {
        if (!printedInitialStatement)
        {
          {
            std::ostringstream oss;
            oss << "[";

            for (size_t iModule = 0;
                 iModule < domainToPlaneByModule.size();
                 iModule++)
            {
              oss << "[";
              oss << vecs(domainToPlaneByModule[iModule][0]) << ",";
              oss << vecs(domainToPlaneByModule[iModule][1]);
              oss << "],";
            }
            oss << "]" << std::endl;
            NTA_INFO << "domainToPlaneByModule:" << std::endl << oss.str();
          }

          {
            std::ostringstream oss;
            oss << "[";
            for (size_t iModule = 0;
                 iModule < latticeBasisByModule.size();
                 iModule++)
            {
              oss << "[";
              oss << vecs(latticeBasisByModule[iModule][0]) << ",";
              oss << vecs(latticeBasisByModule[iModule][1]);
              oss << "],";
            }
            oss << "]" << std::endl;

            NTA_INFO << "latticeBasisByModule:" << std::endl << oss.str();
          }

          NTA_INFO << "readout resolution: " << readoutResolution;

          printedInitialStatement = true;
        }

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
    }
  }

  messages.put(Message::Exiting);
  messageThread.join();

  switch (exitReason.load())
  {
    case ExitReason::Timeout:
      // Python code may check for the precise string "timeout".
      NTA_THROW << "timeout";
    case ExitReason::Interrupt:
      NTA_THROW << "interrupt";
    case ExitReason::Completed:
    default:
      return {state.foundPointBaselineRadius, state.pointWithGridCodeZero};
  }
}

bool tryFindGridCodeZero_noModulo(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  size_t numDims,
  const double x0[],
  const double dims[],
  double rSquared,
  double vertexBuffer[])
{
  for (size_t iDim = 0; iDim < numDims; iDim++)
  {
    vertexBuffer[iDim] = x0[iDim] + (dims[iDim]/2);
  }

  for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
  {
    const pair<double, double> pointOnPlane =
      transformND(domainToPlaneByModule[iModule], vertexBuffer);

    if (pow(pointOnPlane.first, 2) + pow(pointOnPlane.second, 2) > rSquared)
    {
      return false;
    }
  }

  return true;
}

bool tryProveGridCodeZeroImpossible_noModulo_1D(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  size_t numDims,
  const double x0[],
  const double dims[],
  double rSquared)
{
  const double point1 = x0[0];
  const double point2 = x0[0] + dims[0];

  for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
  {
    const pair<double,double> p1 = transformND(domainToPlaneByModule[iModule],
                                               &point1);
    const pair<double,double> p2 = transformND(domainToPlaneByModule[iModule],
                                               &point2);

    if (!lineSegmentIntersectsCircle(p1, p2, {0.0, 0.0}, rSquared))
    {
      // This module never gets near grid code zero for the provided range of
      // locations. So this range can't possibly contain grid code zero.
      return true;
    }
  }

  return false;
}

bool tryProveGridCodeZeroImpossible_noModulo(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  size_t numDims,
  const double x0[],
  const double dims[],
  double r,
  double rSquared,
  double vertexBuffer[],
  vector<vector<vector<pair<double,double>>>>& cachedShadows,
  vector<vector<vector<LineInfo2D>>>& cachedShadowLines,
  size_t frameNumber)
{
  if (numDims == 1)
  {
    return tryProveGridCodeZeroImpossible_noModulo_1D(
      domainToPlaneByModule, numDims, x0, dims, rSquared);
  }

  NTA_ASSERT(frameNumber <= cachedShadows.size());

  if (frameNumber == cachedShadows.size())
  {
    vector<vector<pair<double,double>>> shadowByModule;
    shadowByModule.reserve(domainToPlaneByModule.size());

    vector<vector<LineInfo2D>> linesByModule;
    linesByModule.reserve(domainToPlaneByModule.size());

    for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
    {
      const vector<pair<double, double>> shadow = getShadowConvexHull(
        domainToPlaneByModule[iModule], numDims, dims, vertexBuffer);
      shadowByModule.push_back(shadow);
      linesByModule.push_back(computeShadowLines(shadow));
    }

    cachedShadows.push_back(shadowByModule);
    cachedShadowLines.push_back(linesByModule);
  }

  for (size_t iModule = 0; iModule < domainToPlaneByModule.size(); iModule++)
  {
    const pair<double,double> shift =
      transformND(domainToPlaneByModule[iModule], x0);

    if (!latticePointOverlapsShadow({0.0, 0.0},
                                    cachedShadows[frameNumber][iModule],
                                    cachedShadowLines[frameNumber][iModule],
                                    shift, numDims, rSquared))
    {
      // This module never gets near grid code zero for the provided range of
      // locations. So this range can't possibly contain grid code zero.
      return true;
    }
  }

  return false;
}

bool findGridCodeZeroHelper_noModulo(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  size_t numDims,
  double x0[],
  double dims[],
  double r,
  double rSquaredPositive,
  double rSquaredNegative,
  double vertexBuffer[],
  vector<vector<vector<pair<double,double>>>>& cachedShadows,
  vector<vector<vector<LineInfo2D>>>& cachedShadowLines,
  size_t frameNumber,
  std::atomic<bool>& shouldContinue)
{
  if (!shouldContinue)
  {
    return false;
  }

  if (tryProveGridCodeZeroImpossible_noModulo(
        domainToPlaneByModule, numDims, x0, dims, r, rSquaredNegative,
        vertexBuffer, cachedShadows, cachedShadowLines, frameNumber))
  {
    return false;
  }

  if (tryFindGridCodeZero_noModulo(domainToPlaneByModule, numDims, x0, dims,
                                   rSquaredPositive, vertexBuffer))
  {
    return true;
  }

  size_t iWidestDim = std::distance(dims,
                                    std::max_element(dims, dims + numDims));
  {
    SwapValueRAII swap1(&dims[iWidestDim], dims[iWidestDim] / 2);
    if (findGridCodeZeroHelper_noModulo(
          domainToPlaneByModule, numDims, x0, dims, r, rSquaredPositive,
          rSquaredNegative, vertexBuffer, cachedShadows, cachedShadowLines,
          frameNumber + 1, shouldContinue))
    {
      return true;
    }

    {
      SwapValueRAII swap2(&x0[iWidestDim], x0[iWidestDim] + dims[iWidestDim]);
      return findGridCodeZeroHelper_noModulo(
        domainToPlaneByModule, numDims, x0, dims, r, rSquaredPositive,
        rSquaredNegative, vertexBuffer, cachedShadows, cachedShadowLines,
        frameNumber + 1, shouldContinue);
    }
  }
}

bool findGridCodeZero_noModulo(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  const vector<double>& x0,
  const vector<double>& dims,
  double readoutResolution,
  std::atomic<bool>& shouldContinue,
  vector<double>* pointWithGridCodeZero = nullptr)
{
  // Avoid doing any allocations in each recursion.
  vector<double> x0Copy(x0);
  vector<double> dimsCopy(dims);

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

  vector<vector<vector<pair<double,double>>>> cachedShadows;
  vector<vector<vector<LineInfo2D>>> cachedShadowLines;

  // Add a small epsilon to handle situations where floating point math causes a
  // vertex to be non-zero-overlapping here and zero-overlapping in
  // tryProveGridCodeZeroImpossible. With this addition, anything
  // zero-overlapping in tryProveGridCodeZeroImpossible is guaranteed to be
  // zero-overlapping here, so the program won't get caught in infinite
  // recursion.
  const double rSquaredPositive = pow(readoutResolution/2 + 0.000000001, 2);
  const double rSquaredNegative = pow(readoutResolution/2, 2);

  return findGridCodeZeroHelper_noModulo(
    domainToPlaneByModule, dimsCopy.size(), x0Copy.data(), dimsCopy.data(),
    readoutResolution/2, rSquaredPositive, rSquaredNegative,
    pointWithGridCodeZero->data(), cachedShadows, cachedShadowLines, 0,
    shouldContinue);
}

bool findGridCodeZeroAtRadius(
  double radius,
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  double readoutResolution,
  std::atomic<bool>& shouldContinue)
{
  const size_t numDims = domainToPlaneByModule[0][0].size();

  for (size_t iDim = 0; iDim < numDims; ++iDim)
  {
    // Test the hyperplanes formed by setting this dimension to r and -r.
    vector<double> x0(numDims, -radius);
    vector<double> dims(numDims, 2*radius);

    // Optimization: for the final dimension, don't go negative. Half of the
    // hypercube will be equal-and-opposite phases of the other half, so we
    // ignore the lower half of the final dimension.
    x0[numDims - 1] = 0;
    dims[numDims - 1] = radius;

    dims[iDim] = 0;

    // Test -r
    if (iDim != numDims - 1)
    {
      if (findGridCodeZero_noModulo(domainToPlaneByModule,
                                    x0, dims, readoutResolution,
                                    shouldContinue))
      {
        return true;
      }
    }

    // Test +r
    x0[iDim] = radius;
    if (findGridCodeZero_noModulo(domainToPlaneByModule,
                                  x0, dims, readoutResolution,
                                  shouldContinue))
    {
      return true;
    }
  }

  return false;
}

double
nupic::experimental::grid_uniqueness::computeBinSidelength(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  double readoutResolution,
  double resultPrecision,
  double upperBound,
  double timeout)
{
  typedef std::chrono::steady_clock Clock;

  enum ExitReason {
    Timeout,
    Interrupt,
    Completed
  };

  std::atomic<ExitReason> exitReason(ExitReason::Completed);

  ThreadSafeQueue<Message> messages;
  CaptureInterruptsRAII captureInterrupts(&messages);

  std::atomic<bool> shouldContinue(true);
  std::atomic<bool> interrupted(false);
  std::thread messageThread(
    [&]() {
      while (true)
      {
        switch (messages.take())
        {
          case Message::Interrupt:
            shouldContinue = false;
            exitReason = ExitReason::Interrupt;
            interrupted = true;
            break;
          case Message::Timeout:
            shouldContinue = false;
            exitReason = ExitReason::Timeout;
            break;
          case Message::Exiting:
            return;
        }
      }
    });

  ScheduledTask* scheduledTask = nullptr;
  if (timeout > 0)
  {
    scheduledTask = new ScheduledTask(
      Clock::now() + std::chrono::duration<double>(timeout),
      [&messages](){
        messages.put(Message::Timeout);
      });
  }

  double tested = 0;
  double radius = 0.5;

  while (findGridCodeZeroAtRadius(radius,
                                  domainToPlaneByModule,
                                  readoutResolution,
                                  shouldContinue))
  {
    tested = radius;
    radius *= 2;

    if (radius > upperBound)
    {
      return -1.0;
    }
  }

  // The radius needs to be twice as precise to get the sidelength sufficiently
  // precise.
  const double resultPrecision2 = resultPrecision / 2;

  double dec = (radius - tested) / 2;

  // The possible error is equal to dec*2.
  while (shouldContinue && dec*2 > resultPrecision2)
  {
    const double testRadius = radius - dec;

    if (!findGridCodeZeroAtRadius(testRadius,
                                  domainToPlaneByModule,
                                  readoutResolution,
                                  shouldContinue))
    {
      radius = testRadius;
    }

    dec /= 2;
  }

  if (scheduledTask != nullptr)
  {
    delete scheduledTask;
    scheduledTask = nullptr;
  }

  messages.put(Message::Exiting);
  messageThread.join();

  switch (exitReason.load())
  {
    case ExitReason::Timeout:
      // Python code may check for the precise string "timeout".
      NTA_THROW << "timeout";
    case ExitReason::Interrupt:
      NTA_THROW << "interrupt";
    case ExitReason::Completed:
    default:
      return 2*radius;
  }
}

vector<double>
nupic::experimental::grid_uniqueness::computeBinRectangle(
  const vector<vector<vector<double>>>& domainToPlaneByModule,
  double readoutResolution,
  double resultPrecision,
  double upperBound,
  double timeout)
{
  typedef std::chrono::steady_clock Clock;

  enum ExitReason {
    Timeout,
    Interrupt,
    Completed
  };

  std::atomic<ExitReason> exitReason(ExitReason::Completed);

  ThreadSafeQueue<Message> messages;
  CaptureInterruptsRAII captureInterrupts(&messages);

  std::atomic<bool> shouldContinue(true);
  std::atomic<bool> interrupted(false);
  std::thread messageThread(
    [&]() {
      while (true)
      {
        switch (messages.take())
        {
          case Message::Interrupt:
            shouldContinue = false;
            exitReason = ExitReason::Interrupt;
            interrupted = true;
            break;
          case Message::Timeout:
            shouldContinue = false;
            exitReason = ExitReason::Timeout;
            break;
          case Message::Exiting:
            return;
        }
      }
    });

  ScheduledTask* scheduledTask = nullptr;
  if (timeout > 0)
  {
    scheduledTask = new ScheduledTask(
      Clock::now() + std::chrono::duration<double>(timeout),
      [&messages](){ messages.put(Message::Timeout); });
  }

  const size_t numDims = domainToPlaneByModule[0][0].size();

  double radius = 0.5;

  while (findGridCodeZeroAtRadius(radius,
                                  domainToPlaneByModule,
                                  readoutResolution,
                                  shouldContinue))
  {
    radius *= 2;

    if (radius > upperBound)
    {
      return {};
    }
  }

  // The radius needs to be twice as precise to get the sidelength sufficiently
  // precise.
  const double resultPrecision2 = resultPrecision / 2;

  vector<double> radii(numDims, radius);

  for (size_t iDim = 0; iDim < numDims; ++iDim)
  {
    double dec = radius / 2;

    // The possible error is equal to dec*2.
    while (shouldContinue && dec*2 > resultPrecision2)
    {
      const double testRadius = radii[iDim] - dec;

      vector<double> x0(numDims);
      vector<double> dims(numDims);

      for (size_t iDim2 = 0; iDim2 < numDims; ++iDim2)
      {
        if (iDim2 == numDims - 1)
        {
          // Optimization: for the final dimension, don't go negative. Half of the
          // hypercube will be equal-and-opposite phases of the other half, so we
          // ignore the lower half of the final dimension.
          x0[iDim2] = 0;
          dims[iDim2] = radii[iDim2];
        }
        else
        {
          x0[iDim2] = -radii[iDim2];
          dims[iDim2] = 2*radii[iDim2];
        }
      }

      dims[iDim] = 0;

      bool foundZero = false;
      if (iDim != numDims - 1)
      {
        // Test -r
        x0[iDim] = -testRadius;
        foundZero = findGridCodeZero_noModulo(domainToPlaneByModule,
                                              x0, dims, readoutResolution,
                                              shouldContinue);
      }

      if (!foundZero)
      {
        // Test r
        x0[iDim] = testRadius;;
        foundZero = findGridCodeZero_noModulo(domainToPlaneByModule,
                                              x0, dims, readoutResolution,
                                              shouldContinue);;
      }

      if (!foundZero)
      {
        radii[iDim] = testRadius;
      }
      dec /= 2;
    }
  }

  if (scheduledTask != nullptr)
  {
    delete scheduledTask;
    scheduledTask = nullptr;
  }

  messages.put(Message::Exiting);
  messageThread.join();

  switch (exitReason.load())
  {
    case ExitReason::Timeout:
      // Python code may check for the precise string "timeout".
      NTA_THROW << "timeout";
    case ExitReason::Interrupt:
      NTA_THROW << "interrupt";
    case ExitReason::Completed:
    default:
    {
      vector<double> sidelengths(numDims);
      for (size_t iDim = 0; iDim < numDims; ++iDim)
      {
        sidelengths[iDim] = radii[iDim]*2;
      }

      return sidelengths;
    }
  }
}
