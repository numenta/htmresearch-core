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
 * Functions that analyze the uniqueness of locations in grid cell spaces
 */

#ifndef NTA_GRID_UNIQUENESS_HPP
#define NTA_GRID_UNIQUENESS_HPP

#include <nupic/types/Types.hpp>
#include <vector>
#include <utility>

namespace nupic {
  namespace experimental {
    namespace grid_uniqueness {

      /**
       * Determine whether any points in a k-dimensional rectangle have a grid
       * code equal to the grid code at the origin.
       *
       * If this function says there's a grid code zero at point (42.0, 42.0),
       * that implies that every location has the same grid code as the location
       * that is 42.0 units up and 42.0 units right.
       *
       * This function uses a recursive "divide and conquer" algorithm. It first
       * quickly samples a few points to see if any of them have grid code zero.
       * Then it tries to prove that grid code zero can't possibly happen in
       * this hyperrectangle by showing that at least one module never has grid
       * code zero for any point in this hyperrectangle. Finally, if neither
       * attempt succeeds, it divides the hyperrectangle in half, and then it
       * tries again on both halves.
       *
       * @param x0
       * The lowest corner of the k-dimensional rectangle that will be searched.
       *
       * @param dims
       * The dimensions of the k-dimensional rectangle that will be searched.
       *
       * @param pointWithGridCodeZero
       * Output parameter. A point with grid code zero in this hyperrectangle.
       * Only populated if the function returns true.
       *
       * @return
       * true if grid code zero is found, false otherwise.
       */
      bool findGridCodeZero(
        const std::vector<std::vector<std::vector<Real64>>>& domainToPlaneByModule,
        const std::vector<std::vector<std::vector<Real64>>>& latticeBasisByModule,
        const std::vector<Real64>& x0,
        const std::vector<Real64>& dims,
        Real64 readoutResolution,
        std::vector<Real64>* pointWithGridCodeZero = nullptr);

      /**
       * Given a set of grid cell module parameters, determines the diameter of
       * the k-dimensional cube in which every location has a unique grid cell
       * representation.
       *
       * This function iteratively builds out a larger hypercube out of smaller
       * hyperrectangles, relying on findGridCodeZero to analyze each
       * hyperrectangle. The hypercube expands outward from the origin, forming
       * a larger hypercube that contains positive and negative numbers. After
       * the nearest grid code zero is found, we divide the hypercube in half to
       * get a hypercube in which every location is guaranteed to be unique
       * (rather than just being different from the starting location). The
       * returned diameter is the diameter of this smaller hypercube.
       *
       * Each grid cell module is specified by a pair of matrices. The first one
       * projects a k-dimensional point to a 2D plane, and the second matrix
       * specifies the basis of the grid lattice on this plane, specifying which
       * points on the plane have equivalent grid codes. The distance between
       * two grid codes is equal to the shortest distance between them on the
       * plane. The readout resolution parameter is measured in units of this
       * distance.
       *
       * There's not a strictly "correct" way to configure these matrices and
       * the readout resolution, but a typical way is to use unit vectors as
       * lattice basis vectors, use a plane projection that normalizes distances
       * so that 1 unit is 1 "scale", and then set the readout resolution to
       * some number in the range (0, 1) so that it is measured in units of
       * "scales".
       *
       * @param domainToPlaneByModule
       * A list of 2*k matrices, one per module. The matrix converts from a
       * point in the domain to a point on a plane, normalizing for grid cell
       * scale.
       *
       * @param latticeBasisByModule
       * A list of m 2*2 matrices, one per module. This matrix contains the
       * basis vectors for a lattice, specifying which points on the plane have
       * equivalent location representations in this module.
       *
       * @param readoutResolution
       * The precision of readout of this grid code, measured in distance on the
       * plane. For example, if this is 0.2, then all points on the plane are
       * indistinguishable from those in their surrounding +- 0.1 range.
       *
       * @param ignoredCenterDiameter
       * The diameter of the hypercube at the center which should be ignored,
       * because it contains the *actual* grid code zero. Set this to be
       * sufficiently large to get away from this actual zero.
       *
       * @return
       * - The diameter of the hypercube that contains no collisions.
       * - A point just outside this hypercube that collides with the origin.
       */
      std::pair<Real64,std::vector<Real64>> computeGridUniquenessHypercube(
        const std::vector<std::vector<std::vector<Real64>>>& domainToPlaneByModule,
        const std::vector<std::vector<std::vector<Real64>>>& latticeBasisByModule,
        Real64 readoutResolution,
        Real64 ignoredCenterDiameter);
    }
  }
}

#endif // NTA_GRID_UNIQUENESS_HPP
