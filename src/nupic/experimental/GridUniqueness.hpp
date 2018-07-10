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
       * Determine whether any displacements in a k-dimensional rectangle would
       * result in an unmodified grid cell representation.
       *
       * Rather than looking for collisions between individual grid cell
       * representations, this looks for displacements don't change the grid
       * cell representation -- it looks for "grid displacement zero". For
       * example, if this function says there's a grid displacement zero at
       * displacement (42.0, 42.0), that implies that at any location, if you
       * move up 42.0 units and right 42.0 units, your grid cells will activate
       * the same location representation as they activated where you
       * started. On the other hand, if there's not a grid displacement zero at
       * a particular displacement, then that displacement will always result in
       * a different representation from where you started.
       *
       * This function uses a recursive "divide and conquer" algorithm. It first
       * quickly samples a few displacements to see if any of them have grid
       * displacement zero. Then it tries to prove that grid displacement zero
       * can't possibly happen in this hyperrectangle by showing that at least
       * one module never has displacement zero for any displacement in this
       * hyperrectangle. Finally, if neither attempt succeeds, it divides the
       * hyperrectangle in half, and then it tries again on both halves.
       *
       * @param A
       * A vector of matrices, one per grid cell module. Given that the grid
       * cells are encoding k-dimensional variables, each matrix should have 2
       * rows and k columns. Each column represents a dimension, describing how
       * a movement along that dimension should change the active phase of the
       * 2D module.
       *
       * @param x0
       * The lowest corner of the k-dimensional rectangle that will be searched.
       *
       * @param dims
       * The dimensions of the k-dimensional rectangle that will be searched.
       *
       * @param phaseResolution
       * The precision of phase readout of this grid code. For example, if this
       * is 0.2, then all phases are indistinguishable from those in their
       * surrounding +- 0.1 range.
       *
       * @param displacementWithGridCodeZero
       * Output parameter. A physical displacement with grid code zero in this
       * hyperrectangle. Only populated if the function returns true.
       *
       * @return
       * true if grid displacement zero is found, false otherwise.
       */
      bool findGridDisplacementZero(
        const std::vector<std::vector<std::vector<Real64>>>& A,
        const std::vector<Real64>& x0,
        const std::vector<Real64>& dims,
        Real64 phaseResolution,
        std::vector<Real64>& displacementWithGridCodeZero);

      /**
       * Given a set of grid cell module parameters, determines the diameter of
       * the k-dimensional cube in which every location has a unique grid cell
       * representation.
       *
       * This function iteratively builds out a larger hypercube out of smaller
       * hyperrectangles, relying on findGridDisplacementZero to analyze each
       * hyperrectangle. The hypercube expands outward from the origin, forming
       * a larger hypercube that contains positive and negative numbers. After
       * the nearest grid displacement zero is found, we divide the hypercube in
       * half to get a hypercube in which every location is guaranteed to be
       * unique (rather than just being different from the starting
       * location). The returned diameter is the diameter of this smaller
       * hypercube.
       *
       * @param A
       * A vector of matrices, one per grid cell module. Given that the grid
       * cells are encoding k-dimensional variables, each matrix should have 2
       * rows and k columns. Each column represents a dimension, describing how
       * a movement along that dimension should change the active phase of the
       * 2D module.
       *
       * @param phaseResolution
       * The precision of phase readout of this grid code. For example, if this
       * is 0.2, then all phases are indistinguishable from those in their
       * surrounding +- 0.1 range.
       *
       * @param ignoredCenterDiameter
       * The diameter of the hypercube at the center which should be ignored,
       * because it contains the *actual* displacement zero. Set this to be
       * sufficiently large to get away from this actual zero.
       *
       * @return
       * - The diameter of the hypercube that contains no collisions.
       * - A displacement just outside this hypercube that causes a collision.
       */
      std::pair<Real64,std::vector<Real64>> computeGridUniquenessHypercube(
        const std::vector<std::vector<std::vector<Real64>>>& A,
        Real64 phaseResolution,
        Real64 ignoredCenterDiameter);
    }
  }
}

#endif // NTA_GRID_UNIQUENESS_HPP
