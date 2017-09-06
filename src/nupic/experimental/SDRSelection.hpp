/* ---------------------------------------------------------------------
 * Numenta Platform for Intelligent Computing (NuPIC)
 * Copyright (C) 2017, Numenta, Inc.  Unless you have an agreement
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

#ifndef NTA_SDR_SELECTION_HPP
#define NTA_SDR_SELECTION_HPP

#include <vector>

#include <nupic/types/Types.hpp>

namespace nupic {
  namespace experimental {
    namespace sdr_selection {
      /**
       * Generate a set of distant SDRs. No two SDRs on the set will overlap on
       * "threshold" or more bits.
       *
       * This generates the list brute-force. This approach is really only
       * intended to work with n < 25. If you get much larger than that, it
       * might run for hours / millenia.
       */
      std::vector<std::vector<UInt>> enumerateDistantSDRsBruteForce(
        UInt n,
        UInt w,
        UInt threshold);
    }
  }
}

#endif // NTA_SDR_SELECTION_HPP
