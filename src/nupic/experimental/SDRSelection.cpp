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

#include <vector>

#include <nupic/experimental/SDRSelection.hpp>

using std::vector;
using namespace nupic;

template<typename Iter1, typename Iter2>
UInt overlap(Iter1 aBegin, Iter1 aEnd, Iter2 bBegin, Iter2 bEnd)
{
  UInt count = 0;

  auto it1 = aBegin;
  auto it2 = bBegin;

  if (it1 != aEnd && it2 != bEnd)
  {
    while (true)
    {
      if (*it1 == *it2)
      {
        count++;
        it1++;
        it2++;
        if (it1 == aEnd || it2 == bEnd)
        {
          break;
        }
      }
      else if (*it1 < *it2)
      {
        it1++;
        if (it1 == aEnd)
        {
          break;
        }
      }
      else
      {
        it2++;
        if (it2 == bEnd)
        {
          break;
        }
      }
    }
  }

  return count;
}

void enumerateDistantSDRsHelper(
  vector<vector<UInt>>& results,
  vector<UInt>& current,
  UInt digit,
  UInt begin,
  UInt n,
  UInt w,
  UInt threshold)
{
  for (UInt i = begin; i < n; i++)
  {
    current[digit] = i;

    bool qualified = true;
    for (const vector<UInt>& sdr : results)
    {
      if (overlap(current.begin(), current.begin() + digit + 1,
                  sdr.begin(), sdr.end()) >= threshold)
      {
        qualified = false;
        break;
      }
    }

    if (qualified)
    {
      if (digit == w - 1)
      {
        results.push_back(current);
      }
      else
      {
        enumerateDistantSDRsHelper(results, current, digit + 1, i + 1, n, w,
                                   threshold);
      }
    }
  }
}

vector<vector<UInt>>
nupic::experimental::sdr_selection::enumerateDistantSDRsBruteForce(
  UInt n,
  UInt w,
  UInt threshold)
{
  vector<vector<UInt>> results;
  vector<UInt> current(w, (UInt)-1);
  enumerateDistantSDRsHelper(results, current, 0, 0, n, w, threshold);

  return results;
}
