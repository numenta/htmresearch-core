/* ---------------------------------------------------------------------
 * Numenta Platform for Intelligent Computing (NuPIC)
 * Copyright (C) 2013-2016, Numenta, Inc.  Unless you have an agreement
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
 * Implementation of ExtendedTemporalMemory
 *
 * The functions in this file use the following parameter ordering
 * convention:
 *
 * 1. Output / mutated params
 * 2. Traditional parameters to the function, i.e. the ones that would still
 *    exist if this function were a method on a class
 * 3. Model state (marked const)
 * 4. Model parameters (including "learn")
 */

#include <cstring>
#include <climits>
#include <iostream>
#include <string>
#include <iterator>
#include <vector>

#include <capnp/message.h>
#include <capnp/serialize.h>
#include <kj/std/iostream.h>

#include <nupic/algorithms/Connections.hpp>
#include <nupic/experimental/ExtendedTemporalMemory.hpp>
#include <nupic/utils/GroupBy.hpp>

using namespace std;
using namespace nupic;
using namespace nupic::algorithms::connections;
using namespace nupic::experimental::extended_temporal_memory;

static const Permanence EPSILON = 0.000001;
static const UInt EXTENDED_TM_VERSION = 1;
static const UInt32 MIN_PREDICTIVE_THRESHOLD = 2;



template<typename Iterator>
bool isSortedWithoutDuplicates(Iterator begin, Iterator end)
{
  if (std::distance(begin, end) >= 2)
  {
    Iterator now = begin;
    Iterator next = begin + 1;
    while (next != end)
    {
      if (*now >= *next)
      {
        return false;
      }

      now = next++;
    }
  }

  return true;
}


ExtendedTemporalMemory::ExtendedTemporalMemory()
{
}

ExtendedTemporalMemory::ExtendedTemporalMemory(
  UInt columnCount,
  UInt basalInputSize,
  UInt apicalInputSize,
  UInt cellsPerColumn,
  UInt activationThreshold,
  Permanence initialPermanence,
  Permanence connectedPermanence,
  UInt minThreshold,
  UInt sampleSize,
  Permanence permanenceIncrement,
  Permanence permanenceDecrement,
  Permanence predictedSegmentDecrement,
  bool learnOnOneCell,
  Int seed,
  UInt maxSegmentsPerCell,
  UInt maxSynapsesPerSegment,
  bool checkInputs)
{
  initialize(
    columnCount,
    basalInputSize,
    apicalInputSize,
    cellsPerColumn,
    activationThreshold,
    initialPermanence,
    connectedPermanence,
    minThreshold,
    sampleSize,
    permanenceIncrement,
    permanenceDecrement,
    predictedSegmentDecrement,
    learnOnOneCell,
    seed,
    maxSegmentsPerCell,
    maxSynapsesPerSegment,
    checkInputs);
}

ExtendedTemporalMemory::~ExtendedTemporalMemory()
{
}

void ExtendedTemporalMemory::initialize(
  UInt columnCount,
  UInt basalInputSize,
  UInt apicalInputSize,
  UInt cellsPerColumn,
  UInt activationThreshold,
  Permanence initialPermanence,
  Permanence connectedPermanence,
  UInt minThreshold,
  UInt sampleSize,
  Permanence permanenceIncrement,
  Permanence permanenceDecrement,
  Permanence predictedSegmentDecrement,
  bool learnOnOneCell,
  Int seed,
  UInt maxSegmentsPerCell,
  UInt maxSynapsesPerSegment,
  bool checkInputs)
{
  NTA_CHECK(columnCount > 0);
  NTA_CHECK(cellsPerColumn > 0);
  NTA_CHECK(initialPermanence >= 0.0 && initialPermanence <= 1.0);
  NTA_CHECK(connectedPermanence >= 0.0 && connectedPermanence <= 1.0);
  NTA_CHECK(permanenceIncrement >= 0.0 && permanenceIncrement <= 1.0);
  NTA_CHECK(permanenceDecrement >= 0.0 && permanenceDecrement <= 1.0);
  NTA_CHECK(minThreshold <= activationThreshold);

  columnCount_ = columnCount;
  basalInputSize_ = basalInputSize;
  apicalInputSize_ = apicalInputSize;
  cellsPerColumn_ = cellsPerColumn;
  activationThreshold_ = activationThreshold;
  initialPermanence_ = initialPermanence;
  connectedPermanence_ = connectedPermanence;
  minThreshold_ = minThreshold;
  sampleSize_ = sampleSize;
  learnOnOneCell_ = learnOnOneCell;
  checkInputs_ = checkInputs;
  permanenceIncrement_ = permanenceIncrement;
  permanenceDecrement_ = permanenceDecrement;
  predictedSegmentDecrement_ = predictedSegmentDecrement;
  maxSegmentsPerCell_ = maxSegmentsPerCell;
  maxSynapsesPerSegment_ = maxSynapsesPerSegment;
  iteration_ = 0;

  basalConnections = Connections(numberOfCells());
  apicalConnections = Connections(numberOfCells());

  seed_((UInt64)(seed < 0 ? rand() : seed));

  activeCells_.clear();
  winnerCells_.clear();
  predictedActiveCells_.clear();
  activeBasalSegments_.clear();
  matchingBasalSegments_.clear();
  activeApicalSegments_.clear();
  matchingApicalSegments_.clear();
  chosenCellForColumn_.clear();
}

static UInt32 predictiveScore(
  vector<Segment>::const_iterator cellActiveBasalBegin,
  vector<Segment>::const_iterator cellActiveBasalEnd,
  vector<Segment>::const_iterator cellActiveApicalBegin,
  vector<Segment>::const_iterator cellActiveApicalEnd)
{
  UInt32 score = 0;

  if (cellActiveBasalBegin != cellActiveBasalEnd)
  {
    score += 2;
  }

  if (cellActiveApicalBegin != cellActiveApicalEnd)
  {
    score += 1;
  }

  return score;
}

static tuple<vector<Segment>::const_iterator,
             vector<Segment>::const_iterator>
segmentsForCell(vector<Segment>::const_iterator segmentsStart,
                vector<Segment>::const_iterator segmentsEnd,
                CellIdx cell,
                const Connections& connections)
{
  const auto cellBegin = std::find_if(
    segmentsStart, segmentsEnd,
    [&](Segment segment)
    {
      return connections.cellForSegment(segment) == cell;
    });
  const auto cellEnd = std::find_if(
    cellBegin, segmentsEnd,
    [&](Segment segment)
    {
      return connections.cellForSegment(segment) != cell;
    });
  return std::make_tuple(cellBegin, cellEnd);
}

static CellIdx getLeastUsedCell(
  Random& rng,
  UInt column,
  const Connections& connections,
  UInt cellsPerColumn)
{
  const CellIdx start = column * cellsPerColumn;
  const CellIdx end = start + cellsPerColumn;

  UInt32 minNumSegments = UINT_MAX;
  UInt32 numTiedCells = 0;
  for (CellIdx cell = start; cell < end; cell++)
  {
    const UInt32 numSegments = connections.numSegments(cell);
    if (numSegments < minNumSegments)
    {
      minNumSegments = numSegments;
      numTiedCells = 1;
    }
    else if (numSegments == minNumSegments)
    {
      numTiedCells++;
    }
  }

  const UInt32 tieWinnerIndex = rng.getUInt32(numTiedCells);

  UInt32 tieIndex = 0;
  for (CellIdx cell = start; cell < end; cell++)
  {
    if (connections.numSegments(cell) == minNumSegments)
    {
      if (tieIndex == tieWinnerIndex)
      {
        return cell;
      }
      else
      {
        tieIndex++;
      }
    }
  }

  NTA_THROW << "getLeastUsedCell failed to find a cell";
}

static void adaptSegment(
  Connections& connections,
  Segment segment,
  const vector<bool>& activeInputDense,
  Permanence permanenceIncrement,
  Permanence permanenceDecrement)
{
  const vector<Synapse>& synapses = connections.synapsesForSegment(segment);

  for (SynapseIdx i = 0; i < synapses.size();)
  {
    const SynapseData& synapseData = connections.dataForSynapse(synapses[i]);

    Permanence permanence = synapseData.permanence;
    if (activeInputDense[synapseData.presynapticCell])
    {
      permanence += permanenceIncrement;
    }
    else
    {
      permanence -= permanenceDecrement;
    }

    permanence = min(permanence, (Permanence)1.0);
    permanence = max(permanence, (Permanence)0.0);

    if (permanence < EPSILON)
    {
      connections.destroySynapse(synapses[i]);
      // Synapses vector is modified in-place, so don't update `i`.
    }
    else
    {
      connections.updateSynapsePermanence(synapses[i], permanence);
      i++;
    }
  }

  if (synapses.size() == 0)
  {
    connections.destroySegment(segment);
  }
}

static void destroyMinPermanenceSynapses(
  Connections& connections,
  Random& rng,
  Segment segment,
  Int nDestroy,
  const CellIdx* excludeCellsBegin,
  const CellIdx* excludeCellsEnd)
{
  // Don't destroy any cells that are in excludeCells.
  vector<Synapse> destroyCandidates;
  for (Synapse synapse : connections.synapsesForSegment(segment))
  {
    const CellIdx presynapticCell =
      connections.dataForSynapse(synapse).presynapticCell;

    if (!std::binary_search(excludeCellsBegin, excludeCellsEnd,
                            presynapticCell))
    {
      destroyCandidates.push_back(synapse);
    }
  }

  // Find cells one at a time. This is slow, but this code rarely runs, and it
  // needs to work around floating point differences between environments.
  for (Int32 i = 0; i < nDestroy && !destroyCandidates.empty(); i++)
  {
    Permanence minPermanence = std::numeric_limits<Permanence>::max();
    vector<Synapse>::iterator minSynapse = destroyCandidates.end();

    for (auto synapse = destroyCandidates.begin();
         synapse != destroyCandidates.end();
         synapse++)
    {
      const Permanence permanence =
        connections.dataForSynapse(*synapse).permanence;

      // Use special EPSILON logic to compensate for floating point
      // differences between C++ and other environments.
      if (permanence < minPermanence - EPSILON)
      {
        minSynapse = synapse;
        minPermanence = permanence;
      }
    }

    connections.destroySynapse(*minSynapse);
    destroyCandidates.erase(minSynapse);
  }
}

static void growSynapses(
  Connections& connections,
  Random& rng,
  Segment segment,
  UInt32 nDesiredNewSynapses,
  const CellIdx* growthCandidatesBegin,
  const CellIdx* growthCandidatesEnd,
  Permanence initialPermanence,
  UInt maxSynapsesPerSegment)
{
  // It's possible to optimize this, swapping candidates to the end as
  // they're used. But this is awkward to mimic in other
  // implementations, especially because it requires iterating over
  // the existing synapses in a particular order.

  vector<CellIdx> candidates(growthCandidatesBegin, growthCandidatesEnd);

  NTA_ASSERT(std::is_sorted(candidates.begin(), candidates.end()));

  // Remove cells that are already synapsed on by this segment
  for (Synapse synapse : connections.synapsesForSegment(segment))
  {
    CellIdx presynapticCell =
      connections.dataForSynapse(synapse).presynapticCell;
    auto ineligible = std::lower_bound(candidates.begin(), candidates.end(),
                                       presynapticCell);
    if (ineligible != candidates.end() && *ineligible == presynapticCell)
    {
      candidates.erase(ineligible);
    }
  }

  const UInt32 nActual = std::min(nDesiredNewSynapses,
                                  (UInt32)candidates.size());

  // Check if we're going to surpass the maximum number of synapses.
  const Int32 overrun = (connections.numSynapses(segment) +
                         nActual - maxSynapsesPerSegment);
  if (overrun > 0)
  {
    destroyMinPermanenceSynapses(connections, rng, segment, overrun,
                                 growthCandidatesBegin, growthCandidatesEnd);
  }

  // Recalculate in case we weren't able to destroy as many synapses as needed.
  const UInt32 nActualWithMax = std::min(nActual,
                                         maxSynapsesPerSegment -
                                         connections.numSynapses(segment));

  // Pick nActualWithMax cells randomly.
  for (UInt32 c = 0; c < nActualWithMax; c++)
  {
    size_t i = rng.getUInt32(candidates.size());
    connections.createSynapse(segment, candidates[i], initialPermanence);
    candidates.erase(candidates.begin() + i);
  }
}

static Segment createSegment(
  Connections& connections,
  vector<UInt64>& lastUsedIterationForSegment,
  CellIdx cell,
  UInt64 iteration,
  UInt maxSegmentsPerCell)
{
  while (connections.numSegments(cell) >= maxSegmentsPerCell)
  {
    const vector<Segment>& destroyCandidates =
      connections.segmentsForCell(cell);

    auto leastRecentlyUsedSegment = std::min_element(
      destroyCandidates.begin(), destroyCandidates.end(),
      [&](Segment a, Segment b)
      {
        return (lastUsedIterationForSegment[a] <
                lastUsedIterationForSegment[b]);
      });

    connections.destroySegment(*leastRecentlyUsedSegment);
  }

  const Segment segment = connections.createSegment(cell);
  lastUsedIterationForSegment.resize(connections.segmentFlatListLength());
  lastUsedIterationForSegment[segment] = iteration;

  return segment;
}

static void learnOnCell(
  Connections& connections,
  Random& rng,
  vector<UInt64>& lastUsedIterationForSegment,
  CellIdx cell,
  vector<Segment>::const_iterator cellActiveSegmentsBegin,
  vector<Segment>::const_iterator cellActiveSegmentsEnd,
  vector<Segment>::const_iterator cellMatchingSegmentsBegin,
  vector<Segment>::const_iterator cellMatchingSegmentsEnd,
  const vector<bool>& activeInputDense,
  const CellIdx* growthCandidatesBegin,
  const CellIdx* growthCandidatesEnd,
  const vector<UInt32>& potentialOverlaps,
  UInt64 iteration,
  UInt sampleSize,
  Permanence initialPermanence,
  Permanence permanenceIncrement,
  Permanence permanenceDecrement,
  UInt maxSegmentsPerCell,
  UInt maxSynapsesPerSegment)
{
  if (cellActiveSegmentsBegin != cellActiveSegmentsEnd)
  {
    // Learn on every active segment.

    auto activeSegment = cellActiveSegmentsBegin;
    do
    {
      adaptSegment(connections,
                   *activeSegment,
                   activeInputDense,
                   permanenceIncrement, permanenceDecrement);

      const Int32 nGrowDesired = sampleSize -
        potentialOverlaps[*activeSegment];
      if (nGrowDesired > 0)
      {
        growSynapses(connections, rng,
                     *activeSegment, nGrowDesired,
                     growthCandidatesBegin, growthCandidatesEnd,
                     initialPermanence, maxSynapsesPerSegment);
      }
    } while (++activeSegment != cellActiveSegmentsEnd);
  }
  else if (cellMatchingSegmentsBegin != cellMatchingSegmentsEnd)
  {
    // No active segments.
    // Learn on the best matching segment.

    const Segment bestMatchingSegment = *std::max_element(
      cellMatchingSegmentsBegin, cellMatchingSegmentsEnd,
      [&](Segment a, Segment b)
      {
        return (potentialOverlaps[a] <
                potentialOverlaps[b]);
      });

    adaptSegment(connections,
                 bestMatchingSegment,
                 activeInputDense,
                 permanenceIncrement, permanenceDecrement);

    const Int32 nGrowDesired = sampleSize -
      potentialOverlaps[bestMatchingSegment];
    if (nGrowDesired > 0)
    {
      growSynapses(connections, rng,
                   bestMatchingSegment, nGrowDesired,
                   growthCandidatesBegin, growthCandidatesEnd,
                   initialPermanence, maxSynapsesPerSegment);
    }
  }
  else
  {
    // No matching segments.
    // Grow a new segment and learn on it.

    // Don't grow a segment that will never match.
    const UInt32 nGrowExact = std::min(sampleSize,
                                       (UInt32)std::distance(
                                         growthCandidatesBegin,
                                         growthCandidatesEnd));
    if (nGrowExact > 0)
    {
      const Segment segment = createSegment(connections,
                                            lastUsedIterationForSegment, cell,
                                            iteration, maxSegmentsPerCell);
      growSynapses(connections, rng,
                   segment, nGrowExact,
                   growthCandidatesBegin, growthCandidatesEnd,
                   initialPermanence, maxSynapsesPerSegment);
      NTA_ASSERT(connections.numSynapses(segment) == nGrowExact);
    }
  }
}

static void activatePredictedColumn(
  vector<CellIdx>& activeCells,
  vector<CellIdx>& winnerCells,
  vector<CellIdx>& predictedActiveCells,
  Connections& basalConnections,
  Connections& apicalConnections,
  Random& rng,
  vector<UInt64>& lastUsedIterationForBasalSegment,
  vector<UInt64>& lastUsedIterationForApicalSegment,
  vector<CellIdx>::const_iterator columnPredictedCellsBegin,
  vector<CellIdx>::const_iterator columnPredictedCellsEnd,
  vector<Segment>::const_iterator columnActiveBasalBegin,
  vector<Segment>::const_iterator columnActiveBasalEnd,
  vector<Segment>::const_iterator columnMatchingBasalBegin,
  vector<Segment>::const_iterator columnMatchingBasalEnd,
  vector<Segment>::const_iterator columnActiveApicalBegin,
  vector<Segment>::const_iterator columnActiveApicalEnd,
  vector<Segment>::const_iterator columnMatchingApicalBegin,
  vector<Segment>::const_iterator columnMatchingApicalEnd,
  const vector<bool>& basalInputDense,
  const vector<bool>& apicalInputDense,
  const CellIdx* basalGrowthCandidatesBegin,
  const CellIdx* basalGrowthCandidatesEnd,
  const CellIdx* apicalGrowthCandidatesBegin,
  const CellIdx* apicalGrowthCandidatesEnd,
  const vector<UInt32>& basalPotentialOverlaps,
  const vector<UInt32>& apicalPotentialOverlaps,
  UInt64 iteration,
  UInt sampleSize,
  Permanence initialPermanence,
  Permanence permanenceIncrement,
  Permanence permanenceDecrement,
  UInt maxSegmentsPerCell,
  UInt maxSynapsesPerSegment,
  bool learn)
{
  const auto cellForBasalSegment = [&](Segment segment)
    { return basalConnections.cellForSegment(segment); };
  const auto cellForApicalSegment = [&](Segment segment)
    { return apicalConnections.cellForSegment(segment); };

  for (auto& cellData : iterGroupBy(
         columnPredictedCellsBegin, columnPredictedCellsEnd, identity<CellIdx>,
         columnActiveBasalBegin, columnActiveBasalEnd, cellForBasalSegment,
         columnMatchingBasalBegin, columnMatchingBasalEnd, cellForBasalSegment,
         columnActiveApicalBegin, columnActiveApicalEnd, cellForApicalSegment,
         columnMatchingApicalBegin, columnMatchingApicalEnd, cellForApicalSegment))
  {
    CellIdx cell;
    vector<CellIdx>::const_iterator
      cellPredictedCellsBegin, cellPredictedCellsEnd;
    vector<Segment>::const_iterator
      cellActiveBasalBegin, cellActiveBasalEnd,
      cellMatchingBasalBegin, cellMatchingBasalEnd,
      cellActiveApicalBegin, cellActiveApicalEnd,
      cellMatchingApicalBegin, cellMatchingApicalEnd;
    tie(cell,
        cellPredictedCellsBegin, cellPredictedCellsEnd,
        cellActiveBasalBegin, cellActiveBasalEnd,
        cellMatchingBasalBegin, cellMatchingBasalEnd,
        cellActiveApicalBegin, cellActiveApicalEnd,
        cellMatchingApicalBegin, cellMatchingApicalEnd) = cellData;

    const bool isPredictedCell = (cellPredictedCellsBegin !=
                                  cellPredictedCellsEnd);

    if (isPredictedCell)
    {
      activeCells.push_back(cell);
      winnerCells.push_back(cell);
      predictedActiveCells.push_back(cell);

      if (learn)
      {
        learnOnCell(basalConnections, rng, lastUsedIterationForBasalSegment,
                    cell,
                    cellActiveBasalBegin, cellActiveBasalEnd,
                    cellMatchingBasalBegin, cellMatchingBasalEnd,
                    basalInputDense,
                    basalGrowthCandidatesBegin, basalGrowthCandidatesEnd,
                    basalPotentialOverlaps, iteration,
                    sampleSize, initialPermanence,
                    permanenceIncrement, permanenceDecrement,
                    maxSegmentsPerCell, maxSynapsesPerSegment);

        learnOnCell(apicalConnections, rng, lastUsedIterationForApicalSegment,
                    cell,
                    cellActiveApicalBegin, cellActiveApicalEnd,
                    cellMatchingApicalBegin, cellMatchingApicalEnd,
                    apicalInputDense,
                    apicalGrowthCandidatesBegin, apicalGrowthCandidatesEnd,
                    apicalPotentialOverlaps, iteration,
                    sampleSize, initialPermanence,
                    permanenceIncrement, permanenceDecrement,
                    maxSegmentsPerCell, maxSynapsesPerSegment);
      }
    }
  }
}

static void burstColumn(
  vector<CellIdx>& activeCells,
  vector<CellIdx>& winnerCells,
  Connections& basalConnections,
  Connections& apicalConnections,
  Random& rng,
  vector<UInt64>& lastUsedIterationForBasalSegment,
  vector<UInt64>& lastUsedIterationForApicalSegment,
  map<UInt, CellIdx>& chosenCellForColumn,
  UInt column,
  vector<Segment>::const_iterator columnActiveBasalBegin,
  vector<Segment>::const_iterator columnActiveBasalEnd,
  vector<Segment>::const_iterator columnMatchingBasalBegin,
  vector<Segment>::const_iterator columnMatchingBasalEnd,
  vector<Segment>::const_iterator columnActiveApicalBegin,
  vector<Segment>::const_iterator columnActiveApicalEnd,
  vector<Segment>::const_iterator columnMatchingApicalBegin,
  vector<Segment>::const_iterator columnMatchingApicalEnd,
  const vector<bool>& basalInputDense,
  const vector<bool>& apicalInputDense,
  const CellIdx* basalGrowthCandidatesBegin,
  const CellIdx* basalGrowthCandidatesEnd,
  const CellIdx* apicalGrowthCandidatesBegin,
  const CellIdx* apicalGrowthCandidatesEnd,
  const vector<UInt32>& basalPotentialOverlaps,
  const vector<UInt32>& apicalPotentialOverlaps,
  UInt64 iteration,
  UInt cellsPerColumn,
  UInt sampleSize,
  Permanence initialPermanence,
  Permanence permanenceIncrement,
  Permanence permanenceDecrement,
  UInt maxSegmentsPerCell,
  UInt maxSynapsesPerSegment,
  bool learnOnOneCell,
  bool learn)
{
  // Calculate the active cells.
  const CellIdx start = column * cellsPerColumn;
  const CellIdx end = start + cellsPerColumn;
  for (CellIdx cell = start; cell < end; cell++)
  {
    activeCells.push_back(cell);
  }

  // Mini optimization: don't search for the best basal segment twice.
  auto basalCandidatesBegin = columnMatchingBasalBegin;
  auto basalCandidatesEnd = columnMatchingBasalEnd;

  // Calculate the winner cell.
  CellIdx winnerCell;
  if (learnOnOneCell && chosenCellForColumn.count(column))
  {
    winnerCell = chosenCellForColumn.at(column);
  }
  else
  {
    if (columnMatchingBasalBegin != columnMatchingBasalEnd)
    {
      auto bestBasalSegment = std::max_element(
        columnMatchingBasalBegin, columnMatchingBasalEnd,
        [&](Segment a, Segment b)
        {
          return (basalPotentialOverlaps[a] <
                  basalPotentialOverlaps[b]);
        });

      basalCandidatesBegin = bestBasalSegment;
      basalCandidatesEnd = bestBasalSegment + 1;

      winnerCell = basalConnections.cellForSegment(*bestBasalSegment);
    }
    else
    {
      winnerCell = getLeastUsedCell(rng, column, basalConnections,
                                    cellsPerColumn);
    }

    if (learnOnOneCell)
    {
      chosenCellForColumn[column] = winnerCell;
    }
  }
  winnerCells.push_back(winnerCell);

  // Learn.
  if (learn)
  {
    vector<Segment>::const_iterator
      cellActiveBasalBegin, cellActiveBasalEnd,
      cellMatchingBasalBegin, cellMatchingBasalEnd,
      cellActiveApicalBegin, cellActiveApicalEnd,
      cellMatchingApicalBegin, cellMatchingApicalEnd;
    tie(cellActiveBasalBegin,
        cellActiveBasalEnd) = segmentsForCell(columnActiveBasalBegin,
                                              columnActiveBasalEnd,
                                              winnerCell,
                                              basalConnections);
    tie(cellMatchingBasalBegin,
        cellMatchingBasalEnd) = segmentsForCell(basalCandidatesBegin,
                                                basalCandidatesEnd,
                                                winnerCell,
                                                basalConnections);
    tie(cellActiveApicalBegin,
        cellActiveApicalEnd) = segmentsForCell(columnActiveApicalBegin,
                                               columnActiveApicalEnd,
                                               winnerCell,
                                               apicalConnections);
    tie(cellMatchingApicalBegin,
        cellMatchingApicalEnd) = segmentsForCell(columnMatchingApicalBegin,
                                                 columnMatchingApicalEnd,
                                                 winnerCell,
                                                 apicalConnections);

    learnOnCell(basalConnections, rng, lastUsedIterationForBasalSegment,
                winnerCell,
                cellActiveBasalBegin, cellActiveBasalEnd,
                cellMatchingBasalBegin, cellMatchingBasalEnd,
                basalInputDense,
                basalGrowthCandidatesBegin, basalGrowthCandidatesEnd,
                basalPotentialOverlaps, iteration,
                sampleSize, initialPermanence,
                permanenceIncrement, permanenceDecrement,
                maxSegmentsPerCell, maxSynapsesPerSegment);

    learnOnCell(apicalConnections, rng, lastUsedIterationForApicalSegment,
                winnerCell,
                cellActiveApicalBegin, cellActiveApicalEnd,
                cellMatchingApicalBegin, cellMatchingApicalEnd,
                apicalInputDense,
                apicalGrowthCandidatesBegin, apicalGrowthCandidatesEnd,
                apicalPotentialOverlaps, iteration,
                sampleSize, initialPermanence,
                permanenceIncrement, permanenceDecrement,
                maxSegmentsPerCell, maxSynapsesPerSegment);
  }
}

static void punishPredictedColumn(
  Connections& connections,
  vector<Segment>::const_iterator matchingSegmentsBegin,
  vector<Segment>::const_iterator matchingSegmentsEnd,
  const vector<bool>& activeInputDense,
  Permanence predictedSegmentDecrement)
{
  if (predictedSegmentDecrement > 0.0)
  {
    for (auto matchingSegment = matchingSegmentsBegin;
         matchingSegment != matchingSegmentsEnd; matchingSegment++)
    {
      adaptSegment(connections, *matchingSegment,
                   activeInputDense,
                   -predictedSegmentDecrement, 0.0);
    }
  }
}

void ExtendedTemporalMemory::activateCells_(
  const UInt* activeColumnsBegin,
  const UInt* activeColumnsEnd,
  const CellIdx* basalInputBegin,
  const CellIdx* basalInputEnd,
  const CellIdx* apicalInputBegin,
  const CellIdx* apicalInputEnd,
  const CellIdx* basalGrowthCandidatesBegin,
  const CellIdx* basalGrowthCandidatesEnd,
  const CellIdx* apicalGrowthCandidatesBegin,
  const CellIdx* apicalGrowthCandidatesEnd,
  bool learn)
{
  activeCells_.clear();
  winnerCells_.clear();
  predictedActiveCells_.clear();

  // Perf: Densify these inputs so adaptSegment can quickly check
  // whether a synapse is active.
  vector<bool> basalInputDense(basalInputSize_, false);
  for (auto it = basalInputBegin; it != basalInputEnd; it++)
  {
    basalInputDense[*it] = true;
  }
  vector<bool> apicalInputDense(apicalInputSize_, false);
  for (auto it = apicalInputBegin; it != apicalInputEnd; it++)
  {
    apicalInputDense[*it] = true;
  }

  const auto columnForCellFn = [&](CellIdx cell)
    { return this->columnForCell(cell); };
  const auto columnForBasalSegment = [&](Segment segment)
    { return basalConnections.cellForSegment(segment) / cellsPerColumn_; };
  const auto columnForApicalSegment = [&](Segment segment)
    { return apicalConnections.cellForSegment(segment) / cellsPerColumn_; };

  for (auto& columnData : iterGroupBy(
         activeColumnsBegin, activeColumnsEnd, identity<UInt>,
         predictedCells_.begin(), predictedCells_.end(), columnForCellFn,
         activeBasalSegments_.begin(),
         activeBasalSegments_.end(), columnForBasalSegment,
         matchingBasalSegments_.begin(),
         matchingBasalSegments_.end(), columnForBasalSegment,
         activeApicalSegments_.begin(),
         activeApicalSegments_.end(), columnForApicalSegment,
         matchingApicalSegments_.begin(),
         matchingApicalSegments_.end(), columnForApicalSegment))
  {
    UInt column;
    const UInt
      *columnActiveColumnsBegin, *columnActiveColumnsEnd;
    vector<CellIdx>::const_iterator
      columnPredictedCellsBegin, columnPredictedCellsEnd;
    vector<Segment>::const_iterator
      columnActiveBasalBegin, columnActiveBasalEnd,
      columnMatchingBasalBegin, columnMatchingBasalEnd,
      columnActiveApicalBegin, columnActiveApicalEnd,
      columnMatchingApicalBegin, columnMatchingApicalEnd;
    tie(column,
        columnActiveColumnsBegin, columnActiveColumnsEnd,
        columnPredictedCellsBegin, columnPredictedCellsEnd,
        columnActiveBasalBegin, columnActiveBasalEnd,
        columnMatchingBasalBegin, columnMatchingBasalEnd,
        columnActiveApicalBegin, columnActiveApicalEnd,
        columnMatchingApicalBegin, columnMatchingApicalEnd) = columnData;

    const bool isActiveColumn = (columnActiveColumnsBegin !=
                                 columnActiveColumnsEnd);
    const bool isPredictedColumn = (columnPredictedCellsBegin !=
                                    columnPredictedCellsEnd);

    if (isActiveColumn)
    {
      if (isPredictedColumn)
      {
        activatePredictedColumn(
          activeCells_, winnerCells_, predictedActiveCells_,
          basalConnections, apicalConnections, rng_,
          lastUsedIterationForBasalSegment_, lastUsedIterationForApicalSegment_,
          columnPredictedCellsBegin, columnPredictedCellsEnd,
          columnActiveBasalBegin, columnActiveBasalEnd,
          columnMatchingBasalBegin, columnMatchingBasalEnd,
          columnActiveApicalBegin, columnActiveApicalEnd,
          columnMatchingApicalBegin, columnMatchingApicalEnd,
          basalInputDense, apicalInputDense,
          basalGrowthCandidatesBegin, basalGrowthCandidatesEnd,
          apicalGrowthCandidatesBegin, apicalGrowthCandidatesEnd,
          basalPotentialOverlaps_,
          apicalPotentialOverlaps_, iteration_,
          sampleSize_,
          initialPermanence_, permanenceIncrement_, permanenceDecrement_,
          maxSegmentsPerCell_, maxSynapsesPerSegment_,
          learn);
      }
      else
      {
        burstColumn(
          activeCells_, winnerCells_,
          basalConnections, apicalConnections, rng_,
          lastUsedIterationForBasalSegment_, lastUsedIterationForApicalSegment_,
          chosenCellForColumn_,
          column,
          columnActiveBasalBegin, columnActiveBasalEnd,
          columnMatchingBasalBegin, columnMatchingBasalEnd,
          columnActiveApicalBegin, columnActiveApicalEnd,
          columnMatchingApicalBegin, columnMatchingApicalEnd,
          basalInputDense, apicalInputDense,
          basalGrowthCandidatesBegin, basalGrowthCandidatesEnd,
          apicalGrowthCandidatesBegin, apicalGrowthCandidatesEnd,
          basalPotentialOverlaps_,
          apicalPotentialOverlaps_, iteration_,
          cellsPerColumn_, sampleSize_,
          initialPermanence_, permanenceIncrement_, permanenceDecrement_,
          maxSegmentsPerCell_, maxSynapsesPerSegment_,
          learnOnOneCell_, learn);
      }
    }
    else
    {
      if (learn)
      {
        punishPredictedColumn(
          basalConnections,
          columnMatchingBasalBegin, columnMatchingBasalEnd,
          basalInputDense,
          predictedSegmentDecrement_);

        // Don't punish apical segments.
      }
    }
  }
}

static void calculateExcitation(
  vector<UInt32>& overlaps,
  vector<Segment>& activeSegments,
  vector<UInt32>& potentialOverlaps,
  vector<Segment>& matchingSegments,
  Connections& connections,
  const CellIdx* activeInputBegin,
  const CellIdx* activeInputEnd,
  Permanence connectedPermanence,
  UInt activationThreshold,
  UInt minThreshold)
{
  const UInt32 length = connections.segmentFlatListLength();
  overlaps.assign(length, 0);
  potentialOverlaps.assign(length, 0);

  for (auto cell = activeInputBegin; cell != activeInputEnd; cell++)
  {
    connections.computeActivity(overlaps, potentialOverlaps,
                                *cell, connectedPermanence);
  }

  // Active segments, connected synapses.
  activeSegments.clear();
  for (Segment segment = 0;
       segment < overlaps.size();
       segment++)
  {
    if (overlaps[segment] >= activationThreshold)
    {
      activeSegments.push_back(segment);
    }
  }
  std::sort(activeSegments.begin(), activeSegments.end(),
            [&](Segment a, Segment b)
            {
              return connections.compareSegments(a, b);
            });

  // Matching segments, potential synapses.
  matchingSegments.clear();
  for (Segment segment = 0;
       segment < potentialOverlaps.size();
       segment++)
  {
    if (potentialOverlaps[segment] >= minThreshold)
    {
      matchingSegments.push_back(segment);
    }
  }
  std::sort(matchingSegments.begin(), matchingSegments.end(),
            [&](Segment a, Segment b)
            {
              return connections.compareSegments(a, b);
            });
}

void ExtendedTemporalMemory::depolarizeCells_(
  const CellIdx* basalInputBegin,
  const CellIdx* basalInputEnd,
  const CellIdx* apicalInputBegin,
  const CellIdx* apicalInputEnd,
  bool learn)
{
  predictedCells_.clear();

  calculateExcitation(
    basalOverlaps_, activeBasalSegments_,
    basalPotentialOverlaps_, matchingBasalSegments_,
    basalConnections,
    basalInputBegin, basalInputEnd,
    connectedPermanence_, activationThreshold_, minThreshold_);

  calculateExcitation(
    apicalOverlaps_, activeApicalSegments_,
    apicalPotentialOverlaps_, matchingApicalSegments_,
    apicalConnections,
    apicalInputBegin, apicalInputEnd,
    connectedPermanence_, activationThreshold_, minThreshold_);

  const auto columnForBasalSegment = [&](Segment segment)
    { return basalConnections.cellForSegment(segment) / cellsPerColumn_; };
  const auto columnForApicalSegment = [&](Segment segment)
    { return apicalConnections.cellForSegment(segment) / cellsPerColumn_; };

  for (auto& columnData : groupBy(activeBasalSegments_, columnForBasalSegment,
                                  activeApicalSegments_, columnForApicalSegment))
  {
    UInt column;
    vector<Segment>::const_iterator
      columnActiveBasalBegin, columnActiveBasalEnd,
      columnActiveApicalBegin, columnActiveApicalEnd;
    tie(column,
        columnActiveBasalBegin, columnActiveBasalEnd,
        columnActiveApicalBegin, columnActiveApicalEnd) = columnData;

    const auto cellForBasalSegment = [&](Segment segment)
      { return basalConnections.cellForSegment(segment); };
    const auto cellForApicalSegment = [&](Segment segment)
      { return apicalConnections.cellForSegment(segment); };

    const auto groupedByCell = iterGroupBy(
      columnActiveBasalBegin, columnActiveBasalEnd, cellForBasalSegment,
      columnActiveApicalBegin, columnActiveApicalEnd, cellForApicalSegment);

    UInt32 maxDepolarization = 0;
    for (auto& cellData : groupedByCell)
    {
      vector<Segment>::const_iterator
        cellActiveBasalBegin, cellActiveBasalEnd,
        cellActiveApicalBegin, cellActiveApicalEnd;
      tie(std::ignore,
          cellActiveBasalBegin, cellActiveBasalEnd,
          cellActiveApicalBegin, cellActiveApicalEnd) = cellData;

      maxDepolarization = std::max(maxDepolarization,
                                   predictiveScore(cellActiveBasalBegin,
                                                   cellActiveBasalEnd,
                                                   cellActiveApicalBegin,
                                                   cellActiveApicalEnd));
    }

    if (maxDepolarization >= MIN_PREDICTIVE_THRESHOLD)
    {
      for (auto& cellData : groupedByCell)
      {
        CellIdx cell;
        vector<Segment>::const_iterator
          cellActiveBasalBegin, cellActiveBasalEnd,
          cellActiveApicalBegin, cellActiveApicalEnd;
        tie(cell,
            cellActiveBasalBegin, cellActiveBasalEnd,
            cellActiveApicalBegin, cellActiveApicalEnd) = cellData;

        if (predictiveScore(cellActiveBasalBegin, cellActiveBasalEnd,
                            cellActiveApicalBegin, cellActiveApicalEnd)
            == maxDepolarization)
        {
          predictedCells_.push_back(cell);
        }
      }
    }
  }

  if (learn)
  {
    for (Segment segment : activeBasalSegments_)
    {
      lastUsedIterationForBasalSegment_[segment] = iteration_;
    }
    for (Segment segment : activeApicalSegments_)
    {
      lastUsedIterationForApicalSegment_[segment] = iteration_;
    }

    iteration_++;
  }
}

void ExtendedTemporalMemory::compute(
  const UInt* activeColumnsBegin,
  const UInt* activeColumnsEnd,
  const CellIdx* basalInputBegin,
  const CellIdx* basalInputEnd,
  const CellIdx* apicalInputBegin,
  const CellIdx* apicalInputEnd,
  const CellIdx* basalGrowthCandidatesBegin,
  const CellIdx* basalGrowthCandidatesEnd,
  const CellIdx* apicalGrowthCandidatesBegin,
  const CellIdx* apicalGrowthCandidatesEnd,
  bool learn)
{
  if (checkInputs_)
  {
    NTA_CHECK(isSortedWithoutDuplicates(activeColumnsBegin, activeColumnsEnd))
      << "activeColumns must be sorted without duplicates.";
    NTA_CHECK(isSortedWithoutDuplicates(basalInputBegin, basalInputEnd))
      << "basalInput must be sorted without duplicates.";
    NTA_CHECK(isSortedWithoutDuplicates(apicalInputBegin, apicalInputEnd))
      << "apicalInput must be sorted without duplicates.";
    NTA_CHECK(isSortedWithoutDuplicates(basalGrowthCandidatesBegin,
                                        basalGrowthCandidatesEnd))
      << "basalGrowthCandidates must be sorted without duplicates.";
    NTA_CHECK(isSortedWithoutDuplicates(apicalGrowthCandidatesBegin,
                                        apicalGrowthCandidatesEnd))
      << "apicalGrowthCandidates must be sorted without duplicates.";

    NTA_CHECK(std::all_of(activeColumnsBegin, activeColumnsEnd,
                          [&](UInt c) { return c < columnCount_; }))
      << "Values in activeColumns must be within the range "
      << "[0," << columnCount_ << ").";
    NTA_CHECK(std::all_of(basalInputBegin, basalInputEnd,
                          [&](UInt c) { return c < basalInputSize_; }))
      << "Values in basalInput must be within the range "
      << "[0," << basalInputSize_ << ").";
    NTA_CHECK(std::all_of(apicalInputBegin, apicalInputEnd,
                          [&](UInt c) { return c < apicalInputSize_; }))
      << "Values in apicalInput must be within the range "
      << "[0," << apicalInputSize_ << ").";
    NTA_CHECK(std::all_of(basalGrowthCandidatesBegin, basalGrowthCandidatesEnd,
                          [&](UInt c) { return c < basalInputSize_; }))
      << "Values in basalGrowthCandidates must be within the range " <<
      "[0," << basalInputSize_ << ").";
    NTA_CHECK(std::all_of(apicalGrowthCandidatesBegin, apicalGrowthCandidatesEnd,
                          [&](UInt c) { return c < apicalInputSize_; }))
      << "Values in apicalGrowthCandidates must be within the range "
      << "[0," << apicalInputSize_ << ").";
  }

  depolarizeCells_(basalInputBegin, basalInputEnd,
                   apicalInputBegin, apicalInputEnd,
                   learn);

  activateCells_(activeColumnsBegin, activeColumnsEnd,
                 basalInputBegin, basalInputEnd,
                 apicalInputBegin, apicalInputEnd,
                 basalGrowthCandidatesBegin, basalGrowthCandidatesEnd,
                 apicalGrowthCandidatesBegin, apicalGrowthCandidatesEnd,
                 learn);
}

void ExtendedTemporalMemory::compute(
  const vector<UInt>& activeColumns,
  const vector<CellIdx>& basalInput,
  const vector<CellIdx>& apicalInput,
  const vector<CellIdx>& basalGrowthCandidates,
  const vector<CellIdx>& apicalGrowthCandidates,
  bool learn)
{
  compute(activeColumns.data(),
          activeColumns.data() + activeColumns.size(),
          basalInput.data(), basalInput.data() + basalInput.size(),
          apicalInput.data(), apicalInput.data() + apicalInput.size(),
          basalGrowthCandidates.data(),
          basalGrowthCandidates.data() + basalGrowthCandidates.size(),
          apicalGrowthCandidates.data(),
          apicalGrowthCandidates.data() + apicalGrowthCandidates.size(),
          learn);
}

void ExtendedTemporalMemory::sequenceMemoryCompute(
  const UInt* activeColumnsBegin,
  const UInt* activeColumnsEnd,
  const CellIdx* apicalInputBegin,
  const CellIdx* apicalInputEnd,
  const CellIdx* apicalGrowthCandidatesBegin,
  const CellIdx* apicalGrowthCandidatesEnd,
  bool learn)
{
  const vector<CellIdx> prevActiveCells(activeCells_);
  const vector<CellIdx> prevWinnerCells(winnerCells_);

  compute(activeColumnsBegin, activeColumnsEnd,
          prevActiveCells.data(),
          prevActiveCells.data() + prevActiveCells.size(),
          apicalInputBegin, apicalInputEnd,
          prevWinnerCells.data(),
          prevWinnerCells.data() + prevWinnerCells.size(),
          apicalGrowthCandidatesBegin, apicalGrowthCandidatesEnd,
          learn);
}

void ExtendedTemporalMemory::sequenceMemoryCompute(
  const vector<UInt>& activeColumns,
  const vector<CellIdx>& apicalInput,
  const vector<CellIdx>& apicalGrowthCandidates,
  bool learn)
{
  const vector<CellIdx> prevActiveCells(activeCells_);
  const vector<CellIdx> prevWinnerCells(winnerCells_);

  compute(activeColumns.data(),
          activeColumns.data() + activeColumns.size(),
          prevActiveCells.data(),
          prevActiveCells.data() + prevActiveCells.size(),
          apicalInput.data(),
          apicalInput.data() + apicalInput.size(),
          prevWinnerCells.data(),
          prevWinnerCells.data() + prevWinnerCells.size(),
          apicalGrowthCandidates.data(),
          apicalGrowthCandidates.data() + apicalGrowthCandidates.size(),
          learn);
}

void ExtendedTemporalMemory::reset(void)
{
  activeCells_.clear();
  winnerCells_.clear();
  predictedCells_.clear();
  predictedActiveCells_.clear();
  activeBasalSegments_.clear();
  matchingBasalSegments_.clear();
  activeApicalSegments_.clear();
  matchingApicalSegments_.clear();
  chosenCellForColumn_.clear();
}

// ==============================
//  Helper methods
// ==============================

Segment ExtendedTemporalMemory::createBasalSegment(CellIdx cell)
{
  return ::createSegment(basalConnections, lastUsedIterationForBasalSegment_,
                         cell, iteration_, maxSegmentsPerCell_);
}

Segment ExtendedTemporalMemory::createApicalSegment(CellIdx cell)
{
  return ::createSegment(apicalConnections, lastUsedIterationForApicalSegment_,
                         cell, iteration_, maxSegmentsPerCell_);
}

UInt ExtendedTemporalMemory::columnForCell(CellIdx cell)
{
  _validateCell(cell);

  return cell / cellsPerColumn_;
}

vector<CellIdx> ExtendedTemporalMemory::cellsForColumn(UInt column)
{
  NTA_CHECK(column < numberOfColumns()) << "Invalid column " << column;

  const CellIdx start = cellsPerColumn_ * column;
  const CellIdx end = start + cellsPerColumn_;

  vector<CellIdx> cellsInColumn;
  for (CellIdx i = start; i < end; i++)
  {
    cellsInColumn.push_back(i);
  }

  return cellsInColumn;
}

UInt ExtendedTemporalMemory::numberOfCells(void)
{
  return numberOfColumns() * cellsPerColumn_;
}

vector<CellIdx> ExtendedTemporalMemory::getActiveCells() const
{
  return activeCells_;
}

vector<CellIdx> ExtendedTemporalMemory::getPredictedCells() const
{
  return predictedCells_;
}

vector<CellIdx> ExtendedTemporalMemory::getPredictedActiveCells() const
{
  return predictedActiveCells_;
}

vector<CellIdx> ExtendedTemporalMemory::getWinnerCells() const
{
  return winnerCells_;
}

vector<Segment> ExtendedTemporalMemory::getActiveBasalSegments() const
{
  return activeBasalSegments_;
}

vector<Segment> ExtendedTemporalMemory::getMatchingBasalSegments() const
{
  return matchingBasalSegments_;
}

vector<Segment> ExtendedTemporalMemory::getActiveApicalSegments() const
{
  return activeApicalSegments_;
}

vector<Segment> ExtendedTemporalMemory::getMatchingApicalSegments() const
{
  return matchingApicalSegments_;
}

UInt ExtendedTemporalMemory::numberOfColumns() const
{
  return columnCount_;
}

bool ExtendedTemporalMemory::_validateCell(CellIdx cell)
{
  if (cell < numberOfCells())
    return true;

  NTA_THROW << "Invalid cell " << cell;
  return false;
}

UInt ExtendedTemporalMemory::getBasalInputSize() const
{
  return basalInputSize_;
}

UInt ExtendedTemporalMemory::getApicalInputSize() const
{
  return apicalInputSize_;
}

UInt ExtendedTemporalMemory::getCellsPerColumn() const
{
  return cellsPerColumn_;
}

UInt ExtendedTemporalMemory::getActivationThreshold() const
{
  return activationThreshold_;
}

void ExtendedTemporalMemory::setActivationThreshold(UInt activationThreshold)
{
  activationThreshold_ = activationThreshold;
}

Permanence ExtendedTemporalMemory::getInitialPermanence() const
{
  return initialPermanence_;
}

void ExtendedTemporalMemory::setInitialPermanence(Permanence initialPermanence)
{
  initialPermanence_ = initialPermanence;
}

Permanence ExtendedTemporalMemory::getConnectedPermanence() const
{
  return connectedPermanence_;
}

void ExtendedTemporalMemory::setConnectedPermanence(
  Permanence connectedPermanence)
{
  connectedPermanence_ = connectedPermanence;
}

UInt ExtendedTemporalMemory::getMinThreshold() const
{
  return minThreshold_;
}

void ExtendedTemporalMemory::setMinThreshold(UInt minThreshold)
{
  minThreshold_ = minThreshold;
}

UInt ExtendedTemporalMemory::getSampleSize() const
{
  return sampleSize_;
}

void ExtendedTemporalMemory::setSampleSize(UInt sampleSize)
{
  sampleSize_ = sampleSize;
}

bool ExtendedTemporalMemory::getLearnOnOneCell() const
{
  return learnOnOneCell_;
}

void ExtendedTemporalMemory::setLearnOnOneCell(bool learnOnOneCell)
{
  learnOnOneCell_ = learnOnOneCell;
}

Permanence ExtendedTemporalMemory::getPermanenceIncrement() const
{
  return permanenceIncrement_;
}

void ExtendedTemporalMemory::setPermanenceIncrement(
  Permanence permanenceIncrement)
{
  permanenceIncrement_ = permanenceIncrement;
}

Permanence ExtendedTemporalMemory::getPermanenceDecrement() const
{
  return permanenceDecrement_;
}

void ExtendedTemporalMemory::setPermanenceDecrement(
  Permanence permanenceDecrement)
{
  permanenceDecrement_ = permanenceDecrement;
}

Permanence ExtendedTemporalMemory::getPredictedSegmentDecrement() const
{
  return predictedSegmentDecrement_;
}

void ExtendedTemporalMemory::setPredictedSegmentDecrement(
  Permanence predictedSegmentDecrement)
{
  predictedSegmentDecrement_ = predictedSegmentDecrement;
}

UInt ExtendedTemporalMemory::getMaxSegmentsPerCell() const
{
  return maxSegmentsPerCell_;
}

UInt ExtendedTemporalMemory::getMaxSynapsesPerSegment() const
{
  return maxSynapsesPerSegment_;
}

bool ExtendedTemporalMemory::getCheckInputs() const
{
  return checkInputs_;
}

void ExtendedTemporalMemory::setCheckInputs(bool checkInputs)
{
  checkInputs_ = checkInputs;
}

UInt ExtendedTemporalMemory::version() const
{
  return EXTENDED_TM_VERSION;
}

/**
* Create a RNG with given seed
*/
void ExtendedTemporalMemory::seed_(UInt64 seed)
{
  rng_ = Random(seed);
}

void ExtendedTemporalMemory::write(ExtendedTemporalMemoryProto::Builder& proto) const
{
  proto.setColumnCount(columnCount_);
  proto.setCellsPerColumn(cellsPerColumn_);
  proto.setActivationThreshold(activationThreshold_);
  proto.setInitialPermanence(initialPermanence_);
  proto.setConnectedPermanence(connectedPermanence_);
  proto.setMinThreshold(minThreshold_);
  proto.setSampleSize(sampleSize_);
  proto.setPermanenceIncrement(permanenceIncrement_);
  proto.setPermanenceDecrement(permanenceDecrement_);
  proto.setPredictedSegmentDecrement(predictedSegmentDecrement_);

  proto.setMaxSegmentsPerCell(maxSegmentsPerCell_);
  proto.setMaxSynapsesPerSegment(maxSynapsesPerSegment_);

  auto _basalConnections = proto.initBasalConnections();
  basalConnections.write(_basalConnections);

  auto _apicalConnections = proto.initApicalConnections();
  apicalConnections.write(_apicalConnections);

  auto random = proto.initRandom();
  rng_.write(random);

  auto activeCells = proto.initActiveCells(activeCells_.size());
  UInt i = 0;
  for (CellIdx cell : activeCells_)
  {
    activeCells.set(i++, cell);
  }

  auto predictedCells = proto.initPredictedActiveCells(predictedCells_.size());
  i = 0;
  for (CellIdx cell : predictedCells_)
  {
    predictedCells.set(i++, cell);
  }

  auto predictedActiveCells = proto.initPredictedActiveCells(predictedActiveCells_.size());
  i = 0;
  for (CellIdx cell : predictedActiveCells_)
  {
    predictedActiveCells.set(i++, cell);
  }

  auto winnerCells = proto.initWinnerCells(winnerCells_.size());
  i = 0;
  for (CellIdx cell : winnerCells_)
  {
    winnerCells.set(i++, cell);
  }


  auto activeBasalSegments = proto.initActiveBasalSegments(
    activeBasalSegments_.size());
  for (UInt i = 0; i < activeBasalSegments_.size(); ++i)
  {
    activeBasalSegments[i].setCell(
      basalConnections.cellForSegment(activeBasalSegments_[i]));
    activeBasalSegments[i].setIdxOnCell(
      basalConnections.idxOnCellForSegment(activeBasalSegments_[i]));
  }

  auto matchingBasalSegments = proto.initMatchingBasalSegments(
    matchingBasalSegments_.size());
  for (UInt i = 0; i < matchingBasalSegments_.size(); ++i)
  {
    matchingBasalSegments[i].setCell(
      basalConnections.cellForSegment(matchingBasalSegments_[i]));
    matchingBasalSegments[i].setIdxOnCell(
      basalConnections.idxOnCellForSegment(matchingBasalSegments_[i]));
  }

  auto activeApicalSegments = proto.initActiveApicalSegments(
    activeApicalSegments_.size());
  for (UInt i = 0; i < activeApicalSegments_.size(); ++i)
  {
    activeApicalSegments[i].setCell(
      apicalConnections.cellForSegment(activeApicalSegments_[i]));
    activeApicalSegments[i].setIdxOnCell(
      apicalConnections.idxOnCellForSegment(activeApicalSegments_[i]));
  }

  auto matchingApicalSegments = proto.initMatchingApicalSegments(
    matchingApicalSegments_.size());
  for (UInt i = 0; i < matchingApicalSegments_.size(); ++i)
  {
    matchingApicalSegments[i].setCell(
      apicalConnections.cellForSegment(matchingApicalSegments_[i]));
    matchingApicalSegments[i].setIdxOnCell(
      apicalConnections.idxOnCellForSegment(matchingApicalSegments_[i]));
  }

  auto basalPotentialOverlaps =
    proto.initNumActivePotentialSynapsesForBasalSegment(
      basalPotentialOverlaps_.size());
  for (Segment segment = 0;
       segment < basalPotentialOverlaps_.size();
       segment++)
  {
    basalPotentialOverlaps[segment].setCell(
      basalConnections.cellForSegment(segment));
    basalPotentialOverlaps[segment].setIdxOnCell(
      basalConnections.idxOnCellForSegment(segment));
    basalPotentialOverlaps[segment].setNumber(
      basalPotentialOverlaps_[segment]);
  }

  auto apicalPotentialOverlaps =
    proto.initNumActivePotentialSynapsesForApicalSegment(
      apicalPotentialOverlaps_.size());
  for (Segment segment = 0;
       segment < apicalPotentialOverlaps_.size();
       segment++)
  {
    apicalPotentialOverlaps[segment].setCell(
      apicalConnections.cellForSegment(segment));
    apicalPotentialOverlaps[segment].setIdxOnCell(
      apicalConnections.idxOnCellForSegment(segment));
    apicalPotentialOverlaps[segment].setNumber(
      apicalPotentialOverlaps_[segment]);
  }


  proto.setIteration(iteration_);

  auto lastUsedIterationForBasalSegment =
    proto.initLastUsedIterationForBasalSegment(lastUsedIterationForBasalSegment_.size());
  for (Segment segment = 0;
       segment < lastUsedIterationForBasalSegment_.size();
       ++segment)
  {
    lastUsedIterationForBasalSegment[segment].setCell(
      basalConnections.cellForSegment(segment));
    lastUsedIterationForBasalSegment[segment].setIdxOnCell(
      basalConnections.idxOnCellForSegment(segment));
    lastUsedIterationForBasalSegment[segment].setNumber(
      lastUsedIterationForBasalSegment_[segment]);
  }

  auto lastUsedIterationForApicalSegment =
    proto.initLastUsedIterationForApicalSegment(lastUsedIterationForApicalSegment_.size());
  for (Segment segment = 0;
       segment < lastUsedIterationForApicalSegment_.size();
       ++segment)
  {
    lastUsedIterationForApicalSegment[segment].setCell(
      apicalConnections.cellForSegment(segment));
    lastUsedIterationForApicalSegment[segment].setIdxOnCell(
      apicalConnections.idxOnCellForSegment(segment));
    lastUsedIterationForApicalSegment[segment].setNumber(
      lastUsedIterationForApicalSegment_[segment]);
  }

  proto.setLearnOnOneCell(learnOnOneCell_);
  auto chosenCellsProto = proto.initChosenCellForColumn(chosenCellForColumn_.size());
  UInt32 chosenIdx = 0;
  for (auto pair : chosenCellForColumn_)
  {
    chosenCellsProto[chosenIdx].setColumnIdx(pair.first);
    chosenCellsProto[chosenIdx].setCellIdx(pair.second);
    ++chosenIdx;
  }
}

// Implementation note: this method sets up the instance using data from
// proto. This method does not call initialize. As such we have to be careful
// that everything in initialize is handled properly here.
void ExtendedTemporalMemory::read(ExtendedTemporalMemoryProto::Reader& proto)
{
  columnCount_ = proto.getColumnCount();
  cellsPerColumn_ = proto.getCellsPerColumn();
  activationThreshold_ = proto.getActivationThreshold();
  initialPermanence_ = proto.getInitialPermanence();
  connectedPermanence_ = proto.getConnectedPermanence();
  minThreshold_ = proto.getMinThreshold();
  sampleSize_ = proto.getSampleSize();
  permanenceIncrement_ = proto.getPermanenceIncrement();
  permanenceDecrement_ = proto.getPermanenceDecrement();
  predictedSegmentDecrement_ = proto.getPredictedSegmentDecrement();

  maxSegmentsPerCell_ = proto.getMaxSegmentsPerCell();
  maxSynapsesPerSegment_ = proto.getMaxSynapsesPerSegment();

  auto _basalConnections = proto.getBasalConnections();
  basalConnections.read(_basalConnections);

  auto _apicalConnections = proto.getApicalConnections();
  apicalConnections.read(_apicalConnections);

  basalOverlaps_.assign(
    basalConnections.segmentFlatListLength(), 0);
  basalPotentialOverlaps_.assign(
    basalConnections.segmentFlatListLength(), 0);

  apicalOverlaps_.assign(
    apicalConnections.segmentFlatListLength(), 0);
  apicalPotentialOverlaps_.assign(
    apicalConnections.segmentFlatListLength(), 0);

  auto random = proto.getRandom();
  rng_.read(random);

  activeCells_.clear();
  for (auto cell : proto.getActiveCells())
  {
    activeCells_.push_back(cell);
  }

  predictedCells_.clear();
  for (auto cell : proto.getPredictedCells())
  {
    predictedCells_.push_back(cell);
  }

  predictedActiveCells_.clear();
  for (auto cell : proto.getPredictedActiveCells())
  {
    predictedActiveCells_.push_back(cell);
  }

  winnerCells_.clear();
  for (auto cell : proto.getWinnerCells())
  {
    winnerCells_.push_back(cell);
  }

  activeBasalSegments_.clear();
  for (auto value : proto.getActiveBasalSegments())
  {
    const Segment segment = basalConnections.getSegment(value.getCell(),
                                                        value.getIdxOnCell());
    activeBasalSegments_.push_back(segment);
  }

  matchingBasalSegments_.clear();
  for (auto value : proto.getMatchingBasalSegments())
  {
    const Segment segment = basalConnections.getSegment(value.getCell(),
                                                        value.getIdxOnCell());
    matchingBasalSegments_.push_back(segment);
  }

  activeApicalSegments_.clear();
  for (auto value : proto.getActiveApicalSegments())
  {
    const Segment segment = apicalConnections.getSegment(value.getCell(),
                                                         value.getIdxOnCell());
    activeApicalSegments_.push_back(segment);
  }

  matchingApicalSegments_.clear();
  for (auto value : proto.getMatchingApicalSegments())
  {
    const Segment segment = apicalConnections.getSegment(value.getCell(),
                                                         value.getIdxOnCell());
    matchingApicalSegments_.push_back(segment);
  }

  basalPotentialOverlaps_.clear();
  basalPotentialOverlaps_.resize(
    basalConnections.segmentFlatListLength());
  for (auto segmentNumPair : proto.getNumActivePotentialSynapsesForBasalSegment())
  {
    const Segment segment = basalConnections.getSegment(
      segmentNumPair.getCell(), segmentNumPair.getIdxOnCell());
    basalPotentialOverlaps_[segment] = segmentNumPair.getNumber();
  }

  apicalPotentialOverlaps_.clear();
  apicalPotentialOverlaps_.resize(
    apicalConnections.segmentFlatListLength());
  for (auto segmentNumPair : proto.getNumActivePotentialSynapsesForApicalSegment())
  {
    const Segment segment = apicalConnections.getSegment(
      segmentNumPair.getCell(), segmentNumPair.getIdxOnCell());
    apicalPotentialOverlaps_[segment] = segmentNumPair.getNumber();
  }

  iteration_ = proto.getIteration();

  lastUsedIterationForBasalSegment_.clear();
  lastUsedIterationForBasalSegment_.resize(basalConnections.segmentFlatListLength());
  for (auto segmentIterationPair : proto.getLastUsedIterationForBasalSegment())
  {
    const Segment segment = basalConnections.getSegment(
      segmentIterationPair.getCell(), segmentIterationPair.getIdxOnCell());
    lastUsedIterationForBasalSegment_[segment] = segmentIterationPair.getNumber();
  }

  lastUsedIterationForApicalSegment_.clear();
  lastUsedIterationForApicalSegment_.resize(apicalConnections.segmentFlatListLength());
  for (auto segmentIterationPair : proto.getLastUsedIterationForApicalSegment())
  {
    const Segment segment = apicalConnections.getSegment(
      segmentIterationPair.getCell(), segmentIterationPair.getIdxOnCell());
    lastUsedIterationForApicalSegment_[segment] = segmentIterationPair.getNumber();
  }

  learnOnOneCell_ = proto.getLearnOnOneCell();
  chosenCellForColumn_.clear();
  for (auto chosenCellProto : proto.getChosenCellForColumn())
  {
    chosenCellForColumn_[chosenCellProto.getColumnIdx()] = chosenCellProto.getCellIdx();
  }
}

//----------------------------------------------------------------------
// Debugging helpers
//----------------------------------------------------------------------

// Print the main TM creation parameters
void ExtendedTemporalMemory::printParameters()
{
  std::cout << "------------CPP ExtendedTemporalMemory Parameters ------------------\n";
  std::cout
    << "version                   = " << EXTENDED_TM_VERSION << std::endl
    << "numColumns                = " << numberOfColumns() << std::endl
    << "cellsPerColumn            = " << getCellsPerColumn() << std::endl
    << "activationThreshold       = " << getActivationThreshold() << std::endl
    << "initialPermanence         = " << getInitialPermanence() << std::endl
    << "connectedPermanence       = " << getConnectedPermanence() << std::endl
    << "minThreshold              = " << getMinThreshold() << std::endl
    << "sampleSize                = " << getSampleSize() << std::endl
    << "learnOnOneCell            = " << getLearnOnOneCell() << std::endl
    << "permanenceIncrement       = " << getPermanenceIncrement() << std::endl
    << "permanenceDecrement       = " << getPermanenceDecrement() << std::endl
    << "predictedSegmentDecrement = " << getPredictedSegmentDecrement() << std::endl;
}
