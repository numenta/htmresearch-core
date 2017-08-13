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
 * Definitions for the Extended Temporal Memory in C++
 */

#ifndef NTA_EXTENDED_TEMPORAL_MEMORY_HPP
#define NTA_EXTENDED_TEMPORAL_MEMORY_HPP

#include <vector>
#include <nupic/types/Serializable.hpp>
#include <nupic/types/Types.hpp>
#include <nupic/utils/Random.hpp>
#include <nupic/algorithms/Connections.hpp>
#include <nupic/proto/ExtendedTemporalMemoryProto.capnp.h>

namespace nupic {
  namespace experimental {
    namespace extended_temporal_memory {

      using namespace algorithms::connections;

      struct PredictionData {
        std::vector<CellIdx> predictedCells;
        std::vector<Segment> activeBasalSegments;
        std::vector<Segment> activeApicalSegments;
        std::vector<Segment> matchingBasalSegments;
        std::vector<Segment> matchingApicalSegments;
        std::vector<UInt> basalOverlaps;
        std::vector<UInt> apicalOverlaps;
        std::vector<UInt> basalPotentialOverlaps;
        std::vector<UInt> apicalPotentialOverlaps;
      };

      /**
       * A fast Temporal Memory implementation with apical dendrites and
       * customizable basal input.
       *
       * The public API uses C arrays, not std::vectors, as inputs. C arrays are
       * a good lowest common denominator. You can get a C array from a vector,
       * but you can't get a vector from a C array without copying it. This is
       * important, for example, when using numpy arrays. The only way to
       * convert a numpy array into a std::vector is to copy it, but you can
       * access a numpy array's internal C array directly.
       */
      class ExtendedTemporalMemory :
        public Serializable<ExtendedTemporalMemoryProto> {
      public:
        ExtendedTemporalMemory();

        /**
         * Initialize the temporal memory (TM) using the given parameters.
         *
         * @param columnCount
         * The number of minicolumns
         *
         * @param basalInputSize
         * The number of bits in the basal input
         *
         * @param apicalInputSize
         * The number of bits in the apical input
         *
         * @param cellsPerColumn
         * Number of cells per column
         *
         * @param activationThreshold
         * If the number of active connected synapses on a segment is at least
         * this threshold, the segment is said to be active.
         *
         * @param initialPermanence
         * Initial permanence of a new synapse.
         *
         * @param connectedPermanence
         * If the permanence value for a synapse is greater than this value, it
         * is said to be connected.
         *
         * @param minThreshold
         * If the number of synapses active on a segment is at least this
         * threshold, it is selected as the best matching cell in a bursting
         * column.
         *
         * @param sampleSize
         * How much of the active SDR to sample with synapses.
         *
         * @param permanenceIncrement
         * Amount by which permanences of synapses are incremented during
         * learning.
         *
         * @param permanenceDecrement
         * Amount by which permanences of synapses are decremented during
         * learning.
         *
         * @param predictedSegmentDecrement
         * Amount by which active permanences of synapses of previously
         * predicted but inactive segments are decremented.
         *
         * @param seed
         * Seed for the random number generator.
         *
         * @param maxSegmentsPerCell
         * The maximum number of segments per cell.
         *
         * @param maxSynapsesPerSegment
         * The maximum number of synapses per segment.
         *
         * Notes:
         *
         * predictedSegmentDecrement: A good value is just a bit larger than
         * (the column-level sparsity * permanenceIncrement). So, if column-level
         * sparsity is 2% and permanenceIncrement is 0.01, this parameter should be
         * something like 4% * 0.01 = 0.0004).
         */
        ExtendedTemporalMemory(
          UInt columnCount,
          UInt basalInputSize,
          UInt apicalInputSize,
          UInt cellsPerColumn = 32,
          UInt activationThreshold = 13,
          Permanence initialPermanence = 0.21,
          Permanence connectedPermanence = 0.50,
          UInt minThreshold = 10,
          UInt sampleSize = 20,
          Permanence permanenceIncrement = 0.10,
          Permanence permanenceDecrement = 0.10,
          Permanence predictedSegmentDecrement = 0.0,
          bool learnOnOneCell = false,
          Int seed = 42,
          UInt maxSegmentsPerCell=255,
          UInt maxSynapsesPerSegment=255,
          bool checkInputs = true);

        virtual void initialize(
          UInt columnCount,
          UInt basalInputSize,
          UInt apicalInputSize,
          UInt cellsPerColumn = 32,
          UInt activationThreshold = 13,
          Permanence initialPermanence = 0.21,
          Permanence connectedPermanence = 0.50,
          UInt minThreshold = 10,
          UInt sampleSize = 20,
          Permanence permanenceIncrement = 0.10,
          Permanence permanenceDecrement = 0.10,
          Permanence predictedSegmentDecrement = 0.0,
          bool learnOnOneCell = false,
          Int seed = 42,
          UInt maxSegmentsPerCell=255,
          UInt maxSynapsesPerSegment=255,
          bool checkInputs = true);

        virtual ~ExtendedTemporalMemory();

        //----------------------------------------------------------------------
        //  Main functions
        //----------------------------------------------------------------------

        /**
         * Get the version number of for the TM implementation.
         *
         * @returns Integer version number.
         */
        virtual UInt version() const;

        /**
         * This *only* updates _rng to a new Random using seed.
         *
         * @returns Integer version number.
         */
        void seed_(UInt64 seed);

        /**
         * Indicates the start of a new sequence.
         * Resets sequence state of the TM.
         */
        virtual void reset();

        /**
         * Perform one time step of the Temporal Memory algorithm. Use the
         * inputs to calculate a set of predicted cells, then calculate the set
         * of active cells.
         *
         * @param activeColumns
         * Sorted list of indices of active columns.
         *
         * @param basalInput
         * Sorted list of active input bits for the basal dendrite segments.
         *
         * @param apicalInput
         * Sorted list of active input bits for the apical dendrite segments.
         *
         * @param basalGrowthCandidates
         * List of bits that the active cells may grow new basal synapses to.
         * In traditional TemporalMemory, this would be the prevWinnerCells.
         *
         * @param apicalGrowthCandidates
         * List of bits that the active cells may grow new apical synapses to.
         *
         * @param learn
         * Whether or not learning is enabled.
         */
        void compute(
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
          bool learn = true);

        /**
         * Perform one time step of the Temporal Memory algorithm. Use the
         * inputs to calculate a set of predicted cells, then calculate the set
         * of active cells.
         *
         * @param activeColumns
         * Sorted list of indices of active columns.
         *
         * @param basalInput
         * Sorted list of active input bits for the basal dendrite segments.
         *
         * @param apicalInput
         * Sorted list of active input bits for the apical dendrite segments.
         *
         * @param basalGrowthCandidates
         * List of bits that the active cells may grow new basal synapses to.
         * In traditional TemporalMemory, this would be the prevWinnerCells.
         *
         * @param apicalGrowthCandidates
         * List of bits that the active cells may grow new apical synapses to.
         *
         * @param learn
         * Whether or not learning is enabled.
         */
        void compute(
          const std::vector<UInt>& activeColumns,
          const std::vector<CellIdx>& basalInput,
          const std::vector<CellIdx>& apicalInput,
          const std::vector<CellIdx>& basalGrowthCandidates,
          const std::vector<CellIdx>& apicalGrowthCandidates,
          bool learn = true);

        /**
         * Equivalent to:
         *
         *  etm.compute(activeColumns,
         *              etm.getActiveCells(),
         *              apicalInput,
         *              etm.getWinnerCells(),
         *              apicalGrowthCandidates);
         *
         * @param activeColumns
         * Sorted list of indices of active columns.
         *
         * @param apicalInput
         * Sorted list of active input bits for the apical dendrite segments.
         *
         * @param apicalGrowthCandidates
         * List of bits that the active cells may grow new apical synapses to.
         *
         * @param learn
         * Whether or not learning is enabled.
         */
        void sequenceMemoryCompute(
          const UInt* activeColumnsBegin,
          const UInt* activeColumnsEnd,
          const CellIdx* apicalInputBegin = nullptr,
          const CellIdx* apicalInputEnd = nullptr,
          const CellIdx* apicalGrowthCandidatesBegin = nullptr,
          const CellIdx* apicalGrowthCandidatesEnd = nullptr,
          bool learn = true);

        /**
         * Equivalent to:
         *
         *  etm.compute(activeColumns,
         *              etm.getActiveCells(),
         *              apicalInput,
         *              etm.getWinnerCells(),
         *              apicalGrowthCandidates);
         *
         * @param activeColumns
         * Sorted list of indices of active columns.
         *
         * @param apicalInput
         * Sorted list of active input bits for the apical dendrite segments.
         *
         * @param apicalGrowthCandidates
         * List of bits that the active cells may grow new apical synapses to.
         *
         * @param learn
         * Whether or not learning is enabled.
         */
        void sequenceMemoryCompute(
          const std::vector<UInt>& activeColumns,
          const std::vector<CellIdx>& apicalInput = {},
          const std::vector<CellIdx>& apicalGrowthCandidates = {},
          bool learn = true);

        // ==============================
        //  Helper functions
        // ==============================

        /**
         * Create a segment on the specified cell. This method calls
         * createSegment on the underlying connections, and it does some extra
         * bookkeeping. Unit tests should call this method, and not
         * connections.createSegment().
         *
         * @param cell
         * Cell to add a segment to.
         *
         * @return Segment
         * The created segment.
         */
        Segment createBasalSegment(CellIdx cell);

        /**
         * Create a segment on the specified cell. This method calls
         * createSegment on the underlying connections, and it does some extra
         * bookkeeping. Unit tests should call this method, and not
         * connections.createSegment().
         *
         * @param cell
         * Cell to add a segment to.
         *
         * @return Segment
         * The created segment.
         */
        Segment createApicalSegment(CellIdx cell);

        /**
         * Returns the indices of cells that belong to a column.
         *
         * @param column Column index
         *
         * @return Cell indices
         */
        std::vector<CellIdx> cellsForColumn(UInt column);

        /**
         * Returns the number of cells in this layer.
         *
         * @return Number of cells
         */
        UInt numberOfCells() const;

        /**
         * Calculate the cells that would be predicted by the given basal and
         * apical input.
         *
         * @param basalInput
         * Sorted list of active input bits for the basal dendrite segments.
         *
         * @param apicalInput
         * Sorted list of active input bits for the apical dendrite segments.
         *
         * @returns
         * Cells/segments that would be predicted.
         */
        PredictionData getPredictionsForInput(
          const CellIdx* basalInputBegin,
          const CellIdx* basalInputEnd,
          const CellIdx* apicalInputBegin,
          const CellIdx* apicalInputEnd) const;

        /**
         * Equivalent to:
         *
         *  etm.getPredictionsForInput(etm.getActiveCells(), apicalInput)
         *
         * @param apicalInput
         * Sorted list of active input bits for the apical dendrite segments.
         *
         * @returns
         * Cells/segments that would be predicted by the current active cells.
         */
        PredictionData getSequenceMemoryPredictions(
          const CellIdx* apicalInputBegin,
          const CellIdx* apicalInputEnd) const;

        /**
        * Returns the indices of the active cells.
        *
        * @returns Vector of indices of active cells.
        */
        std::vector<CellIdx> getActiveCells() const;

        /**
        * Returns the indices of the cells that were predicted.
        *
        * @returns Indices of predicted cells.
        */
        std::vector<CellIdx> getPredictedCells() const;

        /**
         * Returns the indices of the active cells that were predicted.
         *
        * @returns Indices of predicted active cells.
         */
        std::vector<CellIdx> getPredictedActiveCells() const;

        /**
        * Returns the indices of the winner cells.
        *
        * @returns (std::vector<CellIdx>) Vector of indices of winner cells.
        */
        std::vector<CellIdx> getWinnerCells() const;

        std::vector<Segment> getActiveBasalSegments() const;
        std::vector<Segment> getMatchingBasalSegments() const;
        std::vector<Segment> getActiveApicalSegments() const;
        std::vector<Segment> getMatchingApicalSegments() const;

        /**
         * @returns Total number of cells in the basal input.
         */
        UInt getBasalInputSize() const;

        /**
         * @returns Total number of cells in the apical input.
         */
        UInt getApicalInputSize() const;

        /**
         * @returns Total number of minicolumns.
         */
        UInt numberOfColumns() const;

        /**
         * Returns the number of cells per column.
         *
         * @returns Integer number of cells per column
         */
        UInt getCellsPerColumn() const;

        /**
         * Returns the activation threshold.
         *
         * @returns Integer number of the activation threshold
         */
        UInt getActivationThreshold() const;
        void setActivationThreshold(UInt);

        /**
         * Returns the initial permanence.
         *
         * @returns Initial permanence
         */
        Permanence getInitialPermanence() const;
        void setInitialPermanence(Permanence);

        /**
         * Returns the connected permanance.
         *
         * @returns Returns the connected permanance
         */
        Permanence getConnectedPermanence() const;
        void setConnectedPermanence(Permanence);

        /**
         * Returns the minimum threshold.
         *
         * @returns Integer number of minimum threshold
         */
        UInt getMinThreshold() const;
        void setMinThreshold(UInt);

        /**
         * @returns Synapse sample size
         */
        UInt getSampleSize() const;
        void setSampleSize(UInt);

        /**
         * Returns whether to always choose the same cell when bursting a column
         * until the next reset occurs.
         *
         * @returns the learnOnOneCell parameter
         */
        bool getLearnOnOneCell() const;
        void setLearnOnOneCell(bool learnOnOneCell);

        /**
         * Returns the permanence increment.
         *
         * @returns Returns the Permanence increment
         */
        Permanence getPermanenceIncrement() const;
        void setPermanenceIncrement(Permanence);

        /**
         * Returns the permanence decrement.
         *
         * @returns Returns the Permanence decrement
         */
        Permanence getPermanenceDecrement() const;
        void setPermanenceDecrement(Permanence);

        /**
         * Returns the predicted Segment decrement.
         *
         * @returns Returns the segment decrement
         */
        Permanence getPredictedSegmentDecrement() const;
        void setPredictedSegmentDecrement(Permanence);

        /**
         * Returns the maxSegmentsPerCell.
         *
         * @returns Max segments per cell
         */
        UInt getMaxSegmentsPerCell() const;

        /**
         * Returns the maxSynapsesPerSegment.
         *
         * @returns Max synapses per segment
         */
        UInt getMaxSynapsesPerSegment() const;

        /**
         * Returns the checkInputs parameter.
         *
         * @returns the checkInputs parameter
         */
        bool getCheckInputs() const;
        void setCheckInputs(bool checkInputs);

        /**
         * Raises an error if cell index is invalid.
         *
         * @param cell Cell index
         */
        bool _validateCell(CellIdx cell);

        using Serializable::write;
        virtual void write(
          ExtendedTemporalMemoryProto::Builder& proto) const override;

        using Serializable::read;
        virtual void read(ExtendedTemporalMemoryProto::Reader& proto) override;

        bool operator==(const ExtendedTemporalMemory& other);
        bool operator!=(const ExtendedTemporalMemory& other);

        //----------------------------------------------------------------------
        // Debugging helpers
        //----------------------------------------------------------------------

        /**
         * Print the main TM creation parameters
         */
        void printParameters();

        /**
         * Returns the index of the column that a cell belongs to.
         *
         * @param cell Cell index
         *
         * @return (int) Column index
         */
        UInt columnForCell(CellIdx cell);

      protected:

        /**
         * Calculate the active cells, using the current active columns and
         * dendrite segments. Grow and reinforce synapses.
         *
         * @param activeColumnsSize
         * Size of activeColumns.
         *
         * @param activeColumns
         * A sorted list of active column indices.
         *
         * @param reinforceCandidatesExternalBasalSize
         * Size of reinforceCandidatesExternalBasal.
         *
         * @param reinforceCandidatesExternalBasal
         * Sorted list of external cells. Any learning basal dendrite segments
         * will use this list to decide which synapses to reinforce and which
         * synapses to punish. Typically this list should be the
         * 'activeCellsExternalBasal' from the prevous time step.
         *
         * @param reinforceCandidatesExternalApical
         * Size of reinforceCandidatesExternalApical.
         *
         * @param reinforceCandidatesExternalApical
         * Sorted list of external cells. Any learning apical dendrite segments will use
         * this list to decide which synapses to reinforce and which synapses to
         * punish. Typically this list should be the 'activeCellsExternalApical' from
         * the prevous time step.
         *
         * @param growthCandidatesExternalBasal
         * Size of growthCandidatesExternalBasal.
         *
         * @param growthCandidatesExternalBasal
         * Sorted list of external cells. Any learning basal dendrite segments can grow
         * synapses to cells in this list. Typically this list should be a subset of
         * the 'activeCellsExternalBasal' from the previous 'activateDendrites'.
         *
         * @param growthCandidatesExternalApical
         * Size of growthCandidatesExternalApical.
         *
         * @param growthCandidatesExternalApical
         * Sorted list of external cells. Any learning apical dendrite segments can grow
         * synapses to cells in this list. Typically this list should be a subset of
         * the 'activeCellsExternalApical' from the previous 'activateDendrites'.
         *
         * @param learn
         * If true, reinforce / punish / grow synapses.
         */
        void activateCells_(
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
          bool learn);

        /**
         * Calculate dendrite segment activity, using the current active cells.
         *
         * @param activeCellsExternalBasalSize
         * Size of activeCellsExternalBasal.
         *
         * @param activeCellsExternalBasal
         * Sorted list of active external cells for activating basal dendrites.
         *
         * @param activeCellsExternalApicalSize
         * Size of activeCellsExternalApical.
         *
         * @param activeCellsExternalApical
         * Sorted list of active external cells for activating apical dendrites.
         *
         * @param learn
         * If true, segment activations will be recorded. This information is
         * used during segment cleanup.
         */
        void depolarizeCells_(
          const CellIdx* basalInputBegin,
          const CellIdx* basalInputEnd,
          const CellIdx* apicalInputBegin,
          const CellIdx* apicalInputEnd,
          bool learn = true);

      protected:

        UInt columnCount_;
        UInt basalInputSize_;
        UInt apicalInputSize_;
        UInt cellsPerColumn_;
        UInt activationThreshold_;
        UInt minThreshold_;
        UInt sampleSize_;
        bool checkInputs_;
        Permanence initialPermanence_;
        Permanence connectedPermanence_;
        Permanence permanenceIncrement_;
        Permanence permanenceDecrement_;
        Permanence predictedSegmentDecrement_;

        std::vector<CellIdx> activeCells_;
        std::vector<CellIdx> predictedCells_;
        std::vector<CellIdx> predictedActiveCells_;
        std::vector<CellIdx> winnerCells_;

        std::vector<Segment> activeBasalSegments_;
        std::vector<Segment> matchingBasalSegments_;
        std::vector<UInt32> basalOverlaps_;
        std::vector<UInt32> basalPotentialOverlaps_;

        std::vector<Segment> activeApicalSegments_;
        std::vector<Segment> matchingApicalSegments_;
        std::vector<UInt32> apicalOverlaps_;
        std::vector<UInt32> apicalPotentialOverlaps_;

        bool learnOnOneCell_;
        std::map<UInt, CellIdx> chosenCellForColumn_;

        UInt maxSegmentsPerCell_;
        UInt maxSynapsesPerSegment_;
        UInt64 iteration_;
        std::vector<UInt64> lastUsedIterationForBasalSegment_;
        std::vector<UInt64> lastUsedIterationForApicalSegment_;

        Random rng_;

      public:
        Connections basalConnections;
        Connections apicalConnections;
      };

    } // end namespace extended_temporal_memory
  } // end namespace algorithms
} // end namespace nupic

#endif // NTA_EXTENDED_TEMPORAL_MEMORY_HPP
