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
 * Declarations for the ApicalTiebreakTemporalMemory
 */

#ifndef NTA_APICAL_TIEBREAK_TM_HPP
#define NTA_APICAL_TIEBREAK_TM_HPP

#include <vector>
#include <nupic/types/Serializable.hpp>
#include <nupic/types/Types.hpp>
#include <nupic/utils/Random.hpp>
#include <nupic/algorithms/Connections.hpp>
#include <nupic/proto/ApicalTiebreakTemporalMemoryProto.capnp.h>

namespace nupic {
  namespace experimental {
    namespace apical_tiebreak_temporal_memory {

      using namespace algorithms::connections;

      /**
       * A fast generalized Temporal Memory implementation with apical dendrites
       * that add a "tiebreak".
       *
       * Basal connections are used to implement traditional Temporal Memory.
       *
       * The apical connections are used for further disambiguation. If multiple
       * cells in a minicolumn have active basal segments, each of those cells
       * is predicted, unless one of them also has an active apical segment, in
       * which case only the cells with active basal and apical segments are
       * predicted.
       *
       * In other words, the apical connections have no effect unless the basal
       * input is a union of SDRs (e.g. from bursting minicolumns).
       *
       * This class is generalized in two ways:
       *
       * - This class does not specify when a 'timestep' begins and ends. It
       *   exposes two main methods: 'depolarizeCells' and 'activateCells', and
       *   callers or subclasses can introduce the notion of a timestep.
       * - This class is unaware of whether its 'basalInput' or 'apicalInput'
       *   are from internal or external cells. They are just cell numbers. The
       *   caller knows what these cell numbers mean, but this class doesn't.
       */
      class ApicalTiebreakTemporalMemory {
      public:
        ApicalTiebreakTemporalMemory();

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
         * @param basalPredictedSegmentDecrement
         * Amount by which active permanences of synapses of previously
         * predicted but inactive segments are decremented.
         *
         * @param apicalPredictedSegmentDecrement
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
        ApicalTiebreakTemporalMemory(
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
          Permanence basalPredictedSegmentDecrement = 0.0,
          Permanence apicalPredictedSegmentDecrement = 0.0,
          bool learnOnOneCell = false,
          Int seed = 42,
          UInt maxSegmentsPerCell=255,
          UInt maxSynapsesPerSegment=255,
          bool checkInputs = true);

        virtual ~ApicalTiebreakTemporalMemory();

        //----------------------------------------------------------------------
        //  Main functions
        //----------------------------------------------------------------------

        /**
         * Reinitialize the random number generator.
         *
         * @returns Integer version number.
         */
        void seed(UInt64 seed);

        /**
         * Clear all active / predicted cells and segments. With
         * 'learnOnOneCell', this also causes the TM to choose a different set
         * of cells to learn on in.
         */
        virtual void reset();

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
        virtual std::vector<CellIdx> getPredictedCells() const;

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
        Permanence getBasalPredictedSegmentDecrement() const;
        void setBasalPredictedSegmentDecrement(Permanence);
        Permanence getApicalPredictedSegmentDecrement() const;
        void setApicalPredictedSegmentDecrement(Permanence);

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

        void write(ApicalTiebreakTemporalMemoryProto::Builder& proto) const;
        void read(ApicalTiebreakTemporalMemoryProto::Reader& proto);

        bool operator==(const ApicalTiebreakTemporalMemory& other);
        bool operator!=(const ApicalTiebreakTemporalMemory& other);

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
         * @param basalReinforceCandidates
         * Sorted list of inputs. Any learning basal dendrite segments will use
         * this list to decide which synapses to reinforce and which synapses to
         * punish. Typically this list should be the 'basalInput' from the
         * prevous time step.
         *
         * @param apicalReinforceCandidates
         * Sorted list of inputs. Any learning apical dendrite segments will use
         * this list to decide which synapses to reinforce and which synapses to
         * punish. Typically this list should be the 'apicalInput' from the
         * prevous time step.
         *
         * @param basalGrowthCandidates
         * Sorted list of inputs. Any learning basal dendrite segments can grow
         * synapses to cells in this list. Typically this list should be a
         * subset of the 'basalInput' from the previous 'depolarizeCells'.
         *
         * @param apicalGrowthCandidates
         * Sorted list of inputs. Any learning apical dendrite segments can grow
         * synapses to cells in this list. Typically this list should be a
         * subset of the 'apicalInput' from the previous 'depolarizeCells'.
         *
         * @param learn
         * If true, reinforce / punish / grow synapses.
         */
        void activateCells(
          const UInt* activeColumnsBegin,
          const UInt* activeColumnsEnd,
          const CellIdx* basalReinforceCandidatesBegin,
          const CellIdx* basalReinforceCandidatesEnd,
          const CellIdx* apicalReinforceCandidatesBegin,
          const CellIdx* apicalReinforceCandidatesEnd,
          const CellIdx* basalGrowthCandidatesBegin,
          const CellIdx* basalGrowthCandidatesEnd,
          const CellIdx* apicalGrowthCandidatesBegin,
          const CellIdx* apicalGrowthCandidatesEnd,
          bool learn);

        /**
         * Calculate dendrite segment activity, using the current active cells.
         *
         * @param basalInput
         * Sorted list of active inputs for activating basal dendrites.
         *
         * @param apicalInput
         * Sorted list of active inputs for activating apical dendrites.
         *
         * @param learn
         * If true, segment activations will be recorded. This information is
         * used during segment cleanup.
         */
        void depolarizeCells(
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
        Permanence basalPredictedSegmentDecrement_;
        Permanence apicalPredictedSegmentDecrement_;

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

      /**
       * Associates the basal input with the active columns as a "pair".
       */
      class ApicalTiebreakPairMemory :
        public ApicalTiebreakTemporalMemory,
        public Serializable<ApicalTiebreakTemporalMemoryProto> {
      public:
        ApicalTiebreakPairMemory(
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
          Permanence basalPredictedSegmentDecrement = 0.0,
          Permanence apicalPredictedSegmentDecrement = 0.0,
          bool learnOnOneCell = false,
          Int seed = 42,
          UInt maxSegmentsPerCell=255,
          UInt maxSynapsesPerSegment=255,
          bool checkInputs = true);

       /**
         * Perform one timestep. Use the basal and apical input to form a set of
         * predictions, then activate the specified columns, then learn.
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
         * Perform one timestep. Use the basal and apical input to form a set of
         * predictions, then activate the specified columns, then learn.
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
         * Returns the cells with active basal segments.
         */
        std::vector<CellIdx> getBasalPredictedCells() const;

        /**
         * Returns the cells with active apical segments.
         */
        std::vector<CellIdx> getApicalPredictedCells() const;

        using Serializable::write;
        virtual void write(
          ApicalTiebreakTemporalMemoryProto::Builder& proto) const override;

        using Serializable::read;
        virtual void read(ApicalTiebreakTemporalMemoryProto::Reader& proto) override;
      };

      /**
       * Traditional TM sequence memory, with apical tiebreak.
       *
       * This exposes a 'getNextPredictedCells' method which predicts the next
       * input.
       */
      class ApicalTiebreakSequenceMemory :
        public ApicalTiebreakTemporalMemory,
        public Serializable<ApicalTiebreakSequenceMemoryProto> {
      public:
        ApicalTiebreakSequenceMemory();
        ApicalTiebreakSequenceMemory(
          UInt columnCount,
          UInt apicalInputSize = 0,
          UInt cellsPerColumn = 32,
          UInt activationThreshold = 13,
          Permanence initialPermanence = 0.21,
          Permanence connectedPermanence = 0.50,
          UInt minThreshold = 10,
          UInt sampleSize = 20,
          Permanence permanenceIncrement = 0.10,
          Permanence permanenceDecrement = 0.10,
          Permanence basalPredictedSegmentDecrement = 0.0,
          Permanence apicalPredictedSegmentDecrement = 0.0,
          bool learnOnOneCell = false,
          Int seed = 42,
          UInt maxSegmentsPerCell=255,
          UInt maxSynapsesPerSegment=255,
          bool checkInputs = true);

        /**
         * Perform one timestep. Activate the specified columns, using the
         * predictions from the previous timestep, then learn. Then form a new
         * set of predictions using the new active cells and the apicalInput.
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
        void compute(
          const UInt* activeColumnsBegin,
          const UInt* activeColumnsEnd,
          const CellIdx* apicalInputBegin = nullptr,
          const CellIdx* apicalInputEnd = nullptr,
          const CellIdx* apicalGrowthCandidatesBegin = nullptr,
          const CellIdx* apicalGrowthCandidatesEnd = nullptr,
          bool learn = true);

        /**
         * Perform one timestep. Activate the specified columns, using the
         * predictions from the previous timestep, then learn. Then form a new
         * set of predictions using the new active cells and the apicalInput.
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
        void compute(
          const std::vector<UInt>& activeColumns,
          const std::vector<CellIdx>& apicalInput = {},
          const std::vector<CellIdx>& apicalGrowthCandidates = {},
          bool learn = true);

        /**
         * Returns the indices of the cells that were predicted.
         *
         * @returns Indices of predicted cells.
         */
        virtual std::vector<CellIdx> getPredictedCells() const override;

        /**
         * Returns the indices of the cells that are predicted for the next
         * timestep.
         *
         * @returns Indices of predicted cells.
         */
        std::vector<CellIdx> getNextPredictedCells() const;

        /**
         * Returns the cells with active basal segments.
         */
        std::vector<CellIdx> getNextBasalPredictedCells() const;

        /**
         * Returns the cells with active apical segments.
         */
        std::vector<CellIdx> getNextApicalPredictedCells() const;

        /**
         * Indicates the start of a new sequence.
         */
        virtual void reset() override;

        using Serializable::write;
        virtual void write(
          ApicalTiebreakSequenceMemoryProto::Builder& proto) const override;

        using Serializable::read;
        virtual void read(ApicalTiebreakSequenceMemoryProto::Reader& proto) override;

      protected:
        std::vector<CellIdx> prevApicalInput_;
        std::vector<CellIdx> prevApicalGrowthCandidates_;
        std::vector<CellIdx> prevPredictedCells_;
      };

    } // end namespace apical_tiebreak_temporal_memory
  } // end namespace experimental
} // end namespace nupic

#endif // NTA_APICAL_TIEBREAK_TM_HPP
