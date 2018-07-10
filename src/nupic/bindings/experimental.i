/* ---------------------------------------------------------------------
 * Numenta Platform for Intelligent Computing (NuPIC)
 * Copyright (C) 2016-2017, Numenta, Inc.  Unless you have an agreement
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
 * ---------------------------------------------------------------------
 */

%module(package="htmresearch_core") experimental
%import <nupic/bindings/algorithms.i>

%pythoncode %{
import os

try:
  # NOTE need to import capnp first to activate the magic necessary for
  # ApicalTiebreakTemporalMemoryProto_capnp, etc.
  import capnp
except ImportError:
  capnp = None
else:
  from htmresearch_core.proto.ApicalTiebreakTemporalMemoryProto_capnp import (
    ApicalTiebreakTemporalMemoryProto, ApicalTiebreakSequenceMemoryProto)


_EXPERIMENTAL = _experimental
%}

%{
/* ---------------------------------------------------------------------
 * Numenta Platform for Intelligent Computing (NuPIC)
 * Copyright (C) 2013-2015, Numenta, Inc.  Unless you have an agreement
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
 * ---------------------------------------------------------------------
 */

#include <Python.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include <nupic/experimental/ApicalTiebreakTemporalMemory.hpp>

#include <nupic/proto/ApicalTiebreakTemporalMemoryProto.capnp.h>

#include <nupic/py_support/NumpyVector.hpp>
#include <nupic/py_support/PyCapnp.hpp>
#include <nupic/py_support/PythonStream.hpp>
#include <nupic/py_support/PyHelpers.hpp>

// Hack to fix SWIGPY_SLICE_ARG not found bug
#if PY_VERSION_HEX >= 0x03020000
# define SWIGPY_SLICE_ARG(obj) ((PyObject*) (obj))
#else
# define SWIGPY_SLICE_ARG(obj) ((PySliceObject*) (obj))
#endif

using namespace nupic::experimental::apical_tiebreak_temporal_memory;
using namespace nupic;

%}

%{
#include <nupic/experimental/SDRSelection.hpp>
%}

%inline {
  PyObject* enumerateDistantSDRsBruteForce(UInt n, UInt w, UInt threshold)
  {
    std::vector<std::vector<UInt> > results =
      nupic::experimental::sdr_selection::enumerateDistantSDRsBruteForce(
      n, w, threshold);

    PyObject* pyResults = PyTuple_New(results.size());
    for (size_t i = 0; i < results.size(); i++)
    {
      PyTuple_SetItem(pyResults, i,
                      nupic::NumpyVectorT<UInt>(results[i].size(),
                                                results[i].data())
                      .forPython());
    }

    return pyResults;
  }
}

//
// Numpy API
//
%{
#include <nupic/py_support/NumpyArrayObject.hpp>
%}
%init %{
  nupic::initializeNumpy();
%}

%naturalvar;


%{
  #include <nupic/experimental/GridUniqueness.hpp>
%}

%pythoncode %{
  def computeGridUniquenessHypercube(A, phaseResolution, ignoredCenterDiameter):
    A = numpy.asarray(A, dtype="float64")

    return _computeGridUniquenessHypercube(A, phaseResolution,
      ignoredCenterDiameter)
%}

%inline {
  PyObject* _computeGridUniquenessHypercube(PyObject* py_A,
                                            Real64 phaseResolution,
                                            Real64 ignoredCenterDiameter)
  {
    PyArrayObject* pyArr_A = (PyArrayObject*)py_A;
    NTA_CHECK(PyArray_NDIM(pyArr_A) == 3);
    npy_intp* npy_dims = PyArray_DIMS(pyArr_A);

    std::vector<std::vector<std::vector<Real64 > > > A;
    for (size_t i = 0; i < npy_dims[0]; i++)
    {
      std::vector<std::vector<Real64> > module;
      for (size_t j = 0; j < npy_dims[1]; j++)
      {
        std::vector<Real64> row;
        for (size_t k = 0; k < npy_dims[2]; k++)
        {
          row.push_back(*(Real64*)PyArray_GETPTR3(pyArr_A, i, j, k));
        }
        module.push_back(row);
      }
      A.push_back(module);
    }

    std::pair<Real64,std::vector<Real64>> result =
      nupic::experimental::grid_uniqueness::computeGridUniquenessHypercube(
        A, phaseResolution, ignoredCenterDiameter);
    PyObject* pyResult = PyTuple_New(2);
    PyTuple_SetItem(pyResult, 0, PyFloat_FromDouble(result.first));
    PyTuple_SetItem(pyResult, 1, nupic::NumpyVectorT<Real64>(result.second.size(),
                                                             result.second.data())
                    .forPython());

    return pyResult;
  }
}

//--------------------------------------------------------------------------------
// Apical Tiebreak Temporal Memory
//--------------------------------------------------------------------------------
%pythoncode %{
  import numpy

  # Without this, Python scripts that haven't imported nupic.bindings.algorithms
  # will get a SwigPyObject rather than a SWIG-wrapped Connections instance
  # when accessing the ApicalTiebreakTemporalMemory's connections.
  import nupic.bindings.algorithms

%}

%extend nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakTemporalMemory
{
  inline PyObject* getActiveCells()
  {
    const std::vector<CellIdx> activeCells = self->getActiveCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      activeCells.size(), activeCells.data()
    ).forPython();
  }


  inline PyObject* getPredictedActiveCells()
  {
    const std::vector<CellIdx> predictedActiveCells =
      self->getPredictedActiveCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      predictedActiveCells.size(), predictedActiveCells.data()
    ).forPython();
  }

  inline PyObject* getPredictedCells()
  {
    const std::vector<CellIdx> predictedCells = self->getPredictedCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      predictedCells.size(), predictedCells.data()
    ).forPython();
  }

  inline PyObject* getWinnerCells()
  {
    const std::vector<CellIdx> winnerCells = self->getWinnerCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      winnerCells.size(), winnerCells.data()
    ).forPython();
  }

  inline PyObject* cellsForColumn(UInt columnIdx)
  {
    const std::vector<CellIdx> cellIdxs = self->cellsForColumn(columnIdx);

    return nupic::NumpyVectorT<nupic::UInt32>(
      cellIdxs.size(), cellIdxs.data()
    ).forPython();
  }
}

%extend nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakPairMemory
{
  %pythoncode %{
    def __init__(self,
                 columnCount=2048,
                 basalInputSize=0,
                 apicalInputSize=0,
                 cellsPerColumn=32,
                 activationThreshold=13,
                 initialPermanence=0.21,
                 connectedPermanence=0.50,
                 minThreshold=10,
                 sampleSize=20,
                 permanenceIncrement=0.10,
                 permanenceDecrement=0.10,
                 basalPredictedSegmentDecrement=0.00,
                 apicalPredictedSegmentDecrement=0.00,
                 learnOnOneCell=False,
                 maxSegmentsPerCell=255,
                 maxSynapsesPerSegment=255,
                 seed=42,
                 checkInputs=True,
                 basalInputPrepend=False):
      """
      @param columnCount (int)
      The number of minicolumns

      @param basalInputSize (sequence)
      The number of bits in the basal input

      @param apicalInputSize (int)
      The number of bits in the apical input

      @param cellsPerColumn (int)
      Number of cells per column

      @param activationThreshold (int)
      If the number of active connected synapses on a segment is at least this
      threshold, the segment is said to be active.

      @param initialPermanence (float)
      Initial permanence of a new synapse

      @param connectedPermanence (float)
      If the permanence value for a synapse is greater than this value, it is said
      to be connected.

      @param minThreshold (int)
      If the number of potential synapses active on a segment is at least this
      threshold, it is said to be "matching" and is eligible for learning.

      @param sampleSize (int)
      How much of the active SDR to sample with synapses.

      @param permanenceIncrement (float)
      Amount by which permanences of synapses are incremented during learning.

      @param permanenceDecrement (float)
      Amount by which permanences of synapses are decremented during learning.

      @param predictedSegmentDecrement (float)
      Amount by which basal segments are punished for incorrect predictions.

      @param learnOnOneCell (bool)
      Whether to always choose the same cell when bursting a column until the
      next reset occurs.

      @param maxSegmentsPerCell (int)
      The maximum number of segments per cell.

      @param maxSynapsesPerSegment (int)
      The maximum number of synapses per segment.

      @param seed (int)
      Seed for the random number generator.

      @param basalInputPrepend (bool)
      If true, this TM will automatically insert its activeCells and winnerCells
      into the basalInput and basalGrowthCandidates, respectively.
      """

      if basalInputPrepend:
        basalInputSize += columnCount * cellsPerColumn

      self.this = _EXPERIMENTAL.new_ApicalTiebreakPairMemory(
        columnCount, basalInputSize, apicalInputSize,
        cellsPerColumn, activationThreshold,
        initialPermanence, connectedPermanence,
        minThreshold, sampleSize, permanenceIncrement,
        permanenceDecrement, basalPredictedSegmentDecrement,
        apicalPredictedSegmentDecrement,
        learnOnOneCell, seed, maxSegmentsPerCell,
        maxSynapsesPerSegment, checkInputs)

      self.basalInputPrepend = basalInputPrepend


    def __getstate__(self):
      # Save the local attributes but override the C++ temporal memory with the
      # string representation.
      d = dict(self.__dict__)
      d["this"] = self.getCState()
      return d


    def __setstate__(self, state):
      # Create an empty C++ temporal memory and populate it from the serialized
      # string.
      self.this = _EXPERIMENTAL.new_ApicalTiebreakPairMemory()
      if isinstance(state, str):
        self.loadFromString(state)
        self.valueToCategory = {}
      else:
        self.loadFromString(state["this"])
        # Use the rest of the state to set local Python attributes.
        del state["this"]
        self.__dict__.update(state)


    def compute(self,
                activeColumns,
                basalInput=(),
                apicalInput=(),
                basalGrowthCandidates=None,
                apicalGrowthCandidates=None,
                learn=True):
      """
      Perform one time step of the Temporal Memory algorithm.

      @param activeColumns (sequence)
      Sorted list of active columns.

      @param basalInput (sequence)
      Sorted list of active input bits for the basal dendrite segments.

      @param apicalInput (sequence)
      Sorted list of active input bits for the apical dendrite segments

      @param basalGrowthCandidates (sequence)
      List of bits that the active cells may grow new basal synapses to.
      If None, the basalInput is assumed to be growth candidates.

      @param apicalGrowthCandidates (sequence)
      List of bits that the active cells may grow new apical synapses to
      If None, the apicalInput is assumed to be growth candidates.

      @param learn (bool)
      Whether or not learning is enabled
      """

      npBasal = numpy.asarray(basalInput, "uint32")
      npApical = numpy.asarray(apicalInput, "uint32")
      npBasalGrowth = (numpy.asarray(basalGrowthCandidates, "uint32")
                       if basalGrowthCandidates is not None
                       else npBasal)
      npApicalGrowth = (numpy.asarray(apicalGrowthCandidates, "uint32")
                        if apicalGrowthCandidates is not None
                        else npApical)

      if self.basalInputPrepend:
        npBasal = numpy.append(self.getActiveCells(),
                               npBasal + self.numberOfCells())
        npBasalGrowth = numpy.append(self.getWinnerCells(),
                                     npBasalGrowth + self.numberOfCells())

      self.convertedCompute(
        numpy.asarray(activeColumns, "uint32"),
        npBasal, npApical, npBasalGrowth, npApicalGrowth,
        learn)


    @classmethod
    def read(cls, proto):
      instance = cls()
      instance.convertedRead(proto)
      return instance

    def write(self, pyBuilder):
      """Serialize the ApicalTiebreakTemporalMemory instance using capnp.

      :param: Destination ApicalTiebreakTemporalMemoryProto message builder
      """
      reader = ApicalTiebreakTemporalMemoryProto.from_bytes(
        self._writeAsCapnpPyBytes()) # copy
      pyBuilder.from_dict(reader.to_dict())  # copy


    def convertedRead(self, proto):
      """Initialize the ApicalTiebreakTemporalMemory instance from the given
      ApicalTiebreakTemporalMemoryProto reader.

      :param proto: ApicalTiebreakTemporalMemoryProto message reader containing data
                    from a previously serialized ApicalTiebreakTemporalMemory
                    instance.

      """
      self._initFromCapnpPyBytes(proto.as_builder().to_bytes()) # copy * 2
  %}

  inline PyObject* _writeAsCapnpPyBytes() const
  {
    return nupic::PyCapnpHelper::writeAsPyBytes(*self);
  }

  inline void _initFromCapnpPyBytes(PyObject* pyBytes)
  {
    nupic::PyCapnpHelper::initFromPyBytes(*self, pyBytes);
  }

  inline void convertedCompute(
    PyObject *py_activeColumns,
    PyObject *py_basalInput,
    PyObject *py_apicalInput,
    PyObject *py_basalGrowthCandidates,
    PyObject *py_apicalGrowthCandidates,
    bool learn)
  {
    nupic::NumpyVectorWeakRefT<nupic::UInt> activeColumns(py_activeColumns);
    nupic::NumpyVectorWeakRefT<nupic::UInt> basalInput(py_basalInput);
    nupic::NumpyVectorWeakRefT<nupic::UInt> apicalInput(py_apicalInput);
    nupic::NumpyVectorWeakRefT<nupic::UInt>
      basalGrowthCandidates(py_basalGrowthCandidates);
    nupic::NumpyVectorWeakRefT<nupic::UInt>
      apicalGrowthCandidates(py_apicalGrowthCandidates);

    self->compute(activeColumns.begin(), activeColumns.end(),
                  basalInput.begin(), basalInput.end(),
                  apicalInput.begin(), apicalInput.end(),
                  basalGrowthCandidates.begin(),
                  basalGrowthCandidates.end(),
                  apicalGrowthCandidates.begin(),
                  apicalGrowthCandidates.end(),
                  learn);
  }

  inline PyObject* getBasalPredictedCells()
  {
    const std::vector<CellIdx> cells = self->getBasalPredictedCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      cells.size(), cells.data()
    ).forPython();
  }

  inline PyObject* getApicalPredictedCells()
  {
    const std::vector<CellIdx> cells = self->getApicalPredictedCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      cells.size(), cells.data()
    ).forPython();
  }
}

%extend nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakSequenceMemory
{
  %pythoncode %{
    def __init__(self,
                 columnCount=2048,
                 apicalInputSize=0,
                 cellsPerColumn=32,
                 activationThreshold=13,
                 initialPermanence=0.21,
                 connectedPermanence=0.50,
                 minThreshold=10,
                 sampleSize=20,
                 permanenceIncrement=0.10,
                 permanenceDecrement=0.10,
                 basalPredictedSegmentDecrement=0.00,
                 apicalPredictedSegmentDecrement=0.00,
                 learnOnOneCell=False,
                 maxSegmentsPerCell=255,
                 maxSynapsesPerSegment=255,
                 seed=42,
                 checkInputs=True,
                 basalInputPrepend=False):
      """
      @param columnCount (int)
      The number of minicolumns

      @param apicalInputSize (int)
      The number of bits in the apical input

      @param cellsPerColumn (int)
      Number of cells per column

      @param activationThreshold (int)
      If the number of active connected synapses on a segment is at least this
      threshold, the segment is said to be active.

      @param initialPermanence (float)
      Initial permanence of a new synapse

      @param connectedPermanence (float)
      If the permanence value for a synapse is greater than this value, it is said
      to be connected.

      @param minThreshold (int)
      If the number of potential synapses active on a segment is at least this
      threshold, it is said to be "matching" and is eligible for learning.

      @param sampleSize (int)
      How much of the active SDR to sample with synapses.

      @param permanenceIncrement (float)
      Amount by which permanences of synapses are incremented during learning.

      @param permanenceDecrement (float)
      Amount by which permanences of synapses are decremented during learning.

      @param predictedSegmentDecrement (float)
      Amount by which basal segments are punished for incorrect predictions.

      @param learnOnOneCell (bool)
      Whether to always choose the same cell when bursting a column until the
      next reset occurs.

      @param maxSegmentsPerCell (int)
      The maximum number of segments per cell.

      @param maxSynapsesPerSegment (int)
      The maximum number of synapses per segment.

      @param seed (int)
      Seed for the random number generator.

      @param basalInputPrepend (bool)
      If true, this TM will automatically insert its activeCells and winnerCells
      into the basalInput and basalGrowthCandidates, respectively.
      """

      self.this = _EXPERIMENTAL.new_ApicalTiebreakSequenceMemory(
        columnCount, apicalInputSize,
        cellsPerColumn, activationThreshold,
        initialPermanence, connectedPermanence,
        minThreshold, sampleSize, permanenceIncrement,
        permanenceDecrement, basalPredictedSegmentDecrement,
        apicalPredictedSegmentDecrement,
        learnOnOneCell, seed, maxSegmentsPerCell,
        maxSynapsesPerSegment, checkInputs)


    def __getstate__(self):
      # Save the local attributes but override the C++ temporal memory with the
      # string representation.
      d = dict(self.__dict__)
      d["this"] = self.getCState()
      return d


    def __setstate__(self, state):
      # Create an empty C++ temporal memory and populate it from the serialized
      # string.
      self.this = _EXPERIMENTAL.new_ApicalTiebreakSequenceMemory()
      if isinstance(state, str):
        self.loadFromString(state)
        self.valueToCategory = {}
      else:
        self.loadFromString(state["this"])
        # Use the rest of the state to set local Python attributes.
        del state["this"]
        self.__dict__.update(state)


    def compute(self,
                activeColumns,
                apicalInput=(),
                apicalGrowthCandidates=None,
                learn=True):
      """
      Perform one time step of the Temporal Memory algorithm.

      @param activeColumns (sequence)
      Sorted list of active columns.

      @param apicalInput (sequence)
      Sorted list of active input bits for the apical dendrite segments

      @param apicalGrowthCandidates (sequence)
      List of bits that the active cells may grow new apical synapses to
      If None, the apicalInput is assumed to be growth candidates.

      @param learn (bool)
      Whether or not learning is enabled
      """

      npApical = numpy.asarray(apicalInput, "uint32")
      npApicalGrowth = (numpy.asarray(apicalGrowthCandidates, "uint32")
                        if apicalGrowthCandidates is not None
                        else npApical)

      self.convertedCompute(
        numpy.asarray(activeColumns, "uint32"),
        npApical, npApicalGrowth,
        learn)


    @classmethod
    def read(cls, proto):
      instance = cls()
      instance.convertedRead(proto)
      return instance

    def write(self, pyBuilder):
      """Serialize the ApicalTiebreakTemporalMemory instance using capnp.

      :param: Destination ApicalTiebreakSequenceMemoryProto message builder
      """
      reader = ApicalTiebreakSequenceMemoryProto.from_bytes(
        self._writeAsCapnpPyBytes()) # copy
      pyBuilder.from_dict(reader.to_dict())  # copy


    def convertedRead(self, proto):
      """Initialize the ApicalTiebreakTemporalMemory instance from the given
      ApicalTiebreakSequenceMemoryProto reader.

      :param proto: ApicalTiebreakSequenceMemoryProto message reader containing data
                    from a previously serialized ApicalTiebreakTemporalMemory
                    instance.

      """
      self._initFromCapnpPyBytes(proto.as_builder().to_bytes()) # copy * 2
  %}

  inline PyObject* _writeAsCapnpPyBytes() const
  {
    return nupic::PyCapnpHelper::writeAsPyBytes(*self);
  }

  inline void _initFromCapnpPyBytes(PyObject* pyBytes)
  {
    nupic::PyCapnpHelper::initFromPyBytes(*self, pyBytes);
  }

  inline void convertedCompute(
    PyObject *py_activeColumns,
    PyObject *py_apicalInput,
    PyObject *py_apicalGrowthCandidates,
    bool learn)
  {
    nupic::NumpyVectorWeakRefT<nupic::UInt> activeColumns(py_activeColumns);
    nupic::NumpyVectorWeakRefT<nupic::UInt> apicalInput(py_apicalInput);
    nupic::NumpyVectorWeakRefT<nupic::UInt>
      apicalGrowthCandidates(py_apicalGrowthCandidates);

    self->compute(activeColumns.begin(), activeColumns.end(),
                  apicalInput.begin(), apicalInput.end(),
                  apicalGrowthCandidates.begin(),
                  apicalGrowthCandidates.end(),
                  learn);
  }

  inline PyObject* getPredictedCells()
  {
    const std::vector<CellIdx> predictedCells = self->getPredictedCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      predictedCells.size(), predictedCells.data()
    ).forPython();
  }

  inline PyObject* getNextPredictedCells()
  {
    const std::vector<CellIdx> predictedCells = self->getNextPredictedCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      predictedCells.size(), predictedCells.data()
    ).forPython();
  }

  inline PyObject* getNextBasalPredictedCells()
  {
    const std::vector<CellIdx> cells = self->getNextBasalPredictedCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      cells.size(), cells.data()
    ).forPython();
  }

  inline PyObject* getNextApicalPredictedCells()
  {
    const std::vector<CellIdx> cells = self->getNextApicalPredictedCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      cells.size(), cells.data()
    ).forPython();
  }
}

%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakTemporalMemory::getActiveCells;
%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakTemporalMemory::getPredictedActiveCells;
%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakTemporalMemory::getPredictedCells;
%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakTemporalMemory::getWinnerCells;
%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakTemporalMemory::cellsForColumn;

%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakSequenceMemory::getPredictedCells;
%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakSequenceMemory::getNextPredictedCells;

%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakPairMemory::getBasalPredictedCells;
%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakPairMemory::getApicalPredictedCells;
%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakSequenceMemory::getNextBasalPredictedCells;
%ignore nupic::experimental::apical_tiebreak_temporal_memory::ApicalTiebreakSequenceMemory::getNextApicalPredictedCells;

%include <nupic/experimental/ApicalTiebreakTemporalMemory.hpp>
