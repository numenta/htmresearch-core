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
 * ---------------------------------------------------------------------
 */

%module(package="htmresearch_core") experimental
%import <nupic/bindings/algorithms.i>

%pythoncode %{
import os

try:
  # NOTE need to import capnp first to activate the magic necessary for
  # ExtendedTemporalMemoryProto_capnp, etc.
  import capnp
except ImportError:
  capnp = None
else:
  from htmresearch_core.proto.ExtendedTemporalMemoryProto_capnp import \
       ExtendedTemporalMemoryProto


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

#include <nupic/experimental/ExtendedTemporalMemory.hpp>

#include <nupic/proto/ExtendedTemporalMemoryProto.capnp.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <nupic/py_support/NumpyVector.hpp>
#if !CAPNP_LITE
#include <nupic/py_support/PyCapnp.hpp>
#endif
#include <nupic/py_support/PythonStream.hpp>
#include <nupic/py_support/PyHelpers.hpp>

// Hack to fix SWIGPY_SLICE_ARG not found bug
#if PY_VERSION_HEX >= 0x03020000
# define SWIGPY_SLICE_ARG(obj) ((PyObject*) (obj))
#else
# define SWIGPY_SLICE_ARG(obj) ((PySliceObject*) (obj))
#endif

using namespace nupic::experimental::extended_temporal_memory;
using namespace nupic;

#define CHECKSIZE(var) \
  NTA_ASSERT(PyArray_DESCR(var)->elsize == 4) << " elsize:" << PyArray_DESCR(var)->elsize

%}

%naturalvar;

//--------------------------------------------------------------------------------
// Extended Temporal Memory
//--------------------------------------------------------------------------------
%inline %{
  template <typename IntType>
  inline PyObject* vectorToList(const vector<IntType> &cellIdxs)
  {
    PyObject *list = PyList_New(cellIdxs.size());
    for (size_t i = 0; i < cellIdxs.size(); i++)
    {
      PyObject *pyIdx = PyInt_FromLong(cellIdxs[i]);
      PyList_SET_ITEM(list, i, pyIdx);
    }
    return list;
  }
%}

%pythoncode %{
  import numpy

  # Without this, Python scripts that haven't imported nupic.bindings.algorithms
  # will get a SwigPyObject rather than a SWIG-wrapped Connections instance
  # when accessing the ExtendedTemporalMemory's connections.
  import nupic.bindings.algorithms

%}

%extend nupic::experimental::extended_temporal_memory::ExtendedTemporalMemory
{
  %pythoncode %{
    def __init__(self,
                 columnDimensions=(2048,),
                 basalInputDimensions=(),
                 apicalInputDimensions=(),
                 cellsPerColumn=32,
                 activationThreshold=13,
                 initialPermanence=0.21,
                 connectedPermanence=0.50,
                 minThreshold=10,
                 maxNewSynapseCount=20,
                 permanenceIncrement=0.10,
                 permanenceDecrement=0.10,
                 predictedSegmentDecrement=0.00,
                 formInternalBasalConnections=True,
                 learnOnOneCell=False,
                 maxSegmentsPerCell=255,
                 maxSynapsesPerSegment=255,
                 seed=42,
                 checkInputs=True):
      self.this = _EXPERIMENTAL.new_ExtendedTemporalMemory(
        columnDimensions, basalInputDimensions, apicalInputDimensions,
        cellsPerColumn, activationThreshold,
        initialPermanence, connectedPermanence,
        minThreshold, maxNewSynapseCount, permanenceIncrement,
        permanenceDecrement, predictedSegmentDecrement,
        formInternalBasalConnections, learnOnOneCell, seed, maxSegmentsPerCell,
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
      self.this = _EXPERIMENTAL.new_ExtendedTemporalMemory()
      if isinstance(state, str):
        self.loadFromString(state)
        self.valueToCategory = {}
      else:
        self.loadFromString(state["this"])
        # Use the rest of the state to set local Python attributes.
        del state["this"]
        self.__dict__.update(state)


    def activateCells(self,
                      activeColumns,
                      reinforceCandidatesExternalBasal=(),
                      reinforceCandidatesExternalApical=(),
                      growthCandidatesExternalBasal=(),
                      growthCandidatesExternalApical=(),
                      learn=True):
      """
      Calculate the active cells, using the current active columns and dendrite
      segments. Grow and reinforce synapses.

      @param activeColumns (sequence)
      A sorted sequence of active column indices.

      @param reinforceCandidatesExternalBasal (sequence)
      Sorted list of external cells. Any learning basal dendrite segments will
      use this list to decide which synapses to reinforce and which synapses to
      punish. Typically this list should be the 'activeCellsExternalBasal' from
      the prevous time step.

      @param reinforceCandidatesExternalApical (sequence)
      Sorted list of external cells. Any learning apical dendrite segments will
      use this list to decide which synapses to reinforce and which synapses to
      punish. Typically this list should be the 'activeCellsExternalApical' from
      the prevous time step.

      @param growthCandidatesExternalBasal (sequence)
      Sorted list of external cells. Any learning basal dendrite segments can
      grow synapses to cells in this list. Typically this list should be a
      subset of the 'activeCellsExternalBasal' from the previous
      'depolarizeCells'.

      @param growthCandidatesExternalApical (sequence)
      Sorted list of external cells. Any learning apical dendrite segments can
      grow synapses to cells in this list. Typically this list should be a
      subset of the 'activeCellsExternalApical' from the previous
      'depolarizeCells'.

      @param learn (boolean)
      Whether to grow / reinforce / punish synapses.
      """
      columnsArray = numpy.array(sorted(activeColumns), "uint32")

      self.convertedActivateCells(
          numpy.asarray(activeColumns, "uint32"),
          numpy.asarray(reinforceCandidatesExternalBasal, "uint32"),
          numpy.asarray(reinforceCandidatesExternalApical, "uint32"),
          numpy.asarray(growthCandidatesExternalBasal, "uint32"),
          numpy.asarray(growthCandidatesExternalApical, "uint32"),
          learn)


    def depolarizeCells(self,
                        activeCellsExternalBasal=(),
                        activeCellsExternalApical=(),
                        learn=True):
      """
      Calculate dendrite segment activity, using the current active cells.

      @param activeCellsExternalBasal (sequence)
      Sorted list of active external cells for activating basal dendrites.

      @param activeCellsExternalApical (sequence)
      Sorted list of active external cells for activating apical dendrites.

      @param learn (bool)
      If true, segment activations will be recorded. This information is used
      during segment cleanup.

      """

      self.convertedDepolarizeCells(
          numpy.asarray(activeCellsExternalBasal, "uint32"),
          numpy.asarray(activeCellsExternalApical, "uint32"),
          learn)


    def compute(self,
                activeColumns,
                activeCellsExternalBasal=(),
                activeCellsExternalApical=(),
                reinforceCandidatesExternalBasal=(),
                reinforceCandidatesExternalApical=(),
                growthCandidatesExternalBasal=(),
                growthCandidatesExternalApical=(),
                learn=True):
      """
      Perform one time step of the Temporal Memory algorithm.

      This method calls activateCells, then calls depolarizeCells. Using the
      TemporalMemory via its compute method ensures that you'll always be able to
      call getPredictiveCells to get predictions for the next time step.

      @param activeColumns (sequence)
      Sorted list of active columns.

      @param activeCellsExternalBasal (sequence)
      Sorted list of active external cells for activating basal dendrites at the
      end of this time step.

      @param activeCellsExternalApical (sequence)
      Sorted list of active external cells for activating apical dendrites at the
      end of this time step.

      @param reinforceCandidatesExternalBasal (sequence)
      Sorted list of external cells. Any learning basal dendrite segments will use
      this list to decide which synapses to reinforce and which synapses to
      punish. Typically this list should be the 'activeCellsExternalBasal' from
      the prevous time step.

      @param reinforceCandidatesExternalApical (sequence)
      Sorted list of external cells. Any learning apical dendrite segments will use
      this list to decide which synapses to reinforce and which synapses to
      punish. Typically this list should be the 'activeCellsExternalApical' from
      the prevous time step.

      @param growthCandidatesExternalBasal (sequence)
      Sorted list of external cells. Any learning basal dendrite segments can grow
      synapses to cells in this list. Typically this list should be a subset of
      the 'activeCellsExternalBasal' from the prevous time step.

      @param growthCandidatesExternalApical (sequence)
      Sorted list of external cells. Any learning apical dendrite segments can grow
      synapses to cells in this list. Typically this list should be a subset of
      the 'activeCellsExternalApical' from the prevous time step.

      @param learn (bool)
      Whether or not learning is enabled
      """

      # Don't call the C++ compute. Implement it in Python. This ensures that we
      # call any Python overrides on the `activateCells` and `depolarizeCells`
      # methods.
      self.activateCells(activeColumns,
                         reinforceCandidatesExternalBasal,
                         reinforceCandidatesExternalApical,
                         growthCandidatesExternalBasal,
                         growthCandidatesExternalApical,
                         learn)
      self.depolarizeCells(activeCellsExternalBasal,
                           activeCellsExternalApical,
                           learn)


    def reset(self):
      _EXPERIMENTAL.ExtendedTemporalMemory_reset(self)


    @classmethod
    def read(cls, proto):
      instance = cls()
      instance.convertedRead(proto)
      return instance

    def write(self, pyBuilder):
      """Serialize the ExtendedTemporalMemory instance using capnp.

      :param: Destination ExtendedTemporalMemoryProto message builder
      """
      reader = ExtendedTemporalMemoryProto.from_bytes(
        self._writeAsCapnpPyBytes()) # copy
      pyBuilder.from_dict(reader.to_dict())  # copy


    def convertedRead(self, proto):
      """Initialize the ExtendedTemporalMemory instance from the given
      ExtendedTemporalMemoryProto reader.

      :param proto: ExtendedTemporalMemoryProto message reader containing data
                    from a previously serialized ExtendedTemporalMemory
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

  inline PyObject* getActiveCells()
  {
    const vector<CellIdx> activeCells = self->getActiveCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      activeCells.size(), activeCells.data()
    ).forPython();
  }

  inline PyObject* getPredictedActiveCells()
  {
    const vector<CellIdx> predictedActiveCells =
      self->getPredictedActiveCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      predictedActiveCells.size(), predictedActiveCells.data()
    ).forPython();
  }

  inline PyObject* getPredictiveCells()
  {
    const vector<CellIdx> predictiveCells = self->getPredictiveCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      predictiveCells.size(), predictiveCells.data()
    ).forPython();
  }

  inline PyObject* getWinnerCells()
  {
    const vector<CellIdx> winnerCells = self->getWinnerCells();

    return nupic::NumpyVectorT<nupic::UInt32>(
      winnerCells.size(), winnerCells.data()
    ).forPython();
  }

  inline PyObject* cellsForColumn(UInt columnIdx)
  {
    const vector<CellIdx> cellIdxs = self->cellsForColumn(columnIdx);
    return vectorToList(cellIdxs);
  }

  inline void convertedActivateCells(
    PyObject *py_activeColumns,
    PyObject *py_reinforceCandidatesExternalBasal,
    PyObject *py_reinforceCandidatesExternalApical,
    PyObject *py_growthCandidatesExternalBasal,
    PyObject *py_growthCandidatesExternalApical,
    bool learn)
  {
    PyArrayObject* _activeColumns =
      (PyArrayObject*) py_activeColumns;
    size_t activeColumnsSize =
      PyArray_DIMS(_activeColumns)[0];
    UInt32* activeColumns =
      (UInt32*)PyArray_DATA(_activeColumns);

    PyArrayObject* _reinforceCandidatesExternalBasal =
      (PyArrayObject*) py_reinforceCandidatesExternalBasal;
    size_t reinforceCandidatesExternalBasalSize =
      PyArray_DIMS(_reinforceCandidatesExternalBasal)[0];
    CellIdx* reinforceCandidatesExternalBasal =
      (CellIdx*)PyArray_DATA(_reinforceCandidatesExternalBasal);

    PyArrayObject* _reinforceCandidatesExternalApical =
      (PyArrayObject*) py_reinforceCandidatesExternalApical;
    size_t reinforceCandidatesExternalApicalSize =
      PyArray_DIMS(_reinforceCandidatesExternalApical)[0];
    CellIdx* reinforceCandidatesExternalApical =
      (CellIdx*)PyArray_DATA(_reinforceCandidatesExternalApical);

    PyArrayObject* _growthCandidatesExternalBasal =
      (PyArrayObject*) py_growthCandidatesExternalBasal;
    size_t growthCandidatesExternalBasalSize =
      PyArray_DIMS(_growthCandidatesExternalBasal)[0];
    CellIdx* growthCandidatesExternalBasal =
      (CellIdx*)PyArray_DATA(_growthCandidatesExternalBasal);

    PyArrayObject* _growthCandidatesExternalApical =
      (PyArrayObject*) py_growthCandidatesExternalApical;
    size_t growthCandidatesExternalApicalSize =
      PyArray_DIMS(_growthCandidatesExternalApical)[0];
    CellIdx* growthCandidatesExternalApical =
      (CellIdx*)PyArray_DATA(_growthCandidatesExternalApical);

    self->activateCells(activeColumnsSize,
                        activeColumns,
                        reinforceCandidatesExternalBasalSize,
                        reinforceCandidatesExternalBasal,
                        reinforceCandidatesExternalApicalSize,
                        reinforceCandidatesExternalApical,
                        growthCandidatesExternalBasalSize,
                        growthCandidatesExternalBasal,
                        growthCandidatesExternalApicalSize,
                        growthCandidatesExternalApical,
                        learn);
  }

  inline void convertedDepolarizeCells(PyObject *py_activeCellsExternalBasal,
                                       PyObject *py_activeCellsExternalApical,
                                       bool learn)
  {
    PyArrayObject* _activeCellsExternalBasal =
      (PyArrayObject*) py_activeCellsExternalBasal;
    size_t activeCellsExternalBasalSize =
      PyArray_DIMS(_activeCellsExternalBasal)[0];
    CellIdx* activeCellsExternalBasal =
      (CellIdx*)PyArray_DATA(_activeCellsExternalBasal);

    PyArrayObject* _activeCellsExternalApical =
      (PyArrayObject*) py_activeCellsExternalApical;
    size_t activeCellsExternalApicalSize =
      PyArray_DIMS(_activeCellsExternalApical)[0];
    CellIdx* activeCellsExternalApical =
      (CellIdx*)PyArray_DATA(_activeCellsExternalApical);

    self->depolarizeCells(activeCellsExternalBasalSize,
                          activeCellsExternalBasal,
                          activeCellsExternalApicalSize,
                          activeCellsExternalApical,
                          learn);
  }
}

%ignore nupic::experimental::extended_temporal_memory::ExtendedTemporalMemory::getActiveCells;
%ignore nupic::experimental::extended_temporal_memory::ExtendedTemporalMemory::getPredictedActiveCells;
%ignore nupic::experimental::extended_temporal_memory::ExtendedTemporalMemory::getPredictiveCells;
%ignore nupic::experimental::extended_temporal_memory::ExtendedTemporalMemory::getWinnerCells;
%ignore nupic::experimental::extended_temporal_memory::ExtendedTemporalMemory::cellsForColumn;

%include <nupic/experimental/ExtendedTemporalMemory.hpp>
