import numpy as np

from typing import List
from enum import IntEnum
from typing import List, Tuple, Union

from dream3d.Filter import Filter, FilterDelegatePy
from dream3d.simpl import NumericTypes
from dream3d.simpl import *


class NormalizeGrayscale(Filter):
    def __init__(self) -> None:
        self.selected_data_array_path: DataArrayPath = DataArrayPath('', '', '')
        self.axes: int = 0
        self.created_data_array_path: DataArrayPath = DataArrayPath('', '', '')

    def _set_selected_data_array_path(self, value: DataArrayPath) -> None:
        self.selected_data_array_path = value

    def _get_selected_data_array_path(self) -> DataArrayPath:
        return self.selected_data_array_path

    def _set_axes(self, value: int) -> None:
        self.axes = value

    def _get_axes(self) -> int:
        return self.axes

    def _set_created_data_array_path(self, value: DataArrayPath) -> None:
        self.created_data_array_path = value

    def _get_created_data_array_path(self) -> DataArrayPath:
        return self.created_data_array_path

    @staticmethod
    def name() -> str:
        return 'Grayscale Normalizer [Python]'

    @staticmethod
    def uuid() -> str:
        return '{b98fa052-4c74-4d25-96ce-5b95074fcecf}'

    @staticmethod
    def group_name() -> str:
        return 'Processing'

    @staticmethod
    def sub_group_name() -> str:
        return 'Normalization'

    @staticmethod
    def human_label() -> str:
        return 'Normalize Grayscale Values [Python]'

    @staticmethod
    def version() -> str:
        return '1.0.0'

    @staticmethod
    def compiled_lib_name() -> str:
        return 'Python'

    def setup_parameters(self) -> List[FilterParameter]:
        # Create a DataArraySelectionFilterParameter Requirement Type that only takes the following:
        # Image Geometry
        # AttributeMatrix must be a CellType
        # DataArray Component dimensions are 1, i.e., a scalar value
        req = DataArraySelectionFilterParameter.RequirementType()
        req.dcGeometryTypes = [IGeometry.Type.Image]
        req.amTypes = [AttributeMatrix.Type.Cell]
        # req.daTypes = ["float", "double", "int32_t"]
        req.componentDimensions = [VectorSizeT([1])]

        return [
            DataArraySelectionFilterParameter('Input Data Array', 'selected_data_array_path',
                                              self.selected_data_array_path, FilterParameter.Category.RequiredArray,
                                              self._set_selected_data_array_path, self._get_selected_data_array_path,
                                              req, -1),
            ChoiceFilterParameter('Slice Axis', 'axes', self.axes, FilterParameter.Category.Parameter, self._set_axes,
                                  self._get_axes, ["x", "y", "z"], False, -1),
            DataArrayCreationFilterParameter('New Array', 'new_array', self.created_data_array_path,
                                             FilterParameter.Category.CreatedArray, self._set_created_data_array_path,
                                             self._get_created_data_array_path,
                                             DataArrayCreationFilterParameter.RequirementType(), -1)
        ]

    def data_check(self, dca: DataContainerArray,
                   status_delegate: Union[FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> Tuple[int, str]:
        # Methodically check for each level of the selected data array path.
        # Ensure the proper geometry is available
        # Ensure the number of components of the input data array is 1
        # Error messages should give as much information as possible. users will rely on the
        #  error message to debug their pipelines, sometimes from the command line so printing
        #  out the DataArrayPath that was passed in is helpful for the user.
        # Error Codes: Should start at a unique value somewhere between 10,000 and 100,000 and then just increment up from there

        dc: DataContainer = dca.getDataContainer(self.selected_data_array_path)
        if dc is None:
            return -61550, 'DataContainer does not exist: ' + self.selected_data_array_path.DataContainerName

        if not dc.Geometry:
            return (-61551, 'Selected DataContainer has no geometry.')

        if dc.Geometry.getGeometryType() != IGeometry.Type.Image:
            return (-61552,
                    'This filter requires a DataContainer with an ImageGeometry. Selected DataContainer has ' + dc.Geometry.getGeometryTypeAsString())

        am: AttributeMatrix = dca.getAttributeMatrix(self.selected_data_array_path)
        if am is None:
            return -61553, 'AttributeMatrix does not exist: ' + self.selected_data_array_path.AttributeMatrixName

        selected_data_array: IDataArray = am.getAttributeArray(self.selected_data_array_path.DataArrayName)
        if selected_data_array is None:
            return -61554, 'DataArray does not exist: ' + self.selected_data_array_path.DataArrayName

        if selected_data_array.getNumberOfComponents() != 1:
            return -61555, 'Selected DataArray must have a single component: ' + str(
                selected_data_array.getNumberOfComponents())

        # If we are creating an array, what other options are there, then we need to create the array (but not allocate the memory during preflight)
        # using the createNonPrereqArrayFromPath should do this for us.
        arr_dtype = selected_data_array.dtype
        if arr_dtype == "int8":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.Int8, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "uint8":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.UInt8, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "int16":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.Int16, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "uint16":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.UInt16, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "int32":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.Int32, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "uint32":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.UInt32, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "int64":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.Int64, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "uint64":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.UInt64, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "float32":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.Float, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "float64":
            dca.createNonPrereqArrayFromPath(status_delegate, NumericTypes.Double, self.created_data_array_path, 0,
                                             VectorSizeT([1]))
        elif arr_dtype == "bool":
            return -61556, 'Bool type arrays are not supported as inputs to this algorithm'

        return (0, 'Success')

    def _execute_impl(self, dca: DataContainerArray,
                      status_delegate: Union[FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> Tuple[
        int, str]:

        dc: DataContainer = dca.getDataContainer(self.selected_data_array_path)
        am: AttributeMatrix = dc.getAttributeMatrix(self.selected_data_array_path)
        selected_data_array: IDataArray = am.getAttributeArray(self.selected_data_array_path)
        created_data_array: IDataArray = am.getAttributeArray(self.created_data_array_path)
        status_delegate.notifyStatusMessage(f" Selected D3D data array is {selected_data_array}")

        # Note that the dimensions come back as X, Y, Z (Fastest to slowest) so we will need to invert the order if using them in a numpy array
        # Note however that the data is actually stored in the SIMPL DataArray in the correct order
        imageGeomDims: SizeVec3 = dc.Geometry.getDimensions()
        shape = [np.int64(imageGeomDims[2]), np.int64(imageGeomDims[1]), np.int64(imageGeomDims[0])]

        npv_selected_data = selected_data_array.npview()
        npv_selected_data = np.reshape(npv_selected_data, shape)
        selected_data_dtype = npv_selected_data.dtype

        # Get a numpy view into the data and then reshape the numpy array
        npv_created_data_array = created_data_array.npview()
        npv_created_data_array = np.reshape(npv_created_data_array, shape)

        # Get the average of the entire data set
        arr_avg = np.average(npv_selected_data)
        slice_axis = self.axes

        # slices are in x-axis
        if slice_axis == 0:
            layer_avg = np.average(npv_selected_data, axis=(0, 1))
            correction_factor = layer_avg / arr_avg
            numpy_created_data_array = npv_selected_data / correction_factor[None, None, :]  # Makes a copy

        # slices are in y-axis
        elif slice_axis == 1:
            layer_avg = np.average(npv_selected_data, axis=(0, 2))
            correction_factor = layer_avg / arr_avg
            numpy_created_data_array = npv_selected_data / correction_factor[None, :, None]  # Makes a copy

        # slices are in z-axis
        elif slice_axis == 2:
            layer_avg = np.average(npv_selected_data, axis=(1, 2))
            correction_factor = layer_avg / arr_avg
            numpy_created_data_array = npv_selected_data / correction_factor[:, None, None]  # Makes a copy

        # Assigns the data from the numpy array numpy_created_data_array to the view into the SIMPL.DataArray
        # through a copy operation.
        # https://numpy.org/doc/stable/user/basics.indexing.html#assigning-values-to-indexed-arrays
        npv_created_data_array[...] = numpy_created_data_array

        return (0, 'Success')


filters = [NormalizeGrayscale]
