"""
AttributeMatrixCreationFilterParameter
AttributeMatrixSelectionFilterParameter
AxisAngleFilterParameter
+ BooleanFilterParameter
CalculatorFilterParameter
ChoiceFilterParameter
ComparisonSelectionAdvancedFilterParameter
ComparisonSelectionFilterParameter
ConstrainedDoubleFilterParameter
ConstrainedIntFilterParameter
DataArrayCreationFilterParameter
+ DataArraySelectionFilterParameter
DataContainerArrayProxyFilterParameter
DataContainerCreationFilterParameter
DataContainerGridSelectionFilterParameter
DataContainerReaderFilterParameter
DataContainerSelectionFilterParameter
DoubleFilterParameter
DynamicChoiceFilterParameter
DynamicTableFilterParameter
FileListInfoFilterParameter
FloatFilterParameter
FloatVec2FilterParameter
FloatVec3FilterParameter
FourthOrderPolynomialFilterParameter
GenerateColorTableFilterParameter
ImportHDF5DatasetFilterParameter
+ InputFileFilterParameter
InputPathFilterParameter
* IntFilterParameter
IntVec2FilterParameter
IntVec3FilterParameter
LinkedBooleanFilterParameter
LinkedChoicesFilterParameter
LinkedDataContainerSelectionFilterParameter
MontageSelectionFilterParameter
MontageStructureSelectionFilterParameter
MultiAttributeMatrixSelectionFilterParameter
MultiDataContainerSelectionFilterParameter
MultiDataArraySelectionFilterParameter
NumericTypeFilterParameter
OutputFileFilterParameter
OutputPathFilterParameter
ParagraphFilterParameter
PreflightUpdatedValueFilterParameter
RangeFilterParameter
ReadASCIIDataFilterParameter
ScalarTypeFilterParameter
SecondOrderPolynomialFilterParameter
SeparatorFilterParameter
ShapeTypeSelectionFilterParameter
+ StringFilterParameter
ThirdOrderPolynomialFilterParameter
UInt64FilterParameter
UnknownFilterParameter
MultiInputFileFilterParameter"""


from typing import List, Tuple, Union

from dream3d.Filter import Filter as SIMPLFilter
from dream3d.Filter import FilterDelegatePy as SIMPLFilterDelegatePy
from dream3d.simpl import *
import dream3d.simpl as simpl


class D3DReviewTestFilter(SIMPLFilter):
  def __init__(self) -> None:
    self.int_param: int = 5
    self.dap_param: simpl.DataArrayPath = simpl.DataArrayPath('', '', '')
    self.str_param: str = "Something"
    self.bool_param:bool = False
    self.input_file_param:str = "/No/Path/anywhere.txt"
    self.input_path_param:str = "/Some/Folder/Somewhere"
    self.float_param:float = 33.44

  def _set_float_param(self, value: float) -> None:
      self.float_param = value

  def _get_float_param(self) -> float:
      return self.float_param

  def _set_int(self, value: int) -> None:
    self.int_param = value

  def _get_int(self) -> int:
    return self.int_param

  def _set_dap(self, value: simpl.DataArrayPath) -> None:
    self.dap_param = value

  def _get_dap(self) -> simpl.DataArrayPath:
    return self.dap_param

  def _get_str(self) -> str:
      return self.str_param

  def _set_str(self, value: str) -> None:
      self.str_param = value

  def _get_bool(self) -> bool:
      return self.bool_param

  def _set_bool(self, value ) -> None:
      self.bool_param = value

  def _get_input_file(self) -> str:
      return self.input_file_param

  def _set_input_file(self, value: str ) -> None:
      self.input_file_param = value
  
  def _get_input_path_param(self) -> str:
      return self.input_path_param  

  def _set_input_path_param(self, value: str) -> None:
      self.input_path_param = value

  @staticmethod
  def name() -> str:
    return 'D3DReviewTestFilter'

  @staticmethod
  def uuid() -> str:
    return '{6128bd81-8c49-54c4-b4c7-06fab6ae456a}'

  @staticmethod
  def group_name() -> str:
    return 'Test'

  @staticmethod
  def sub_group_name() -> str:
    return 'Test'

  @staticmethod
  def human_label() -> str:
    return '_Test Filter [Python]'

  @staticmethod
  def version() -> str:
    return '0.1.0'

  @staticmethod
  def compiled_lib_name() -> str:
    return 'DREAM3DReview [Python]'


  def setup_parameters(self) -> List[simpl.FilterParameter]:

    """
    inline const QString Bool("bool");
    inline const QString Float("float");
    inline const QString Double("double");
    inline const QString Int8("int8_t");
    inline const QString UInt8("uint8_t");
    inline const QString Int16("int16_t");
    inline const QString UInt16("uint16_t");
    inline const QString Int32("int32_t");
    inline const QString UInt32("uint32_t");
    inline const QString Int64("int64_t");
    inline const QString UInt64("uint64_t");
    """
    # Create a DataArraySelectionFilterParameter Requirement Type that only takes the following:
    # Image Geometry
    # AttributeMatrix must be a CellType
    # DataArrayTypes are only float, double, int32_t  # See the list from above
    # DataArray Component dimensions are 1 or 3, i.e., a scalar or vector of size 3 only
    req = simpl.DataArraySelectionFilterParameter.RequirementType()
    req.dcGeometryTypes = [simpl.IGeometry.Type.Image]
    req.amTypes = [simpl.AttributeMatrix.Type.Cell]
    req.daTypes = ["float", "double", "int32_t"]
    req.componentDimensions = [VectorSizeT([1]), VectorSizeT([3])]

    return [
        simpl.IntFilterParameter('Integer', 'int_param', self.int_param, simpl.FilterParameter.Category.Parameter, self._set_int, self._get_int, -1),
        simpl.DataArraySelectionFilterParameter('Data Array Path Selection', 'dap_param', self.dap_param, simpl.FilterParameter.Category.RequiredArray, self._set_dap, self._get_dap, req, -1),
        simpl.StringFilterParameter('String', 'str_param', self.str_param, simpl.FilterParameter.Category.Parameter, self._set_str, self._get_str, -1),
        simpl.BooleanFilterParameter('Boolean', 'bool_param', self.bool_param, simpl.FilterParameter.Category.Parameter, self._set_bool, self._get_bool, -1),
        simpl.InputFileFilterParameter('Input File', 'input_file_param', self.input_file_param,
                                    simpl.FilterParameter.Category.Parameter, self._set_input_file, self._get_input_file,'*.ang', 'EDAX Ang', -1),
        simpl.InputPathFilterParameter('Input Directory', 'input_path_param', self.input_path_param,
                                    simpl.FilterParameter.Category.Parameter, self._set_input_path_param, self._get_input_path_param, -1)
    ]

  def data_check(self, dca: simpl.DataContainerArray, delegate: Union[simpl.FilterDelegateCpp, SIMPLFilterDelegatePy] = SIMPLFilterDelegatePy()) -> Tuple[int, str]:
    am = dca.getAttributeMatrix(self.dap_param)

    if am is None:
      return -1, 'AttributeMatrix is None'

    da = am.getAttributeArray(self.dap_param)
    if da is None:
      return -2, 'DataArray is None'

    delegate.notifyStatusMessage('data_check finished!')

    return 0, 'Success'

  def _execute_impl(self, dca: simpl.DataContainerArray, delegate: Union[simpl.FilterDelegateCpp, SIMPLFilterDelegatePy] = SIMPLFilterDelegatePy()) -> Tuple[int, str]:
    delegate.notifyStatusMessage(f'int_param = {self.int_param}')

    da = dca.getAttributeMatrix(self.dap_param).getAttributeArray(self.dap_param)

    data = da.npview()
    delegate.notifyStatusMessage(f'before = {data}')
    data += self.int_param
    delegate.notifyStatusMessage(f'after = {data}')

    delegate.notifyStatusMessage('execute finished!')
    return 0, 'Success'

filters = [D3DReviewTestFilter]