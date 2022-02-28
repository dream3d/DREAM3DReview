#This filter takes a series of TDMS files and converts them to individual part files, where each part is stored in an HDF5 File.

import numpy as np
import os
import h5py
import nptdms
import re

from contextlib import ExitStack
from pathlib import Path
from typing import Any, AnyStr, Dict, Generator, List, Match, Pattern
from enum import IntEnum
from typing import List, Tuple, Union

from dream3d.Filter import Filter, FilterDelegatePy
from dream3d.simpl import *

VERSION_KEY: str = 'Version'
SLICES_KEY: str = 'TDMSData'
INDEX_KEY: str = 'Index'
LAYER_START_TIME_KEY: str = 'LayerStartTime'
LAYER_END_TIME_KEY: str = 'LayerEndTime'
PART_START_TIME_KEY: str = 'PartStartTime'
PART_END_TIME_KEY: str = 'PartEndTime'
TDMS_GROUP_NAME_KEY: str = 'TDMS_GroupName'
VERTICES_KEY: str = 'Vertices'
  
def get_file_list(end_index: int, file_extension: str, file_prefix: str, file_suffix: str, increment: int, input_path: str, ordering: int, padding_digits: int, start_index: int, status_delegate: Union[FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> None:
  file_list = []

  index = 0
  i = 0
  while i < (end_index - start_index) + 1:

    if(ordering):
      index = end_index - i
    else:
      index = start_index + i

    file = input_path + '/' + file_prefix + str(index).zfill(padding_digits) + file_suffix + '.' + file_extension
    file_list.append(file)
    i = i + increment

  return file_list

def tdms2h5(file_info: StackFileListInfo, file_list: List, output_dir: Path, area_offset: int, intensity_offset: int, laser_offset: int, groups: List[str] = [], status_delegate: Union[FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> None:

  largest_offset = max(area_offset, intensity_offset)

  with ExitStack() as exitStack:
    h5_files: Dict[str, h5py.File] = {}
    slice_indices: List[int] = []

    regex_name: Pattern[AnyStr] = re.compile(fr'{file_info.FilePrefix}(\d+)')

    path: Path

    paths = []
    for file in file_list:
      paths.append(Path(file))

    paths_size = len(paths)
    path_ind = 1
    for path in paths:
      status_delegate.notifyStatusMessage(f'Converting {path_ind} of {paths_size}: \"{path}\"')
      path_ind = path_ind + 1

      match: Match[AnyStr] = regex_name.search(path.stem)
      slice_index = int(match.group(1))
      slice_indices.append(slice_index)

      with nptdms.TdmsFile(path) as tdmsFile:
        build_firmware_version: int = 3
        try:
          status_delegate.notifyStatusMessage(f'Trying to find old style Bitgain property..')
          bitgain_os_1: float = tdmsFile.properties['Bitgain OS 1']
          bitgain_os_2: float = tdmsFile.properties['Bitgain OS 2']
        except:
          try:
            status_delegate.notifyStatusMessage(f'Trying to find new style Bitgain property..')
            build_firmware_version = 4
            bitgain_os_1: float = tdmsFile.properties['Scanner.Bitgain.OS1']
            bitgain_os_2: float = tdmsFile.properties['Scanner.Bitgain.OS2']
          except:
            status_delegate.notifyStatusMessage(f'Bitgain Property not found. Unable to determine firmware version. Exiting Function.')
            return

        group: nptdms.TdmsGroup
        for ind, group in enumerate(tdmsFile.groups()):
          if groups and not any(re.match(pattern, group.name) for pattern in groups):
            if ind == len(tdmsFile.groups())-1:
              status_delegate.notifyStatusMessage(f' Group(s) not located')
            continue

          output_file_path = output_dir / f'{group.name}.h5'
          if group.name not in h5_files:
            h5_files[group.name] = exitStack.enter_context(h5py.File(output_file_path, 'w'))
            h5_file = h5_files[group.name]
            h5_file.attrs[VERSION_KEY] = build_firmware_version
            h5_group = h5_file.create_group(SLICES_KEY)
            h5_group.attrs[TDMS_GROUP_NAME_KEY] = group.name
          h5_file = h5_files[group.name]
          h5_group: h5py.Group = h5_file[SLICES_KEY].create_group(str(slice_index))

          layer_replacements = {
            'StartTime': LAYER_START_TIME_KEY,
            'EndTime': LAYER_END_TIME_KEY
          }
          _write_tdms_properties(h5_group, tdmsFile.properties, layer_replacements)

          part_replacements = {
            'StartTime' : PART_START_TIME_KEY,
            'EndTime' : PART_END_TIME_KEY,
            'Part.StartTime': PART_START_TIME_KEY,
            'Part.EndTime' : PART_END_TIME_KEY
          }
          _write_tdms_properties(h5_group, group.properties, part_replacements)

          # LaserTTL only uses laser_offset. The end has to be adjusted to make the resulting array consistent
          laser_channel: nptdms.TdmsChannel = group['LaserTTL']
          laser_end_index: int = len(laser_channel) - (laser_offset)
          h5_group.create_dataset(laser_channel.name, data=laser_channel[laser_offset: laser_end_index])

          # Intensity and Area use laser_offset
          area_channel: nptdms.TdmsChannel = group['Area']
          # At this point for illustrative purposes, since the laser_offset is always the largest
          end_index: int = len(area_channel) - (area_offset)
          h5_group.create_dataset(area_channel.name, data=area_channel[area_offset: end_index])

          intensity_channel: nptdms.TdmsChannel = group['Intensity']
          end_index: int = len(intensity_channel) - (intensity_offset)
          h5_group.create_dataset(intensity_channel.name, data=intensity_channel[intensity_offset: end_index])

          # Have not figured out how to correlate parameter to the actual parameter used, just use the same as Laser TTL since it is a machine setting
          parameter_channel: nptdms.TdmsChannel = group['Parameter']
          h5_group.create_dataset(parameter_channel.name, data=parameter_channel[:])

          # X and Y channels just adjust the maximum
          x_channel: nptdms.TdmsChannel = group['X-Axis']
          x_dataset = h5_group.create_dataset(x_channel.name, data=(x_channel[:] / bitgain_os_1), dtype=np.float32)
          x_dataset.attrs['Units'] = 'μm'

          y_channel: nptdms.TdmsChannel = group['Y-Axis']
          y_dataset = h5_group.create_dataset(y_channel.name, data=(y_channel[:] / bitgain_os_2), dtype=np.float32)
          y_dataset.attrs['Units'] = 'μm'

          # Resulting slices will be aligned with the same number of data points for each channel
          
    if not h5_files:
      status_delegate.notifyStatusMessage(f' No TDMS files located')
      return -1

    slice_indices = sorted(slice_indices)

    for h5_file in h5_files.values():
      index_dataset = np.zeros((len(slice_indices), 3), dtype=np.int64)
      for i, index in enumerate(slice_indices):
        index_dataset[i][0] = index
        index_dataset[i][1] = h5_file[SLICES_KEY][str(index)].attrs['layerThickness']
        index_dataset[i][2] = h5_file[SLICES_KEY][str(index)]['X-Axis'].size
      dataset: h5py.Dataset = h5_file.create_dataset(INDEX_KEY, data=index_dataset)
      dataset.attrs['Column0'] = 'SliceIndex'
      dataset.attrs['Column1'] = 'LayerThickness (μm)'
      dataset.attrs['Column2'] = 'NumVertices'
      
    status_delegate.notifyStatusMessage('\nWrote files:')
    h5_file: h5py.File
    for h5_file in h5_files.values():
      status_delegate.notifyStatusMessage(f' \"{h5_file.filename}\"')

def _write_tdms_properties(h5_group: h5py.Group, tdms_dict: Dict[str, Any], replacements: Dict[str, str]) -> None:
  key: str
  value: Any
  for key, value in tdms_dict.items():
    if key in replacements:
      key = replacements[key]
    if isinstance(value, np.datetime64):
      h5_group.attrs[key] = str(np.datetime_as_string(value, unit='us', timezone='UTC'))
    else:
      h5_group.attrs[key] = value

class TDMStoH5(Filter):
  def __init__(self) -> None:
    self.area_offset: int = 0
    self.intensity_offset: int = 0
    self.laser_offset: int = 0
    self.output_folder: str = ''
    self.group: str = ''
    self.stack_info_param: StackFileListInfo = StackFileListInfo(0, 0, 0, 0, 1, '', 'Prefix_', '_Suffix', 'tdms')

  def _set_area_offset(self, value: int) -> None:
    self.area_offset = value

  def _get_area_offset(self) -> int:
    return self.area_offset

  def _set_intensity_offset(self, value: int) -> None:
    self.intensity_offset = value

  def _get_intensity_offset(self) -> int:
    return self.intensity_offset

  def _set_laser_offset(self, value: int) -> None:
    self.laser_offset = value

  def _get_laser_offset(self) -> int:
    return self.laser_offset

  def _set_output_folder(self, value: str) -> None:
    self.output_folder = value

  def _get_output_folder(self) -> str:
    return self.output_folder

  def _set_group(self, value: str) -> None:
    self.group = value

  def _get_group(self) -> str:
    return self.group
  
  def _get_stack_info_param(self) -> StackFileListInfo:
    return self.stack_info_param

  def _set_stack_info_param(self, value) -> None:
    self.stack_info_param = value

  @staticmethod
  def name() -> str:
    return 'TDMS to H5'

  @staticmethod
  def uuid() -> str:
    return '{8b069c55-6d94-4db1-b012-cdfb7dd08ad6}'

  @staticmethod
  def group_name() -> str:
    return 'Example'

  @staticmethod
  def sub_group_name() -> str:
    return 'Sub Example'

  @staticmethod
  def human_label() -> str:
    return 'Convert TDMS to HDF5'

  @staticmethod
  def version() -> str:
    return '1.0.0'

  @staticmethod
  def compiled_lib_name() -> str:
    return 'Python'

  def setup_parameters(self) -> List[FilterParameter]:
    req = MultiDataArraySelectionFilterParameter.RequirementType([IGeometry.Type.Image], [AttributeMatrix.Type.Cell], [], [])
    return [
      IntFilterParameter('Area Offset', 'area_offset', self.area_offset, FilterParameter.Category.Parameter, self._set_area_offset, self._get_area_offset, -1),
      IntFilterParameter('Intensity Offset', 'intensity_offset', self.intensity_offset, FilterParameter.Category.Parameter, self._set_intensity_offset, self._get_intensity_offset, -1),
      IntFilterParameter('Laser Offset', 'laser_offset', self.laser_offset, FilterParameter.Category.Parameter, self._set_laser_offset, self._get_laser_offset, -1),
      OutputPathFilterParameter('Output Folder', 'Output Folder', '', FilterParameter.Category.Parameter,
                               self._set_output_folder, self._get_output_folder, -1),
      StringFilterParameter('Group', 'group', self.group, FilterParameter.Category.Parameter, self._set_group,
                           self._get_group, -1),
      FileListInfoFilterParameter('File List', 'filelist', self.stack_info_param, FilterParameter.Category.Parameter,
                           self._set_stack_info_param, self._get_stack_info_param)
    ]

  def data_check(self, dca: DataContainerArray, status_delegate: Union[FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> Tuple[int, str]:
    if not self.stack_info_param.InputPath:
      return (-5550, 'An input folder must be selected')

    if not self.output_folder:
      return (-301, 'An output folder must be selected')

    if not os.path.exists(self.stack_info_param.InputPath):
      return (-302, f' Input path {self.input_folder} does not exist')

    if not os.path.exists(self.output_folder):
      status_delegate.setWarningCondition(1, f' Output folder {self.output_folder} does not exist; creating a new folder')

    return (0, 'Success')

  def _execute_impl(self, dca: DataContainerArray, status_delegate: Union[FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> Tuple[int, str]:
    
    file_list = get_file_list(self.stack_info_param.EndIndex, self.stack_info_param.FileExtension, self.stack_info_param.FilePrefix, self.stack_info_param.FileSuffix, self.stack_info_param.IncrementIndex, self.stack_info_param.InputPath, self.stack_info_param.Ordering, self.stack_info_param.PaddingDigits, self.stack_info_param.StartIndex)

    if not os.path.exists(self.output_folder):
      os.mkdir(self.output_folder)

    if tdms2h5(self.stack_info_param, file_list, Path(self.output_folder), self.area_offset, self.intensity_offset, self.laser_offset, self.group, status_delegate) == -1:
      return (-1, 'could not create HDF5 files')

    return (0, 'Success')


filters = [TDMStoH5]
