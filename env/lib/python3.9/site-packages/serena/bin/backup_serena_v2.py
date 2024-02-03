"""
File that defines the main RNA sequence data
"""


from enum import Enum
from attrs import define, field
from collections import namedtuple
from typing import List, Dict, Any,TypeVar, Type
from pathlib import Path

from data_squirrel.config.dynamic_data_nut import (
	Nut,
	Value,
	GenericAttribute,
	AtrClass,
	CustomAttribute
)


from serena.utilities.ensemble_structures import Sara2SecondaryStructure
from serena.utilities.ensemble_structures import Sara2StructureList
from serena.interfaces.Sara2_API_Python3 import DesignPerformanceData
from serena.interfaces.Sara2_API_Python3 import DesignInformation
from serena.interfaces.Sara2_API_Python3 import WetlabData
from serena.interfaces.nupack4_0_28_wsl2_interface import NupackSettings
from serena.interfaces.nupack4_0_28_wsl2_interface import MaterialParameter

class Nut_Attributes(Enum):
	Data = "data_db"


class PNASData(Nut):

	def __init__(self, working_folder:Path, var_name:str, use_db:bool = False) -> None:
		super().__init__(enum_list=Nut_Attributes,
			use_db=True,
			db=None,
			var_name=var_name,
			working_folder=working_folder)


		self.data_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="design_info_db",
			atr_type=['DesignPerformanceData', 'DesignInformation', 'WetlabData']))

		self.data_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="nupack_settings_db",
			atr_type=['NupackSettings', 'MaterialParameter']))

		self.data_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="structs_db",
			atr_type=['Sara2StructureList', 'Sara2SecondaryStructure']))

		self.data_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="fmn_folded_mfe_db",
			atr_type=['Sara2SecondaryStructure']))

		self.data_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="fmn_folded_weighted_db",
			atr_type=['Sara2SecondaryStructure']))

class Data(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def design_info(self)->DesignPerformanceData:
		self.parent.nut_filter.yaml_operations.yaml.register_class(DesignPerformanceData)
		self.parent.nut_filter.yaml_operations.yaml.register_class(DesignInformation)
		self.parent.nut_filter.yaml_operations.yaml.register_class(WetlabData)
		return self.parent.design_info_db

	@design_info.setter
	def design_info(self, value:DesignPerformanceData):
		if isinstance(value, DesignPerformanceData) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(DesignPerformanceData)
		self.parent.nut_filter.yaml_operations.yaml.register_class(DesignInformation)
		self.parent.nut_filter.yaml_operations.yaml.register_class(WetlabData)
		self.parent.design_info_db = value


	@property
	def nupack_settings(self)->NupackSettings:
		self.parent.nut_filter.yaml_operations.yaml.register_class(NupackSettings)
		self.parent.nut_filter.yaml_operations.yaml.register_class(MaterialParameter)
		return self.parent.nupack_settings_db

	@nupack_settings.setter
	def nupack_settings(self, value:NupackSettings):
		if isinstance(value, NupackSettings) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(NupackSettings)
		self.parent.nut_filter.yaml_operations.yaml.register_class(MaterialParameter)
		self.parent.nupack_settings_db = value


	@property
	def structs(self)->Sara2StructureList:
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2StructureList)
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2SecondaryStructure)
		return self.parent.structs_db

	@structs.setter
	def structs(self, value:Sara2StructureList):
		if isinstance(value, Sara2StructureList) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2StructureList)
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2SecondaryStructure)
		self.parent.structs_db = value


	@property
	def fmn_folded_mfe(self)->Sara2SecondaryStructure:
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2SecondaryStructure)
		return self.parent.fmn_folded_mfe_db

	@fmn_folded_mfe.setter
	def fmn_folded_mfe(self, value:Sara2SecondaryStructure):
		if isinstance(value, Sara2SecondaryStructure) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2SecondaryStructure)
		self.parent.fmn_folded_mfe_db = value


	@property
	def fmn_folded_weighted(self)->Sara2SecondaryStructure:
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2SecondaryStructure)
		return self.parent.fmn_folded_weighted_db

	@fmn_folded_weighted.setter
	def fmn_folded_weighted(self, value:Sara2SecondaryStructure):
		if isinstance(value, Sara2SecondaryStructure) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2SecondaryStructure)
		self.parent.fmn_folded_weighted_db = value


class ArchiveSecondaryStructureList(PNASData):

	def __init__(self, working_folder:str, var_name:str, use_db:bool = False) -> None:
		super().__init__(use_db=use_db,
			var_name=var_name,
			working_folder=Path(working_folder))


		self._data: Data = Data(save_value=True,
			current=None,
			parent=self.data_db)

	@property
	def data(self)->Data:
		return self._data

	@data.setter
	def data(self, struct:Data):
		if isinstance(struct, Data) == False:
			raise ValueError("Invalid value assignment")
		self._data = struct


