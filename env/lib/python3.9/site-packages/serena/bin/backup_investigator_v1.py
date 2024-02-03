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
from serena.analysis.investigator import ComparisonEvalResults
from serena.analysis.investigator import RatioResults
from serena.utilities.comparison_structures import ComparisonNucResults
from serena.utilities.comparison_structures import ComparisonNucCounts
from serena.utilities.local_minima_variation import ComparisonLMVResponse
from serena.utilities.local_minima_variation import ComparisonLMV
from serena.analysis.investigator import LMVAssertionResult
from serena.utilities.ensemble_variation import EV
from serena.analysis.scoring import BasicScoreResults
from serena.analysis.scoring import AdvancedScoreResults
from serena.interfaces.Sara2_API_Python3 import DesignInformation
from serena.interfaces.Sara2_API_Python3 import WetlabData
from serena.utilities.weighted_structures import WeightedEnsembleResult
from serena.analysis.ensemble_analysis import ReferenceStructures

class Nut_Attributes(Enum):
	Investigator = "investigator_db"
	Scores = "scores_db"
	DesignParameters = "design_info_db"


class PNASData(Nut):

	def __init__(self, working_folder:Path, var_name:str, use_db:bool = False) -> None:
		super().__init__(enum_list=Nut_Attributes,
			use_db=True,
			db=None,
			var_name=var_name,
			working_folder=working_folder)


		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.PARENT,
			attribute="comparison_eval_result_db",
			atr_type=None))

		self.investigator_db.comparison_eval_result_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="ratios_db",
			atr_type=['RatioResults', 'CLASS']))

		self.investigator_db.comparison_eval_result_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="BRaise_list_db",
			atr_type=float))

		self.investigator_db.comparison_eval_result_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="BUratio_list_db",
			atr_type=float))

		self.investigator_db.comparison_eval_result_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="bound_total_list_db",
			atr_type=float))

		self.investigator_db.comparison_eval_result_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="unbound_total_list_db",
			atr_type=float))

		self.investigator_db.comparison_eval_result_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="nuc_penatly_count_db",
			atr_type=int))

		self.investigator_db.comparison_eval_result_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="first_BUratio_db",
			atr_type=float))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.PARENT,
			attribute="lmv_values_db",
			atr_type=None))

		self.investigator_db.lmv_values_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="lmv_comps_db",
			atr_type=['ComparisonLMV', 'EV', 'CLASS']))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.PARENT,
			attribute="lmv_assertions_db",
			atr_type=None))

		self.investigator_db.lmv_assertions_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="comp_compare_to_mfe_db",
			atr_type=str))

		self.investigator_db.lmv_assertions_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="unbouund_pronounced_db",
			atr_type=bool))

		self.investigator_db.lmv_assertions_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="bound_pronounced_db",
			atr_type=bool))

		self.investigator_db.lmv_assertions_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="is_on_off_switch_db",
			atr_type=bool))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="comp_nuc_counts_db",
			atr_type=['ComparisonNucResults', 'ComparisonNucCounts']))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="num_groups_db",
			atr_type=int))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="total_structures_ensemble_db",
			atr_type=int))

		self.investigator_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="lmv_references_db",
			atr_type=['ReferenceStructures', 'WeightedEnsembleResult', 'Sara2SecondaryStructure']))

		self.scores_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="basic_scores_db",
			atr_type=['BasicScoreResults']))

		self.scores_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="advanced_scores_db",
			atr_type=['AdvancedScoreResults']))

		self.scores_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="number_structures_db",
			atr_type=int))

		self.design_info_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="design_info_db",
			atr_type=['DesignInformation']))

		self.design_info_db.new_attr(GenericAttribute(atr_class=AtrClass.CHILD,
			attribute="wetlab_data_db",
			atr_type=['WetlabData']))

class ComparisonEvalResults(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def ratios(self)->List[RatioResults]:
		self.parent.nut_filter.yaml_operations.yaml.register_class(RatioResults)
		return self.parent.ratios_db

	@ratios.setter
	def ratios(self, value:List[RatioResults]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, RatioResults) == False:
				raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(RatioResults)
		self.parent.ratios_db = value


	@property
	def BRaise_list(self)->List[float]:
		return self.parent.BRaise_list_db

	@BRaise_list.setter
	def BRaise_list(self, value:List[float]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, float) == False:
				raise ValueError("Invalid value assignment")
		self.parent.BRaise_list_db = value


	@property
	def BUratio_list(self)->List[float]:
		return self.parent.BUratio_list_db

	@BUratio_list.setter
	def BUratio_list(self, value:List[float]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, float) == False:
				raise ValueError("Invalid value assignment")
		self.parent.BUratio_list_db = value


	@property
	def bound_total_list(self)->List[float]:
		return self.parent.bound_total_list_db

	@bound_total_list.setter
	def bound_total_list(self, value:List[float]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, float) == False:
				raise ValueError("Invalid value assignment")
		self.parent.bound_total_list_db = value


	@property
	def unbound_total_list(self)->List[float]:
		return self.parent.unbound_total_list_db

	@unbound_total_list.setter
	def unbound_total_list(self, value:List[float]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, float) == False:
				raise ValueError("Invalid value assignment")
		self.parent.unbound_total_list_db = value


	@property
	def nuc_penatly_count(self)->int:
		return self.parent.nuc_penatly_count_db

	@nuc_penatly_count.setter
	def nuc_penatly_count(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nuc_penatly_count_db = value


	@property
	def first_BUratio(self)->float:
		return self.parent.first_BUratio_db

	@first_BUratio.setter
	def first_BUratio(self, value:float):
		if isinstance(value, float) == False:
			raise ValueError("Invalid value assignment")
		self.parent.first_BUratio_db = value


class ComparisonLMVResponse(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def lmv_comps(self)->List[ComparisonLMV]:
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonLMV)
		self.parent.nut_filter.yaml_operations.yaml.register_class(EV)
		return self.parent.lmv_comps_db

	@lmv_comps.setter
	def lmv_comps(self, value:List[ComparisonLMV]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, ComparisonLMV) == False:
				raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonLMV)
		self.parent.nut_filter.yaml_operations.yaml.register_class(EV)
		self.parent.lmv_comps_db = value


class LMVAssertionResult(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def comp_compare_to_mfe(self)->List[str]:
		return self.parent.comp_compare_to_mfe_db

	@comp_compare_to_mfe.setter
	def comp_compare_to_mfe(self, value:List[str]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, str) == False:
				raise ValueError("Invalid value assignment")
		self.parent.comp_compare_to_mfe_db = value


	@property
	def unbouund_pronounced(self)->List[bool]:
		return self.parent.unbouund_pronounced_db

	@unbouund_pronounced.setter
	def unbouund_pronounced(self, value:List[bool]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, bool) == False:
				raise ValueError("Invalid value assignment")
		self.parent.unbouund_pronounced_db = value


	@property
	def bound_pronounced(self)->List[bool]:
		return self.parent.bound_pronounced_db

	@bound_pronounced.setter
	def bound_pronounced(self, value:List[bool]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, bool) == False:
				raise ValueError("Invalid value assignment")
		self.parent.bound_pronounced_db = value


	@property
	def is_on_off_switch(self)->List[bool]:
		return self.parent.is_on_off_switch_db

	@is_on_off_switch.setter
	def is_on_off_switch(self, value:List[bool]):
		if isinstance(value, list) == False:
			raise ValueError("Invalid value assignment")
		if len(value) < 1:
			raise Exception("Empty lists not allowed")

		for item in value:
			if isinstance(item, bool) == False:
				raise ValueError("Invalid value assignment")
		self.parent.is_on_off_switch_db = value


class DesignParameters(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def design_info(self)->DesignInformation:
		self.parent.nut_filter.yaml_operations.yaml.register_class(DesignInformation)
		return self.parent.design_info_db

	@design_info.setter
	def design_info(self, value:DesignInformation):
		if isinstance(value, DesignInformation) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(DesignInformation)
		self.parent.design_info_db = value


	@property
	def wetlab_data(self)->WetlabData:
		self.parent.nut_filter.yaml_operations.yaml.register_class(WetlabData)
		return self.parent.wetlab_data_db

	@wetlab_data.setter
	def wetlab_data(self, value:WetlabData):
		if isinstance(value, WetlabData) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(WetlabData)
		self.parent.wetlab_data_db = value


class Investigator(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value
		self._comparison_eval_result: ComparisonEvalResults = ComparisonEvalResults(save_value=True,
			current=None,
			parent=self.parent.comparison_eval_result_db)

		self._lmv_values: ComparisonLMVResponse = ComparisonLMVResponse(save_value=True,
			current=None,
			parent=self.parent.lmv_values_db)

		self._lmv_assertions: LMVAssertionResult = LMVAssertionResult(save_value=True,
			current=None,
			parent=self.parent.lmv_assertions_db)


	@property
	def comparison_eval_result(self)->ComparisonEvalResults:
		return self._comparison_eval_result

	@comparison_eval_result.setter
	def comparison_eval_result(self, value:ComparisonEvalResults):
		if isinstance(value, ComparisonEvalResults) == False:
			raise ValueError("Invalid value assignment")
		self._comparison_eval_result = value


	@property
	def lmv_values(self)->ComparisonLMVResponse:
		return self._lmv_values

	@lmv_values.setter
	def lmv_values(self, value:ComparisonLMVResponse):
		if isinstance(value, ComparisonLMVResponse) == False:
			raise ValueError("Invalid value assignment")
		self._lmv_values = value


	@property
	def lmv_assertions(self)->LMVAssertionResult:
		return self._lmv_assertions

	@lmv_assertions.setter
	def lmv_assertions(self, value:LMVAssertionResult):
		if isinstance(value, LMVAssertionResult) == False:
			raise ValueError("Invalid value assignment")
		self._lmv_assertions = value


	@property
	def comp_nuc_counts(self)->ComparisonNucResults:
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonNucResults)
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonNucCounts)
		return self.parent.comp_nuc_counts_db

	@comp_nuc_counts.setter
	def comp_nuc_counts(self, value:ComparisonNucResults):
		if isinstance(value, ComparisonNucResults) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonNucResults)
		self.parent.nut_filter.yaml_operations.yaml.register_class(ComparisonNucCounts)
		self.parent.comp_nuc_counts_db = value


	@property
	def num_groups(self)->int:
		return self.parent.num_groups_db

	@num_groups.setter
	def num_groups(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.num_groups_db = value


	@property
	def total_structures_ensemble(self)->int:
		return self.parent.total_structures_ensemble_db

	@total_structures_ensemble.setter
	def total_structures_ensemble(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.total_structures_ensemble_db = value


	@property
	def lmv_references(self)->ReferenceStructures:
		self.parent.nut_filter.yaml_operations.yaml.register_class(ReferenceStructures)
		self.parent.nut_filter.yaml_operations.yaml.register_class(WeightedEnsembleResult)
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2SecondaryStructure)
		return self.parent.lmv_references_db

	@lmv_references.setter
	def lmv_references(self, value:ReferenceStructures):
		if isinstance(value, ReferenceStructures) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(ReferenceStructures)
		self.parent.nut_filter.yaml_operations.yaml.register_class(WeightedEnsembleResult)
		self.parent.nut_filter.yaml_operations.yaml.register_class(Sara2SecondaryStructure)
		self.parent.lmv_references_db = value


class Scores(CustomAttribute):
	def __init__(self, parent: Any, current:Any, save_value:bool) -> None:
		self.parent = parent
		self.current = current
		self.do_save = save_value

	@property
	def basic_scores(self)->BasicScoreResults:
		self.parent.nut_filter.yaml_operations.yaml.register_class(BasicScoreResults)
		return self.parent.basic_scores_db

	@basic_scores.setter
	def basic_scores(self, value:BasicScoreResults):
		if isinstance(value, BasicScoreResults) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(BasicScoreResults)
		self.parent.basic_scores_db = value


	@property
	def advanced_scores(self)->AdvancedScoreResults:
		self.parent.nut_filter.yaml_operations.yaml.register_class(AdvancedScoreResults)
		return self.parent.advanced_scores_db

	@advanced_scores.setter
	def advanced_scores(self, value:AdvancedScoreResults):
		if isinstance(value, AdvancedScoreResults) == False:
			raise ValueError("Invalid value assignment")
		self.parent.nut_filter.yaml_operations.yaml.register_class(AdvancedScoreResults)
		self.parent.advanced_scores_db = value


	@property
	def number_structures(self)->int:
		return self.parent.number_structures_db

	@number_structures.setter
	def number_structures(self, value:int):
		if isinstance(value, int) == False:
			raise ValueError("Invalid value assignment")
		self.parent.number_structures_db = value


class ArchiveInvestigator(PNASData):

	def __init__(self, working_folder:str, var_name:str, use_db:bool = False) -> None:
		super().__init__(use_db=use_db,
			var_name=var_name,
			working_folder=Path(working_folder))


		self._investigator: Investigator = Investigator(save_value=True,
			current=None,
			parent=self.investigator_db)

		self._scores: Scores = Scores(save_value=True,
			current=None,
			parent=self.scores_db)

		self._design_info: DesignParameters = DesignParameters(save_value=True,
			current=None,
			parent=self.design_info_db)

	@property
	def investigator(self)->Investigator:
		return self._investigator

	@investigator.setter
	def investigator(self, struct:Investigator):
		if isinstance(struct, Investigator) == False:
			raise ValueError("Invalid value assignment")
		self._investigator = struct


	@property
	def scores(self)->Scores:
		return self._scores

	@scores.setter
	def scores(self, struct:Scores):
		if isinstance(struct, Scores) == False:
			raise ValueError("Invalid value assignment")
		self._scores = struct


	@property
	def design_info(self)->DesignParameters:
		return self._design_info

	@design_info.setter
	def design_info(self, struct:DesignParameters):
		if isinstance(struct, DesignParameters) == False:
			raise ValueError("Invalid value assignment")
		self._design_info = struct


