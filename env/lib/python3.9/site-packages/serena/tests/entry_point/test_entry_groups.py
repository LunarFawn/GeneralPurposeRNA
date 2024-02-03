import pytest
from typing import List

#from serena.generate_structures import SecondaryStructures, EnsembleGroups
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs, MakeEnsembleGroups



def test_make_single_ensemble_group(secondary_structures_list_2_item:Sara2StructureList, ensemble_state_mfe_structs:EnsembleSwitchStateMFEStructs):
    ensemble_group:MakeEnsembleGroups = MakeEnsembleGroups()
    group:SingleEnsembleGroup = ensemble_group.make_singel_ensemble_group(ensemble_structures=secondary_structures_list_2_item,
                                                                            mfe_switch_structures=ensemble_state_mfe_structs,                                                                           
                                                                            kcal_start=-50,
                                                                            kcal_end=-30)
    assert group.group == secondary_structures_list_2_item
    assert group.kcal_start == -50
    assert group.kcal_end == -30
    assert group.kcal_span == 20

def test_make_multiple_ensemble_groups(single_ensemble_group:SingleEnsembleGroup, single_ensemble_group_2:SingleEnsembleGroup, ensemble_state_mfe_structs:EnsembleSwitchStateMFEStructs):
    ensemble_groups_list:List[SingleEnsembleGroup] = []
    ensemble_groups_list.append(single_ensemble_group)
    ensemble_groups_list.append(single_ensemble_group_2)

    ensemble_groups:MakeEnsembleGroups = MakeEnsembleGroups()
    groups:MultipleEnsembleGroups = ensemble_groups.make_multiple_ensemple_groups(ensemble_groups=ensemble_groups_list,
                                                                                  mfe_switch_structures=ensemble_state_mfe_structs)
    assert groups.groups == [single_ensemble_group, single_ensemble_group_2]
    assert groups.num_groups == 2
    assert groups.non_switch_state_structure == ensemble_state_mfe_structs.non_switch_mfe_struct
