import pytest

from serena.ensemble_variation import RunEnsembleVariation
from serena.utilities.ensemble_structures import Sara2StructureList,Sara2SecondaryStructure

def get_ev_from_structures_list(secondary_structures_list_2_item:Sara2StructureList,  secondary_structure_5:Sara2SecondaryStructure):
    ensemble_variation:RunEnsembleVariation = RunEnsembleVariation()
    ev:float = ensemble_variation.ev_from_structures_list(structures_list=secondary_structures_list_2_item,
                                                            mfe_structure=secondary_structure_5)
    assert ev == 3.0