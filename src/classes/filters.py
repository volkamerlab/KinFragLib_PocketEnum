from kinfraglib import filters
from pathlib import Path


class Filter:
    """
    Represents all custom kinfraglib filters

    Attributes
    ----------
    name: str
        Name of the filter
    params: list
        List of filterspecific parameters
    """

    def __init__(self, name, params: dict):
        self.name = name.lower()
        self.params = params

    def get_param(self, param_name):
        """
        Get a filter-parameter by it's name

        Returns
        ----------
        Paramter "param_name" if exists else None

        Parameters
        ----------
        param_name: str
            Name of the parameter
        """
        if param_name in self.params.keys():
            return self.params[param_name]
        return None

    def apply_filter(self, fragment_library):
        """
        Applies the filter to the given fragment_library

        Returns
        ----------
        Filtered library

        Parameters
        ----------
        fragment_library: Dict
            Library containing all fragments where the index should match to the fragment ids
        """
        if self.name == "pains":
            fragment_library, _ = filters.unwanted_substructures.get_pains(fragment_library)
        elif self.name == "brenk":
            fragment_library, _ = filters.unwanted_substructures.get_brenk(
                fragment_library, Path(self.get_param("path_data"))
            )
        elif self.name == "ro3":
            fragment_library = filters.drug_likeness.get_ro3_frags(fragment_library)
        elif self.name == "qed":
            fragment_library = filters.drug_likeness.get_qed(
                fragment_library, cutoff_val=self.get_param("cutoff_val")
            )
        elif self.name == "bb":
            fragment_library = filters.synthesizability.check_building_blocks(
                fragment_library,
                str(str(self.get_param("path_data")) + "/Enamine_Building_Blocks.sdf"),
            )
        elif self.name == "syba":
            fragment_library = filters.synthesizability.calc_syba(
                fragment_library, cutoff=self.get_param("cutoff_val")
            )

        # remove all filtered fragments from the fragment library
        for sp in fragment_library.keys():
            t = fragment_library[sp]["bool_" + self.name]
            fragment_library[sp] = fragment_library[sp][
                fragment_library[sp]["bool_" + self.name] == 1
            ]

        return fragment_library
