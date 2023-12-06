import subprocess
import os
from rdkit import Chem
from pathlib import Path

from classes.ligand import Ligand
from classes.config import Config
from classes.recombination import Recombination

def hyde_scoring(path_docking_results: Path, path_config: Path, path_output: Path, path_hyde: Path, ligand: Ligand, print_output=False) -> list:
    """
    runs hydescoring

    Returns
    ----------
    List of all poses (docking result) as Mol (Properties: BIOSOLVEIT.DOCKING_SCORE, pose)

    Parameters
    ----------
    path_fragment: pathlib.path
        Path to fragment sdf-file
    path_config: pathlib.path
        Path to hyde-config file.
    path_output: pathlib.path
        Path to output file
    path_hyde: pathlib.path
        Path to Hyde
    ligand: Ligand
        Ligand object of molecule to score 
    """
    output_text = subprocess.run(
        [
            str(Path('.') / path_hyde),
            "-i",
            str(path_docking_results),
            "--binding-site-definition",
            str(path_config),
            "-o",
            str(path_output),
        ],
        capture_output=True
    )
    if print_output:
        print(output_text.stderr)

    # read results from sdf
    opt_fragments = []
    if os.stat(str(path_output)).st_size and os.stat(str(path_docking_results)).st_size:  # only acces reult file if at least one pose was generated
        for molecule_docking, molecule_opt in zip(Chem.SDMolSupplier(str(path_docking_results)), Chem.SDMolSupplier(str(path_output))):
            # clear all hyde properties, that are not needed
            superfluos_props = ['BIOSOLVEIT.HYDE_ATOM_SCORES [kJ/mol]', 'BIOSOLVEIT.HYDE_LIGAND_EFFICIENCY range: ++, +, 0, -, --', 'BIOSOLVEIT.HYDE_LIGAND_LIPOPHILIC_EFFICIENCY range: ++, +, 0, -, --', 'BIOSOLVEIT.INTER_CLASH range: red, yellow, green',
                                'BIOSOLVEIT.INTRA_CLASH range: red, yellow, green', 'BIOSOLVEIT.INTRA_CLASH range: red, yellow, green', 'BIOSOLVEIT.MOLECULE_CHECKSUM', 'BIOSOLVEIT.TORSION_QUALITY range: red, yellow, green, not rotatable']
            for prop in superfluos_props:
                molecule_opt.ClearProp(prop)
            # set proporties
            molecule_opt.SetProp('fragment_ids', str(ligand.fragment_ids))
            molecule_opt.SetProp('smiles_ligand', Chem.MolToSmiles(molecule_opt))
            molecule_opt.SetProp('smiles_fragments_dummy', str(ligand.smiles_dummy))
            molecule_opt.SetProp('smiles_fragments', str(ligand.smiles))
            # copy docking from docked molecule to optimzed molecule
            molecule_opt.SetProp('BIOSOLVEIT.DOCKING_SCORE', molecule_docking.GetProp('BIOSOLVEIT.DOCKING_SCORE'))
            opt_fragments.append(molecule_opt)
    return opt_fragments

def prepare_core_fragments(fragment_library: dict, config: Config) -> list:
    """
    Initializes all core fragments from the fragment library

    Returns
    ----------
    List of all core fragments as Ligand

    Parameters
    ----------
    fragment_library: dict
        KinFragLib fragment library
    config: Config
        Config object, storing program configurations
    """

    core_fragments = []

    # prepare all core fragments
    for i in fragment_library[config.core_subpocket].index:
        smiles = fragment_library[config.core_subpocket]['smiles'][i]
        smiles_dummy = fragment_library[config.core_subpocket]['smiles_dummy'][i]
        fragment_recombination = Recombination([config.core_subpocket + "_" + str(i)], [], {config.core_subpocket: smiles}, {config.core_subpocket: smiles_dummy})
        core_fragments.append(Ligand(fragment_library[config.core_subpocket]['ROMol'][i], {config.core_subpocket: i}, 
                                                   fragment_recombination, 
                                                   {config.core_subpocket: smiles_dummy}, {config.core_subpocket: smiles}))
    
    return core_fragments
