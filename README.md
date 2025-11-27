# KinFragLib_PocketEnum
This pipeline automatically places and expands fragments in a binding pocket of a given kinase structure – based on a user-defined subpocket path – to generate numerous possible kinase inhibitors. It allows the user to customize this process by e.g. setting the starting subpocket and the subpocket path for fragment growing and defining thresholds.

**Note:** This repository is currently a pre-release.
## Usage

### Requirements
BioSolveIT's [SeeSAR](https://www.biosolveit.de/SeeSAR) (3D desktop modeling platform to prepare the protein input files) and the following [SeeSAR command-line tools](https://www.biosolveit.de/download/) need to be installed and attached with a valid license:
* **FlexX** - for docking
* **HYDE** - for scoring and optimization 
The SeeSAR command-line tools should be be preferable placed within the root directory of the project.

**Note:** since we employ the FlexX docking and HYDE scoring tool, you need a valid [SeeSAR license](https://www.biosolveit.de/license/) to succesfully run this pipeline.

### Installation & Dependencies
Create a conda environment containing all required packages:
```bash
conda env create -f environment.yml
```
```bash
# When using a MacBook with an M1 chip you may need:
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
```
Activate the new environment:
```bash
conda activate kinfraglib_pocket_enum
```

Download and install `KinFragLib` package:
```bash
git clone https://github.com/volkamerlab/KinFragLib.git

# now we can continue installing the package
pip install -e KinFragLib
```

Install the `KinFragLib_PocketEnum` package itself (required to run the notebooks):
```bash
cd ..
# install package
pip install -e KinFragLib_PocketEnum
```

### Input
* `.flexx` and `.hydescorer` binding pocket, and pharmacophore (i.e., subpocket) definition files need to be prepared and downloaded for evey subpocket using the [SeeSAR GUI](https://www.biosolveit.de/SeeSAR).
These files All these files need to be placed in one directoy that is named **only** by the pdb ID  of the structure (e.g., `5n1f`) and should be named according to the following scheme: `<subpocket>.flexx` (`<subpocket>.hydescorer`), `<subpocket>` needs to be replaced with the subpocket of this file. E.g. the FlexX file for the AP subpocket needs to be stored as `AP.flexx`. By default, the program will search for this directory within the `config` directory, however, this can be changed wihtin the `settings.json` file. For an example see the [`config/5n1f`](config/5n1f) directory.
* Create a JSON configuration file, such that the structure pdb id, the core subpocket, path of subpockets for fragment growing, the path to the fragment library, the path to the FlexX and HYDE executeable, and the path to the folder conatining folder (named with the PDB ID) with the `.flexx` and `.hydescorer` files is defined. A template file is given in [`config/templates/settings.json`](config/templates/settings.json), where all **required** arguments are labeled with a `TODO` (`TODO` needs to be replaced by the argument) and all *optional* arguments are set by their deafult value. [`config/5n1f/settings.json`](config/5n1f/settings.json) provides an example configuration file.

**Note:** We provide a script ([`pocket_definition/pocket_definition.py`](pocket_definition)) that determines the 85 KLIFS binding pocket and subpocket-specific anchor residues which might be handy if the kinase structure of interest is not listed in [KLIFS](https://klifs.net/).

### Run subpocket based docking programm
```bash
python3 src/fragment_docking.py -s <JSON_settings_file> -r <path_to_results_folder> 
```
Here, `<JSON_settings_file>` should be replaced by the path to the JSON configuration file (e.g., `conf/5n1f/settings.json`), and `<path_to_results_folder>` by a path to a result folder, by default this will be set to `results`.

For help run:
```bash
python3 src/fragment_docking.py -h
```

**Note:** The process can be tracked on https://wandb.ai/home:
    To enable this, just run the subpocket based docking porgram (as instructed) and follow the instruction (either select (1) to create an account or (2) to login if an account exists already)
If this is not wanted choose the option (3)

### Output
All output files are located in a folder named with the PDB ID followed by the submitting date within the specified result directory (`<path_to_results_folder> / <PDB ID>_<submit_date>`, e.g. `results/5n1f_2024_01_31`):
* `program_statistics.json`: contains program statistics, such as the runtime, or quantities of the ligand generation process
* `results.sdf`: comprising the generated ligands (for each ligand only the docking conformation of the highest scoring pose is given)
* `SPx.sdf`: comprising *all* docking conformations of the ligands that have been docking in the `x`-th iteration

## Citation.
Buchthal K, Kramer PL, Hubach D, Bach G, Wagner N, Krieger J, et al. Novel Kinase Ligand Generation Using Subpocket-Based Docking. *ChemRxiv.* **2025**; [doi:10.26434/chemrxiv-2025-f9ctg](https://chemrxiv.org/engage/chemrxiv/article-details/68cffe6af4163037709b289f). **This content is a preprint and has not been peer-reviewed.**
