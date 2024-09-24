# KinFragLib_PocketEnum
This pipeline automatically places and expands fragments in a binding pocket of a given kinase structure – based on a user-defined subpocket path – to generate numerous possible kinase inhibitors. The pipeline allows the user to customize this process by e.g. setting the starting subpocket and the subpocket path for fragment growing and defining thresholds

## Usage

### Requirements
BioSolveIT's [SeeSAR](https://www.biosolveit.de/SeeSAR) (3D desktop modeling platform to prepare the protein input files) and the following [SeeSAR command-line tools](https://www.biosolveit.de/download/) need to be installed and attached with a valid license:
* **FlexX** - for docking
* **HYDE** - for scoring and optimization 
The SeeSAR command-line tools should be be preferable placed within the root directory of the project.

### Installation & Dependencies
Create a conda environment containing all required packages:
```bash
conda env create -f environment.yml
# When using a MacBook with an M1 chip you may need:
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
```
Activate the new environment:
```bash
conda activate king_frag_lib_pocket_enum
```

Download and install `KinFragLib` package:
```bash
git clone https://github.com/volkamerlab/KinFragLib.git
# we need to change the branch since Custom-KinFragLib is not on the main branch yet:
cd KinFragLib
git checkout custom-kinfraglib
cd ..
# now we can continue installing the package
pip install -e KinFragLib
```

### Input
* `.flexx` and `.hydescorer` files need to be prepared and dowloaded for evey subpocket using the [SeeSAR GUI](https://www.biosolveit.de/SeeSAR). All these files need to be placed in one directoy that is named **only** by the pdb ID  of the structure (e.g., `5n1f`) and should be named according to the following scheme: `<subpocket>.flexx` (`<subpocket>.hydescorer`), `<subpocket>` needs to be replaced with the subpocket of this file. E.g. the FlexX file for the AP subpocket needs to be stored as `AP.flexx`. By default, the program will search for this directory within the `config` directory, however, this can be changed wihtin the `settings.json` file. For an example see the [`config/5n1f`](config/5n1f) directory.
* Create a JSON configuration file, such that the structure pdb id, the core subpocket, path of subpockets for fragment growing, the path to the fragment library, the path to the FlexX and HYDE executeable, and the path to the folder conatining folder (named with the PDB ID) with the `.flexx` and `.hydescorer` files is defined. A template file is given in [`config/templates/settings.json`](config/templates/settings.json), where all **required** arguments are labeled with a `TODO` (`TODO` needs to be replaced by the argument) and all *optional* arguments are set by their deafult value. [`config/5n1f/settings.json`](config/5n1f/settings.json) provides an example configuration file.

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
All output files are located in a folder named with the PDB ID within the specified result directory (`<path_to_results_folder> / <PDB ID>`, e.g. `results/5n1f`). 

