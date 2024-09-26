# KinFragLib_PocketEnum
This pipeline automatically places and expands fragments in a binding pocket of a given kinase structure – based on a user-defined subpocket path – to generate numerous possible kinase inhibitors. The pipeline allows the user to customize this process by e.g. setting the starting subpocket and the subpocket path for fragment growing and defining thresholds

## Usage

### Requirements
The following BioSolveIT command-line tools need to be installed and attached with a valid license:
* FlexX for docking
* HYDE for scoring and optimization 
Preferable they should be placed within the root directory of the project.

### Installation & dependencies
Create a conda environment containing all required packages:
```bash
conda env create -f environment.yml
```
When using a MacBook with an M1 chip, you may need
```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
```
Activate the new environment:
```bash
# activate the environment
conda activate kinfraglib_pocket_enum
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

Install the `KinFragLib_PocketEnum` package itself (required to run the notebooks):
```bash
cd ..
# install package
pip install -e KinFragLib_PocketEnum
```

### Input
* `.flexx` and `.hydescorer` files need to be prepared and dowloaded for evey subpocket using the SeeSAR GUI. All these files need to be placed in one folder named with the pdb id of the sturture and should be named according to the following scheme: `<subpocket>.flexx` (`<subpocket>.hydescorer`), `<subpocket>` needs to be replaced with the subpocket of this file. E.g. the FlexX file for the AP subpocket needs to be stored as `AP.flexx`. For an example see the folder config_smarts_exclude. 
* Change to provided JSON configuration file (`definition.json`) or create your own, such that the structure pdb id, the core subpocket, path of subpockets for fragment growing, the path to the fragment library, the path to the FlexX and HYDE executeable, and the path to the folder conatining folder (named with the PDB ID) with the `.flexx` and `.hydescorer` files is defined.

### Run subpocket based docking programm
```bash
python3 src/fragment_docking.py -d <JSON_definition_file> -o <JSON_outputfile> -r <path_to_results_folder>
```
for help run:
```bash
python3 src/fragment_docking.py -h
```

Note: The process can be tracked on https://wandb.ai/home:
    To enable this, just run the subpocket based docking porgram (as instructed) and follow the instruction (either select (1) to create an account or (2) to login if an account exists already)
If this is not wanted choose the option (3)
