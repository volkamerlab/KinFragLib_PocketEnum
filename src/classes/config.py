import json
import os
from pathlib import Path

from classes.filters import Filter


class Config:
    """
    Samples all configurations for a subpocket-based docking run

    Attributes
    ----------
    pdb_code: str
        PDB code of strcuture of interest
    core_subpocket: str
        Name of the subpocket to start the subpocket-based docking procedure
    subpockets: list(str)
        List of supockets, dictating the order of fragment growing
    path_kinfraglib: Path
        Path to the fragment library
    path_structure_config: Path
        Path to folder containing the .flexx and .hydescorer configuration file of the structure of interest
    fragments_per_iteration: int
        Amount of fragments to choose per docking iteration (according to docking score)
    poses_per_fragment: int
        Amount of conformers to choose per docked fragment  (according to docking score and diversity)
    filters: list
        List of all custom kinfraglib filters, that should be applied on the fragment library
    cluster_based: bool
        If true, a cluster based approach is applied to select docking poses
    cluster_threshold: float = 1.5
        Distance threshold that should be used for pose  (if cluster_based == true)
    self.use_hyde: bool
        If true, hyde scoring is performed after each docking run
    self.path_hyde: Path
        Path to hyde (only required if use_hyde == true)
    self.hyde_displacement_cutoff: float = 2.5
        RMSD cutoff dictating when a optimized/displaced pose should be droped (if hyde should be applied)
    self.path_flexx: Path
        Path to FlexX
    self.num_thread: int
        Number of threads that should be used
    self.path_temp: Path
        Path to folder of temp-files
    path_results: Path
        Path to folder where reults should be placed
    seed: Int
        Seed used for random algorithm
    """

    def __init__(self) -> None:
        self.pdb_code: str = None
        self.core_subpocket: str = None
        self.subpockets: list = None
        self.path_kinfraglib = None
        self.path_structure_config = None
        self.fragments_per_iteration: int = None
        self.poses_per_fragment: int = None
        self.filters: list = None
        self.cluster_based_pose_selection: bool = None
        self.cluster_threshold: float = None
        self.cluster_based_fragment_selection: bool = None
        self.P: int = None
        self.hyde_displacement_cutoff: float = None
        self.use_hyde: bool = None
        self.path_hyde = None
        self.path_flexx = None
        self.num_threads: int = None
        self.path_temp = None
        self.path_results = None
        self._result = None
        self.seed = None

    def parse(self, config_file) -> None:
        """
        Parses the configuration file

        ----------
        config_file: Path
            Path to JSON-config file
        path_results: Path
            Path to folder where reults should be placed
        """

        # load program definitions
        with open(config_file, "rt") as json_file:
            definitions = json.load(json_file)

        self.pdb_code: str = definitions["pdbCode"]
        self.core_subpocket: str = definitions["CoreSubpocket"]
        self.subpockets: list = definitions["Subpockets"]
        self.fragments_per_iteration: int = definitions["NumberFragmentsPerIterations"]
        self.poses_per_fragment: int = definitions["NumberPosesPerFragment"]
        self.cluster_based_pose_selection: bool = definitions[
            "UseClusterBasedPoseFiltering"
        ]

        self.cluster_based_fragment_selection: bool = definitions[
            "UseClusterBasedFragmentFiltering"
        ]

        if self.cluster_based_pose_selection:
            self.cluster_threshold = definitions.get("HydeDisplacementCutoff") or 2.5

        if self.cluster_based_fragment_selection:
            self.P = definitions.get("PSoftMin") or 1

        self.use_hyde = definitions["UseHyde"]

        if self.use_hyde:
            self.hyde_displacement_cutoff = (
                definitions.get("HydeDisplacementCutoff") or 2.5
            )

        self.num_threads = definitions.get("NumberThreads")

        self.filters = [
            Filter(name, values) for name, values in definitions["Filters"].items()
        ]

        self.seed = definitions.get("Seed") or 42

        # Paths
        self.path_kinfraglib = Path(definitions["KinFragLib"])
        self.path_flexx = Path(definitions["FlexX"])
        self.path_structure_config = Path(definitions["Config"]) / self.pdb_code
        self.path_hyde = Path(definitions["Hyde"]) if self.use_hyde else None

    def initialize_folders(self, path_results) -> None:
        """
        Checks if results folder exists, creates the temp folder and the pdb-code specific folders

        Parameters
        ----------
        path_results: Path
            Path to folder where reults should be placed
        """
        HERE = Path().resolve()

        # create temp folder if it does not exists

        if not os.path.exists(HERE / "temp"):
            os.mkdir(HERE / "temp")
        self.path_temp = HERE / "temp" / self.pdb_code

        if not os.path.exists(self.path_temp):
            os.mkdir(self.path_temp)

        if not os.path.exists(HERE / path_results):
            raise OSError(f"Result folder {str(HERE / path_results)} does not exists")

        self.path_results = HERE / path_results / self.pdb_code

        if not os.path.exists(self.path_results):
            os.mkdir(self.path_results)
        else:
            # clear file, if it already exists
            with open(self.path_results / "results.sdf", "wt") as file:
                pass
