import json
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
        self.cluster_based: bool = None
        self.cluster_threshold: float = None
        self.hyde_displacement_cutoff: float = None
        self.use_hyde: bool = None
        self.path_hyde = None
        self.path_flexx = None
        self.num_threads: int = None
        self.path_temp = None
        self.path_results = None

    def parse(self, config_file, path_results) -> None:
        """
        Parses the configuration file

        Parameterspath_results:
            Path to folder where reults should be placed
        ----------
        config_file: Path
            Path to JSON-config file
        path_results: Path
            Path to folder where reults should be placed
        """

        # load program definitions
        with open(config_file, 'rt') as json_file:
            definitions = json.load(json_file)

        self.pdb_code: str = definitions['pdbCode']
        self.core_subpocket: str = definitions['CoreSubpocket']
        self.subpockets: list = definitions['Subpockets']
        self.fragments_per_iteration: int = definitions['NumberFragmentsPerIterations']
        self.poses_per_fragment: int = definitions['NumberPosesPerFragment']
        self.cluster_based: bool = definitions["UseClusterBasedPoseFiltering"]
        
        if self.cluster_based:
            self.cluster_threshold = definitions.get('HydeDisplacementCutoff') or 2.5
        
        self.use_hyde = definitions['UseHyde']

        if self.use_hyde:
            self.hyde_displacement_cutoff = definitions.get('HydeDisplacementCutoff') or 2.5

        self.num_threads = definitions.get('NumberThreads')

        self.filters = [Filter(name, values) for name, values in definitions['Filters'].items()]

        # Paths
        HERE = Path().resolve()
        self.path_temp = HERE / 'temp' / definitions['pdbCode']
        self.path_kinfraglib = Path(definitions['KinFragLib'])
        self.path_flexx = Path(definitions['FlexX'])
        self.path_structure_config = Path(definitions['Config']) / definitions['pdbCode']
        self.path_results = HERE / path_results / definitions['pdbCode']
        self.path_hyde = Path(definitions['Hyde']) if self.use_hyde else None
