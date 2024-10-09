# imports
from rdkit.Chem import PandasTools
from rdkit import Chem


# define constants
PATH_COMPOUNDS = "results_5n1f_25_02/5n1f/results_post_filtered.sdf"
PATH_CHEMBL = "evaluation/chembl_33.sdf"
PATH_OUTPUT = "results_5n1f_25_02/5n1f/results_chembl.sdf"

if __name__ == "__main__":

    # load compounds to dataframe
    compounds = PandasTools.LoadSDF(PATH_COMPOUNDS, embedProps=True)

    print(f"compounds loaded")

    # load chembl
    chembl_compounds = PandasTools.LoadSDF(PATH_CHEMBL, embedProps=True)

    print(f"chembl loaded")

    # add murcko scaffolds (smiles are canonical)
    PandasTools.AddMurckoToFrame(compounds)
    print(f"added murckos to compounds")

    PandasTools.AddMurckoToFrame(chembl_compounds)
    print(f"added murckos to chembl")

    # find matchings
    compounds["matches"] = compounds["Murcko_SMILES"].apply(
        lambda murcko_comp: [
            chembl_comp["chembl_id"]
            for _, chembl_comp in chembl_compounds.iterrows()
            if chembl_comp["Murcko_SMILES"] == murcko_comp
        ]
    )

    PandasTools.WriteSDF(compounds, PATH_OUTPUT, properties=list(compounds.columns))
