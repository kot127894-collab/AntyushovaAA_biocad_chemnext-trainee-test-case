from pathlib import Path
from rdkit import Chem
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
import parmed

INPUT_SDF = "example.sdf"

OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

ff = ForceField("openff_unconstrained-2.1.0.offxml")


def clean(name, i):
    if not name:
        return f"ligand_{i}"
    return "".join(c if c.isalnum() else "_" for c in name)


suppl = Chem.SDMolSupplier(INPUT_SDF, removeHs=False)

for i, mol in enumerate(suppl):

    try:

        if mol is None:
            print(f"Skipping molecule {i}: RDKit failed")
            continue

        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        name = clean(name, i)

        print(f"\nProcessing: {name}")

        out = OUTPUT_DIR / name
        out.mkdir(exist_ok=True)

        # RDKit -> OpenFF
        off_mol = Molecule.from_rdkit(
            mol,
            allow_undefined_stereo=True
        )

        # charges
        off_mol.assign_partial_charges("gasteiger")

        # PDB
        off_mol.to_file(
            str(out / "ligand.pdb"),
            file_format="PDB"
        )

        topology = off_mol.to_topology()

        system = ff.create_openmm_system(
            topology,
            charge_from_molecules=[off_mol]
        )

        structure = parmed.openmm.load_topology(
            topology.to_openmm(),
            system,
            xyz=off_mol.conformers[0].to_openmm()
        )

        # topology
        structure.save(
            str(out / "topol.top"),
            overwrite=True
        )

        # gro
        structure.save(
            str(out / "conf.gro"),
            overwrite=True
        )

        print(f"Saved: {name}")

    except Exception as e:

        print(f"ERROR for molecule {i}: {e}")
        continue

print("\nDONE")