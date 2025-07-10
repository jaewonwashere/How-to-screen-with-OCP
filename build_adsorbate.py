from ase.data.pubchem import pubchem_atoms_search, pubchem_atoms_conformer_search
from ase.visualize import view
from ase.io import write

# This code reads in custom adsorbates (lactic/pyruvic acid) here and writes them as .xyz files to be read into the mlrelax_uma.py execution.

# Fetch the main structure and conformers for lactic acid
lacid = pubchem_atoms_search(name='lactic%20acid')
lacid_c = pubchem_atoms_conformer_search(name='lactic%20acid')

# Export lactic acid conformers as .xyz files
for i, conf in enumerate(lacid_c):
    filename = f"lactic_acid_conf_{i+1}.xyz"
    write(filename, conf)
    print(f"Saved {filename}")

# Fetch the main structure and conformers for pyruvic acid
pacid = pubchem_atoms_search(name='pyruvic%20acid')
pacid_c = pubchem_atoms_conformer_search(name='pyruvic%20acid')

# Export pyruvic acid conformers as .xyz files
for i, conf in enumerate(pacid_c):
    filename = f"pyruvic_acid_conf_{i+1}.xyz"
    write(filename, conf)
    print(f"Saved {filename}")

# Optionally visualize pyruvic acid conformers
for i, conf in enumerate(pacid_c):
    print(f"Showing conformer {i + 1}")
    view(conf, viewer='ase')
