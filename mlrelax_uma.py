#!/usr/bin/env python3

# Last modified: 2025/07/07
# Updated to use UMA model (uma-s-1p1) with AdsorbML pipeline via import.
# Supports adsorbate 'ipxch' (*O intermediate) and integrates anomaly detection.
# Uses LBFGS optimizer. Fixed Adsorbate handling for XYZ-based inputs.
# Takes four argparse flags: --adsorbates, --username, --slab_type, --mp_id.

import os
import argparse
import time
from ase.io import read
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from fairchem.data.oc.core import Adsorbate, Bulk, Slab
from fairchem.core.components.calculate.recipes.adsorbml import run_adsorbml
from ase.build import molecule
import torch

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate adsorption energies using UMA and AdsorbML pipeline.")
    parser.add_argument('--adsorbates', type=str, nargs='+', required=True, help="Specify the adsorbates (e.g., ch3, oh, h, ipxo, acetone, acetonec)")
    parser.add_argument('--username', type=str, required=True, help="Specify the username for output directory")
    parser.add_argument('--slab_type', type=str, choices=['monometallic', 'alloy'], required=True, help="Specify the type of slab: 'monometallic' or 'alloy'.")
    parser.add_argument('--mp_id', type=str, nargs='*', help="Optional list of mp_id's to process. If not provided, will use stable mp_ids.")
    return parser.parse_args()

def create_adsorbate_kwargs(adsorbates):
    adsorbate_kwargs_list = []
    for adsorbate in adsorbates:
        adsorbate = adsorbate.lower()
        kwargs = {}

        if adsorbate in ['ch3', 'oh', 'h', 'oh2']:
            kwargs['adsorbate_smiles_from_db'] = '*' + adsorbate.upper()

        elif adsorbate in ['ipx', 'ipxo', 'ipxc']:
            file_path = f'/home/{args.username}/ocp/utils/structures/ipxoh.xyz'
            adsorbate_atoms = read(file_path)
            binding_index = 2 if adsorbate == 'ipxo' else 1
            kwargs['adsorbate'] = Adsorbate(adsorbate_atoms=adsorbate_atoms, adsorbate_binding_indices=[binding_index])

        elif adsorbate in ['ipa']:
            file_path = f'/home/{args.username}/ocp/utils/structures/ipa.xyz'
            adsorbate_atoms = read(file_path)
            kwargs['adsorbate'] = Adsorbate(adsorbate_atoms=adsorbate_atoms, adsorbate_binding_indices=[0])
            
        elif adsorbate in ['ipxch']:
            file_path = f'/home/{args.username}/ocp/utils/structures/ipxch.xyz'
            adsorbate_atoms = read(file_path)
            kwargs['adsorbate'] = Adsorbate(adsorbate_atoms=adsorbate_atoms, adsorbate_binding_indices=[3])  # O binding, C hydrogenation

        elif adsorbate in ['acetone', 'acetoneo', 'acetonec']:
            adsorbate_atoms = molecule("CH3COCH3")
            binding_index = 0 if adsorbate == 'acetoneo' else 1
            kwargs['adsorbate'] = Adsorbate(adsorbate_atoms=adsorbate_atoms, adsorbate_binding_indices=[binding_index])

        else:
            raise ValueError(f"Unsupported adsorbate: {adsorbate}")

        adsorbate_kwargs_list.append(kwargs)

    return adsorbate_kwargs_list

def run_adsorbml_wrapper(slab, adsorbate_obj, calculator, output_dir, num_placements=3, fmax=0.02, steps=100):
    # Wrapper to adapt run_adsorbml for Adsorbate objects
    adsorbate_kwargs = {
        'adsorbate_smiles_from_db': adsorbate_obj.smiles if adsorbate_obj.smiles else None,
        'adsorbate_atoms': adsorbate_obj.atoms if not adsorbate_obj.smiles else None,
        'adsorbate_binding_indices': adsorbate_obj.binding_indices
    }
    results = run_adsorbml(
        slab=slab,
        adsorbate=adsorbate_kwargs['adsorbate_smiles_from_db'] or adsorbate_obj,
        calculator=calculator,
        optimizer_cls=LBFGS,
        num_placements=num_placements,
        fmax=fmax,
        steps=steps,
        reference_ml_energies=True,
    )
    return results

def main():
    global args # Make args accessible in create_adsorbate_kwargs
    tinit = time.time()
    args = parse_arguments()

    adsorbates = args.adsorbates
    username = args.username
    slab_type = args.slab_type
    mp_id_list = args.mp_id

    # Define stable mp_ids based on slab_type
    stable_monometallic_mp_ids = ['mp-124', 'mp-81', 'mp-94', 'mp-102', 'mp-90', 'mp-30', 'mp-13', 'mp-101', 
                                  'mp-129', 'mp-75', 'mp-23', 'mp-49', 'mp-20483', 'mp-2', 'mp-126', 
                                  'mp-8', 'mp-74', 'mp-33', 'mp-67', 'mp-117', 'mp-84', 'mp-91', 'mp-112', 
                                  'mp-79', 'mp-131']
    
    stable_alloy_mp_ids = ['mp-2369', 'mp-20536', 'mp-1851', 'mp-317', 'mp-45']

    # If mp_id_list is provided, override the stable mp_ids
    mp_ids = mp_id_list if mp_id_list else (stable_monometallic_mp_ids if slab_type == 'monometallic' else stable_alloy_mp_ids)

    # Initialize UMA model
    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda" if torch.cuda.is_available() else "cpu")
    calculator = FAIRChemCalculator(predictor, task_name="oc20")

    adsorbate_kwargs_list = create_adsorbate_kwargs(adsorbates)

    for bulk_mp_id in mp_ids:
        for adsorbate, kwargs in zip(adsorbates, adsorbate_kwargs_list):
            bulk = Bulk(bulk_src_id_from_db=bulk_mp_id)
            slabs = bulk.get_slabs(max_miller=1)

            bulk_element_str = bulk.atoms.get_chemical_formula()
            for slab in slabs:
                millers_str = str(slab.millers).replace(', ', '').replace(',', '')
                output_dir = f"/home/{username}/ocp/data/{slab_type}/{adsorbate.lower()}/{bulk_element_str}_{bulk_mp_id}_{millers_str}_{adsorbate.lower()}"
                os.makedirs(output_dir, exist_ok=True)

                try:
                    results = run_adsorbml_wrapper(
                        slab=slab,
                        adsorbate_obj=kwargs['adsorbate'],
                        calculator=calculator,
                        output_dir=output_dir,
                        num_placements=3,
                        fmax=0.02,
                        steps=100,
                    )
                    # Save trajectories for top candidates
                    for idx, adslab_result in enumerate(results['adslabs']):
                        atoms = adslab_result['atoms']
                        atoms.write(f"{output_dir}/{idx}.traj")
                    print(f"Completed: {output_dir}, Top candidates: {len(results['adslabs'])}")
                    if results['adslab_anomalies']:
                        print(f"Anomalies detected in {output_dir}: {results['adslab_anomalies']}")
                except Exception as e:
                    print(f"Failed: {output_dir}, Error: {str(e)}")

    print(f'Elapsed time: {time.time() - tinit:1.1f} seconds')

if __name__ == "__main__":
    main()