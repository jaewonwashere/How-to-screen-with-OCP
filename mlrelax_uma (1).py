#!/usr/bin/env python3

# Last modified: 2025/07/08
# Updated to use UMA model (uma-s-1p1) with AdsorbML pipeline via import.
# Supports adsorbate 'ipxch' (*O intermediate) and integrates anomaly detection.
# Uses LBFGS optimizer. Fixed Adsorbate handling, global args syntax, UnboundLocalError, idx error, anomaly message, CSV anomaly consistency, and list index error.
# Added enhanced debug prints and ensured adslab_anomalies aligns with adslabs.
# Generates CSV summary in ~/ocp/data/summary.csv with adsorbate, slab type, element, mp-id, Miller index, steps, energy, force, and anomalies.
# Appends to existing summary.csv instead of overwriting.
# Takes four argparse flags: --adsorbates, --username, --slab_type, --mp_id.

import os
import argparse
import time
import numpy as np
import pandas as pd
from ase.io import read
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from fairchem.data.oc.core import Adsorbate, Bulk, Slab
from fairchem.core.components.calculate.recipes.adsorbml import adsorb_ml_pipeline
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
            kwargs['adsorbate'] = Adsorbate(adsorbate_smiles_from_db=kwargs['adsorbate_smiles_from_db'])

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
    print("Adsorbate kwargs:", adsorbate_kwargs_list)
    return adsorbate_kwargs_list

def run_ml_relax_job(initial_atoms, calc, optimizer_cls, fmax, steps, output_dir, idx=0):
    atoms = initial_atoms.copy()
    atoms.calc = calc
    dyn = optimizer_cls(atoms, trajectory=f"{output_dir}/{idx}.traj")
    dyn.run(fmax=fmax, steps=steps)
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    max_force = np.max(np.abs(forces))
    nsteps = dyn.nsteps
    atoms.calc = None
    return {
        "input_atoms": {"atoms": initial_atoms},
        "atoms": atoms,
        "results": {"energy": energy, "forces": forces},
        "nsteps": nsteps,
        "max_force": max_force
    }

def run_adsorbml_wrapper(slab, adsorbate_obj, calculator, output_dir, num_placements=1, fmax=0.02, steps=10):
    # Wrapper to adapt run_adsorbml for Adsorbate objects
    adsorbate_kwargs = {
        'adsorbate_atoms': adsorbate_obj.atoms,
        'adsorbate_binding_indices': adsorbate_obj.binding_indices
    } if not adsorbate_obj.smiles else {
        'adsorbate_smiles_from_db': adsorbate_obj.smiles
    }
    print(f"Running adsorb_ml_pipeline with adsorbate_kwargs: {adsorbate_kwargs}")
    
    # Replicate run_adsorbml logic to handle adsorbate_kwargs directly
    atomic_reference_energies = {
        "H": -3.477, "O": -7.204, "C": -7.282, "N": -8.083
    }
    ml_relax_job = lambda atoms: run_ml_relax_job(
        initial_atoms=atoms,
        calc=calculator,
        optimizer_cls=LBFGS,
        fmax=fmax,
        steps=steps,
        output_dir=output_dir,
        idx=0  # Default idx=0, as num_placements=1
    )
    
    try:
        results = adsorb_ml_pipeline(
            slab=slab,
            adsorbates_kwargs=[adsorbate_kwargs],
            multiple_adsorbate_slab_config_kwargs={"num_configurations": num_placements},
            ml_slab_adslab_relax_job=ml_relax_job,
            reference_ml_energies=True,
            atomic_reference_energies=atomic_reference_energies
        )
        # Ensure adslab_anomalies has one element per adslab
        adslabs = results.get('adslabs', [])
        adslab_anomalies = results.get('adslab_anomalies', [[] for _ in range(len(adslabs))])
        print(f"Debug: results keys: {list(results.keys())}")
        print(f"Debug: adslabs length: {len(adslabs)}")
        print(f"Debug: adslab_anomalies length: {len(adslab_anomalies)}, content: {adslab_anomalies}")
        results['adslab_anomalies'] = adslab_anomalies
        return results
    except Exception as e:
        print(f"Debug: adsorb_ml_pipeline failed with error: {str(e)}")
        return {'adslabs': [], 'adslab_anomalies': []}

def write_summary_csv(summary_data, username):
    csv_path = f"/home/{username}/ocp/data/summary.csv"
    columns = [
        "Adsorbate_Type", "Slab_Type", "Element", "MP_ID", "Miller_Index",
        "Steps_Taken", "Final_Energy", "Final_Force", "Anomaly"
    ]
    new_df = pd.DataFrame(summary_data, columns=columns)
    
    # Check if CSV exists and append, otherwise create new
    if os.path.exists(csv_path):
        existing_df = pd.read_csv(csv_path)
        updated_df = pd.concat([existing_df, new_df], ignore_index=True)
    else:
        updated_df = new_df
    
    updated_df.to_csv(csv_path, index=False)
    print(f"Summary appended to {csv_path}")

def main():
    global args
    tinit = time.time()
    args = parse_arguments()
    print("Parsed arguments:", vars(args))

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
    print("Initializing UMA model...")
    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda" if torch.cuda.is_available() else "cpu")
    calculator = FAIRChemCalculator(predictor, task_name="oc20")
    print("UMA model initialized")

    # Initialize summary data list
    summary_data = []

    adsorbate_kwargs_list = create_adsorbate_kwargs(adsorbates)

    for bulk_mp_id in mp_ids:
        print(f"Processing bulk mp_id: {bulk_mp_id}")
        bulk = Bulk(bulk_src_id_from_db=bulk_mp_id)
        slabs = bulk.get_slabs(max_miller=1)
        print(f"Generated {len(slabs)} slabs")

        bulk_element_str = bulk.atoms.get_chemical_formula()
        for slab in slabs:
            millers_str = str(slab.millers).replace(', ', '').replace(',', '')
            # Relax bare slab
            output_dir = f"/home/{username}/ocp/data/{slab_type}/slab/{bulk_element_str}_{bulk_mp_id}_{millers_str}_slab"
            os.makedirs(output_dir, exist_ok=True)
            print(f"Output directory: {output_dir}")
            try:
                slab_atoms = slab.atoms.copy()
                slab_atoms.calc = calculator
                dyn = LBFGS(slab_atoms, trajectory=f"{output_dir}/slab.traj")
                dyn.run(fmax=0.02, steps=10)
                energy = slab_atoms.get_potential_energy()
                forces = slab_atoms.get_forces()
                max_force = np.max(np.abs(forces))
                nsteps = dyn.nsteps
                summary_data.append([
                    "slab", slab_type, bulk_element_str, bulk_mp_id, millers_str,
                    nsteps, energy, max_force, "None"
                ])
                print(f"Completed bare slab: {output_dir}")
            except Exception as e:
                print(f"Failed bare slab: {output_dir}, Error: {str(e)}")
                summary_data.append([
                    "slab", slab_type, bulk_element_str, bulk_mp_id, millers_str,
                    0, 0.0, 0.0, f"Error: {str(e)}"
                ])

            for adsorbate, kwargs in zip(adsorbates, adsorbate_kwargs_list):
                output_dir = f"/home/{username}/ocp/data/{slab_type}/{adsorbate.lower()}/{bulk_element_str}_{bulk_mp_id}_{millers_str}_{adsorbate.lower()}"
                os.makedirs(output_dir, exist_ok=True)
                print(f"Output directory: {output_dir}")

                try:
                    results = run_adsorbml_wrapper(
                        slab=slab,
                        adsorbate_obj=kwargs['adsorbate'],
                        calculator=calculator,
                        output_dir=output_dir,
                        num_placements=100,
                        fmax=0.02,
                        steps=100,
                    )
                    print(f"Debug: Processing results for {output_dir}")
                    for idx, adslab_result in enumerate(results.get('adslabs', [])):
                        atoms = adslab_result['atoms']
                        atoms.write(f"{output_dir}/{idx}.traj")
                        nsteps = adslab_result.get('nsteps', 0)
                        energy = adslab_result['results']['energy']
                        max_force = adslab_result.get('max_force', np.max(np.abs(adslab_result['results']['forces'])))
                        anomalies = results.get('adslab_anomalies', [[]])[idx] if idx < len(results.get('adslab_anomalies', [[]])) else []
                        anomaly = "None" if not anomalies or anomalies == [] else str(anomalies)
                        print(f"Debug: idx={idx}, energy={energy}, max_force={max_force}, anomaly={anomaly}")
                        summary_data.append([
                            adsorbate, slab_type, bulk_element_str, bulk_mp_id, millers_str,
                            nsteps, energy, max_force, anomaly
                        ])
                    print(f"Completed: {output_dir}, Top candidates: {len(results.get('adslabs', []))}")
                    if any(results.get('adslab_anomalies', [])):
                        print(f"Anomalies detected in {output_dir}: {results['adslab_anomalies']}")
                except Exception as e:
                    print(f"Failed: {output_dir}, Error: {str(e)}")
                    summary_data.append([
                        adsorbate, slab_type, bulk_element_str, bulk_mp_id, millers_str,
                        0, 0.0, 0.0, f"Error: {str(e)}"
                    ])

    # Append summary to CSV
    write_summary_csv(summary_data, username)
    print(f'Elapsed time: {time.time() - tinit:1.1f} seconds')

if __name__ == "__main__":
    main()