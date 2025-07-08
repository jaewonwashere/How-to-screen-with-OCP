# A Beginner's Guide to Using `mlrelax_uma.py` on the Expanse Cluster

Welcome! This document explains the `mlrelax_uma.py` Python script, which is used to study how molecules stick to metal surfaces using a computer model called the Universal Molecular Assistant (UMA). We’ll focus on testing one metal, Rhodium (Rh), with a molecule called `ipxo`. This guide is written for beginners who have never used Linux, Conda, or installed Python packages. We’ll walk you through what the script does, how it works, and how to run it on the Expanse supercomputer.

## What Does the Script Do?

The `mlrelax_uma.py` script helps scientists understand how molecules (called adsorbates) bind to metal surfaces (called slabs). This is important in fields like chemistry and materials science, for example, in designing catalysts that speed up chemical reactions.

Here’s what the script does in simple terms:
1. **Loads a Metal Structure**: It starts with Rhodium (Rh), identified by `mp-74`, a code from a materials database.
2. **Adds a Molecule**: It places the `ipxo` molecule (defined in a file called `ipxoh.xyz`) on the metal surface.
3. **Simulates Relaxation**: It uses a machine learning model (UMA) to “relax” the system, meaning it adjusts the positions of atoms to find the most stable arrangement.
4. **Calculates Energies**: It computes the energy of the molecule sticking to the surface (adsorption energy).
5. **Checks for Problems**: It looks for issues, like the molecule breaking apart or not sticking properly.
6. **Saves Results**: It saves the final atom positions in files called trajectory files (`.traj`) for later analysis.

The script runs on a supercomputer called Expanse, which has powerful GPUs (graphics processing units) to handle these calculations quickly.

## How Does the Script Work?

Let’s break down the main parts of the script (`mlrelax_uma.py`) and what each part does. Don’t worry about the code details yet—we’ll focus on the big picture.

### 1. **Imports and Setup**
The script starts by importing tools (libraries) it needs:
- **Python Libraries**: Tools like `os` (for handling files), `argparse` (for reading user inputs), and `time` (to track how long it runs).
- **ASE (Atomic Simulation Environment)**: A library for working with atoms and molecules, used to load and manipulate structures.
- **FairChem**: A library from Meta AI that provides the UMA model and tools for adsorption calculations.
- **PyTorch**: A machine learning library that helps UMA run calculations, especially on GPUs.

### 2. **Reading User Inputs (`parse_arguments`)**
The script expects you to provide four pieces of information when you run it:
- `--adsorbates`: The molecule to study (e.g., `ipxo`).
- `--username`: Your Expanse username (e.g., `jlee36`).
- `--slab_type`: The type of metal surface (`monometallic` for single metals like Rh, or `alloy` for mixed metals).
- `--mp_id`: The ID of the metal (e.g., `mp-74` for Rhodium).

These inputs are processed by the `parse_arguments` function, which makes sure you’ve provided valid information.

### 3. **Creating the Molecule (`create_adsorbate_kwargs`)**
The `create_adsorbate_kwargs` function prepares the `ipxo` molecule:
- It loads the molecule’s structure from a file (`/home/jlee36/ocp/utils/structures/ipxoh.xyz`), which describes the positions of atoms in `ipxo`.
- It creates an `Adsorbate` object, a special format used by FairChem, specifying which atom in `ipxo` (index 2) will bind to the metal surface.
- For other molecules (like `ch3` or `acetone`), it uses predefined formats or ASE’s `molecule` function.

### 4. **Loading the Metal Surface**
The script uses the `Bulk` and `Slab` classes from FairChem to:
- Load Rhodium (`mp-74`) as a bulk material (a 3D crystal structure).
- Create a surface (slab) from this bulk, focusing on specific crystal planes (Miller indices up to 1, e.g., (100) surface).

### 5. **Initializing the UMA Model**
The script loads the UMA model (`uma-s-1p1`), a machine learning model trained to predict how atoms interact:
- It uses `pretrained_mlip.get_predict_unit` to load the model.
- It sets up a `FAIRChemCalculator` to use UMA for energy and force calculations.
- It checks if a GPU is available (using `torch.cuda.is_available()`) to run faster; otherwise, it uses the CPU.

### 6. **Running the Adsorption Simulation (`run_adsorbml_wrapper`)**
The `run_adsorbml_wrapper` function is the core of the calculation:
- It uses FairChem’s `adsorb_ml_pipeline` to:
  1. Relax the metal surface (adjust atom positions to minimize energy).
  2. Place the `ipxo` molecule on the surface in one position (`num_placements=1` for testing).
  3. Relax the combined system (metal + molecule) using the LBFGS optimizer, which adjusts atoms until forces are small (`fmax=0.02`) or 100 steps are reached.
  4. Calculate the adsorption energy by comparing the energy of the combined system to the separate metal and molecule.
  5. Check for problems (anomalies) like the molecule breaking apart or not sticking.
- It saves the relaxed atom positions as `.traj` files in a directory like `/home/jlee36/ocp/data/monometallic/ipxo/Rh_mp-74_100_ipxo/`.

### 7. **Main Loop (`main`)**
The `main` function ties everything together:
- It reads your inputs (e.g., `ipxo`, `jlee36`, `monometallic`, `mp-74`).
- It loads the Rhodium surface and `ipxo` molecule.
- It runs the simulation for each surface plane of Rhodium.
- It prints results, including the number of valid configurations and any issues (anomalies).
- It tracks and prints the total time taken.

## Setting Up and Running the Script on Expanse

To run the script, you’ll use the Expanse supercomputer’s login node for testing. The login node is like a control center where you can prepare and test the script, but heavy calculations should later use the GPU partition. Follow these steps carefully.

### Step 1: Log into Expanse
1. Open a terminal on your computer (e.g., Terminal on Mac, Command Prompt or PowerShell on Windows).
2. Type the following command, replacing `jlee36` with your Expanse username:
   ```bash
   ssh jlee36@login.expanse.sdsc.edu
   ```
3. Enter your password when prompted (you won’t see it as you type).
4. You’re now on the Expanse login node, where you’ll run commands.

### Step 2: Set Up the Environment
Expanse uses a system called Conda to manage Python and its libraries. We’ll use a pre-installed environment to save time.

1. Load/install the Conda module. You can use miniconda.

2. Create a conda environment. This will be your workspace, where you will have all your dependencies installed and not clutter up your computer.
   ```bash
   conda create -n <your_environment_name>
   ```
   This creates your environment. I named mine 'uma'.


3. Activate the pre-installed PyTorch environment:
   ```bash
   conda activate uma
   ```
   You’ll see `(uma)` before your prompt, indicating the environment is active.

4. Install required dependencies:
  ```bash
  pip install fairchem-core fairchem-data-oc fairchem-applications-cattsunami
  ```

  Make sure you have a Hugging Face account, have already applied for model access to the UMA model repository, and have logged in to Hugging Face using an access token. You can use the following to save an auth token,
  ```bash
  huggingface-cli login
  ```

5. Check if required libraries are installed:
   ```bash
   pip list | grep -E "fairchem|torch|ase"
   ```
   You should see something like:
   ```
   fairchem-core   <version>
   torch           <version>
   ase             <version>
   ```
   If any are missing, ask your mentor to install them, as this requires special permissions. For reference, the command would be:
   ```bash
   pip install fairchem-core torch ase
   ```

### Step 3: Verify the Molecule File
The script needs the `ipxoh.xyz` file to define the `ipxo` molecule. Check if it exists:
```bash
ls /home/jlee36/ocp/utils/structures/ipxoh.xyz
```
If you see the file name, it’s there. If not, ask your mentor for the correct file.

### Step 4: Save the Script
1. Create or edit the script file:
   ```bash
   nano mlrelax_uma.py
   ```
   Nano is a simple text editor. It will open a blank file or the existing script.

2. Copy and paste the entire script from your mentor’s provided code (the code block above, starting with `#!/usr/bin/env python3` and ending with `main()`).

3. Save and exit:
   - Press `Ctrl + O`, then `Enter` to save.
   - Press `Ctrl + X` to exit.

### Step 5: Run the Script
Run the script with the specific inputs for Rhodium and `ipxo`:
```bash
python3 mlrelax_uma.py --adsorbates ipxo --username jlee36 --slab_type monometallic --mp_id mp-74
```
- This tells the script to:
  - Use `ipxo` as the molecule.
  - Use your username (`jlee36`).
  - Use a single metal (`monometallic`).
  - Use Rhodium (`mp-74`).

- You should see output like:
  ```
  Parsed arguments: {'adsorbates': ['ipxo'], 'username': 'jlee36', 'slab_type': 'monometallic', 'mp_id': ['mp-74']}
  Initializing UMA model...
  UMA model initialized
  Adsorbate kwargs: [{'adsorbate': <Adsorbate object>}]
  Processing bulk mp_id: mp-74
  Generated <number> slabs
  Output directory: /home/jlee36/ocp/data/monometallic/ipxo/Rh_mp-74_100_ipxo
  Running adsorb_ml_pipeline with adsorbate_kwargs: {'adsorbate_atoms': <Atoms object>, 'adsorbate_binding_indices': [2]}
  Completed: /home/jlee36/ocp/data/monometallic/ipxo/Rh_mp-74_100_ipxo, Top candidates: <number>
  Elapsed time: <time> seconds
  ```

### Step 6: Check Results
- Look for trajectory files:
  ```bash
  ls /home/jlee36/ocp/data/monometallic/ipxo/Rh_mp-74_100_ipxo/*.traj
  ```
  These files store the final atom positions.
- If you see messages about “Anomalies detected,” it means the molecule didn’t stick properly (e.g., it broke apart). Ask your mentor to review.

### Step 7: Troubleshooting
If you see errors:
- **“FileNotFoundError”**: The `ipxoh.xyz` file is missing or in the wrong place. Double-check the path or recreate it.
- **“ModuleNotFoundError”**: A library (e.g., `fairchem`) is missing. Ask your mentor to install it.
- **Script Stops Early**: The login node may be too slow. Tell your mentor it needs to run on the GPU partition.

## What’s Next?
This script is set up for testing with one molecule placement (`num_placements=1`) to avoid overloading the login node. For a full run with more placements (e.g., 100), your mentor will help you create a SLURM script to use Expanse’s GPU partition, which is faster and designed for heavy calculations.

## Key Terms to Know
- **Adsorbate**: The molecule (e.g., `ipxo`) that sticks to the metal surface.
- **Slab**: The flat metal surface (e.g., Rhodium) where the molecule attaches.
- **UMA**: A machine learning model that predicts how atoms move and interact.
- **LBFGS**: A method to adjust atom positions to find the most stable structure.
- **Trajectory File (.traj)**: A file that saves the final positions of atoms.
- **Expanse**: A supercomputer at UC San Diego for running big calculations.
- **Login Node**: The part of Expanse where you set up and test scripts.
- **GPU Partition**: The powerful part of Expanse for running heavy calculations.

## Questions?
If anything is unclear or the script doesn’t work, ask your mentor (Jaewon) for help. Share any error messages you see in the terminal—they’ll help figure out what went wrong. Once this test works, you can try other metals or molecules or move to the GPU partition for faster runs.

Happy learning, and welcome to computational chemistry!
