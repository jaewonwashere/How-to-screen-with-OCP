# A Beginner's Guide to Using `mlrelax_uma.py` on the UCSB CSC Pod Cluster

Welcome! This document explains the `mlrelax_uma.py` Python script, which is used to study how molecules stick to metal surfaces using a computer model called the Universal Molecular Assistant (UMA). We’ll focus on testing one metal, Rhodium (Rh), with a molecule called `ipxo`. This guide is written for beginners who have never used Linux, Conda, or installed Python packages. We’ll walk you through what the script does, how it works, and how to run it on the UCSB CSC Pod cluster, a supercomputer at UC Santa Barbara.

## What Does the Script Do?

The `mlrelax_uma.py` script helps scientists understand how molecules (called adsorbates) bind to metal surfaces (called slabs). This is important in fields like chemistry and materials science, for example, in designing catalysts that speed up chemical reactions.

Here’s what the script does in simple terms:
1. **Loads a Metal Structure**: It starts with Rhodium (Rh), identified by `mp-74`, a code from a materials database.
2. **Adds a Molecule**: It places the `ipxo` molecule (defined in a file called `ipxoh.xyz`) on the metal surface.
3. **Simulates Relaxation**: It uses a machine learning model (UMA) to “relax” the system, meaning it adjusts the positions of atoms to find the most stable arrangement.
4. **Calculates Energies**: It computes the energy of the molecule sticking to the surface (adsorption energy).
5. **Checks for Problems**: It looks for issues, like the molecule breaking apart or not sticking properly.
6. **Saves Results**: It saves the final atom positions in files called trajectory files (`.traj`) for later analysis.

The script runs on the CSC Pod cluster, which has powerful GPUs (graphics processing units) to handle these calculations quickly.

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
- `--username`: Your CSC Pod username (e.g., `jaewon_lee`).
- `--slab_type`: The type of metal surface (`monometallic` for single metals like Rh, or `alloy` for mixed metals).
- `--mp_id`: The ID of the metal (e.g., `mp-74` for Rhodium).

These inputs are processed by the `parse_arguments` function, which makes sure you’ve provided valid information.

### 3. **Creating the Molecule (`create_adsorbate_kwargs`)**
The `create_adsorbate_kwargs` function prepares the `ipxo` molecule:
- It loads the molecule’s structure from a file (`/home/jaewon_lee/ocp/utils/structures/ipxoh.xyz`), which describes the positions of atoms in `ipxo`.
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
- It reads your inputs (e.g., `ipxo`, `jaewon_lee`, `monometallic`, `mp-74`).
- It loads the Rhodium surface and `ipxo` molecule.
- It runs the simulation for each surface plane of Rhodium.
- It prints results, including the number of valid configurations and any issues (anomalies).
- It tracks and prints the total time taken.

## Setting Up and Running the Script on the CSC Pod Cluster

To run the script, you’ll use the CSC Pod cluster’s login node for testing. The login node is like a control center where you can prepare and test the script, but heavy calculations should later use the GPU partition. Follow these steps carefully.

### Step 1: Log into the CSC Pod Cluster
1. Open a terminal on your computer (e.g., Terminal on Mac, Command Prompt or PowerShell on Windows).
2. Type the following command, replacing `jlee36` with your CSC Pod username:
   ```bash
   ssh jaewon_lee@pod-login.csc.ucsb.edu
   ```
3. Enter your password when prompted (you won’t see it as you type). If you’re using two-factor authentication, follow the prompts (e.g., enter a code from your phone).
4. You’re now on the CSC Pod login node, where you’ll run commands.

### Step 2: Set Up the Environment
The CSC Pod cluster uses a system called Conda to manage Python and its libraries. Since there’s no pre-installed PyTorch environment like on Expanse, we’ll create a new Conda environment to install the required tools.

1. **Load the Conda Module**:
   The CSC Pod cluster typically has Conda or Miniconda available. Load it with:
   ```bash
   module load miniconda3
   ```
   If this command fails, ask your mentor to confirm Miniconda is installed or use:
   ```bash
   source ~/miniconda3/bin/activate
   ```

2. **Create a Conda Environment**:
   Create a new environment named `uma` for your workspace:
   ```bash
   conda create -n uma python=3.10
   ```
   This sets up a clean environment with Python 3.10, which is compatible with FairChem.

3. **Activate the Environment**:
   Activate the `uma` environment:
   ```bash
   conda activate uma
   ```
   You’ll see `(uma)` before your prompt, indicating the environment is active.

4. **Install Required Dependencies**:
   Install the necessary Python packages:
   ```bash
   pip install fairchem-core fairchem-data-oc fairchem-applications-cattsunami torch ase
   ```
   This installs:
   - `fairchem-core`: For the UMA model and adsorption tools.
   - `fairchem-data-oc`: For handling materials data.
   - `fairchem-applications-cattsunami`: For additional FairChem tools.
   - `torch`: For machine learning calculations.
   - `ase`: For working with atoms and molecules.

5. **Set Up Hugging Face Access**:
   The UMA model (`uma-s-1p1`) is hosted on Hugging Face, and you need an account with access to download it.
   - If you don’t have a Hugging Face account, create one at [huggingface.co](https://huggingface.co).
   - Request access to the UMA model repository (ask your mentor for the exact repository name, likely `meta-ai/uma-s-1p1`).
   - Log in to Hugging Face in the terminal:
     ```bash
     huggingface-cli login
     ```
   - Enter your Hugging Face access token (provided by Hugging Face) when prompted.

6. **Verify Installation**:
   Check that the libraries are installed:
   ```bash
   pip list | grep -E "fairchem|torch|ase"
   ```
   You should see:
   ```
   fairchem-core               <version>
   fairchem-data-oc            <version>
   fairchem-applications-cattsunami <version>
   torch                       <version>
   ase                         <version>
   ```
   If any are missing, repeat the `pip install` command or ask your mentor for help.

### Step 3: Verify the Molecule File
The script needs the `ipxoh.xyz` file to define the `ipxo` molecule. Check if it exists:
```bash
ls /home/jaewon_lee/ocp/utils/structures/ipxoh.xyz
```
If the file is listed, it’s ready. If not, ask your mentor.

### Step 4: Save the Script
1. Create or edit the script file:
   ```bash
   nano mlrelax_uma.py
   ```
   Nano is a simple text editor that opens in the terminal.

2. Copy and paste the entire script provided by your mentor (starting with `#!/usr/bin/env python3` and ending with `main()`). This is the code that runs the UMA simulation.

3. Save and exit:
   - Press `Ctrl + O`, then `Enter` to save.
   - Press `Ctrl + X` to exit.

### Step 5: Run the Script
Run the script on the login node with the specific inputs for Rhodium and `ipxo`:
```bash
python3 mlrelax_uma.py --adsorbates ipxo --username jaewon_lee --slab_type monometallic --mp_id mp-74
```
- This tells the script to:
  - Use `ipxo` as the molecule.
  - Use your username (`jaewon_lee`).
  - Use a single metal (`monometallic`).
  - Use Rhodium (`mp-74`).

- You should see output like:
  ```
  Parsed arguments: {'adsorbates': ['ipxo'], 'username': 'jaewon_lee', 'slab_type': 'monometallic', 'mp_id': ['mp-74']}
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
  ls /home/jaewon_lee/ocp/data/monometallic/ipxo/Rh_mp-74_100_ipxo/*.traj
  ```
  These files store the final atom positions from the simulation.
- If you see messages about “Anomalies detected,” it means the molecule didn’t stick properly (e.g., it broke apart). Show these to your mentor.

### Step 7: Troubleshooting
If you see errors:
- **“FileNotFoundError”**: The `ipxoh.xyz` file is missing or in the wrong place. Check the path or ask your mentor for the file.
- **“ModuleNotFoundError”**: A library (e.g., `fairchem-core`) is missing. Repeat the `pip install` command or ask your mentor.
- **“HuggingFace Hub: HTTPError”**: The UMA model couldn’t download. Verify your Hugging Face login and access token.
- **Script Stops Early**: The login node may be too slow for calculations. Tell your mentor you need to run on the GPU partition.

## What’s Next?
This script is set up for testing with one molecule placement (`num_placements=1`) to avoid overloading the login node. For a full run with more placements (e.g., 100), your mentor will help you create a SLURM script to use the CSC Pod’s GPU partition, which is faster and designed for heavy calculations. A SLURM script tells the cluster to reserve a GPU and run your script efficiently.

## Key Terms to Know
- **Adsorbate**: The molecule (e.g., `ipxo`) that sticks to the metal surface.
- **Slab**: The flat metal surface (e.g., Rhodium) where the molecule attaches.
- **UMA**: A machine learning model that predicts how atoms move and interact.
- **LBFGS**: A method to adjust atom positions to find the most stable structure.
- **Trajectory File (.traj)**: A file that saves the final positions of atoms.
- **CSC Pod Cluster**: A supercomputer at UC Santa Barbara for running big calculations.
- **Login Node**: The part of the cluster where you set up and test scripts.
- **GPU Partition**: The powerful part of the cluster for running heavy calculations.
- **Conda**: A tool to manage Python environments and libraries.
- **Hugging Face**: A website hosting the UMA model, requiring an account to access.

## Questions?
If anything is unclear or the script doesn’t work, ask your mentor (Jaewon) for help. Share any error messages you see in the terminal—they’ll help figure out what went wrong. Once this test works, you can try other metals or molecules or move to the GPU partition for faster runs.

Happy learning, and welcome to computational chemistry!
