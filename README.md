# OCP_screening
Codebase for screening catalysts using the Open Catalyst Project. Largely based off of the Open Catalyst Project tutorials available on the fair-chem repository.
Includes scripts for determining relevant bulk structures, cleaving surfaces, setting adsorption sites, configuring adslabs, relaxation, and post-processing.

## Setup
The first thing we need to do is get you set up on the UCSB cluster. Here is the UCSB CNSI CSC's guideline on HPC (High Perfomance Computing): [https://csc.cnsi.ucsb.edu/sites/default/files/2023-01/HPC_Workshop_Winter_23.pdf]

1. Request an account on the UCSB Center for Scientific Computing: [https://csc.cnsi.ucsb.edu/]
2. Referring to the document above, set up the client that suits your local operating system. You will access the campus cluster (Pod) via _ssh_, and each operating system (Windows, MacOS, Linux) has different ways of using _ssh_.
   * I use _wsl_ on my Windows system and the default terminal on my Mac.
4. You'll get an email that has your login information (with a temporary password). _ssh_ in to your account by typing in:
   ```bash
   ssh your-login@pod-login1.cnsi.ucsb.edu
   ```
   This will log you into Pod's login node, which is where you set up your environments, scripts, and send your jobs to other nodes in the cluster.
   Then, follow the instructions on the email to change your temporary password to the one that you want!
6. 
