# OCP_screening
Codebase for screening catalysts using the [Open Catalyst Project](https://opencatalystproject.org/index.html). Largely based off of the Open Catalyst Project tutorials available on Meta AI Research's [fair-chem repository](https://github.com/FAIR-Chem/fairchem).
Includes scripts for determining relevant bulk structures, cleaving surfaces, setting adsorption sites, configuring adslabs, relaxation, and post-processing.

## Setup
The first thing we need to do is get you set up on the UCSB cluster.
Here is the UCSB CNSI CSC's [overview](https://csc.cnsi.ucsb.edu/sites/default/files/2023-01/HPC_Workshop_Winter_23.pdf) on HPC (High Perfomance Computing):

1. [**Request an account**](https://csc.cnsi.ucsb.edu/) on the UCSB Center for Scientific Computing. It will take a couple of days to get it set up.
   
2. Referring to the document above, **set up the client that suits your local operating system**. You will access the campus cluster (Pod) via _ssh_, and each operating system (Windows, MacOS, Linux) has different ways of using _ssh_.
   * I use [https://learn.microsoft.com/en-us/windows/wsl/install](_wsl_) on my Windows system (you can also use Putty) and the default terminal on my Mac.
     
3. You'll get an email that has your login information (with a temporary password). _ssh_ in to your account by typing in:
   ```bash
   ssh your-login@pod-login1.cnsi.ucsb.edu
   ```
   It will ask your password. Just type in the temporary password in your email, and press Enter.
   This will **log you into Pod's login node**, which is where you set up your environments, scripts, and send your jobs to other nodes in the cluster.
   Then, follow the instructions on the email to **change your temporary password** to the one that you want!
   Review the information on CSC's guidlines to set up your file transfer system, and VPN if you plan on logging in to the cluster off-campus.
  
4. If you are logged in properly, you should see something like:
   ```bash
   ----------------------------
   
   Welcome to Pod
   For basic documentation to get started please see
   http://csc.cnsi.ucsb.edu/docs/pod-cluster
   ```
   
5. Now, let's get ready to install stuff! We will use **[_anaconda_](https://docs.conda.io/projects/conda/en/stable/)** to manage our packages and environments.
   Specifically, I like to use [_Miniconda_](https://docs.anaconda.com/miniconda/), which is a minimal version of Anaconda.
   Follow the installation guides from the _Miniconda_ link above. You can use the **"Quick Command Line Install"** part at the bottom of the page.
   It will take a while to install. Sit back and relax!
   
6. If everything went okay, you should be able to see a 'miniconda3' folder in your home directory. This is a good time to go over **basic linux commands**.
    |list files/folders| change directory | pring working directory |
    |------------------|------------------|-------------------------|
    |`ls`              | `cd {directory_where_you_want_to_move_to}` | `pwd`|
    
    (The brackets {} are a notation for placeholders which you can replace with your relevant string. Don't include the brackets!)
    (FYI: the 'home directory' for you account can be accessed by `cd ~`. This is the same directory as `/home/{your_user_name}/`)

    CSC's [HPC overview](https://csc.cnsi.ucsb.edu/sites/default/files/2023-01/HPC_Workshop_Winter_23.pdf), as well as resources such as this [cheat sheet](https://www.stationx.net/linux-command-line-cheat-sheet/) will guide you through more essential Linux CLI(command line interface) commands.

8. Now, let's install some packages. An easy way to install (almost) everything you need to run ocp calculations is by following the guidelines [here](https://fair-chem.github.io/core/install.html). However, the default environment uses cuda 11.8, which is not an available module on Pod. Therefore, we will set it up so that it uses cuda 12.1. I have set up my local environment as such, and the environment file is available on Pod as: `fairchem_Pod.yml`.
   Copy the `.yml` file into your home directory:
   `cp /home/jaewon_lee/fairchem_Pod.yml ./'
   
9. The `.yml` file acts like the seed of a tree. We can use it to reproduce the environment that I have on my system:
   `conda env create -f fairchem_Pod.yml`
   (If you run into problems, reach out to `jaewon_lee@ucsb.edu`)
   
10. Now you should have the proper `fair-chem` environment on your login node!  
   

   '
   
