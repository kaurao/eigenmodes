This repository contains the code for reconstructing brain activity patterns as described in the paper [Pang et al. 2023, Nature](https://www.nature.com/articles/s41586-023-06098-1).

## Requirements
- MATLAB 2021a or later (we use statistical abd parallel toolboxes)
- [BrainSpace](https://brainspace.readthedocs.io/en/latest/index.html)
  - git clone -b v0.1.10 https://github.com/MICA-MNI/BrainSpace.git
- [matlab_GIfTI](https://github.com/nno/matlab_GIfTI)
  - git clone https://github.com/nno/matlab_GIfTI.git
- [Rand_index](https://github.com/cmccomb/rand_index)
  - git clone https://github.com/cmccomb/rand_index.git
- [Connectome Workbench](https://www.humanconnectome.org/software/get-connectome-workbench)
  - [matlab] curpath = getenv('PATH');
  - [matlab] setenv('PATH',sprintf('%s:%s','/path/to/.../workbench/bin_macosx64',curpath));
  - [matlab] clear curpath  
  
## Data
We use the data as provided by Pang et al.
- First clone the GitHub repo: https://github.com/NSBLab/BrainEigenmodes


- copy files from https://osf.io/xczmp/
  - To BrainEigenmodes/data/empirical
    - S255_tfMRI_ALLTASKS_raw_lh.mat
    - S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat

  - To BrainEigenmodes/data/results
    - basis_connectome_density_matched_midthickness-lh_evec_200.mat  
    - basis_connectome_EDR_midthickness-lh_evec_200.mat              
    - basis_connectome_midthickness-lh_evec_200.mat
    - basis_geometric_midthickness-lh_evec_200.txt
        
## Usage
- [terminal] git clone https://github.com/NSBLab/BrainEigenmodes.git (add data as described above)
- [terminal] git clone -b v0.1.10 https://github.com/MICA-MNI/BrainSpace.git
- [terminal] git clone https://github.com/nno/matlab_GIfTI.git
- [terminal] git clone https://github.com/cmccomb/rand_index.git
- [terminal] git clone https://github.com/kaurao/eigenmodes.git
- [terminal] mkdir -p ./eigenmodes/results
- [terminal] cd ./eigenmodes
- [matlab] generate_figure;
