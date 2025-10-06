
# REWEIGHTING AND ANALYSES REPO FOR PaaA2 CG MODEL SYSTEM

Code for all analyses and Maximum Entropy Reweighting of CG PaaA2 ensemble



## Screenshots

![App Screenshot](https://via.placeholder.com/468x300?text=App+Screenshot+Here)


## Authors


Kaushik Borthakur, Thomas Sisk
## Methods and programs

Details of Maximum Entropy Reweighting method can be found in the manuscript ()

- Chemical shifts from the trajectories were calculated using SPARTA+ (https://spin.niddk.nih.gov/bax/software/SPARTA+/)
- Maximum Entropy Reweighting approach used in this project is from https://doi.org/10.1101/2024.10.04.616700
- SAXS profiles were calculated using Pepsi-SAXS (https://team.inria.fr/nano-d/software/pepsi-saxs/)
- RDCs were calculated using PALES (https://www3.mpibpc.mpg.de/groups/zweckstetter/_links/software_pales.htm) [will require a SINGULARITY BUILD/APPTAINER to run PALES on a newer Linux environment]
- writhe_tools package (pip install writhe-tools) was used to calculate writhe features for the analyses [More information on the method can be found in the preprint https://doi.org/10.1101/2025.04.26.650781
]
## Acknowledgements

The authors acknowledge Korey Reid, Thomas Sisk, Emanuele Scalone for their discussions.


![Logo](https://dev-to-uploads.s3.amazonaws.com/uploads/articles/th5xamgrr6se0x5ro4g6.png)

