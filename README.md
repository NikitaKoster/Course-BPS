# Course-BPS
Project of Chara, Marit and Nikita 

# In silico design of a novel ligand for p38

Nikita Koster (s2213370)
Marit van Arkel (s2509997)
Chara Spyropoulou (s3594971)

# Project description
*What our application does and why we used these technologies*
In this project, we designed a novel ligand for p38 MAPK in silico. We wanted to design a ligand that eventually had a higher docking score than the co-crystallized ligand that was used as a reference.

In the first part (Introduction to Bioinformatics), we first did a literature study to obtain more knowledge on the target protein and the off-target proteins. We looked up the target in the Protein Knowledgebase (UniprotKB) and used BLAST to compare the similarity between the proteins. Subsequently, we compared the 3D orientations of the target and off target proteins. Then, we started off the coding with visualizing the co-crystallized ligand in the target protein.  Next, the protein and ligand were split and saved seperately. The same was done with the off-target protein and this section was finished by aligning the target protein with the off-target protein. This was done to visualize the differences (so conclusions can be drawn on the potential binding of the novel ligand to the off-target protein).

In the second part (Machine Learning), we started off with looking for similar compounds on Zinc and PubChem. We started off the coding with loading in the ChEMBL dataset (in our case CHEMBL260). We filtered the data on assay type 'B'and removed the columns that were not necessary for this part of the project. Subsequently, we added an activity classifier which classified compounds with a pCHEMBL (pIC50) value of higher than 6.5 as active, and classified the other compounds as inactive. Then, we started off with the molecule encoding. We transferred SMILES to MACCS fingerprints and added these in the columns. Thereafter, we used RF, SVM and ANN algorithms to classify the compounds. The model was devided in a 80/20 training and test set and ROC curves were plotted for all three algorithms. Subsequently, we added Morgan fingerprints and we did a cross-validation using the Morgan fingerprints with N_FOLDS = 3 to test the ability of the model to predict new data. Next, we trained a classification model which estimated an MAE of 0.56, which is an acceptable error range. Thereafter, the 15 newly designed molecules (based on scaffolds of molecules that got to stage 2 or 3 in clinical trials) were entered in the RF model to create prediction values. Lig_01-Lig_05, Lig_10 and Lig_15 were based on the same scaffold (the scaffold of reference ligand Lig_00) and showed the best prediction values. Therefore, these ligands will be used in the third part of this application.

In the last part (Molecular Docking), we started with the visualization of the ligand outside of the binding pocket. Then, we prepared the protein and ligands for docking by converting it to the PDBQT file format (which stores information on the protein and ligand). Then, we calculated the box size by calculating the radius of gyration and the center of geometry to maximize the accuracy of finding the proper binding pose. Thereafter, we performed molecular docking with AutoDock Vina and created a mol object from the results. This resulted in a calculated affinity and a calculated pCHEMBL value. These docking scores were visualized in 3D and the key interactions were showed. Based on this, conclusions were drawn. We performed the molecular docking on all ligands based on the Lig_00 scaffold (Lig_01-Lig_05, Lig_10 and Lig_15), and found the best docking scores for Lig_15. However, this predicted pCHEMBL was still lower than the co-crystallized ligand.

*Challenges we faced and future prospect*
One of the challenges we faced during the execution of this project was that we often experiences a lack of time. Since the course only lasted for two weeks and we had to learn how to code, we did not manage to optimize the accuracy of the models we used. In the future we would like to execute the following things:

* Change the random split with a temporal split to make the model less overoptimistic
* Use a VSM model in the Machine Learning section (since this showed slightly better results)
* Generate novel molecules using DrugEx instead of drawing them ourselves based on promising ligands in the CHEMBL database
* Utilize interaction fingerprints besides the Morgan fingerprints to improve the accuracy of the docking
* Use a more accurate docking algorithm than AutoDockVina.

# How to install and run the project
To install and run the project, the following URL needs to be copied: https://github.com/NikitaKoster/Course-BPS. Subsequently, on Jupyter Notebook a +GitRepo should be added and the URL should be entered to copy the file. After this, each code block can be run and should work since everything that needs to be imported is described in the code blocks.

Packages/programs that were imported in order to make these codes work:
*Part 1*
import nglview
import os
import shutil
from Bio.PDB import PDBParser, PDBIO, Select,  PDBList, MMCIFParser, StructureAlignment
import Bio.Align
import os
from pathlib import Path
import rdkit

local scripts
from scripts import viewer
from scripts import bio_align

*Part 2*
from pathlib import Path
from warnings import filterwarnings
import time

import pandas as pd
import numpy as np
from sklearn import svm, metrics, clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import auc, accuracy_score, recall_score
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
Silence some expected warnings
filterwarnings("ignore")
Fix seed for reproducible results
SEED = 22
#seed_everything(SEED)

*Part 3*
install  open drug discovery toolkit (ODDT) and vina
import py3Dmol

import pandas as pd
import glob
import sys

from vina import Vina
#import pybel

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

#from meeko import MoleculePreparation
#from meeko import obutils

import MDAnalysis as mda
from MDAnalysis.coordinates import PDB

#import prolif
#from prolif.plotting.network import LigNetwork

import nglview
from scripts import viewer

import sys, os, shutil
sys.path.insert(1, '/project/jhllei001/JHL_data/Jupyter_Dock/utilities')
#from utils import fix_protein, getbox, generate_ledock_file, pdbqt_to_sdf, dok_to_sdf

import warnings
warnings.filterwarnings("ignore")
%config Completer.use_jedi = False

# Credits
We want to thank Willem Jespers for the interesting lectures and guidance during the project. We feel like we learned a lot! Also, we would like to thank Marina Gorostiola Gonzalez for answering our questions. Lastly, we would like to thank Leiden University for eventually adding us to the course :D
