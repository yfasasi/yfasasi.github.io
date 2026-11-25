---
title: "Variable-Cell Relaxation (vc-relax) of Bulk Copper Using ASE and Quantum ESPRESSO"
excerpt: "Learn the difference between ionic relaxation and variable-cell relaxation (vc-relax), and how to optimize bulk FCC copper's lattice constant using both ASE optimizers and Quantum ESPRESSO."
collection: portfolio
permalink: /portfolio/2026-11-25-vc-relax-bulk-copper-ase-qe
date: 2026-11-25
tags:
  - Quantum ESPRESSO
  - ASE
  - DFT
  - vc-relax
  - Python
---

Now that we know the `ecutwfc` and `ecutrho` values to use from the previous calculation, and we also know the optimal kpts, we can go ahead and do a vc-relax calculation

## Periodic DFT codes distinguish between two main types of optimizations:

### 1. Ionic relaxation

This is also called a `relax calculation` in Quantum Espresso (QE). And I guess it should stay the same way for VASP. In this type of optimization, you fix the cell volume and only allow the atomic positions to change to find the best atomic arrangement for the current configuration. This is activated with the ` "calculation":"relax" ` flag in QE. 

### 2. Cell relaxation

Also called a `variable-cell` relaxation in QE. This is where you allow both the atomic positions and the cell vectors (lattice vector) to find their true optimal configuration. This is activated with the ` "calculation":"vc-relax" ` flag in QE. This can be used for a new material where we dont know the lattice parameters beforehand or any type of calculations where finding the true optimal cell shape is essential. 

We will be using the Atomic simulation enviroment together with Quantum Espresso that we have already downloaded to carry out vc-relax calculation on `Bulk Copper.` 

First We will see how to use ASE to do everything. Then later we will see how to use ASE to only do the vc-relax while we handoff the scf part to QE.

> Heads up : Each cell in a Jupyter notebook must be executed in order. Otherwise, what you imported in the upper cell will not be carried into the bottom cell if you have not run the upper cell. Or you will have to copy and paste what is in the upper cell to the bottom cell.


First we are going to import all that we need from the ASE library. Importing the neccessary modules and functions is essential to do the calculations that we want to do


```python
from ase.optimize import BFGS #This import the BFGS algorithm that move the atoms around changing their coordinates

from ase.filters import ExpCellFilter,FrechetCellFilter #This import the cell vector relaxer. The FrechetCellFilter is new while the ExpCellFilter
#may be deprecated soon

import os,numpy as np,glob

from ase.build import bulk

from pymatgen.io.ase import AseAtomsAdaptor as converter

from mp_api.client import MPRester

from ase.calculators.espresso import EspressoProfile,Espresso #import the espresso calculator from ASE

#I figure it is faster loading the one atom copper from the materials projects rather than building our own from sratch especailly for those who
#have the material project API key and have installed the mp-api with `pip install mp-api` otherwise we can build it with the crystal()




```

# Then we build the bulk Copper as shown below


```python
with MPRester() as mpr:
    structure = mpr.get_structure_by_material_id("mp-30") #two things to note here: I didnt put my API keys inside the MPRester because I already
    # have it saved globally using one of the steps in the tutorial I published earlier. Alos, I had to search the mp-id for bulk copper to know it
    #is mp-30

bulk_Cu = converter.get_atoms(structure) #converting the obtained pymatgen atoms object to ASE atoms object so we can use an ASE calculator

#if you dont have an API key, use this function below. uncomment it by deleting the #

#bulk_Copper = bulk("Cu","fcc")

print(bulk_Cu.cell.cellpar())

#You can also confirm whether bulk_Cu and bulk_Copper print the same atoms object since we get them from different source.

print(np.isclose(bulk_Cu.cell.cellpar(),bulk_Copper.cell.cellpar(),atol=1e-2,rtol=1e-2)) #This prints True for all 6 (abc alpha beta gamma) 

#You must uncomment the bulk_Copper to compare. I was able to run this properly because I uncommented it earlier


   
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #00bdad; text-decoration-color: #00bdad">✓ </span>Retrieved 1 MaterialsDoc documents<span style="color: #7a7a7a; text-decoration-color: #7a7a7a"> in 0.3s</span>
</pre>



    [ 2.5296253   2.52962585  2.52962495 60.00001167 60.00000445 60.00000264]
    [ True  True  True  True  True  True]
    

# Setting up the Quantum Espresso Calculator
Now we will set up a calculator. To do this, we need to know the actual locations on your PC where you downloaded and extracted the Quantum Espresso pw.x into. The `pw.x` program is able to run most calculations such as vc-relax,scf,nscf. Then copy the file path by right clicking onm pw.x on the file menu, then select `copy as path` You will paste that path inside the EspressoProfile. Also, you need to know where you donwloaded the pseudopotentials. we will need the one for copper from there. copy the path as well for the pseudopotential folder. 

For me when I check, this is where I saved my downloaded pw.x on my PC (r"C:\Users\Test_PC\Downloads\qe-7.5-win-oneapi-msmpi\bin\pw.exe").  Also this is where my donwloaded and extracted pseudopotential is located (r"C:\Users\Test_PC\Downloads\all_pbe_UPF_v1.5"). You must preface the path with an `r` so window doesnt confuse the `\` for an escape. The pseudopotential is a PBE pseudopotential and we have already determine the optimal cutoff for the wavefunction to include in the planewave in the earlier tutorial with an scf calculation. 

1. You need to set up two things for the calculator to work: profile (which is set up with an EspressoProfile module) This is where you paste the path to your pseudopotential and the pw.x executable. You can use two main flags here: 1. `command`: paste the path to your pw.x executable here. and 2. `pseudo_dir`: paste the path to your pseudopotentials here.
2. Then you set up the calculator with the `Espresso` flag and put the profile you have built for the calculator inside the flag as `profile` = whatever you named the profile as.
3. There are specific cards and flags that goes inside a typical calculation where you put the name of the calculation that you want to run (whether it is vc-relax, relax, scf etc), and also specify the cutoff for the wavefunction,charge density,energy and density convergence criteria (so the calculation will know to stop once the change in energy between consecutive iteration meet those criteria. These values are typically set in Rydberg unit in QE. These flags and cards are set up with the `input_data' keyword inside the Espresso(). You set them up like a python dictionary. See below
>For the input_data dictionary, You dont have to nest it. For instance, `calculation` flag usually belong to the &CONTROL keyword in QE. But you dont have to put that there. ASE will rearrange every flags in the QE_calculation_flags below into their correct keyword. 


```python
#We have already imported all modules including the calculator from the ealier cell. And I have ran that cell so it is active. 
#Otherwise, just copy those imported modules and paste them here again

profile = EspressoProfile(command = r"C:/Users/Test_PC/Downloads/qe-7.5-win-oneapi-msmpi/bin/pw.exe", #note the backslash in command's path
                          pseudo_dir = r"C:\Users\Test_PC\Downloads\all_pbe_UPF_v1.5")

QE_calculation_flags = {"calculation":"vc-relax", #this set up the calculation type
                        "prefix":"Cu_vc_relax", #set up the name of the folder that store the charge density and wavefunctions. you can delete these files omce your calculation finishes.
                        "outdir":"./", #where do you want to store the temp folder? This says "I will save it in the current working dir the prefix folder will be located inside this folder
                       "etot_conv_thr":1e-6,
                       "forc_conv_thr":2e-5,
                       "conv_thr":1e-7,
                       "input_dft":"pbe", #This sets up the DFT functional (method) that you want to use. 
                       "occupations":"smearing", #smearing aid convergence in metals. electrons near the top of the valence band tend to behave weirdly
                        # we need to give them fractional occupations. the smearing method i am using here is gaussian and degauss value set the width
                       "smearing":"gaussian",
                       "degauss":0.01,
                       "ecutwfc":60, # There was a tutorial where we determine the best value for this, and also for the ecutrho, with an scf calculation
                       "ecutrho":600}
from copy import deepcopy

QE_vc_flags = deepcopy(QE_calculation_flags)

QE_vc_flags.update({"calculation":"scf"}) #And only pass this new input data

QE_vc_flags.update({"tstress":True,
                    "tprnfor":True}) #I got an error when i didnt turn on explicit force and stress printing. Because it it scf QE is doing, forc
# and stress arent automatically turn on for scf calculation. 

calculator = Espresso(profile = profile,
                      input_data = QE_vc_flags,
                      directory = os.path.join(os.getcwd(),"Bulk_Cu_vc_relax"),
                      kpts = (6,6,6), # we also have a tutorial where we converged the kpts. And we have sen that kpts around 10x10x10 works best.
                      #but I am using 6,6,6 here for faster calculations. 
                     pseudopotentials={"Cu":"cu_pbe_v1.2.uspp.F.UPF"}) 
                      
```

# Attaching the calculator

1. If you want ASE to handle the cell relaxation with its own BFGS and ExpCellFilter, then you will need to set the `calculation` to `scf`. so that Quantum Espresso only update the forces and stress on each atoms and hand that over to ASE optimizer
2. If you want Quantum Espresso to handle eveyrything, thats fine. Just leave calculation as vc-relax since we are doing a vc-relax calucltaion
   > I got an error when i didnt turn on explicit force and stress printing. Because QE is handling the SCF calculation, force and stress arent automatically turn on for scf calculation.


```python
#Now attach the calculator to the system you are working on. bulk_Cu

bulk_Cu.calc = calculator #This is how you attach the caluclator profile you set up

#And you can change the calculation to scf since you are letting ASE handle the calculation. Just update the dictionary





```

# 1. Letting ASE handle the calculation , and only allow QE to return the forces at each step


```python
optimizer = ExpCellFilter(bulk_Cu)

atom_relaxer = BFGS(optimizer,trajectory = "bulk_Cu_relax.traj")


atom_relaxer.run(fmax=0.01)
```

    C:\Users\Test_PC\AppData\Local\Temp\ipykernel_67880\213379754.py:1: DeprecationWarning: Use FrechetCellFilter for better convergence w.r.t. cell variables.
      optimizer = ExpCellFilter(bulk_Cu)
    

          Step     Time          Energy          fmax
    BFGS:    0 07:43:40    -5493.701156        0.505610
    BFGS:    1 07:44:05    -5493.711733        0.254795
    BFGS:    2 07:44:29    -5493.716923        0.131404
    BFGS:    3 07:45:05    -5493.716841        0.184742
    BFGS:    4 07:45:27    -5493.717556        0.067515
    BFGS:    5 07:45:44    -5493.717112        0.071204
    BFGS:    6 07:46:09    -5493.717428        0.065130
    BFGS:    7 07:46:48    -5493.710475        0.477950
    BFGS:    8 07:47:18    -5493.717562        0.019951
    BFGS:    9 07:47:50    -5493.717554        0.056131
    BFGS:   10 07:48:16    -5493.717559        0.154890
    BFGS:   11 07:48:40    -5493.717557        0.018750
    BFGS:   12 07:49:07    -5493.717556        0.072884
    BFGS:   13 07:49:32    -5493.717557        0.018816
    BFGS:   14 07:50:07    -5493.717557        0.142464
    BFGS:   15 07:50:31    -5493.717557        0.037424
    BFGS:   16 07:51:02    -5493.717557        0.124618
    BFGS:   17 07:51:26    -5493.717557        0.073584
    BFGS:   18 07:51:45    -5493.717555        0.117456
    BFGS:   19 07:52:09    -5493.717557        0.029600
    BFGS:   20 07:52:31    -5493.717557        0.000266
    




    np.True_



Since we have run this calculation inside the current working directory and you have specified, inside the espresso caluclator setup that the folder should be saved as `Bulk_Cu_vc_relax`, you either need to go to your file directory on windows to look for your output files by following the steps below,  or you can just press `ctrl + b` on your keyboard right within jupyter notebook to open the file menu. This will load the files and directories inside the current working directory where your jupyter notebook is being hosted. 

For me, I typed ctrl + b, and I then look at the left, I can see a folder named "Bulk_Cu_vc_relax", all the input and output files will be saved here. The output file might be named as "espresso.pwo", and the input file could appear as "espresso.pwi", and the trajectory file will show up as `bulk_Cu_relax.traj.` since this is what we specified inside the BFGS setup above. 

The trajectory file collect each step of the iteration into a trajectory that monitors what happens to the system from trajectory = 0 (which is essntially the input file as it is) to the final trajectory where all our set criteria have been met and the system is now fully relaxed and done. 

You can open the output file "espresso.pwo", scroll towards the end and check the last result to see the new atomic positions and the lattice vectors if you have explicitly specified the calculation type as `vc-relax` inside the calculation flag in input_data, but rather we did an scf. So when you open the output file, it will only contain electronic eigenvalues and stress tensors. If you did an vc-relax inside the input flag, be aware that Quantum Espresso doesnt print the lattice parameters as scalar quantity. It print the lattice as 3x3 vector where row 1 = lattice a, row2 = lattice b and row3 = lattice c. You will need to get the norm of these vector to get the scalar. 

To get the cell vector (and lattice parameters) and the atomic coordinates, you will need to open the trajectory file saved as `bulk_Cu_relax.traj`. The trajectory file isnt human readable. But you can open it with `ase.io` module. 


```python
import os,glob
from ase.io import read

#I want to use glob to find the traj file to save me the stress
final_optimized_cu = read(glob.glob(os.path.join(os.getcwd(),"*traj*"))[0],index = -1,format="traj")

#glob function returns a list even if it only finds one item which is true in our case. [0] get the first item in the list.
print("\nFinal relaxed atoms object:")

#This returns the final optimized atoms object since we passed index=-1 to grab the last one.If you pass index = 0 it grabe the intial atoms object

print(final_optimized_cu.cell.cellpar()) #This helps us get the scalar lattice vectors (a b c alpha beta gamma) from the cell vector.

unrelaxed_cell_at_index_0 = read(glob.glob(os.path.join(os.getcwd(),"*traj*"))[0],index = 0,format="traj")

print("\ninitial un-relaxed atoms object from the material project:")
print(unrelaxed_cell_at_index_0.cell.cellpar()) #comparing the imported atoms objectv from material project with the new optiomized atoms object.
```

    
    Final relaxed atoms object:
    [ 2.57476212  2.57476274  2.57476075 60.00002544 60.00001755 59.99980469]
    
    initial un-relaxed atoms object from the material project:
    [ 2.5296253   2.52962585  2.52962495 60.00001167 60.00000445 60.00000264]
    





**2. The second method is to let Quantum Espresso drive the entire relaxation where you specify the calculation type as `vc-relax` in `calculation`. This way the output file will contain the relaxed atomic coordinates and the cell vector like We have understood above. This can be done as follow:**


```python
#Since we have already configured the Quantum Espresso input flags in the earlier cells, and we have already ran those cells,
# we dont have to copy those again. We can use python dictionary methods to update and delete. 

QE_calculation_flags.update({"prefix":"Cu_QE_vc_relax"}) # setting a different name for the prefix to differentiate it from the earlier calculation

#print(QE_calculation_flags)

# we can reattach the calculator to bulk. And the reason i am resetting the calculator is to be able to set a different folder for "directory" thand we did earlier. 

calculator = Espresso(profile = profile,
                      input_data = QE_calculation_flags,
                      directory = os.path.join(os.getcwd(),"QE_Bulk_Cu_vc_relax"),
                      kpts = (6,6,6), # we also have a tutorial where we converged the kpts
                     pseudopotentials={"Cu":"cu_pbe_v1.2.uspp.F.UPF"}) 

bulk_Cu.calc = calculator

#Now that the calculation card has the "vc-relax" keyword, we dont have to use ExpCellFilter and ASE BFGS. QE will set the calculation up on its own

print("Relax. Your calculation is getting set up. I will print the energy in a bit")

print("\nThe final optimized energy for single element Copper is",bulk_Cu.get_potential_energy()) #This activates the calculation here

#Note that the calculation will take a bit more time to get done than a typical scf calculation. This is because vc-relax moves the atomic
#coordinates and the lattice parameters around to find the best configuration and arrangement where the forces on all atoms are minimal, almost zero.
```

    Relax. Your calculation is getting set up. I will print the energy in a bit
    The final optimized energy for single element Copper is -5493.716899017709
    

#Now the calculation is done and the energy has been printed. You can check the folders named `QE_Bulk_vc_relax.` click on it, you will see all the wavefunction files, the input files and also the output file named `espresso.pwo.` click on the output file, scroll down to the latest optimized cell logs, and check for the word `End of BFGS Geometry Optimization,` and you can see the cell vectors and the new atomic coordinates which you can copy and paste in a new input file for a different calculation.

* Alternatively, you can read the optimized parameters by using a few line of codes as shown below: 


```python
_file = [file for file in os.listdir(os.getcwd()) if "QE_Bulk" in file]
#print(output_file)
output_file = []
for root,folders,files in os.walk(_file[0]):
    for file in files:
        if file.endswith("pwo") and "checkpoint" not in file:
            output_file.append(os.path.join(root,file))
            
print(read(output_file[0],index=-1,format="espresso-out").get_potential_energy())

#The lines above can be condensed to this if you are able to copy the full path:

total_electronic_energy = read("QE_Bulk_Cu_vc_relax/espresso.pwo", index = -1, format="espresso-out") #This saves the last fully relaxed atoms object

# You can right click on the name of the output file, and select "copy as path", then paste the full path just like it is shown above"

print("The total electronic_energy (potential_energy) is",total_electronic_energy.get_potential_energy())

#Note that the kinetic energy will be zero because in Born-Oppenheimer approximation,
#the atomic nuclei are treated as stationary (zero nuclear kinetic energy) while the electrons relax into their lowest energy state around them.

print("The kinetic energy is",total_electronic_energy.get_kinetic_energy())
```

    -5493.716899017709
    The total electronic_energy (potential_energy) is -5493.716899017709
    The kinetic energy is 0.0
    


```python

```
