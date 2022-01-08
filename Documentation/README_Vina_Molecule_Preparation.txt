
Gordon Chalmers 10/21

This file explains the preparation of the protein.pdbqt file from the protein.pdb file and 
the centering in docking.

--

Vina requires an .pdbqt file.  This is found from the .pdb by using the commands (in Matlab):

prepare_receptor4_path="/home/gordon/mgltools/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py";

An example is, from COX-2 PDB_ID 4ph9 with ligands and solvent deleted : 
protein_pdb_in="/home/gordon/Ligand_GA/version2/Aspirin/4ph9_no_IBP_HOH.pdb"
protein_pdb_out="/home/gordon/Ligand_GA/version2/Aspirin/4ph9_no_IBP_HOH.pdbqt"

protein_molecule_format_conversion=prepare_receptor4_path+" -r " + protein_pdb_in + "-A hydrogens -o " + protein_pdb_out;
system(char(protein_molecule_format_conversion));

--

The output file has the molecules at each iteration.  This can be used in the function Ligand_GA or 
used independently, 

clear sorted_unique_population;
clear sorted_unique_fitness;

[sorted_unique_population,sorted_unique_fitness] = Ligand_GA_Output_Function(OutFileName);
number_pop=size(sorted_unique_fitness,1);
chm_fittest=sorted_unique_population(number_pop);
[molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm ...
    ,adj_atom,avail_heavy_bond,chiral,num_bonds_left,num_bonds_right] = MoleculeStructure(chm_fittest);

The sorted unique lists have the unique molecules and fitnesses.  The OutFileName is a .mat file generated from 
the software having the complete set of molecules and fitnesses.  The molecules are in non-isomeric SMILES format.  

--

Note that these molecules are in a non-isomeric representation.  This means that stereoisomer information is 
not contained in the information.  To find this, use the Ligand_GA_Molecule_Fitness_Stereoisomer.m function.  The 
output has the stereosimoer fitnesses for the input molecule chm.

stereoisomer_fitness = Ligand_GA_Molecule_Fitness_Stereoisomer(chm)

where chm is the input molecule and stereoisomer_fitness is an array of fitnesses.  Keep in mind that the Vina 
docking algorithm uses random initial starting points from the molecule and may have to be re-run to find the 
most bound pose of the molecule - this is done by running Ligand_GA_Molecule_Fitness_Stereoisomer and examining 
the results each time.  

-- 

Viewing the superimposed protein-ligand comples 

Load the protein pdbqt file that was used in the config file.  Load a ligand pdb file, or prefered, 
the ligand pdbqt file.  The pdbqt file has different poses if there are more than 1 and the default 
is a maximum of 9.  

Translate the ligand to the center coordinates from the config file.  In PyMol, type 

translate [center_x, center_y, center_], name_of_ligand_file  

This will put the ligand at the center point that was used to find the binding.  The view 
can be changed from there.

Type: 

center xxx/  

on the residue of the protein that the location is at to center the view there.  



