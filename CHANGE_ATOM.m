%% Gordon Chalmers 10/21

function new_chm = CHANGE_ATOM(chm)

%% CHANGE_ATOM.m in branch or ring

global atom_type; 
global atom_val;
global total_atom_types; 
global percent_atom; 
global alphabet;
global max_tries;
global mutation_type_probability;
global cutoff_bond_percentage;
global inter_bond_distance;
global min_heavy_atoms;
global corina_path;
global ligand_dir;

chm=char(chm);
chm=RING_RENUMBER_CHECK(chm);

attempt=0;
success=0;

%% load molecule information
[molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ... 
    chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);

if num_heavy_atoms<min_heavy_atoms
    attempt=max_tries;
    success=0;
end

%% attempts until max_tries is reached due to random numbers of heavy atoms
while attempt<max_tries && success==0
    
    %% generally not possible to always change atom types
    %%  while statement until a new_chm~=chm
    new_chm='';
    
    %% random heavy atom in chm
    random_atom=ceil(rand*num_heavy_atoms);
    
    %%  new atom type is random from abundance percentage
    total_percent=0;
    prob_atom=rand;
    for i=1:total_atom_types
        total_percent=total_percent+percent_atom(i);
        if prob_atom<=total_percent
            if prob_atom>total_percent-percent_atom(i)
                atom_type_insert=atom_type(i);
                new_val_insert=atom_val(i); 
            end
        end
    end
    
    %% determine valence of atom heavy_idx_chm(heavy_atom)
    for i=1:total_atom_types
        if chm(heavy_idx_chm(random_atom))==atom_type(i)
            val_init=atom_val(i);
        end
    end
    
    insert_idx=heavy_idx_chm(random_atom);
    if new_val_insert>=adj_heavy(random_atom)
        new_chm(1:insert_idx-1)=chm(1:insert_idx-1);
        new_chm(insert_idx)=atom_type_insert;
        new_chm(insert_idx+1:chm_len)=chm(insert_idx+1:chm_len);
    end
   
    %% create pdb file
    system("rm "+ligand_dir+"/molecule/*.*");
    fileID=fopen(ligand_dir+"/molecule/molecule.smi",'w');
    fprintf(fileID,'%s',char(new_chm));
    fclose(fileID);
    
    system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,axchir,msi=50,msc=10,names,preserve "+ligand_dir+"/molecule/molecule.smi "+ligand_dir+"/molecule/molecule.pdb");
    %% if ok the file test.pdb will exist
    success=0;
    if exist(ligand_dir+"/molecule/molecule.001.pdb")>0
        if dir(ligand_dir+"/molecule/molecule.001.pdb").bytes>0
            success=1;
        end
    end
    
    %% check if a bond clash between non-bonded atoms in the geometry from the pdb
    %%  if pdb exists from previous
    if success==1
        no_clash=CLASH_CHECK(new_chm);
        if no_clash==0 || strcmp(new_chm,chm)==1
            success=0;
        end
    end
   
    attempt=attempt+1;
    
end  %% while attempt

if success==0
    new_chm=chm;
end

%% attempt
%% success

end 

