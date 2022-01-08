%% Gordon Chalmers 10/21

function new_chm = OPEN_RING(chm)

%% OPEN_RING at ring points

global atom_type; 
global atom_val;
global total_atom_types; 
global percent_atom; 
global alphabet;
global max_tries;
global min_heavy_atoms;
global mutation_type_probability;
global cutoff_bond_percentage;
global inter_bond_distance;
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

while attempt<max_tries && success==0
    
    new_chm='';
    %% random ring in sequence
    random_ring=ceil(rand*num_rings);
    %% remove ring number
    idx=0;
    for atom=1:chm_len
        if chm(atom)~=int2str(random_ring)
            idx=idx+1;
            new_chm(idx)=chm(atom);
        end
    end
    
%% renumber_check rings 
    new_chm=RING_RENUMBER_CHECK(new_chm);
    
    attempt=attempt+1;

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
        if no_clash==0
            success=0;
        end
    end
    
end

if success==0
    new_chm=chm;
end

%% attempt
%% success

end

%% all chiral/nonchiral cases
%% with multiple rings to 9
%% includes all bond types
