%% Gordon Chalmers 10/21

function new_chm = ADD_ATOM(chm)

%% ADD_ATOM to a branch or ring from l.h.s. of idx(random_atom)

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

attempt=0;
success=0;

chm=char(chm);
chm=RING_RENUMBER_CHECK(chm);

%% load molecule information
[molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ...
    chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);

if num_heavy_atoms<min_heavy_atoms 
    attempt=max_tries;
    success=0;
end

while attempt<max_tries && success==0
    
    new_chm='';   
    %% inserted atom type is random from abundance percentage
    total_percent=0;
    prob_atom=rand;
    for i=1:total_atom_types
        total_percent=total_percent+percent_atom(i);
        if prob_atom<=total_percent
            if prob_atom>total_percent-percent_atom(i)
                %% insert atom type: e.g., C,N,O,P,S, ... (F,Cl,Br,I different and only at ends of branches)
                atom_type_insert=atom_type(i);
                new_val_insert=atom_val(i);  %% available electrons of nonbonded atom
            end
        end
    end
    
    random_atom=ceil(rand*num_heavy_atoms);    
    add_atom=heavy_idx_chm(random_atom);
    
    %% single bond not first heavy atom
    if random_atom>1 && chm(add_atom-1)~="=" && chm(add_atom-1)~="#"
        new_chm(1:add_atom-1)=chm(1:add_atom-1);
        new_chm(add_atom:add_atom)=atom_type_insert;
        new_chm(add_atom+1:chm_len+1)=chm(add_atom:chm_len);
    end
     
    %% first atom - check avail bond
    if random_atom==1 && avail_heavy_bond(random_atom)>0
        new_chm(1:1)=atom_type_insert;
        new_chm(2:chm_len+1)=chm(1:chm_len);
    end
    
    %% double bond 
    if random_atom>1 && chm(add_atom-1)=="=" && new_val_insert>=3
        test=rand; 
        if test>.5
        %% right of =
            new_chm(1:add_atom-1)=chm(1:add_atom-1);
            new_chm(add_atom:add_atom)=atom_type_insert;
            new_chm(add_atom+1:chm_len+1)=chm(add_atom:chm_len);
        end
        if test<=.5 
        %% left of =
            new_chm(1:add_atom-2)=chm(1:add_atom-2);
            new_chm(add_atom-1:add_atom-1)=atom_type_insert;
            new_chm(add_atom:chm_len+1)=chm(add_atom-1:chm_len);
        end
    end
    
    %% triple bond 
    if random_atom>1 && chm(add_atom-1)=="#" && new_val_insert>=4
        test=rand; 
        if test>.5
        %% right of #
            new_chm(1:add_atom-1)=chm(1:add_atom-1);
            new_chm(add_atom:add_atom)=atom_type_insert;
            new_chm(add_atom+1:chm_len+1)=chm(add_atom:chm_len);
        end
        if test<=.5 
        %% left of #
            new_chm(1:add_atom-2)=chm(1:add_atom-2);
            new_chm(add_atom-1:add_atom-1)=atom_type_insert;
            new_chm(add_atom:chm_len+1)=chm(add_atom-1:chm_len);
        end
    end
    
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
