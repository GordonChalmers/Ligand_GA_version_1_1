%% Gordon Chalmers 10/21

function new_chm = SINGLE_DOUBLE_BOND(chm)

%% SINGLE_DOUBLE_BOND  changes bond type from single to double
%%				on l.h.s. of heavy atom

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
    
    %% random atom from sequence
    random_atom=ceil(rand*num_heavy_atoms);
    atom_bond_change_idx=heavy_idx_chm(random_atom);
    
    if random_atom>1
        %% check both sides of possible = for avail bond ok
        %%  if branch then the connect is the minimum of the connectivity of random_atom to other atoms
        min_conect_atom=1000;
        for bond=1:6
            if adj(random_atom,bond)<min_conect_atom && adj(random_atom,bond)~=0
                min_conect_atom=adj(random_atom,bond);
            end
        end
        
        if avail_heavy_bond(random_atom)>=1 && avail_heavy_bond(min_conect_atom)>=1 && chm(atom_bond_change_idx-1)~="=" && chm(atom_bond_change_idx-1)~="#"
            new_chm(1:atom_bond_change_idx-1)=chm(1:atom_bond_change_idx-1);
            new_chm(atom_bond_change_idx)="=";
            new_chm(atom_bond_change_idx+1:chm_len+1)=chm(atom_bond_change_idx:chm_len);
        end
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
            %% check if a bond clash between non-bonded atoms in the geometry from the pdb
            no_clash=CLASH_CHECK(new_chm);
            if no_clash==0
                success=0;
            end
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

%% bonds, rings, and branches are taken into account
