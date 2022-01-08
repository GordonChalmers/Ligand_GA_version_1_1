%% Gordon Chalmers 10/21

function new_chm = DELETE_ATOM(chm)

%%% atom_delete, no ring or branch points

global atom_type;
global atom_val;
global total_atom_types;
global percent_atom;
global alphabet;
global max_tries;
global mutation_type_probability;
global min_heavy_atoms;
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


%% branch and ring points not deleted
while attempt<max_tries && success==0
    
    new_chm='';
    
    %% random atom from sequence
    random_atom=ceil(rand*num_heavy_atoms);
    atom_delete_idx=heavy_idx_chm(random_atom);
    branch=0;
    
    %% remove atom
    %% single atom branch
    if atom_delete_idx-1>=1 && chm(atom_delete_idx-1)=="(" && atom_delete_idx+1<=chm_len && chm(atom_delete_idx+1)==")"
        if atom_delete_idx-2>=1
            new_chm(1:atom_delete_idx-2)=chm(1:atom_delete_idx-2);
        end
        if atom_delete_idx+2<=chm_len
            new_chm(atom_delete_idx-1:chm_len-3)=chm(atom_delete_idx+2:chm_len);
        end
        branch=1;
    end
    if atom_delete_idx-2>=1 && chm(atom_delete_idx-2:atom_delete_idx-1)=="(=" && atom_delete_idx+1<=chm_len && chm(atom_delete_idx+1)==")"
        if atom_delete_idx-3>=1
            new_chm(1:atom_delete_idx-3)=chm(1:atom_delete_idx-3);
        end
        if atom_delete_idx+2<=chm_len
            new_chm(atom_delete_idx-2:chm_len-4)=chm(atom_delete_idx+2:chm_len);
        end
        branch=1;
    end
    if atom_delete_idx-2>=1 && chm(atom_delete_idx-2:atom_delete_idx-1)=="(#" && atom_delete_idx+1<=chm_len && chm(atom_delete_idx+1)==")"
        if atom_delete_idx-3>=1
            new_chm(1:atom_delete_idx-3)=chm(1:atom_delete_idx-3);
        end
        if atom_delete_idx+2<=chm_len
            new_chm(atom_delete_idx-2:chm_len-4)=chm(atom_delete_idx+2:chm_len);
        end
        branch=1;
    end
    
    %% general atom, not ring point
    delete_left=0;
    delete_right=0;
    if branch==0
        if atom_delete_idx+1<=chm_len && chm(atom_delete_idx+1)~="(" && isnan(str2double(chm(atom_delete_idx+1)))==1   %% not a branch or ring point atom
            if atom_delete_idx+1<=chm_len && chm(atom_delete_idx+1)=="=" || chm(atom_delete_idx+1)=="#"  %% take away the double or triple bond on right
                delete_right=1;
            end
            if atom_delete_idx-1>=1 && (chm(atom_delete_idx-1)=="=" || chm(atom_delete_idx-1)=="#")  %% take away the double or triple bond on left
                delete_left=1;
            end
            if atom_delete_idx-delete_left-1>=1
                new_chm(1:atom_delete_idx-delete_left-1)=chm(1:atom_delete_idx-delete_left-1);
            end
            new_chm(atom_delete_idx-delete_left:chm_len-delete_right-delete_left-1)=chm(atom_delete_idx+delete_right+1:chm_len);
        end
        
        if atom_delete_idx==chm_len  %% end of chm with no branch or ring point
            if atom_delete_idx-1>=1 && (chm(atom_delete_idx-1)=="=" || chm(atom_delete_idx-1)=="#")  %% take away the double or triple bond on left
                delete_left=1;
            end
            if atom_delete_idx-delete_left-1
                new_chm(1:atom_delete_idx-delete_left-1)=chm(1:atom_delete_idx-delete_left-1);
            end
            new_chm(atom_delete_idx-delete_left:chm_len-delete_right-delete_left-1)=chm(atom_delete_idx+delete_right+1:chm_len);
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
    
    %% check minimum number of heavy atoms   
    if success==1
        
        [new_molecule,new_chm_len,new_adj,new_adj_heavy,new_num_heavy_atoms,new_heavy_atom_list,new_heavy_idx_chm, ...
            new_num_rings,new_ring_idx_chm,new_adj_atom,new_avail_heavy_bond,new_chiral,new_num_bonds_left,new_num_bonds_right]=MoleculeStructure(new_chm);
        
        if new_num_heavy_atoms<min_heavy_atoms
            success=0;
        end        
    end
    
    attempt=attempt+1;
    
end  %% attempt

if success==0
    new_chm=chm;
end

%% attempt
%% success

end

