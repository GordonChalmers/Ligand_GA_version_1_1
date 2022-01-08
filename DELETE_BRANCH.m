%% Gordon Chalmers 10/21

function new_chm = DELETE_BRANCH(chm)

%% DELETE_BRANCH
%% delete a branch without any ring numbers open,  without any ( open

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

while attempt<max_tries && success==0
    new_chm='';
    %% find location of ('s and )'s
    %% branch is defined by an opening and closing of ( and ) without any ring or ( open
    num_left_pars=0;
    num_right_pars=0;
    left_pars_idx_chm=zeros(1,1);
    right_pars_idx_chm=zeros(1,1);
    for i=1:chm_len
        if chm(i)=="("
            num_left_pars=num_left_pars+1;
            left_pars_idx_chm(num_left_pars)=i;
        end
        if chm(i)==")"
            num_right_pars=num_right_pars+1;
            right_pars_idx_chm(num_right_pars)=i;
        end
    end
    total_pars=num_left_pars+num_right_pars;
    
    rand_left_pars=ceil(rand*num_left_pars);
    if rand_left_pars>0
        
        %% find the nearest closing )
        %% take into account the closing of ring numbers inside the ( and )
        add_par=0;
        ring_total=zeros(num_rings,1);
        done=0;
        initial_idx_branch=left_pars_idx_chm(rand_left_pars);
        idx_branch=left_pars_idx_chm(rand_left_pars);
        closing_branch=0;
        
        %% find nearest closing right pars
        for idx_branch=initial_idx_branch:chm_len
            %% for idx_branch=rand_left_pars:max(right_pars_idx_chm) - not used
            if chm(idx_branch)=="("
                add_par=add_par+1;
            end
            if chm(idx_branch)==")"
                add_par=add_par-1;
            end
            %% add the ring points %2 at ring_number means closed
            %%  if norm(mod(ring_total),2)==0 then no open rings to point idx_branch
            if isnan(str2double(chm(idx_branch)))==0
                ring_total(str2double(chm(idx_branch)))=ring_total(str2double(chm(idx_branch)))+1;
            end
            if add_par==0 && norm(mod(ring_total,2))==0
                closing_branch=idx_branch;  %% break the loop
                break;
            end
        end
        
        %% deleting branch reduces bond number in the branch point and chirality
        %% if closing_branch==0 then no branch (e.g. morphine except for the 1st in sequence before any rings and branches)
        if closing_branch~=0
            new_chm(1:initial_idx_branch-1)=chm(1:initial_idx_branch-1);
            new_chm(initial_idx_branch:chm_len-closing_branch-1+initial_idx_branch)=chm(closing_branch+1:chm_len);  %% deleted (...) w/o open rings
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
        new_chm=RING_RENUMBER_CHECK(new_chm);
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
    
end  %% while attempt

if success==0
    new_chm=chm;
end

%% attempt
%% success

end

%% could add a constraint on limit of number of atoms in deleted branch
