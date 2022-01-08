%% Gordon Chalmers 10/21

function new_chm = CLOSE_RING(chm)

%% CLOSE_RING
%%  creates a ring from 2 atoms in the molecule

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

while attempt<max_tries && success==0 && num_rings<9
    
    new_chm='';
    
    %% closing and adding a ring
    random_atom_one=ceil(rand*num_heavy_atoms);
    random_atom_two=ceil(rand*num_heavy_atoms);
    
    if random_atom_two<random_atom_one
        temp=random_atom_two;
        random_atom_two=random_atom_one;
        random_atom_one=temp;
    end
    
    one_idx_chm=heavy_idx_chm(random_atom_one);
    two_idx_chm=heavy_idx_chm(random_atom_two);
    
    if avail_heavy_bond(random_atom_one)>0 && avail_heavy_bond(random_atom_two)>0 && abs(random_atom_two-random_atom_one)>2
        
        new_chm(1:one_idx_chm)=chm(1:one_idx_chm);   %% add first segment
        atom=one_idx_chm;
        num1=0;
        while atom+num1+1<chm_len  %% add previous rings
            if isnan(str2double(chm(atom+num1+1)))==0
                new_chm(atom+num1+1)=chm(one_idx_chm+num1+1);
                num1=num1+1;
            else
                break;
            end
        end
        new_chm(atom+num1+1)=int2str(num_rings+1);  %% add new ring number
        insert=atom+num1+1;
        
        new_chm(one_idx_chm+num1+2:two_idx_chm+1)=chm(one_idx_chm+num1+1:two_idx_chm);  %% add middle segment
        insert=two_idx_chm+1;
        atom=two_idx_chm;
        num2=0;
        while atom+num2+1<chm_len  %% add previous rings
            if isnan(str2double(chm(atom+num2+1)))==0
                new_chm(insert+num2+1)=chm(two_idx_chm+num2+1);
                num2=num2+1;
            else
                break;
            end
        end
        insert=insert+num2;
        new_chm(insert+1)=int2str(num_rings+1);   %% add new ring number
        insert=insert+1;
        new_chm(insert+1:chm_len-two_idx_chm-num2+insert)=chm(two_idx_chm+num2+1:chm_len);  %% add final segment
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
        %% renumber_check rings - chm
        new_chm=RING_RENUMBER_CHECK(new_chm);
        no_clash=CLASH_CHECK(new_chm);
        if no_clash==0
            success=0;
        end
    end
    
    attempt=attempt+1;
    
end %% possible

if success==0
    new_chm=chm;
end

%% attempt
%% success

end

%% no need to renumber as largest ring number is added
%% limit of 9 rings
%% 3 integers maximum at a heavy atom
%% chirality taken into account
%% double/triple bonds taken into account
