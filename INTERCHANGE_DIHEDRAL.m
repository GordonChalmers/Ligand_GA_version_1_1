%% Gordon Chalmers 10/21

function [new_chm1,new_chm2] = INTERCHANGE_DIHEDRAl(chm1,chm2)

chm1=char(chm1);
chm2=char(chm2);

%% INTERCHANGE_BRANCH of 2 chms
%% this function is the single-point crossover function of 2 chromosomes

global atom_type;
global atom_val;
global total_atom_types;
global percent_atom;
global alphabet;
global max_tries;
global mutation_type_probability;
global success;
global min_heavy_atoms;
global cutoff_bond_percentage;
global inter_bond_distance;
global prepare_ligand4_path;
global corina_path;
global ligand_dir;

%% chm1
%% chm2

%% display('INTERCHANGE_BRANCH and 2 molecules')
chm1=RING_RENUMBER_CHECK(chm1);
chm2=RING_RENUMBER_CHECK(chm2);

attempt=max_tries;
success=0;

[molecule1,chm_len1,adj1,adj_heavy1,num_heavy_atoms1,heavy_atom_list1,heavy_idx_chm1,num_rings1,ring_idx_chm1,adj_atom1,avail_heavy_bond1, ...
    chiral1,num_bonds_left1,num_bonds_right1]=MoleculeStructure(chm1);

[molecule2,chm_len2,adj2,adj_heavy2,num_heavy_atoms2,heavy_atom_list2,heavy_idx_chm2,num_rings2,ring_idx_chm2,adj_atom2,avail_heavy_bond2, ...
    chiral2,num_bonds_left2,num_bonds_right2]=MoleculeStructure(chm2);


%% first prepare the pdbqt files and the dihedral information

%% chm -> molecule.pdb -> molecule.pdbqt

    %% create pdb file
    system("rm "+ligand_dir+"/molecule/*.*");    
    fileID1=fopen(ligand_dir+"/molecule/molecule1.smi",'w');
    fprintf(fileID1,'%s',char(chm1));
    fclose(fileID1);

    system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,axchir,msi=50,msc=10,names,preserve "+ligand_dir+"/molecule/molecule1.smi "+ligand_dir+"/molecule/molecule1.pdb");
    %% create pdb file
    fileID2=fopen(ligand_dir+"/molecule/molecule2.smi",'w');
    fprintf(fileID2,'%s',char(chm2));
    fclose(fileID2);

    system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,axchir,msi=50,msc=10,names,preserve "+ligand_dir+"/molecule/molecule2.smi "+ligand_dir+"/molecule/molecule2.pdb");
    %% if ok the pdb.001 file pdb will exist

%% call to to prepare_ligand4.py to make molecule.fileno -> ligand.pdbqt
%% ../mgltools/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -r ligand.pdb -o ligand.pdbqt
%% e.g., prepare_ligand4_path="/../mgltools/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
molecule_format_conversion1=prepare_ligand4_path+" -l "+ligand_dir+"/molecule/molecule1.001.pdb -o "+ligand_dir+"/molecule/molecule1.pdbqt";
system(char(molecule_format_conversion1));

molecule_format_conversion2=prepare_ligand4_path+" -l "+ligand_dir+"/molecule/molecule2.001.pdb -o "+ligand_dir+"/molecule/molecule2.pdbqt";
system(char(molecule_format_conversion2));

%% check if files exist - default is attempt=max_tries unless files are there

if exist(ligand_dir+"/molecule/molecule1.pdbqt")>0 && dir(ligand_dir+"/molecule/molecule1.pdbqt").bytes>0 && exist(ligand_dir+"/molecule/molecule2.pdbqt")>0 && dir(ligand_dir+"/molecule/molecule2.pdbqt").bytes>0
     
    attempt=0;
    success=0;
    
    %% find the dihedral atom pairs
    fid=fopen(ligand_dir+"/molecule/molecule1.pdbqt");
    frewind(fid);
    
    fgetl(fid);
    fgetl(fid);
    text=fgetl(fid);
    atom1_char_molecule1='';
    atom2_char_molecule1='';
    dihedral_atom1_molecule1=zeros(1,100);
    dihedral_atom2_molecule1=zeros(1,100);
    atom_idx=0;
    while(text(1:4)=='REMA')
        split_text=strsplit(text);
        if isnan(str2double(split_text{2}))~=1  %% not an inactive dihedral bond
            atom_idx=atom_idx+1;
            atom1_char_molecule1(atom_idx,1:length(split_text{6}))=split_text{6};
            atom2_char_molecule1(atom_idx,1:length(split_text{8}))=split_text{8};
            atom1_num_char_molecule1=strsplit(atom1_char_molecule1(atom_idx,:),"_");
            atom2_num_char_molecule1=strsplit(atom2_char_molecule1(atom_idx,:),"_");
            dihedral_atom1_molecule1(atom_idx)=str2double(atom1_num_char_molecule1{2});
            dihedral_atom2_molecule1(atom_idx)=str2double(atom2_num_char_molecule1{2});
        end
        text=fgetl(fid);
    end
    fclose(fid);
    number_dihedrals_molecule1=atom_idx;
    
    %% translate dihedral_atom numbers into heavy atoms
    %% this is from heavy_atom_list(heavy_atom)=pdb atom number
    dihedral_heavy_atom1_molecule1=zeros(1,number_dihedrals_molecule1);
    dihedral_heavy_atom2_molecule1=zeros(1,number_dihedrals_molecule1);
    for atom_idx=1:number_dihedrals_molecule1
        for heavy_atom_idx=1:num_heavy_atoms1
            if dihedral_atom1_molecule1(atom_idx)==heavy_atom_list1(heavy_atom_idx)
                dihedral_heavy_atom1_molecule1(atom_idx)=heavy_atom_list1(heavy_atom_idx);
                dihedral_heavy_atom2_molecule1(atom_idx)=heavy_atom_idx;
            end
        end
    end
    
    %% find the dihedral atom pairs
    fid=fopen(ligand_dir+"/molecule/molecule2.pdbqt");
    frewind(fid);
    
    fgetl(fid);
    fgetl(fid);
    text=fgetl(fid);
    atom1_char_molecule2='';
    atom2_char_molecule2='';
    dihedral_atom1_molecule2=zeros(1,100);
    dihedral_atom2_molecule2=zeros(1,100);
    atom_idx=0;
    while(text(1:4)=='REMA')
        split_text=strsplit(text);
        if isnan(str2double(split_text{2}))~=1
            atom_idx=atom_idx+1;
            atom1_char_molecule2(atom_idx,1:length(split_text{6}))=split_text{6};
            atom2_char_molecule2(atom_idx,1:length(split_text{8}))=split_text{8};
            atom1_num_char_molecule2=strsplit(atom1_char_molecule2(atom_idx,:),"_");
            atom2_num_char_molecule2=strsplit(atom2_char_molecule2(atom_idx,:),"_");
            dihedral_atom1_molecule2(atom_idx)=str2double(atom1_num_char_molecule2{2});
            dihedral_atom2_molecule2(atom_idx)=str2double(atom2_num_char_molecule2{2});
        end
        text=fgetl(fid);
    end
    fclose(fid);
    number_dihedrals_molecule2=atom_idx;
    
    %% translate dihedral_atom numbers into heavy atoms
    %% this is from heavy_atom_list(heavy_atom)=pdb atom number
    dihedral_heavy_atom1_molecule2=zeros(1,number_dihedrals_molecule2);
    dihedral_heavy_atom2_molecule2=zeros(1,number_dihedrals_molecule2);
    for atom_idx=1:number_dihedrals_molecule2
        for heavy_atom_idx=1:num_heavy_atoms2
            if dihedral_atom1_molecule2(atom_idx)==heavy_atom_list2(heavy_atom_idx)
                dihedral_heavy_atom1_molecule2(atom_idx)=heavy_atom_list2(heavy_atom_idx);
                dihedral_heavy_atom2_molecule2(atom_idx)=heavy_atom_idx;
            end
        end
    end
    
    %% check to make sure that there is at least one atom to cross in each molecule
    if number_dihedrals_molecule1==0 || number_dihedrals_molecule2==0
        attempt=max_tries;
        success=0;
    end
    
end  %% pdbqt exists for both molecules


%% molecule check - shouldn't happen - was before file check 
if num_heavy_atoms1<min_heavy_atoms || num_heavy_atoms2<min_heavy_atoms
    attempt=max_tries;
    success=0;
end


while attempt<max_tries && success==0 && num_rings1+num_rings2<=9
    
    %% renumbering rings placed inside of success test due to ring limit
    
    %% renumber rings of molecule2 by adding the num_rings1 to it
    %% this will avoid ring numbering conflicts in multiple rings
    %% then at end renumber new_chm1, new_chm2 from 1 to num_rings1, num_rings2
    total_rings=num_rings1+num_rings2;
    renumbered_chm2=chm2;
    for ring=1:total_rings-num_rings1
        for i=1:chm_len2
            if chm2(i)==num2str(ring)
                renumbered_chm2(i)=num2str(ring+num_rings1);
            end
        end
    end
    chm2=renumbered_chm2;
    
    [molecule2,chm_len2,adj2,adj_heavy2,num_heavy_atoms2,heavy_atom_list2,heavy_idx_chm2,num_rings2,ring_idx_chm2,adj_atom2,avail_heavy_bond2, ...
        chiral2,num_bonds_left2,num_bonds_right2]=MoleculeStructure(chm2);
    
    %% until numbering issue is fixed
    if num_heavy_atoms2>=min_heavy_atoms && num_heavy_atoms1>=min_heavy_atoms
        
        new_chm1='';
        new_chm2='';
        
        %% segments can be either the rest of the molecule or the last part of a branch
        %% not_branch=0 or not_branch=1
        
        %% random_heavy_atom1 is start of 2nd side of molecule 1
        %% same for 2nd side of molecule 2
        
        %% check if branch or side of molecule
        %% if not_branch==-1 then no segment (could be in a ring)
        not_branch1=-1;
        not_branch2=-1;
        
        %% the 2 random heavy atoms for the segment interchange are found from these
        %%   from molecule1 and molecule2
        %% at one dihedral point the atom on either side could be chosen
        %% valences checked
        
        %% the interchange is from dihedral_atom2 from 1st,2nd molecule, which is rotatably bonded to dihedral_atom1 (single bond)
        %% no problem with valence as no bond is changed and there is already a single bond
        
        rand_atom1=ceil(number_dihedrals_molecule1*rand);
        random_heavy_atom1=dihedral_heavy_atom2_molecule1(rand_atom1);
        
        rand_atom2=ceil(number_dihedrals_molecule2*rand);
        random_heavy_atom2=dihedral_heavy_atom2_molecule2(rand_atom2);
        
        %% the rest of the function can be changed because there are no double or triple bonds to the left of the random_heavy_atom
            
        %% first molecule
        no_pars=0;
        ring_segment=zeros(9,1);
        for idx1=heavy_idx_chm1(random_heavy_atom1):chm_len1
            if chm1(idx1)=="("
                no_pars=no_pars+1;
            end
            if chm1(idx1)==")"
                no_pars=no_pars-1;
            end
            if isnan(str2double(chm1(idx1)))==0
                ring_segment(str2double(chm1(idx1)))=ring_segment(str2double(chm1(idx1)))+1;
            end
        end
        %% if completely broken into 2 segments and not a branch
        if no_pars==0 && norm(mod(ring_segment(:),2))==0
            not_branch1=1;
        end
        
        ring_segment=zeros(9,1);
        no_pars=0;
        for idx1=heavy_idx_chm1(random_heavy_atom1):chm_len1
            if chm1(idx1)=="("
                no_pars=no_pars+1;
            end
            if chm1(idx1)==")"
                no_pars=no_pars-1;
            end
            if isnan(str2double(chm1(idx1)))==0
                ring_segment(str2double(chm1(idx1)))=ring_segment(str2double(chm1(idx1)))+1;
            end
            %% branch end at idx1 - closed at )
            if no_pars==-1 && norm(mod(ring_segment(:),2))==0
                not_branch1=0;
                idx1_branch_end=idx1;
                break;
            end
        end
        
        %% second molecule
        no_pars=0;
        ring_segment=zeros(9,1);
        for idx2=heavy_idx_chm2(random_heavy_atom2):chm_len2
            if chm2(idx2)=="("
                no_pars=no_pars+1;
            end
            if chm2(idx2)==")"
                no_pars=no_pars-1;
            end
            if isnan(str2double(chm2(idx2)))==0
                ring_segment(str2double(chm2(idx2)))=ring_segment(str2double(chm2(idx2)))+1;
            end
        end
        %% if completely broken into 2 segments and not a branch
        if no_pars==0 && norm(mod(ring_segment(:),2))==0
            not_branch2=1;
        end
        
        ring_segment=zeros(9,1);
        no_pars=0;
        for idx2=heavy_idx_chm2(random_heavy_atom2):chm_len2
            if chm2(idx2)=="("
                no_pars=no_pars+1;
            end
            if chm2(idx2)==")"
                no_pars=no_pars-1;
            end
            if isnan(str2double(chm2(idx2)))==0
                ring_segment(str2double(chm2(idx2)))=ring_segment(str2double(chm2(idx2)))+1;
            end
            %% branch end at idx1 - closed at )
            if no_pars==-1 && norm(mod(ring_segment(:),2))==0
                not_branch2=0;
                idx2_branch_end=idx2;
                break;
            end
        end
        
        %% exchange the 2 segments of the molecules
        
        %% four cases: not_branch1=0,1 and not_branch2=0,1
        %% exchanging part of branch with side of molecule etc...
        
        delete_one=0;
        delete_two=0;
        
        start_one=heavy_idx_chm1(random_heavy_atom1);  %% random_heavy_atom1 on right side of exchange
        if start_one-1>=1 && (chm1(start_one-1)=="=" || chm1(start_one-1)=="#")  %% non-chiral also
            start_one=start_one-1;
            delete_one=1;  %%  the =,# is dropped after exchange
        end
        
        start_two=heavy_idx_chm2(random_heavy_atom2);  %% random_heavy_atom2 on right side of exchange
        if start_two-1>=1 && (chm2(start_two-1)=="=" || chm2(start_two-1)=="#")  %% non-chiral also
            start_two=start_two-1;
            delete_two=1;  %%  the =,# is dropped after exchange
        end
        
        if start_one-1>=1 && start_two-1>=1 && chm1(start_one-1)=="=" && chm2(start_two-1)=="="
            delete_one=0;
            delete_two=0;
        end
        
        %% initial
        if not_branch1>-1 && not_branch2>-1
            %% unchanged part - left side of expressions
            %%    if a branch and not left side then end ) and rest of molecule has to be inserted
            if start_one-1>=1
                new_chm1(1:start_one-1)=chm1(1:start_one-1);   %% all on left side of start_one
            end
            if start_two-1>=1
                new_chm2(1:start_two-1)=chm2(1:start_two-1);
            end
        end
        
        %% not_branch1==1, not_branch2==1  both right sides
        if not_branch1==1 && not_branch2==1
            new_chm1(start_one:chm_len2-start_two-delete_two+start_one)=chm2(start_two+delete_two:chm_len2);
            new_chm2(start_two:chm_len1-start_one-delete_one+start_two)=chm1(start_one+delete_one:chm_len1);
        end
        
        %% not_branch1==1, not_branch2==0   first is right side and second is part of branch
        if not_branch1==1 && not_branch2==0
            %% in new chm1 insert branch as right side
            new_chm1(start_one:idx2_branch_end-1-start_two-delete_two+start_one)=chm2(start_two+delete_two:idx2_branch_end-1);  %% not including )
            %% in new_chm2 insert right side of first into branch of second, include the )
            new_chm2(start_two:chm_len1-start_one+delete_one+start_two)=chm1(start_one-delete_one:chm_len1);
            total=chm_len1-start_one+delete_one+start_two;
            new_chm2(total+1)=")";
            total=total+1;
            new_chm2(total+1:chm_len2-idx2_branch_end-1+total+1)=chm2(idx2_branch_end+1:chm_len2);
        end
        
        %% not_branch1==0, not_branch2==1   first is part of branch and second is right side
        if not_branch1==0 && not_branch2==1
            %% in new_chm1 insert right side of second into branch of second, include the )
            new_chm1(start_one:chm_len2-start_two+delete_two+start_one)=chm2(start_two-delete_two:chm_len2);
            total=chm_len2-start_two+delete_two+start_one;
            new_chm1(total+1)=")";
            total=total+1;
            new_chm1(total+1:chm_len1-idx1_branch_end-1+total+1)=chm1(idx1_branch_end+1:chm_len1);
            %% in new chm2 insert branch as right side
            new_chm2(start_two:idx1_branch_end-1-start_one-delete_one+start_two)=chm1(start_one+delete_one:idx1_branch_end-1);  %% not including )
        end
        
        %% not_branch1==0, not_branch2==0   first and second are part of branches
        if not_branch1==0 && not_branch2==0
            %% in new_chm1 insert the branch of first into branch of second, the ) is idx1_branch_end
            new_chm1(start_one:idx2_branch_end-start_two-delete_two+start_one)=chm2(start_two+delete_two:idx2_branch_end);
            total=idx2_branch_end-start_two-delete_two+start_one;
            new_chm1(total+1:total+1+chm_len1-idx1_branch_end-1)=chm1(idx1_branch_end+1:chm_len1);
            %% in new_chm1 insert the branch of first into branch of second, the ) is idx1_branch_end
            new_chm2(start_two:idx1_branch_end-start_one-delete_one+start_two)=chm1(start_one+delete_one:idx1_branch_end);
            total=idx1_branch_end-start_one-delete_one+start_two;
            new_chm2(total+1:total+1+chm_len2-idx2_branch_end-1)=chm2(idx2_branch_end+1:chm_len2);
        end
        
        %% create pdb file
        
     %% create pdb file
        system("rm "+ligand_dir+"/molecule/*.*");
    	fileID=fopen(ligand_dir+"/molecule/molecule.smi",'w');
    	fprintf(fileID,'%s',char(new_chm1));
    	fclose(fileID);

    	system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,axchir,msi=50,msc=10,names,preserve "+ligand_dir+"/molecule/molecule.smi "+ligand_dir+"/molecule/molecule.pdb");
    	%% if ok the file test.pdb will exist
    	success1=0;
    	if exist(ligand_dir+"/molecule/molecule.001.pdb")>0
           if dir(ligand_dir+"/molecule/molecule.001.pdb").bytes>0
              success1=1;
           end
        end
  
	  %% create pdb file
    system("rm "+ligand_dir+"/molecule/*.*");
    fileID=fopen(ligand_dir+"/molecule/molecule.smi",'w');
    fprintf(fileID,'%s',char(new_chm2));
    fclose(fileID);

    system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,axchir,msi=50,msc=10,names,preserve "+ligand_dir+"/molecule/molecule.smi "+ligand_dir+"/molecule/molecule.pdb");
    %% if ok the file test.pdb will exist
    success2=0;
    if exist(ligand_dir+"/molecule/molecule.001.pdb")>0
        if dir(ligand_dir+"/molecule/molecule.001.pdb").bytes>0
            success2=1;
        end
    end 
        
        %% if both pdb files are non-empty success
        if success1>0 && success2>0
            success=1;
            %% renumber rings - chm1
            new_chm1=RING_RENUMBER_CHECK(new_chm1);
            %% renumber rings - chm2
            new_chm2=RING_RENUMBER_CHECK(new_chm2);          
            
            %% check minimum number of atoms
            [new_molecule1,new_chm_len1,new_adj1,new_adj_heavy1,new_num_heavy_atoms1,new_heavy_atom_list1,new_heavy_idx_chm1, ...
                new_num_rings1,new_ring_idx_chm1,new_adj_atom1,new_avail_heavy_bond1,new_chiral1,new_num_bonds_left1,new_num_bonds_right1]=MoleculeStructure(new_chm1);
            
            [new_molecule2,new_chm_len2,new_adj2,new_adj_heavy2,new_num_heavy_atoms2,new_heavy_atom_list2,new_heavy_idx_chm2, ...
                new_num_rings2,new_ring_idx_chm2,new_adj_atom2,new_avail_heavy_bond2,new_chiral2,new_num_bonds_left2,new_num_bonds_right2]=MoleculeStructure(new_chm2);
            
            if new_num_heavy_atoms1<min_heavy_atoms || new_num_heavy_atoms2<min_heavy_atoms
                success=0;
            end
            
            %% check clashes in both molecules
            if CLASH_CHECK(new_chm1)==0 && CLASH_CHECK(new_chm2)
                success=0;
            end
            
        end
        
    end %% num_heavy_atoms2>0
    
    attempt=attempt+1;
    
end  %% while attempt

%% attempt
%% success

if success==0
    new_chm1=chm1;
    new_chm2=chm2;
end

end


