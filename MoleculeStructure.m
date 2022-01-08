%% Gordon Chalmers 10/21

function [molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm ...
    ,adj_atom,avail_heavy_bond,chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm)

%% this function returns the various items characterizing the chromosome of the molecule

%% display('MolecularStructure')
%% chm

%% default 
molecule=0;
chm_len=0;
adj=0;
adj_heavy=0;
num_heavy_atoms=0;
heavy_atom_list=0;
heavy_idx_chm=0;
num_rings=0;
ring_idx_chm=0;
adj_atom=0;
avail_heavy_bond=0;
chiral=0;
num_bonds_left=0;
num_bonds_right=0;

chm=char(chm);

global percent_atom;
global total_atom_types;
global atom_type;
global atom_val;
global max_tries;
global min_heavy_atoms;
global corina_path;
global ligand_dir;

%% variables:
%% atom_type(i) is atom type
%% atom_val(i) is valence of unbonded atom
%% percent(i) is percentage of population (e.g., abundance) of atoms
%%      - used as a 0th order approximate not to overload molecules with low or trace occurances of atom types


%% from input chromosome
%%	adj(i,j) is adjacency matrix
%% 	adj_heavy(atom) is connectivity of heavy atoms in molecule - used in chirality calculations
%%  molecule_atom_type(:) molecule heavy atom type from chm
%% 	adj_hydrogen(atom) is that of hydrogens
%% 	ring_points(:) are ring points
%%  ring_atoms(:) are atoms in rings
%%  lccations are chm_idx(ring_atom(:))

%% length of chm
chm_len=length(chm);

%% this is a common program - pause before being used again
pause(.02);

%% create pdb file
    system("rm "+ligand_dir+"/molecule/*.*");
    fileID=fopen(ligand_dir+"/molecule/molecule.smi",'w');
    fprintf(fileID,'%s',char(chm));
    fclose(fileID);

    system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,axchir,msi=50,msc=10,names,preserve "+ligand_dir+"/molecule/molecule.smi "+ligand_dir+"/molecule/molecule.pdb");
    %% if ok the file test.pdb will exist
    success=0;
    if exist(ligand_dir+"/molecule/molecule.001.pdb")>0
        if dir(ligand_dir+"/molecule/molecule.001.pdb").bytes>0
            success=1;
        end
    end



%% create pdb file

if success>0
    %% read into structure
    molecule=pdbread(ligand_dir+"/molecule/molecule.001.pdb");
    %% adj_heavy  num_heavy_atoms
    num_heavy_atoms=0;
    heavy_atom_list=zeros(1,1);
    for i=1:size(molecule.Connectivity,2)
        if molecule.Model.HeterogenAtom(i).element~='H'
            num_heavy_atoms=num_heavy_atoms+1;
            heavy_atom_list(num_heavy_atoms)=i;
        end
    end
    
    adj=zeros(num_heavy_atoms,6);  %% 4 is max, but 6 is default due to hypervalence of S,P,...
    adj_heavy=zeros(1,num_heavy_atoms);
    for heavy_atom=1:num_heavy_atoms
        for bond=1:size(molecule.Connectivity(heavy_atom_list(heavy_atom)).BondAtomList,1)
            %%		if molecule.Connectivity(heavy_atom_list(heavy_atom)).BondAtomList(bond)<=num_heavy_atoms
            if molecule.Model.HeterogenAtom(molecule.Connectivity(heavy_atom_list(heavy_atom)).BondAtomList(bond)).element~='H'
                adj_heavy(heavy_atom)=adj_heavy(heavy_atom)+1;
                %% rewrite to heavy_atom from heavy_atom_list
                for atom=1:num_heavy_atoms
                    if molecule.Connectivity(heavy_atom_list(heavy_atom)).BondAtomList(bond)==heavy_atom_list(atom)
                        adj(heavy_atom,bond)=atom;
                    end
                end
            end
        end
    end
    %% atom_type
    %%molecule_atom_type=zeros(1,1);
    for heavy_atom=1:num_heavy_atoms
        molecule_atom_type(heavy_atom)=molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom)).element;
    end
    
    %% chm location of heavy atoms
    heavy_atom=0;
    heavy_idx_chm=zeros(num_heavy_atoms,1);
    for i=1:chm_len
        for type=1:total_atom_types
            if char(chm(i))==atom_type(type)
                heavy_atom=heavy_atom+1;
                heavy_idx_chm(heavy_atom)=i;
            end
        end
    end
    %%num_heavy_atoms=heavy_atom;
    
    %% chm location of ring atoms - ring_atom(chm location, 1 2 ...  rings)=ring number at 1 2 ... rings
    %% location of ring integers and values - used in open bond
    ring_idx_chm=zeros(1,1);
    ring_num=zeros(1,1);
    num_rings=0;
    for idx=1:chm_len
        if isnan(str2double(chm(idx)))==0
            num_rings=num_rings+1;
            ring_idx_chm(num_rings)=idx;
            ring_num(num_rings)=str2double(chm(idx));
        end
    end
    num_rings=num_rings/2;
        
    %% adjacency of all atoms
    num_atoms=size(molecule.Model.HeterogenAtom,2);
    adj_atom=zeros(1,num_atoms);
    for atom=1:num_atoms
        for bond=1:size(molecule.Connectivity(atom).BondAtomList,1)
            adj_atom(atom)=adj_atom(atom)+1;
        end
    end
    
    %% CHIRALITY is not used in version2
    %% chirality of heavy atoms -
    %%      chiral=1 is @H]
    %%      chiral=2 is @@H]
    %%      chiral=3 is @]
    %%      chiral=4 is @@]
    chiral=zeros(1,num_heavy_atoms);
    %{
for atom=1:num_heavy_atoms
    if AtomCharVal(molecule.Model.HeterogenAtom(heavy_atom_list(atom)).element)==4
        if adj_heavy(atom)==3 && avail_heavy_bond(atom)==1
            if chm(heavy_idx_chm(atom)+1:heavy_idx_chm(atom)+3)=="@H]"
                chiral(atom)=1;
            end
            if chm(heavy_idx_chm(atom)+1:heavy_idx_chm(atom)+4)=="@@H]"
                chiral(atom)=2;
            end
        end
        if adj_heavy(atom)==4 && avail_heavy_bond(atom)==0
            if chm(heavy_idx_chm(atom)+1:heavy_idx_chm(atom)+2)=="@]"
                chiral(atom)=3;
            end
            if chm(heavy_idx_chm(atom)+1:heavy_idx_chm(atom)+3)=="@@]"
                chiral(atom)=4;
            end
        end
    end
end
    %}
    
    %% number of bonds to left and to right
    %% num_bonds_left, num_bonds_right in the expression
    num_bonds_left=zeros(num_heavy_atoms,1);
    num_bonds_right=zeros(num_heavy_atoms,1);
    %% num_bonds_left
    for atom=2:num_heavy_atoms
        for bond=1:6
            if adj(atom,bond)>0
                if adj(atom,bond)>0 && adj(atom,bond)<atom
                    num_bonds_left(atom)=num_bonds_left(atom)+1;
                end
            end
        end
        %% multiple bonds on lhs is always immediately before
        if heavy_idx_chm(atom)-1>=1 && chm(heavy_idx_chm(atom)-1)=="="
            num_bonds_left(atom)=num_bonds_left(atom)+1;
        end
        if heavy_idx_chm(atom)-1>=1 && chm(heavy_idx_chm(atom)-1)=="#"
            num_bonds_left(atom)=num_bonds_left(atom)+1;
        end
    end
    %% num_bonds_right
    %% potential issue if atom is double ring atom - only 1st atom in chromosome can have this  C12= if valence 4 is maximum
    for atom=1:num_heavy_atoms
        for bond=1:6
            if adj(atom,bond)>0 && adj(atom,bond)>atom
                %% count bonds of non-ring atoms
                if heavy_idx_chm(adj(atom,bond))<chm_len && isnan(str2double(chm(heavy_idx_chm(atom)+1)))~=0                    
                    num_bonds_right(atom)=num_bonds_right(atom)+1;
                    %% atom is not a ring to right
                    if heavy_idx_chm(adj(atom,bond))<chm_len && chm(heavy_idx_chm(adj(atom,bond))-1)=="=" && isnan(str2double(chm(heavy_idx_chm(atom)+1)))~=0
                        num_bonds_right(atom)=num_bonds_right(atom)+1;
                    end
                    if heavy_idx_chm(adj(atom,bond))<chm_len && chm(heavy_idx_chm(adj(atom,bond)-1))=="#" && isnan(str2double(chm(heavy_idx_chm(atom)+1)))~=0
                        num_bonds_right(atom)=num_bonds_right(atom)+2;
                    end                   
                end
                %% count bonds of ring atoms
                if heavy_idx_chm(adj(atom,bond))<chm_len && isnan(str2double(chm(heavy_idx_chm(atom)+1)))==0
                    num_bonds_right(atom)=num_bonds_right(atom)+1;                    
                    if heavy_idx_chm(adj(atom,bond))<chm_len && chm(heavy_idx_chm(adj(atom,bond))-1)=="=" && chm(heavy_idx_chm(atom)+1)~=chm(heavy_idx_chm(adj(atom,bond))+1)
                        num_bonds_right(atom)=num_bonds_right(atom)+1;
                    end
                    if heavy_idx_chm(adj(atom,bond))<chm_len && chm(heavy_idx_chm(adj(atom,bond))-1)=="#" && chm(heavy_idx_chm(atom)+1)~=chm(heavy_idx_chm(adj(atom,bond))+1)
                        num_bonds_right(atom)=num_bonds_right(atom)+2;
                    end                   
                end
            end
        end
    end
       
    %% available heavy atom bonds (instead of H)
    avail_heavy_bond=zeros(1,num_heavy_atoms);
    for atom=1:num_heavy_atoms
        avail_heavy_bond(atom)=AtomCharVal(molecule.Model.HeterogenAtom(heavy_atom_list(atom)).element)-num_bonds_left(atom)-num_bonds_right(atom);
    end
    
end  %% success>0

end


