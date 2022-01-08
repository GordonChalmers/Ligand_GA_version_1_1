%% Gordon Chalmers 10/21

%% ADME Penalty
%%  Soft Lipinski
%%  Hard Lipinski
%%  None

function penalty=ADME_Penalty(molecule_idx,chm,heavy_idx_chm,avail_heavy_bond)

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
global max_heavy_atoms;
global max_stereoisomers;
global pop_size;
global atomic_daltons;

global prepare_ligand4_path;
global Vina_path;
global config_file_path;
global corina_path;
global MoleculeDirectories;

global ADME_Penalty_Type;
global max_mass_ADME;
global max_acceptors_ADME;
global max_donors_ADME;
global max_rotatable_dihedrals_ADME;
global penalty_mass_ADME;
global penalty_acceptors_ADME;
global penalty_donors_ADME;
global penalty_rotatable_dihedrals_ADME;
global hard_violation_penalty_ADME;

global GOLD_path;

global gold_config_file;
global ligand_dir;
global OutFileName;

penalty=0;

%% molecular mass in amu of heavy atoms
total_daltons=0;
for chm_idx=1:length(chm)
    for atom_type_idx=1:total_atom_types
        if strcmp(chm(chm_idx),atom_type(atom_type_idx))==1
            total_daltons=total_daltons+atomic_daltons(atom_type_idx);
        end
    end
end

%% hydrogen bond acceptors (all N or O)
total_acceptors=0;
for chm_idx=1:length(chm)
    if strcmp(chm(chm_idx),'N')==1 || strcmp(chm(chm_idx),'O')==1
        total_acceptors=total_acceptors+1;
    end
end

%% hydrogen bond donors (available NH or OH)
total_donors=0;
for atom=1:size(heavy_idx_chm,1)
    if (strcmp(chm(heavy_idx_chm(atom)),'N')==1 || strcmp(chm(heavy_idx_chm(atom)),'O')==1) && avail_heavy_bond(atom)>0
        total_donors=total_donors+1;
    end
end

%% find the dihedral atom pairs
fid=fopen(ligand_dir+"/molecule"+molecule_idx+"/molecule.001.pdbqt");
frewind(fid);

fgetl(fid);
fgetl(fid);
text=fgetl(fid);
total_rotatable_dihedrals=0;
while(text(1:4)=='REMA')
    split_text=strsplit(text);
    if isnan(str2double(split_text{2}))~=1  %% not an inactive dihedral bond
        total_rotatable_dihedrals=total_rotatable_dihedrals+1;
    end
    text=fgetl(fid);
end
fclose(fid);

if strcmp(ADME_Penalty_Type,'SoftLipinski')==1
    
%%    display("SoftLipinski Penalty")
    
    penalty_Lipinski_mass=0;
    penalty_Lipinski_acceptors=0;
    penalty_Lipinski_donors=0;
    penalty_Lipinski_rotatable_dihedrals=0;
    
    %% Lipinski molecular weight
    if total_daltons>max_mass_ADME
        penalty_Lipinski_mass=penalty_mass_ADME*(total_daltons-max_mass_ADME);
    end
   
    %% Lipinski hydrogen bond acceptors (N or O) 
    if total_acceptors>max_acceptors_ADME
        penalty_Lipinski_acceptors=penalty_acceptors_ADME*(total_acceptors-max_acceptors_ADME);
    end
    
    %% Lipinski hydrogen bond donors (N or O) 
    if total_donors>max_donors_ADME
        penalty_Lipinski_donors=penalty_donors_ADME*(total_donors-max_donors_ADME);
    end
    
    %% max dihedrals penalty
    if total_rotatable_dihedrals>max_rotatable_dihedrals_ADME
        penalty_Lipinski_rotatable_dihedrals=penalty_rotatable_dihedrals_ADME*(total_rotatable_dihedrals-max_rotatable_dihedrals_ADME);
    end
    
    penalty=penalty_Lipinski_mass+penalty_Lipinski_donors+penalty_Lipinski_acceptors+penalty_Lipinski_rotatable_dihedrals;
 
    %% total_daltons
    %% total_acceptors
    %% total_donors
    
    %% penalty_Lipinski_mass
    %% penalty_Lipinski_acceptors
    %% penalty_Lipinski_donors
    %% penalty_Lipinski_rotatable_dihedrals
    
end  %% soft Lipinksi penalty

if strcmp(ADME_Penalty_Type,'HardLipinski')==1
    
    violations=0;
    
    %% Lipinski molecular weight
    if total_daltons>max_mass_ADME
        violations=violations+1;
    end
    
    %% Lipinski hydrogen bond acceptors (N or O)
    if total_acceptors>max_acceptors_ADME
        violations=violations+1;
    end
    
    %% Lipinski hydrogen bond donors (N or O) 
    if total_donors>max_donors_ADME
        violations=violations+1;
    end
    
    if total_rotatable_dihedrals>max_rotatable_dihedrals_ADME
        violations=violations+1;
    end
    
    penalty=hard_violation_penalty_ADME*violations;
    
end  %% hard Lipinksi penalty

end
