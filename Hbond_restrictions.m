%% Gordon Chalmers 10/31

%% Hbond restrictions
%%  this checks the highest scoring pose only 
%%  BA3_BO2_1KYJ is hard-coded and the penalty for O3 or O4 of the 3 glycans

function [total_hbond_penalty,total_hbonds]=Hbond_restrictions(molecule_idx,fileno)

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

total_hbond_penalty=0;
hbond_penalty=-20;

%% parse log.out file
log_file=ligand_dir+"/molecule"+molecule_idx+"/stereoisomer"+fileno+"/docked_molecule_file"+fileno+".mol2";

%% check for Hbonds
%%  the hydrogen bonds are in the output gold file
fid=fopen(char(log_file));
total_hbonds=0;  %% this depends on stereoisomer - hydrogen bonds
while ~feof(fid)
    text=fgetl(fid);
    if text=="> <Gold.Chemscore.Hbonds>"
        text=fgetl(fid);
        text=fgetl(fid);
        total_donors=0;
        while strcmp(text,"> <Gold.PLP.Chemscore.Internal.Correction>")==0
                if strcmp(text,"")==1 
                    break;
                end
                hbond_test=strsplit(text);
                if strcmp(hbond_test(2),"P1")==1 && (str2double(hbond_test(4))==74 || str2double(hbond_test(4))==76)
                    total_hbond_penalty=total_hbond_penalty+hbond_penalty;
                end
                if strcmp(hbond_test(2),"P1")==1 && (str2double(hbond_test(4))==99 || str2double(hbond_test(4))==101)
                    total_hbond_penalty=total_hbond_penalty+hbond_penalty;
                end
                if strcmp(hbond_test(2),"P1")==1 && (str2double(hbond_test(4))==124 || str2double(hbond_test(4))==126)
                    total_hbond_penalty=total_hbond_penalty+hbond_penalty;
                end
                
                if strcmp(hbond_test(6),"P1")==1 && (str2double(hbond_test(7))==74 || str2double(hbond_test(7))==76)
                    total_hbond_penalty=total_hbond_penalty+hbond_penalty;
                end
                if strcmp(hbond_test(6),"P1")==1 && (str2double(hbond_test(7))==99 || str2double(hbond_test(7))==101)
                    total_hbond_penalty=total_hbond_penalty+hbond_penalty;
                end
                if strcmp(hbond_test(6),"P1")==1 && (str2double(hbond_test(7))==124 || str2double(hbond_test(7))==126)
                    total_hbond_penalty=total_hbond_penalty+hbond_penalty;
                end
                total_hbonds=total_hbonds+1;
                text=fgetl(fid);
        end
        break;
    end
end
fclose(fid);
    
end
