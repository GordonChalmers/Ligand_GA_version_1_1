%% Gordon Chalmers 10/21

%% fitness function from GOLD PLP docking scores
%%
%%  Inputs: molecules, i.e. chm
%%  Output: fitnesses to drive the GA

function FitnessFcn=Ligand_GA_Fitness_Function_GOLD_ADME(x)

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
global set_max;
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
global max_rotatable_dihedrals;
global penalty_mass_ADME;
global penalty_acceptors_ADME;
global penalty_donors_ADME;
global penalty_rotatable_dihedrals_ADME;
global hard_violation_penalty_ADME;

global GOLD_path;

global gold_config_file;
global ligand_dir;
global OutFileName;

FitnessFcn=zeros(size(x,1),1);
number_files=zeros(size(x,1),1);

%% create the molecules at this iteration
corina_success=zeros(size(x,1),1);
for molecule_idx=1:size(x,1)
    
    %% clean the molecule directory
    system("rm -r "+ligand_dir+"/molecule"+molecule_idx+"/*");
    
    %% create mol2 file to check existence
    %%  each molecule_idx accesses a different directory
    chm=char(x{molecule_idx});
    file_open=ligand_dir+"/molecule"+molecule_idx+"/molecule.smi";
    fileID=fopen(file_open,'w');
    fprintf(fileID,'%s',chm);
    
    %% create mol2 ligand file
    system(corina_path + " -i t=smiles -o t=mol2,xlabel,split -d wh -d stergen,msi="+max_stereoisomers+",msc=10,names,preserve "+ligand_dir+"/molecule" + molecule_idx + "/molecule.smi " + ligand_dir + "/molecule" + molecule_idx+"/molecule.mol2 &");
    %%    display("corina mol2 + "+molecule_idx);
    
    %% create pdb ligand file
    system(corina_path + " -i t=smiles -o t=pdb,xlabel,split -d wh -d stergen,msi="+max_stereoisomers+",msc=10,names,preserve "+ligand_dir+"/molecule" + molecule_idx + "/molecule.smi " + ligand_dir + "/molecule" + molecule_idx+"/molecule.pdb &");
    %%    display("corina pdb + "+molecule_idx);
    
end

%% check if corina jobs are done
completed=0;
while completed==0  %% every 5 seconds check
    system("ps -aux > "+ligand_dir+"/ps_aux_out.txt");
    system("grep corina "+ligand_dir+"/ps_aux_out.txt > "+ligand_dir+"/grep_out.txt");
    if dir(ligand_dir+"/grep_out.txt").bytes==0
        completed=1;
    end
    pause(5);
end

%% corina output check here
for molecule_idx=1:size(x,1)
    %% if ok the file pdb will exist
    corina_success(molecule_idx)=0;
    dir_idx=ligand_dir+"/molecule"+molecule_idx;
    dir_test=dir(dir_idx);
    if size(dir_test,1)>3
        if dir_test(3).bytes>0
            corina_success(molecule_idx)=1;
        end
    end
end

%%  should use own code in matlab and not mgltools
%% create pdbqt files for bond analysis
for molecule_idx=1:size(x,1)
    %% total torsion angles
    %% call to to prepare_ligand4.py to make molecule.fileno -> ligand.pdbqt
    %% ../mgltools/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -r ligand.pdb -o ligand.pdbqt
    %% e.g., prepare_ligand4_path="/../mgltools/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    if corina_success(molecule_idx)>0
        molecule_format_conversion1=prepare_ligand4_path+" -l "+ligand_dir+"/molecule"+molecule_idx+"/molecule.001.pdb -o "+ligand_dir+"/molecule"+molecule_idx+"/molecule.001.pdbqt";
        system(char(molecule_format_conversion1));
        %%        display("prepare_ligand4 + "+molecule_idx);
    end
end

%% wait for jobs to complete
completed=0;
while completed==0  %% every 5 seconds check
    system("ps -aux > "+ligand_dir+"/ps_aux_out.txt");
    system("grep mgltools "+ligand_dir+"/ps_aux_out.txt > "+ligand_dir+"/grep_out.txt");
    if dir(ligand_dir+"/grep_out.txt").bytes==0
        completed=1;
    end
    pause(5);
end

%% pdbqt output check here
pdbqt_success=zeros(1,size(x,1));
for molecule_idx=1:size(x,1)
    %% if ok the file pdb will exist
    pdbqt_success(molecule_idx)=0;
    %%    if exist("molecule"+molecule_idx+"/molecule.001.pdb")>0
    dir_idx=ligand_dir+"/molecule.001.pdbqt";
    dir_test=dir(dir_idx);
    if isfile(ligand_dir+"/molecule1/molecule.001.pdbqt")==1
        pdbqt_success(molecule_idx)=1;
    end
end

%% now use CSD GOLD docking
%% -- set_max at a time to avoid memory issue
for molecule_idx_set=1:ceil(size(x,1)/set_max)
    if size(x,1)>=set_max*molecule_idx_set
        upper=set_max*(molecule_idx_set);
    end
    if size(x,1)<set_max*molecule_idx_set
        upper=size(x,1);
    end
    for molecule_idx=(set_max*(molecule_idx_set-1)+1):upper        
        if corina_success(molecule_idx)>0 && pdbqt_success(molecule_idx)>0
            %% not required -
            %% load molecule information
            %% [molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ...
            %%    chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);
            
            %% structure with directory information
            directory=dir(ligand_dir+"/molecule"+molecule_idx+"/*.*");
            number_files(molecule_idx)=floor((size(directory,1)-4)/2);
            
            if number_files(molecule_idx)>0
                %% calculate fitnesses of all stereoisomers
                for fileno=1:number_files(molecule_idx)
                    %% each stereoisomer gets its own directory
                    mkdir_stereoisomer_idx_fileno="mkdir "+ligand_dir+"/molecule"+molecule_idx+"/stereoisomer"+fileno;
                    system(mkdir_stereoisomer_idx_fileno);
                    
                    display("GOLD docking in progress: molecule "+molecule_idx+", stereoisomer "+fileno);
                    
                    %% change the gold.conf file for molecule_idx and stereoisomer
                    %% could with sed command
                    
                    if fileno<10
                        pdb_file="molecule.00"+fileno+".mol2";
                    end
                    if fileno>=10 && fileno<100
                        pdb_file="molecule.0"+fileno+".mol2";
                    end
                    if fileno>=100
                        pdb_file="molecule."+fileno+".mol2";
                    end
                    
                    %% gold config file
                    gold_config_file_read=gold_config_file;
                    gold_config_file_idx_stereoisomer=ligand_dir+"/molecule"+molecule_idx+"/stereoisomer"+fileno+"/gold.conf";
                    %% copy the previous file_read to the new one
                    scp_gold_conf="scp "+gold_config_file_read+" "+gold_config_file_idx_stereoisomer;
                    system(scp_gold_conf);
                    
                    %% then change the paths in the gold conf file
                    %%   directory, ligand_data_file, concatenated_output
                    change_directory="sed -i 's#"+ligand_dir+"/molecule1/#"+ligand_dir+"/molecule"+molecule_idx+"/stereoisomer"+fileno+"/"+"#g' "+gold_config_file_idx_stereoisomer;
                    system(change_directory);
                    
                    change_ligand_data_file="sed -i 's#molecule.001.mol2#"+pdb_file+"#g' "+gold_config_file_idx_stereoisomer;
                    system(change_ligand_data_file);
                    
                    change_output="sed -i 's#docked_molecule_file.001.mol2#docked_molecule_file"+fileno+".mol2"+"#g' "+gold_config_file_idx_stereoisomer;
                    system(change_output);
                    
                    %% copy the pdb (mol2) file to the directory
                    scp_stereoisomer_idx_fileno="scp "+ligand_dir+"/molecule"+molecule_idx+"/"+pdb_file+" "+ ...
                        ligand_dir+"/molecule"+molecule_idx+"/stereoisomer"+fileno+"/"+pdb_file;
                    system(scp_stereoisomer_idx_fileno);
                    
                    %% GOLD docking
                    %% example /home/gordon/CCDC/Discovery_2021/bin/gold_auto /home/gordon/Vector_Ligand_GA_GOLD/version2/Simeprevir/gold.conf
                    dock_GOLD=GOLD_path + " " + gold_config_file_idx_stereoisomer+" &";
                    system(dock_GOLD);
                    
                end  %% fileno
            end  %% number_files
            
        end  %% success
    end  %% molecule_idx
    
    %% after the jobs are done - check if docking programs still going
    completed=0;
    while completed==0  %% every 5 seconds check
        system("ps -aux > "+ligand_dir+"/ps_aux_out.txt");
        system("grep /gold_linux_64 "+ligand_dir+"/ps_aux_out.txt > "+ligand_dir+"/grep_out.txt");
        dir(ligand_dir+"/grep_out.text").bytes
        if dir(ligand_dir+"/grep_out.txt").bytes==0
            completed=1;
        end
        pause(3);
    end
end %% molecule_idx_set

%% next part is to parse the results for the fitnesses

display('evaluate log files')

for molecule_idx=1:size(x,1)   
    %% minimum PLP score from GOLD for each non chiral centric molecule

    if strcmp(ADME_Penalty_Type,"None")==0
        chm=char(x{molecule_idx});
        [molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ...
            chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);
    end
    
    PLP_fitness_fileno=zeros(number_files(molecule_idx),1);
    if corina_success(molecule_idx) && pdbqt_success(molecule_idx)
        chm=char(x{molecule_idx});
        for fileno=1:number_files(molecule_idx)
            log_file=ligand_dir+"/molecule"+molecule_idx+"/stereoisomer"+fileno+"/bestranking.lst";
            %% parse log.out file for the docking score
            %% gold problem with docking could produce an almost enpty bestranking.lst file 
            if exist(log_file,'file')>0 && dir(log_file).bytes>289
                %% checks molecule_idx and stereoisomer fileno
                fid=fopen(char(log_file));
                for i=1:8
                    text=fgetl(fid);
                end
                text2=strsplit(text);
                PLP_fitness=str2double(text2(2));
                fclose(fid);
                
                %% ADME penalties
                if strcmp(ADME_Penalty_Type,"None")==0
                    PLP_fitness_fileno(fileno)=PLP_fitness-ADME_Penalty(molecule_idx,chm,heavy_idx_chm,avail_heavy_bond);
                end
                
                %% Hbond restrictions
                %% PLP_fitness=PLP_fitness-Hbond_restrictions(molecule_idx,fileno);
                
            end %% exist log file
        end  %% fileno
        
        %% display
        %%        most_PLP_fitness
        %%        most_negative_stereoisomer
        
        %% FitnessFcn(molecule_idx) is most_negative_stereoisomer from the molecule
        %%      normalized per heavy atom, including penalty
        %%        FitnessFcn(molecule_idx)=-most_PLP_fitness/num_heavy_atoms*1.2;  %% 1.4 is a parameter now
        [most_PLP_fitness,stereoisomer]=max(PLP_fitness_fileno(:));
        FitnessFcn(molecule_idx)=-most_PLP_fitness;
        
        display("stereoismer: "+stereoisomer)
        display("score: "+most_PLP_fitness)
        
        %% write the most_PLP_fitness to a file
        
    end %% if files exist
end  %% molecule_idx

end


