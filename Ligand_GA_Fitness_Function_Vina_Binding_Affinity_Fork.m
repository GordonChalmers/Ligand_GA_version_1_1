
%% deprecated

%% Gordon Chalmers 10/21

%% fitness function from Vina docking
%%
%%  Inputs: molecules, i.e. chm
%%  Output: fitnesses to drive the GA
%%
%% no check on clash or min_heavy_atoms as this should be done before

function FitnessFcn=Ligand_GA_Fitness_Function_Vina_Binding_Affinity(x)

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

global prepare_ligand4_path;
global Vina_path;
global config_file_path;
global corina_path;
global MoleculeDirectories;

FitnessFcn=zeros(size(x,1),1);

number_files=zeros(size(x,1),1);

%% the first part will fork the vina calculation of all the chromosomes in the population
for molecule_idx=1:size(x,1)
    
    binding_affinity=0;
    most_negative_binding_affinity=0;
    
    chm=x{molecule_idx};
    chm=char(chm);
    
    %% create pdb file
    %% each molecule_idx accesses a different directory
    rm_temp_dir_idx="rm "+"molecule"+molecule_idx+"/*.*"; 
    system(char(rm_temp_dir_idx));
    file_open="molecule"+molecule_idx+"/molecule.smi";
    fileID=fopen(file_open,'w');
    fprintf(fileID,'%s',chm);
    fclose(fileID);
    
    system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,msi="+max_stereoisomers+",msc=10,names,preserve molecule" + molecule_idx + "/molecule.smi molecule" + molecule_idx+"/molecule.pdb");
    %% if ok the file test.pdb will exist
    success=0;
    %%    if exist("molecule"+molecule_idx+"/molecule.001.pdb")>0
    dir_idx="molecule"+molecule_idx;
    dir_test=dir(dir_idx);
    if size(dir_test,1)>3
        if dir_test(3).bytes>0
            success=1;
        end
    end
    %%    end
    %% print warning if pdb wasn't created - shouldn't happen
    %%    if exist("molecule"+molecule_idx+"/molecule.001.pdb")==0
    %%        success=0;
    %%    end
    
    if success>0
        
        %% load molecule information
        [molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ...
            chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);
        
        %% structure with directory information
        directory=dir("molecule"+molecule_idx+"/*.*");
        number_files(molecule_idx)=size(directory,1)-3;
        
        if number_files(molecule_idx)>0
            
            most_negative_binding_affinity=0;
            %% calculate fitnesses of all stereoisomers
            for fileno=1:number_files(molecule_idx)
                molecule_idx
                fileno
                
                %% no need for Matlab structure - pdb files are in molecule/
                
                %% create inputs to AutoDock Vina and call it to produce the log file
                
                %% call to to prepare_ligand4.py to make molecule.fileno -> ligand.pdbqt
                %% ../mgltools/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -r ligand.pdb -o ligand.pdbqt
                %% e.g., prepare_ligand4_path="/../mgltools/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
                molecule_format_conversion_fileno=prepare_ligand4_path+" -l" + " molecule" + molecule_idx+"/"+directory(fileno+2).name + " -o molecule" + molecule_idx+"/"+directory(fileno+2).name + "_ligand.pdbqt";
                system(char(molecule_format_conversion_fileno));
                
                %% call to use AutoDock Vina on molecule.fileno
                %% Vina_path="../AutoDockVina/autodock_vina_1_1_2_linux_x86/bin/vina"
                autodock_vina=Vina_path+"vina" + " --config "+config_file_path + " --ligand molecule" + molecule_idx+"/"+directory(fileno+2).name + "_ligand.pdbqt" + " --log molecule" + molecule_idx+"/"+directory(fileno+2).name + "_log.out &";
                system(char(autodock_vina));
                
            end  %% fileno
        end  %% number_files
        
    end  %% success
end  %% molecule_idx

%% after the jobs are done - check if vina still going   
vina_completed=0;
while vina_completed==0  %% every 5 seconds check
    system("ps -aux > ps_aux_out.txt");
    system("grep vina ps_aux_out.txt > grep_out.txt"); 
    if dir("grep_out.txt").bytes==0
        vina_completed=1;
    end    
    pause(5);
end

%% next part is to parse the results for the fitnesses

display('evaluate log files')

for molecule_idx=1:size(x,1)
    for fileno=1:number_files(molecule_idx)
        
        %% parse log.out file for binding affinity out of all of the different
        if fileno<10
            log_file="molecule"+molecule_idx+"/molecule.00"+fileno+".pdb_log.out";
        end
        if fileno>=10 && fileno<100
            log_file="molecule"+molecule_idx+"/molecule.0"+fileno+".pdb_log.out";
        end
        if fileno>=100
            log_file="molecule"+molecule_idx+"/molecule."+fileno+".pdb_log.out";
        end
        
%%        fileno
%%        log_file

        fid=fopen(char(log_file));
        while ~feof(fid)
            text=fgetl(fid);
            split_text=strsplit(text);  %% cell array
            %% check for mode 1 - highest binding affinity of this stereoisomer
            %% this is from the default log file
            if size(split_text,2)==5
                if str2double(split_text(2))==1
                    binding_affinity=str2double(split_text(3));
                end
            end
        end
        fclose(fid);
        
        if binding_affinity<most_negative_binding_affinity
            most_negative_binding_affinity=binding_affinity;
            most_negative_stereoisomoer=fileno;
        end
        
    end  %% fileno
    
    most_negative_binding_affinity
    
    %%  FitnessFcn(molecule_idx) is most_negative_binding_affinity from all stereosiomers of the molecule
    FitnessFcn(molecule_idx)=most_negative_binding_affinity;
    
end  %% molecule_idx

end


