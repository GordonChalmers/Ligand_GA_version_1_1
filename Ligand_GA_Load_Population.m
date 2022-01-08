%% Gordon Chalmers 10/21

function pop=Ligand_GA_Load_Population(nvars,FitnessFcn,options)

%% GA uses a custom population - InitialPopFile
%%  this function loads the txt file into a Matlab cell array 

global atom_type; 
global atom_val;
global total_atom_types; 
global percent_atom; 
global alphabet;
global max_tries;
global mutation_type_probability;
global pop_size;
global InitialPopFile;

%% pop_cell_array=readcell(InitialPopFile);  
%% input is single column with SMILES expressed molecules

index=0;
pop_cell_array={};
fid=fopen(InitialPopFile);
while ~feof(fid)
    index=index+1;
    pop_cell_array{index}=fgetl(fid);
end
fclose(fid);
pop_cell_array=pop_cell_array';

pop_size=size(pop_cell_array,1);
pop=pop_cell_array;

end
