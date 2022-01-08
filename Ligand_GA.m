%% Gordon Chalmers 10/21 

%% program to find ideal molecules of given fitness function and constraints
%% 
%% input: Ligand GA Configuration file
%% output: list of molecule sequences
%%
%% molecules are in non-isomeric SMILES notation
%% 
%% external softwares - Corina, GOLD, MGLTools 

function Ligand_GA(Ligand_GA_Config_File)

%% new random seed
rng shuffle;

%% load the global variables - parameters for Ligand_GA and can be changed in Ligand_GA_Conf.m file
run(Ligand_GA_Config_File);

%% clear the directories and create /molecule tmp used 
system("rm -r "+ligand_dir+"/molecule*");
system("mkdir "+ligand_dir+"/molecule");

%% fitness function can be changed to any function
ObjectiveFunction = @(x)Ligand_GA_Fitness_Function_GOLD_ADME(x);
%% ObjectiveFunction = @(x)Ligand_GA_Fitness_Function_Atom_Density(x);
CrossoverFunction = @(x)Ligand_GA_Crossover(x); 
MutationFunction = @(x)Ligand_GA_Mutation(x); 

%% The plot function is used to show the fitness of the population
%%  as a function of the iteration.
Plot_Fitness= @(options,state,flag)Ligand_GA_Fitness_Figure(options,state,flag); 

nvars = 1;    %% Number of variables in chromosome, SMILES expression is nvars=1

%% check if the output file exists in case of overwriting 
%%  user has the option to ctrl-C and rename the output
if exist(OutFileName,'file')==2
    display('Output Population File.mat file already exists');
    display('the previous file will be deleted - when you press a key');
    display('  else, type Ctrl-C');
    pause;
    delete(OutFileName);
end

%% UB: Upper bound - none
%% LB: Lower bound - none

%% set iter_pop=0, total_pop=0 for a fresh run 
iter_pop=0;
total_pop=0;

%% this finds the population size
pop=Load_Population(InitialPopFile);
pop_size=size(pop,1);
%% create the /molecule_idx directories 
for idx=1:pop_size
    system("mkdir "+ligand_dir+"/molecule"+idx);
    system("touch "+ligand_dir+"/molecule"+idx+"/empty.txt");
end
    
%% scan over these rates to avoid partial tuning, but computationally epensive
for mutation_fraction=mutation_minus:mutation_increment:mutation_plus
    for crossover_fraction=crossover_minus:crossover_increment:crossover_plus
        
opts = optimoptions( 'ga', 'UseParallel',false, 'UseVectorized',true, ...
    'CreationFcn',@Ligand_GA_Load_Population, ...
    'InitialScoresMatrix',InitialScores, ...
    'PopulationSize',pop_size, ...
    'Generations', generations, ...
    'TolFun', tolerance_fitness, ...
    'EliteCount', ceil(.15*pop_size), ...
    'StallGenLimit',stall_gen_limit, ...
    'TolCon',tolerance_constraints, ...
    'CrossoverFcn', {@Ligand_GA_Crossover,crossover_fraction}, ...
    'MutationFcn', {@Ligand_GA_Mutation, mutation_fraction}, ...
    'OutputFcns',@Ligand_GA_Pop_Save, ...
    'Display','iter');

%% GA in Matlab Global Optimization Toolbox
ga(ObjectiveFunction, nvars, [], [], [], [], [], [],[],[], opts);

    end    
end

%% generate the unique molecules and fitnesses from a Ligand GA run

[sorted_unique_population,sorted_unique_fitness] = Ligand_GA_Output_Function(OutFileName);

end


%% the blocked sites have not been not included yet due to lack of computing resources
%%  these would be penalties in the fitness function 
%% global blocked_pdb(blocked_num);
%% global blocked_region(blocked_num,blocked_region);



