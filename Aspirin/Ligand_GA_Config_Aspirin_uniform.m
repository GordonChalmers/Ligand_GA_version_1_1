

%% This file is the input to Ligand_GA(Config_File).

%% There are 4 different sets of variables.

%% External software paths
%% Input files and directories to the molecule run
%% Parameters to the mutations and crossover in the molecule evolution
%% GA parameters sent to the the Matlab global optimization toolbox GA in 
%%  Ligand_GA main function

%% clear from previous run in same Matlab session
clear global;

%% PATHS TO EXTERNAL SOFTWARE

%% Corina 
corina_path="/home/gordon/Vector_Ligand_GA/corina/corina";
%% CCDC GOLD
GOLD_path="/home/gordon/CCDC/Discovery_2021/bin/gold_auto";
%% MGLTools, prepare_ligand.py
prepare_ligand4_path="/home/gordon/mgltools/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py";

global corina_path;
global GOLD_path;
global prepare_ligand4_path;

%% PATH NAMES relevant to specific run of Ligand_GA

OutFileName = "Aspirin/Aspirin_no_ADME_uniform.mat";
ligand_dir = "Aspirin";  
InitialPopFile = "Aspirin/Aspirin_Pop_100.txt";
gold_config_file = "Aspirin/gold.conf";

global OutFileName;
global ligand_dir;
global InitialPopFile;
global gold_config_file;


%% MOLECULAR PARAMETERS 

%% atom_types used - single characters only, not dependent on order 
atom_type='CNOPSF'; 
total_atom_types=size(atom_type,2);

%% the associated variables depend on the order in atom_type

%% atom valences 
atom_val=[4,3,2,3,2,1];

%% atomic mass in Daltons 
atomic_daltons=[12.011,14.007,15.999,30.974,32.065,18.998];

%% percent_atom is used to weight the selection of atoms in 
%%  CHANGE_ATOM OR ADD_BRANCH
%% this could be based on natural abundance of elements in organic molecules
%% could be from molecular goal in Ligand GA search

%% uniform probability:
%% unnorm_percent_atom=[1,1,1,1,1,1];

%% no chance of P or S:
unnorm_percent_atom=[1,1,1,1,1,1];
%% normalize it
percent_atom=unnorm_percent_atom/sum(unnorm_percent_atom);

%% The mutation and crossover procedure requires random choices of locations 
%%  in the molecule, choices of atoms, crossover points, etc...
%% If no success after max_tries then no change in the molecule mutation or 
%%   crossover at that iteratio
max_tries=20;

%% The different mutation and crossover functions called from Ligand_GA_Crossover.m 
%%  and Ligand_GA_Mutation.m are called with normalized probabilities [0,1].  
%% 	The unnormalized probabilities, i.e. weights, are set by an array.

%% point mutation operations in Ligand_GA 
%% CHANGE_ATOM(chm)
%% ADD_ATOM(chm)
%% DELETE_ATOM(chm);
%% ADD_BRANCH(chm);
%% DELETE_BRANCH(chm);
%% CLOSE_RING(chm);
%% OPEN_RING(chm);
%% OPEN_BOND(chm);  not included yet

%% SINGLE_DOUBLE_BOND(chm);
%% DOUBLE_SINGLE_BOND(chm);
%% SINGLE_TRIPLE_BOND(chm);
%% TRIPLE_SINGLE_BOND(chm);
%% DOUBLE_TRIPLE_BOND(chm);
%% TRIPLE_DOUBLE_BOND(chm);

%% mutation probabilities of mutation type
clear unnorm_mut_probability;
%% uniform probability - not good, e.g. most atoms in molecule can not change to a triple bond
%% unnorm_mut_probability=[1,1,1,1,1,1,1,1,1,1,1,1,1];

%% OPEN_BOND is not included in this version
%% Tuned for keeping rings and less triple bond chance of mutations: 
unnorm_mut_probability=[1,1,1,1,1,1,1,1,1,1,1,1,1]; 
num_mutations=size(unnorm_mut_probability,2);
%% normalize these 
clear norm_mut_probabiltiy;
norm_mut_probability=zeros(num_mutations,1);
norm_mut_probability(:)=unnorm_mut_probability(:)/sum(unnorm_mut_probability);
%% then break it up into slices for use in Ligand_GA_Mutation.m
clear mutation_type_probability;
mutation_type_probability=zeros(num_mutations,1);
total_mutation_type_probability=0;
mutation_type_probability(1)=norm_mut_probability(1);
%% mutation_type_probability select the [] slots of random numbers to select mutations
total_mutation_type_probability=total_mutation_type_probability+norm_mut_probability(1);
for type=2:num_mutations
    mutation_type_probability(type)=total_mutation_type_probability+norm_mut_probability(type);
    total_mutation_type_probability=total_mutation_type_probability+norm_mut_probability(type);
end

%% No short molecules
%% Minimum wanted atoms
min_heavy_atoms=10;

%% New run - tells OutFileName that this is the 1st generation 
iter_pop=0;

%% Bond length test (in Angstroms) of atom types - used for clash test
%% Percentage used to check for a non CONECT bonded pair of atoms to not have a bond in terms of bond distance
%%  maximum distance of nearby non-bonded atoms can't be less than cutoff_bond_percentage*inter_bond_distance
cutoff_bond_percentage=1.10;

%% single bond lengths - taken from 2-heavy_atom molecules
%% atom_type='CNOPSF' and the 2-d array is in the order of atom_type
inter_bond_distance=zeros(total_atom_types,total_atom_types);

inter_bond_distance(1,1)=1.53;  %% from CC
inter_bond_distance(1,2)=1.49;  %% from CN
inter_bond_distance(1,3)=1.43;  %% from CO
inter_bond_distance(1,4)=1.82;  %% from CP
inter_bond_distance(1,5)=1.81;  %% from CS
inter_bond_distance(1,6)=1.40;  %% from CF

inter_bond_distance(2,2)=1.45;  %% from NN  
inter_bond_distance(2,3)=1.46;  %% from NO
inter_bond_distance(2,4)=1.68;  %% from NP
inter_bond_distance(2,5)=1.77;  %% from NS
inter_bond_distance(2,6)=1.30;  %% from NF

inter_bond_distance(3,3)=1.47;  %% from OO
inter_bond_distance(3,4)=1.61;  %% from OP
inter_bond_distance(3,5)=1.52;  %% from OS
inter_bond_distance(3,6)=1.21;  %% from OF

inter_bond_distance(4,4)=2.18;  %% from PP
inter_bond_distance(4,5)=2.12;  %% from PS 
inter_bond_distance(4,6)=1.67;  %% from PF 

inter_bond_distance(5,5)=2.05;  %% from SS
inter_bond_distance(5,6)=1.61;  %% from SF 

for type1=1:total_atom_types
    for type2=type1+1:total_atom_types
        inter_bond_distance(type2,type1)=inter_bond_distance(type1,type2);
    end
end

%% This parameter is the maximal number of stereoisomers generated by Corina 
%% 64 is good but is set to 512 - this number does exponentiate 2^{N_{chiral}}
%%  chiral centers of C and cis-/trans- of double bonds
max_stereoisomers=512;

%% number of moleceules sent at once to the cluster
set_max=1;

%% ADME restriction - 
%%  choices are: SoftLipinski, HardLipinski, None (anything but the previous 2)
%% penalty function for ADME restriction 
ADME_Penalty_Type="None";
global ADME_Penalty_Type;

%% ADME limits and penalties
max_mass_ADME=500;  global max_mass_ADME;
max_acceptors_ADME=10;  global max_acceptors_ADME;
max_donors_ADME=5;  global max_donors_ADME;
max_rotatable_dihedrals_ADME=10;  global max_rotatable_dihedrals_ADME;
penalty_mass_ADME=4;  global penalty_mass_ADME;
penalty_acceptors_ADME=10;  global penalty_acceptors_ADME;
penalty_donors_ADME=10;  global penalty_donors_ADME;
penalty_rotatable_dihedrals_ADME=10;  global penalty_rotatable_dihedrals_ADME;
hard_violation_penalty_ADME=20;  global hard_violation_penalty_ADME;


%% GA PARAMETERS
%% Genetic algorithm parameters used by the Matlab Global Optimization Toolbox GA algorithm
%%  defaults

%% generations is large to due to an indefinite run
generations=10000;  global generations; 
%% elite number is important, 1/20 to 1/1
elite_fraction=.1;  global elite_fraction;

%% scan over rates 
mutation_minus=.9;  global mutation_minus;
mutation_increment=.1;  global mutation_increment;
mutation_plus=.9;  global mutation_plus;

crossover_minus=.9;  global crossover_minus;
crossover_increment=.1;  global crossover_increment;
crossover_plus=.9;  global crossover_plus;

%% these 3 - useful to end GA search generally, but not in this and are 
%%  defined, but have no real effect due to fitness function landscape (high dimensional)
tolerance_fitness=.02;  global tolerance_fitness;
stall_gen_limit=100;  global stall_gen_limit;
tolerance_constraints=.01;  global tolerance_constraints;  

%% Make the molecular parameters global in the functions 
global atom_type; 
global atom_val;
global total_atom_types; 
global percent_atom; 
global alphabet;
global max_tries;
global mutation_type_probability;
global min_heavy_atoms;
global iter_pop;
global total_pop;
global pop_size;
global inter_bond_distance;
global cutoff_bond_percentage;
global corina_path;
global max_stereoisomers;
global set_max;
global pop_size;
global atomic_daltons;

