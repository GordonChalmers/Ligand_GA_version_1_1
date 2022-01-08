%% Gordon Chalmers 10/21

function mutationChildren = Ligand_GA_Mutation(parents,options,nvars, ...
    FitnessFcn, state, thisScore, thisPopulation, mutationRate)

%% point mutation operations in Ligand_GA
%%   CHANGE_ATOM, ADD_ATOM, DELETE_ATOM
%%   ADD_BRANCH, DELETE_BRANCH
%%   CLOSE_RING, OPEN_RING, OPEN_BOND (not ring point)
%%   bond changes, e.g. SINGLE_DOUBLE_BOND, ...

global atom_type;
global atom_val;
global total_atom_types;
global percent_atom;
global alphabet;
global max_tries;
global mutation_type_probability;
global min_heavy_atoms;
global corina_path;
global ligand_dir;

mutationChildren = cell(length(parents),1);
for i=1:length(parents)
    chm=thisPopulation{parents(i)};
    
    %% this function tries max_tries to make a mutation
    %% each called function also has max_tries
    %%  max total attempts of max_tries^2 in case chromosome gets stuck
    attempt=0;
    success=0;
    new_chm=chm;
    %% set attempt
    while attempt<max_tries && success==0
        attempt=attempt+1;
        new_chm=MUTATE_CHM_GA(chm);
        if strcmp(char(new_chm),char(chm))==0
%%            display("Mutation -- successful")
            success=1;
        end
    end
    mutationChildren{i}=new_chm;
end

end

