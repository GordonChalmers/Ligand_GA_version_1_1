%% Gordon Chalmers 10/21

function xoverKids = Ligand_GA_Crossover(parents,options,nvars,FitnessFcn,thisScore,thisPopulation,state)

%% this function crosses over the parents 2 at a time to produce
%%  a new set of molecules

global atom_type;
global atom_val;
global total_atom_types;
global percent_atom;
global alphabet;
global max_tries;
global mutation_type_probability;
global min_heavy_atoms;

nKids=size(parents,2)/2;
xoverKids=cell(nKids,1);

for index=1:nKids
       
    chm1=char(thisPopulation(parents(2*index-1),:));
    chm2=char(thisPopulation(parents(2*index),:));
    
    new_chm1=chm1;
    new_chm2=chm2;
    
    %% this will try max_tries times the crossover until one is made
    %% some molecules don't have crossover capability due to absence of branches
    attempt=0;
    %% attempt=max_tries-1;  %% try 1 time at max_tries in INTERCHANGE_DIHEDRAL
    success1=0;
    success2=0;
    while success1==0 && success2==0 && attempt<max_tries
	%% INTERCHANGE_DIHEDRAL is crossover
        [new_chm1,new_chm2]=INTERCHANGE_DIHEDRAL(chm1,chm2);
        if strcmp(char(new_chm1),char(chm1))==0
            success1=1;
        end
        if strcmp(char(new_chm2),char(chm2))==0
            success2=1;
        end
        attempt=attempt+1;
    end
    
    xoverKids{index}=chm2;
    
    if success1==1 && success2==1
        choice=rand;
        if choice<=.5
            xoverKids{index}=new_chm1;
        end
        if choice>.5
            xoverKids{index}=new_chm2;
        end
%%        display("Crossover - Successful");
    end
    
end  %% index

end  %% function
