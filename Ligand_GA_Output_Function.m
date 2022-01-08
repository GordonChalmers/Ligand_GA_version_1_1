%% Gordon Chalmers 10/21

function [sorted_unique_population,sorted_unique_fitness] = Ligand_GA_Output_Function(OutFileName)

%% parses OutFileName 
%%  unique list of non-isomeric molecules 
%%  sorted according to fitness

global parastring
global iter;
global pop_size;
global max_fitness;
global total;
global OutFileName;
global Protein_Name;

load(OutFileName);

%% GenData - Population, Fitness, Generation, total

%% there could be duplicates of the chromosomes with differing fitnesses
%% take into account

[temp_unique_population,index_unique]=unique(GenData.Population(:));
number_temp_unique_pop=size(temp_unique_population,1);

pop_size=size(GenData.Population,2);
number_fitnesses=size(Protein_Name,2);

%% find fitnesses, the most negative if repeats in GenData
temp_unique_fitness=zeros(number_temp_unique_pop,number_fitnesses);
for unique_pop_idx=1:number_temp_unique_pop
    temp_fitness=0;
    best_fitness=0;
    for pop_idx=1:pop_size
        if strcmp(GenData.Population(pop_idx),temp_unique_population(unique_pop_idx))==1
            temp_fitness=GenData.Fitness(pop_idx);
            if temp_fitness<best_fitness
                best_fitness=GenData.Fitness(pop_idx);
                temp_unique_fitness(unique_pop_idx)=GenData.Fitness(pop_idx);
            end
        end
    end
end

%% eliminate the zero fitness chromosomes  
unique_fitness=zeros(number_fitnesses,1);
unique_population="";
idx=0;
for pop_idx=1:number_temp_unique_pop
    if norm(temp_unique_fitness(pop_idx))~=0
        idx=idx+1;
        unique_fitness(idx)=temp_unique_fitness(pop_idx);
        unique_population(idx)=temp_unique_population(pop_idx);
    end
end
    
%% sort the unique chromosomes in fitness
    
number_pop=size(unique_population,2);
sorted_unique_population="";
sorted_unique_fitness=zeros(number_pop,number_fitnesses);

[empty,sorted_index_unique]=sort(unique_fitness(:),'descend');

for idx=1:number_pop
    sorted_unique_fitness(idx)=unique_fitness(sorted_index_unique(idx));
    sorted_unique_population(idx)=unique_population(sorted_index_unique(idx));
end

sorted_unique_fitness=sorted_unique_fitness;
sorted_unique_population=sorted_unique_population';

end
