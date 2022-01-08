%% Gordon Chalmers 10/21

%% saves Matlab file of population and fitnesses at each iteration

function [state, options, optchanged] = Ligand_GA_Pop_Save(options, state, flag)

global parastring
global iter_pop;
global pop_size;
global max_fitness;
global total_pop;
global OutFileName;

if (exist(OutFileName)) == 2
    load(OutFileName);
    iter_pop=iter_pop+1;
else
    GenData.Generation=1;
    iter_pop = 1;
end

for i=1:pop_size
        total_pop=total_pop+1;
        GenData.Population(total_pop) = state.Population(i);
        GenData.Fitness(total_pop,:)=state.Score(i,:);
        GenData.total=total_pop;
end

save(OutFileName, 'GenData');
% standard template code follows
optchanged = false;
switch flag
    case 'init'
        disp('Starting the algorithm');
    case {'iter','interrupt'}
        disp('Iterating ...')
    case 'done'
        disp('Performing final task');
end    % Ligand_GA_Pop_Save()
