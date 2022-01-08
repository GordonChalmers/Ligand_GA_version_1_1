
%% test is a vector of GenData structures 
%%   e.g., test(1).GenData.Fitness

test=[file1.mat file2.mat ... ];

index=0

for i=1:length(test)
    for individual=1:size(test(i).GenData.Fitness)
	index=index+1;   
	total.Fitness(index)=test(i).GenData.Fitness(individual);
 	total.Population(index)=test(i).GenData.Population(individual);
    end
end


