%% Gordon Chalmers 10/21

function new_chm = MUTATE_CHM_GA(chm)

%% point mutation operations in Ligand_GA
%%   CHANGE_ATOM, ADD_ATOM, DELETE_ATOM
%%	 ADD_BRANCH, DELETE_BRANCH
%%   CLOSE_RING, OPEN_RING, OPEN_BOND (not ring point) - not included yet
%%   SINGLE_DOUBLE_BOND, DOUBLE_SINGLE_BOND
%%   SINGLE_TRIPLE_BOND, TRIPLE_SINGLE_BOND
%%   DOUBLE_TRIBPLE_BOND, TRIPLE_DOUBLE_BOND   

%% all of these are site non-specific and atom_type non-specific

global atom_type;
global atom_val;
global total_atom_types;
global percent_atom;
global alphabet;
global max_tries;
global mutation_type_probability;
global min_heavy_atoms;
global corina_path;

%% mutation type choice is determined by a random number
prob=rand;

%% these individual functions will return the the input if after max_tries no success in modification
if prob<=mutation_type_probability(1)
%%    display('mutation type CHANGE_ATOM')
    new_chm=CHANGE_ATOM(chm);
end
if prob > mutation_type_probability(1) && prob <=mutation_type_probability(2)
%%    display('mutation type ADD_ATOM')
    new_chm=ADD_ATOM(chm);
end
if prob>mutation_type_probability(2) && prob <=mutation_type_probability(3)
%%    display('mutation type DELETE_ATOM')
    new_chm=DELETE_ATOM(chm);
end
if prob>mutation_type_probability(3) && prob <=mutation_type_probability(4)
%%    display('mutation type ADD_BRANCH')
    new_chm=ADD_BRANCH(chm);
end
if prob>mutation_type_probability(4) && prob <=mutation_type_probability(5)
%%    display('mutation type DELETE_BRANCH')
    new_chm=DELETE_BRANCH(chm);
end
if prob>mutation_type_probability(5) && prob <=mutation_type_probability(6)
%%    display('mutation type CLOSE_RING')
    new_chm=CLOSE_RING(chm);
end
if prob>mutation_type_probability(6) && prob <=mutation_type_probability(7)
%%    display('mutation type OPEN_RING')
    new_chm=OPEN_RING(chm);
end
if prob>mutation_type_probability(7) && prob <=mutation_type_probability(8)
%%    display('mutation type SINGLE_DOUBLE_BOND')
    new_chm=SINGLE_DOUBLE_BOND(chm);
end
if prob>mutation_type_probability(8) && prob <=mutation_type_probability(9)
%%    display('mutation type DOUBLE_SINGLE_BOND')
    new_chm=DOUBLE_SINGLE_BOND(chm);
end
if prob>mutation_type_probability(9) && prob <=mutation_type_probability(10)
%%    display('mutation type SINGLE_TRIPLE_BOND')
    new_chm=SINGLE_TRIPLE_BOND(chm);
end
if prob>mutation_type_probability(10) && prob <=mutation_type_probability(11)
%%    display('mutation type TRIPLE_SINGLE_BOND')
    new_chm=TRIPLE_SINGLE_BOND(chm);
end
if prob>mutation_type_probability(11) && prob <=mutation_type_probability(12)
%%    display('mutation type DOUBLE_TRIPLE_BOND')
    new_chm=DOUBLE_TRIPLE_BOND(chm);
end
if prob>mutation_type_probability(12) && prob <=mutation_type_probability(13)
%%    display('mutation type TRIPLE_DOUBLE_BOND')
    new_chm=TRIPLE_DOUBLE_BOND(chm);
end

end





