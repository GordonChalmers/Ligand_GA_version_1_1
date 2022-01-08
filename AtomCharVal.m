%% Gordon Chalmers 10/21

function val=Atom_Char_Val(atom)

%% function that gives the numerical valence from the atom types

global atom_type;
global atom_val;
global total_atom_types; 

	for idx_atom=1:total_atom_types
		if atom==atom_type(idx_atom)
			val=atom_val(idx_atom);
		end
	end
end
