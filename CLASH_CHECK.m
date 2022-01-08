%% Gordon Chalmers 10/21

function no_clash = CLASH_CHECK(chm)

%% CLASH_CHECK checks single bond lengths against distances in the pdb file
%%  pairs of heavy atoms not bonded in the adjacency matrix are considered and have to be greater than a single bond length*cutoff_bond_percentage
%%  success=1 means that there are no conflicts

global atom_type;
global atom_val;
global total_atom_types;
global percent_atom;
global alphabet;
global max_tries;
global mutation_type_probability;
global min_heavy_atoms;
global cutoff_bond_percentage;
global inter_bond_distance;

no_clash=0;
if length(chm)>0  
    
    %% check if molecule can be constructed from chm in the functions
    %%  this is checked before calling
    molecule_exist_smiles=1;
    
    if molecule_exist_smiles>0
        %% load molecule information
        [molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ...
            chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);
        
        %% check the bond lengths
        no_clash=1;
        for heavy_atom1=1:num_heavy_atoms
            for heavy_atom2=heavy_atom1+1:num_heavy_atoms
                
                %% check if not in adjacency matrix  adj
                in_adj=0;  %% not in adj
                for adj1=1:num_heavy_atoms
                    for adj2=1:6
                        if heavy_atom2==adj(heavy_atom1,adj2)
                            in_adj=1;  %% in adj
                        end
                    end
                end
                if in_adj==0
                    %% pair distance
                    inter_pair_distance=0;
                    inter_pair_distance=(molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom2)).X-molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom1)).X)^2 ...
                        + (molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom2)).Y-molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom1)).Y)^2 ...
                        + (molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom2)).Z-molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom1)).Z)^2;
                    inter_pair_distance=sqrt(inter_pair_distance);
                    
                    %% atom_type in int from list of atom_type
                    atom_int_type1=0;
                    atom_int_type2=0;
                    for atom_type_idx=1:total_atom_types
                        if molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom1)).element==atom_type(atom_type_idx)
                            atom_int_type1=atom_type_idx;
                        end
                        if molecule.Model.HeterogenAtom(heavy_atom_list(heavy_atom2)).element==atom_type(atom_type_idx)
                            atom_int_type2=atom_type_idx;
                        end
                    end
                    
                    %% test single bond length for minimum separation of the 2 unbonded heavy atoms
                    if inter_pair_distance<cutoff_bond_percentage*inter_bond_distance(atom_int_type1,atom_int_type2)
                        no_clash=0;
                    end
                end  %% in_adj
                
            end  %% heavy_atom1
        end  %% heavy_atom2 
        
        for atom=1:num_heavy_atoms 
            if avail_heavy_bond(atom)<0 && chm(heavy_idx_chm(atom))~='S' && chm(heavy_idx_chm(atom))~='P'
                no_clash=0;
            end
            for idx_chm=1:chm_len-1 
                if chm(idx_chm:idx_chm+1)=="=#" || chm(idx_chm:idx_chm+1)=="#=" 
                    no_clash=0;
                end
            end
        end
        
    end  %% molecule_exist_smiles
end  %% null chm check

end







