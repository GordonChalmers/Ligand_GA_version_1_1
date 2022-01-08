
function new_chm = RING_RENUMBER_CHECK(chm)

%% renumber rings to be consecutive, no missing ring numbers

%% load molecule information
[molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ...
    chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);

%% renumber rings
    new_chm=chm;
    
    %% renumber due to possible missing ring numbers
    current_ring=1;
    now_open=0;
    for ring=1:9
        open_ring=0;
        for atom=1:chm_len
            if isnan(str2double(chm(atom)))==0 && str2double(chm(atom))==ring 
                    new_chm(atom)=int2str(current_ring);
                    if mod(open_ring,2)==1
                        current_ring=current_ring+1;
                    end
                    open_ring=open_ring+1;
            end
        end
    end
  
end
