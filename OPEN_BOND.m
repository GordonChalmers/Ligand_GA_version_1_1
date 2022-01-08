%% Gordon Chalmers 10/21

function [new_chm,success]= OPEN_BOND(chm)

%% OPEN_BOND at arbitrary bond in set of ring atoms

%% this function opens a bond (unlink) within a ring and creates the new expression

%% randomly select bond - check if it is in a ring and not in a pure branch that would
%%   not break the molecule into 2 pieces
%% find the maximum ring number of an open ring containing the bond
%% break chm into 4 parts
%% the components are treated differently in reassembling the expression
%% renumber the rings in the final new_chm

%% this algorithm won't fail if both ring points are in the same branch
%%  but will generate chm with ring points in different branches

global atom_type; 
global atom_val;
global total_atom_types; 
global percent_atom; 
global alphabet;
global max_tries;
global mutation_type_probability;
global cutoff_bond_percentage;
global inter_bond_distance;
global min_heavy_atoms;
global corina_path;

attempt=0;
success=0;

chm=char(chm);
chm=RING_RENUMBER_CHECK(chm);

%% load molecule information
[molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ...
    chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);

if num_heavy_atoms<min_heavy_atoms 
    attempt=max_tries;
    success=0;
end

while attempt<max_tries && success==0
      
    new_chm='';

    first_atom_bond=ceil(rand*(num_heavy_atoms-1));
    atom_bond_idx=heavy_idx_chm(random_atom_bond);
    %% bond is broken to the r.h.s. of the atom and last heavy atom is not allowed - random_atom_bond is limited to num_heavy_atoms-1
    %% atom to which random_atom_bond is bond broken is
    second_atom_bond=max(adj(first_atom_bond,:));
    second_atom_bond_idx=heavy_idx_chm(second_atom_bond);
    
    %% find the chm_idx of max open ring idx, closure, and ring_num
    %% maximum open ring at atom_bond_idx
    ring_num=zeros(num_rings,1);
    ring_num_idx=zeros(num_rings,1);
    for ring=1:num_rings
        for idx=1:atom_bond_idx
            if isnan(str2double(chm(idx)))==0
                ring_num(ring)=ring_num(ring)+1;
                ring_num_idx(ring)=idx;
            end
        end
    end
    %% ring_num(ring)==2 before atom_bond_idx then bond is not in that ring
    for ring=1:num_rings
        if ring_num(ring)==1 && atom_bond_idx>ring_num_idx(ring)
            max_open_ring=ring;   %% ring number to be removed
            max_open_ring_idx=ring_num_idx(ring);  %% idx of ring number
        end
    end
    %% max_open_ring_closure_idx
    for idx=1:chm_len
        if isnan(str2dbouble(chm(idx)))==0 && str2double(chm(idx))==max_open_ring
            max_open_ring_closure_idx=idx;
        end
    end
    %% heavy atom before first ring point
    for atom=1:num_heavy_atoms
        if heavy_atom_idx(atom)<max_open_ring_idx
            first_heavy_atom=atom;
            first_heavy_atom_idx=heavy_idx_chm(atom);
        end
    end
    %% heavy atom before closing ring point
    for atom=1:num_heavy_atoms
        if heavy_atom_idx(atom)<max_open_ring_closure_idx
            second_heavy_atom=atom;
            second_heavy_atom_idx=heavy_idx_chm(atom);
        end
    end
    
    %% find the components
    
    %% the first atom_bond atom could be a ring point or it couldn't
    %% the second atom_bond atom could be a ring point or it couldn't
    
    %% first non ring point atom - the ring point atom will keep its chirality as the atom isn't selected
    
    %% find X1minus, X2plus(1)+X1plus+X2plus(all) - transposing X2plus(1) and X2plus(all) in chirality and X2plus(all) in order
    
    %% both not ring points
    if random_atom_bond~=first_heavy_atom
        %% multiple ring points
        
        %% X1-
        %% X1-  no chirality change as no bonds deleted - rearranged in X2 transpose
        %% add any other ring points then one deleted
        while isnan(str2double(chm(total_X1minus+1)))~=0 && str2double(chm(total_X1minus+1))~=max_open_ring
            X1minus(total_X1minus+1)=chm(total_X1minus+1);
            total_X1minus=total_X1minus+1;
        end
        
        
        %% X1+
        %% next is X1plus - coming from chm after the second heavy atom and ring numbers, max_open_ring
        %% no transpose or chirality change
        total_X1plus=max_open_ring_closure_idx;
        while isnan(str2double(chm(total_X1plus+1)))~=0  %% find 1st non-numeric character - after ring points of second heavy atom
            total_X1plus=total_X1plus+1;
        end
        %% including branch notation
        X1plus(1:chm_len-total_X1plus+3)="("+chm(total_X1plus:chm_len)+")";
        
        
        %% X2-
        %% coming from first_heavy_atom to first in broken bond, atom_bond_idx  (random_atom_bond)
        %%  first heavy atom doesn't lose bonds
        %%  first_atom_bond does
        %%  C5CCCCC with  C5CCCC C broken
        %%  C5CCCC(O)C with C5CCCC(O) C broken
        %% take into account chirality and rings
        %% one bond is broken (single possibly double/triple)
        %%  shouldn't be ring point at first atom of broken bond - not done
        if chm(first_atom_bond_idx+1)~="@"  %% non-chiral initially and non-chiral after
            %% any double,triple bonds are dropped
            X2minus(1:first_atom_bond_idx)=chm(total_X1minus+1:first_atom_bond_idx);
            total_X2minus=first_atom_bond_idx;
            while isnan(str2double(chm(total_X2minus+1)))~=0  %% and ring numbers following 2nd atom
                X2minus(total_X2minus+1)=str2double(chm(total_X2minus+1);
                total_X2minus=total_X2minus+1;
            end
        end
        
        %% following can happen with C(X) X, C(X)(X) X  broken bonds for example
        if chm(first_atom_bond_idx+1:first_atom_bond_idx+3)~="@H]"
            X2minus(1:first_atom_bond_idx-2-total_X1minus)=chm(total_X1minus+1:first_atom_bond_idx-2);
            total_X2minus=first_atom_bond_idx-2-total_X1minus;
            X2minus(total_X2minus+1)=chm(first_bond_atom_idx);  %% lost chirality
            total_X2minus=total_X2minus+1;
        end
        if chm(first_atom_bond_idx+1:first_atom_bond_idx+4)~="@@H]"
            X2minus(1:first_atom_bond_idx-2-total_X1minus)=chm(total_X1minus+1:first_atom_bond_idx-2);
            total_X2minus=first_atom_bond_idx-2-total_X1minus;
            X2minus(total_X2minus+1)=chm(first_bond_atom_idx);  %% lost chirality
            total_X2minus=total_X2minus+1;
        end
        if chm(first_atom_bond_idx+1:first_atom_bond_idx+2)~="@]"
            X2minus(1:first_atom_bond_idx-2-total_X1minus)=chm(total_X1minus+1:first_atom_bond_idx-2);
            total_X2minus=first_atom_bond_idx-2-total_X1minus;
            X2minus(total_X2minus+1:total_X2minus+5)="["+chm(first_bond_atom_idx)+"@H]";   %% kept chirality
            total_X2minus=total_X2minus+5;
        end
        if chm(first_atom_bond_idx+1:first_atom_bond_idx+3)~="@@]"
            X2minus(1:first_atom_bond_idx-2-total_X1minus)=chm(total_X1minus+1:first_atom_bond_idx-2);
            total_X2minus=first_atom_bond_idx-2-total_X1minus;
            X2minus(total_X2minus+1:total_X2minus+6)="["+chm(first_bond_atom_idx)+"@@H]";   %% kept chirality
            total_X2minus=total_X2minus+6;
        end
        
        
        %% X2+
        %% add 1st atom of X2+' (transpose) - no bond change
        %% second_heavy_atom - no chirality change and all chirality cases including multiple ring points
        %% no first ring point added - max_open_ring
%%        if chm(second_heavy_atom_idx+1)~="@"   %% no chirality
%%            X2plus(1)=chm(second_heavy_atom_idx);  %% no ring number
%%            total_X1plus=1;
%%        end
%%        if chm(second_heavy_atom_idx+1:second_heavy_atom_idx+3)~="@H]"
%%            X2plus(1:6)="["+chm(second_heavy_atom_idx)+"@@H]";
%%            total_X1plus=6;
%%        end
%%        if chm(second_heavy_atom_idx+1:second_heavy_atom_idx+4)~="@@H]"
%%            X2plus(1:5)="["+chm(second_heavy_atom_idx)+"@H]";
%%            total_X1plus=5;
%%        end
%%        if chm(second_heavy_atom_idx+1:second_heavy_atom_idx+2)~="@]"
%%            X2plus(1:5)="["+chm(second_heavy_atom_idx)+"@@]";
%%            total_X1plus=5;
%%        end
%%        if chm(second_heavy_atom_idx+1:second_heavy_atom_idx+4)~="@@]"
%%            X2plus(1:4)="["+chm(second_heavy_atom_idx)+"@]";
%%            total_X1plus=4;
%%        end
        %% add any additional ring points
        while isnan(str2double(chm(total_X1minus+1)))~=0 && str2double(chm(total_X1minus+1))~=max_open_ring
            total_X2plus=total_X2plus+1;
            X2plus(total_X2plus+1)=chm(total_X2plus);
        end
        
        %% 2nd part of X2plus 
        %% starts from second_heavy_atom-1 and goes until second_bond_atom 
        %%  e.g., C(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5 -> [C@H]1[C@H]3C=C[C@@H](O)[C@H]2OC4=C(O) 
        %% linking in chain is from adj min bond numbers - can be used to include the branches 
        
%{        
        transpose_atom=second_heavy_atom-1;
        while transpose_atom>=second_bond_atom
            %% link in reverse transpose order next is (and bypasses branches)
            next_atom=-max(adj(transpose_atom,:));
            subtract_t=0;
            subtract_n=0;
            if chm(transpose_atom-1)=="["
                subtract_t=1;
            end
            if chm(next_atom-1)=="[" 
                subtract_n=1; 
            end 
            X2plusSecond(total_X2plusSecond+1:total_X2plusSecond+heavy_idx_chm(next_atom)-subtract_t-heavy_idx_chm(transpose_atom)+subtract_n+1)=chm(heavy_idx_chm(transpose_atom)-subtract_t:heavy_idx_chm(next_atom)-subtract_n-1);
            total_X2plusSecond=total_X2plusSecond+heavy_idx_chm(next_atom)-subtract_t-heavy_idx_chm(transpose_atom)+subtract_n+1;
            
            %% preceeding atom 
            min=1000; 
            for bond=1:6
                if adj(transpose_atom,bond)<min && adj(transpose_atom,bond)~=0 
                    min=adj(transpose_atom,bond);
                end
            end
            transpose_atom=min;
        end
%}      
        
        %% transpose and chirality change
        
        transpose_atom=second_heavy_atom-1;
        while transpose_atom>=second_bond_atom
            %% link in reverse transpose order next is (and bypasses branches)
            next_atom=-max(adj(transpose_atom,:));
            subtract_t=0;
            subtract_n=0;
            if chm(transpose_atom-1)=="["
                subtract_t=1;
            end
            if chm(next_atom-1)=="[" 
                subtract_n=1; 
            end 
            
            X2plusSecond(total_X2plusSecond+1:total_X2plusSecond+heavy_idx_chm(next_atom)-subtract_t-heavy_idx_chm(transpose_atom)+subtract_n+1)=chm(heavy_idx_chm(transpose_atom)-subtract_t:heavy_idx_chm(next_atom)-subtract_n-1);
            total_X2plusSecond=total_X2plusSecond+heavy_idx_chm(next_atom)-subtract_t-heavy_idx_chm(transpose_atom)+subtract_n+1;
            
            %% preceeding atom 
            min=1000; 
            for bond=1:6
                if adj(transpose_atom,bond)<min && adj(transpose_atom,bond)~=0 
                    min=adj(transpose_atom,bond);
                end
            end
            transpose_atom=min;
        end
        
        
            
        %% transpose X2+ and insert X1+
        %%  use the set of ring_atoms (backbone of the ring) to transpose order
        %%  change [@->@@] and [@@->@] and after atom
        
        %% X2transpose
        
        %% CN1CC[C@]23C4=C5C=C C(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5
        %%
        %%               first second  ring atom
        %%
        %%  break CC:
        %%
        %% CN1CC[C@]23C4=C5  C=C  C(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5   XXXXX     C C are in rings 1,2,3,4,5
        %% ----------------  ---  ----------------------------------------
        %%
        %%    X1-            X2-            X2+  (no X1+)                    XXXXX     X...X at end is X1+
        %%
        %%    X1- (X2+'(first) (X1+) X2+'(rest) ) X2-
        %%
        %%
        %%   CN1CC[C@]23C4=C ( C (XXXXX) [C@H]1[C@H]3C=C[C@@H](O)[C@H]2OC4=C(O) ) C=C
        %%
        %% because of the direction of branch at C5, reverse the chirality in the X2+
        %%
        %%   CN1CC[C@]23C4=C ( C (XXXXX) [C@@H]1[C@@H]3C=C[C@H](O)[C@@H]2OC4=C(O) ) C=C
        %%
        %% to avoid large branch notation interchange
        %%
        %%    CN1CC[C@]23C4=C (C=C) C (XXXXX) [C@@H]1[C@@H]3C=C[C@H](O)[C@@H]2OC4=C(O)
        
        
        
        %% now break the C[CH] in the large branch
        %%
        %%   CN1CC[C@]23 C4=C(C(XXXXX)[C@@H]1[C@@H]3C=C   [C@H](O)[C@@H]2  OC4=C(O))C=C  branch needs to be taken into account
        
        
        %%       X1-         X2-                                X2+           X1+
        
        %%   CN1CC[C@H]3 ( C (OC4=C(O))C=C) C[@@H](O) )  C4=C(C(XXXXX)[C@@H]1[C@@H]3C=C
        
        %% algorithm can be modified
        
        
        %% calculate min_bb_ring_atom, max_bb_ring_atom of X2+ - not done
        %%  from chm sequence
        
        %% X2plus_X1plus
        %% 1st atom of X2+
        total=0;
        if chm(ring_atom(max_bb_ring_atom)-1)=='['
            delete_one=1;
        end
        X2plusX1plus(total+1:total+1+last_chm_point_X2plus-ring(max_bb_ring_atom)-delete_one-1-1)=chm(ring_atom(max_bb_ring_atom)-delete_one:last_chm_point_X2plus-1);
        total=total+last_chm_point_X2plus-ring(max_bb_ring_atom)-delete_one+1;  %% minus the int
        
        %% add X1plus
        X2plusX1plus(total+1)='(';
        total=total+1;
        X2plusX1plus(total+1:total+1+size(X1plus,2))=X1plus(1:size(X1plus,2));
        total=total+size(X1plus,2);
        X2plus1plus(total+1)=')';
        total=total+1;
        
        %% 2nd to end of X2+ atoms
        for atom=max_bb_ring_atom:-1:min_bb_ring_atom-1
            %%  char array index (location in chm) is ring_atoms(atom)
            %%  takes into account the [ in all ring atom parts - inclusion is from heavy atom -> [ heavy atom
            X2plus(total+1:total+1+ring_atom(atom+1)-ring(atom)-1)=chm(ring_atom(atom):ring_atom(atom+1)-1);
            total=total+ring_atom(atom+1)-ring(atom);
        end
        X2plus(total+1:total+1+ring_atom(atom+1)-ring(atom)-1)=chm(ring_atom(atom):ring_atom(atom+1)-1);
        total=total+ring_atom(atom+1)-ring(atom);
        
    end
    
    
    attempt=attempt+1;
    
    %% create pdb file
    chm
    system("rm molecule/*.*");
    fileID=fopen('molecule.smi','w');
    fprintf(fileID,'%s',char(new_chm));
    fclose(fileID);
    
    system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,axchir,msi=50,msc=10,names,preserve molecule.smi molecule/molecule.pdb");
    %% if ok the file test.pdb will exist
    success1=0;
    if exist("molecule/molecule.001.pdb")>0
        if dir("molecule/molecule.001.pdb").bytes>0
            success=1;
        end
    end
    
    %% print warning if pdb wasn't created - shouldn't happen
    if exist("molecule/molecule.001.pdb")==0
        new_chm=chm;
        success=0;
    end
     
end

if success==0
    new_chm=chm;
end

end

%%   break [C@]23 C4
%%
%% CN1CC[C@]23C4=C5C=CC(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5
%%
%%  CN1CC[C@]23  C4=C5C=CC(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3  [C@H]1C5    ring 3
%%
%%    X1-           X2+                                      X2-
%%
%%  CN1CC[C@]2(was 3)  ( [C@H](was 3) ([C@H]1C5) C=C[C@@H](O)[C@H]2OC4=C(O)C=CC5=C4 )
%%
%%  renumber rings


%%   break C(O)=C
%%
%% CN1CC[C@]23C4=C5C=CC(O) =C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5   X1plus       ring 5
%%
%%  CN1CC[C@]23C4=C5  C=CC(O)
%%
%%    X1-             X2-                 X2+                       no X1+, after 2nd ring atom
%%
%%  CN1CC[C@]23C4=C ( C (was 5) (X1plus) [C@H]1[C@H]3C=C[C@@H](O)[C@H]2OC4=C(O)C=CC5=C4 ) C=CC(O)
%%
%%    X1-                                  X2+'                                   X2-  (then chirality change in X2+ but not X1+)
%%
%%        X1- (X2+'(first) (X1+) X2+'(rest) ) X2-



