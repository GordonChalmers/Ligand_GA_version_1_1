%% Gordon Chalmers 10/21

%% objective function for mass density
%%
%%  Inputs: molecules, i.e. chm
%%  Output: fitnesses to drive the GA

function FitnessFcn=Ligand_GA_Fitness_Function_Atom_Density(x)

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

global ligand_dir;

fitness=zeros(size(x,1),1);
FitnessFcn=zeros(size(x,1),1);

for molecule_idx=1:size(x,1)
    
    largest_mass_density=0;
    
    molecule_idx;
    chm=x{molecule_idx};
    chm=char(chm);
    
    %% create pdb file
    %% each molecule_idx accesses a different directory
    rm_temp_dir_idx="rm -r "+ligand_dir+"/molecule"+molecule_idx+"/*.*";
    system(rm_temp_dir_idx);
    file_open=ligand_dir+"/molecule"+molecule_idx+"/molecule.smi";
    fileID=fopen(file_open,'w');
    fprintf(fileID,'%s',char(chm));
    fclose(fileID);
    
    system(corina_path + " -i t=smiles -o t=pdb,xlabel,pdbelement,split -d wh -d stergen,msi=50,msc=10,names,preserve "+ligand_dir+"/molecule"+molecule_idx+"/molecule.smi "+ligand_dir+"/molecule"+molecule_idx+"/molecule.pdb");
    %%    if exist("molecule"+molecule_idx+"/molecule.001.pdb")>0
    success=0;
    dir_idx=ligand_dir+"/molecule"+molecule_idx;
    dir_test=dir(dir_idx);
    if size(dir_test,1)>3
        if dir_test(3).bytes>0
            success=1;
        end
    end
    
    if success>0
        
        %% load molecule information
        [molecule,chm_len,adj,adj_heavy,num_heavy_atoms,heavy_atom_list,heavy_idx_chm,num_rings,ring_idx_chm,adj_atom,avail_heavy_bond, ...
            chiral,num_bonds_left,num_bonds_right]=MoleculeStructure(chm);
        
        if num_heavy_atoms>min_heavy_atoms
            
            %% structure with directory information
            directory=dir(ligand_dir+"/molecule"+molecule_idx+"/*.*");
            number_files=size(directory,1)-3;
            
            FitnessFcn(molecule_idx)=0;
            
            if number_files>0
                
                %% calculate fitnesses of all stereoisomers
                largest_mass_density=0;
                
                for fileno=1:number_files
                    
                    %% load molecule information - has coordinates
                    if fileno<10
                        molecule_file=ligand_dir+"/molecule"+molecule_idx+"/molecule.00"+fileno+".pdb";
                    end
                    if fileno>=10 && fileno<100
                        molecule_file=ligand_dir+"/molecule"+molecule_idx+"/molecule.0"+fileno+".pdb";
                    end
                    if fileno>=100
                        molecule_file=ligand_dir+"/molecule"+molecule_idx+"/molecule."+fileno+".pdb";
                    end
		    %% molecule_file
                    molecule=pdbread(molecule_file);
                    
                    %% coordinates  molecule.Model.HeterogenAtom(atom).X , Y, Z
                    %% atom type    molecule.Model.HeterogenAtom(atom).element
                    %% coordinates
                    points=zeros(num_heavy_atoms,3);
                    for atom=1:num_heavy_atoms
                        points(atom,1)=molecule.Model.HeterogenAtom(heavy_atom_list(atom)).X;
                        points(atom,2)=molecule.Model.HeterogenAtom(heavy_atom_list(atom)).Y;
                        points(atom,3)=molecule.Model.HeterogenAtom(heavy_atom_list(atom)).Z;
                    end
                    
                    %% center of mass
                    com=zeros(1,3);
                    total_mass=0;
                    for atom=1:num_heavy_atoms
                        if molecule.Model.HeterogenAtom(atom).element=="C"
                            com(1)=points(atom,1)*12.0107+com(1);
                            com(2)=points(atom,2)*12.0107+com(2);
                            com(3)=points(atom,3)*12.0107+com(3);
                            total_mass=total_mass+12.0107;
                        end
                        if molecule.Model.HeterogenAtom(atom).element=="N"
                            com(1)=points(atom,1)*14.0067+com(1);
                            com(2)=points(atom,2)*14.0067+com(2);
                            com(3)=points(atom,3)*14.0067+com(3);
                            total_mass=total_mass+14.0067;
                        end
                        if molecule.Model.HeterogenAtom(atom).element=="O"
                            com(1)=points(atom,1)*15.9994+com(1);
                            com(2)=points(atom,2)*15.9994+com(2);
                            com(3)=points(atom,3)*15.9994+com(3);
                            total_mass=total_mass+15.9994;
                        end
                        if molecule.Model.HeterogenAtom(atom).element=="S"
                            com(1)=points(atom,1)*32.065+com(1);
                            com(2)=points(atom,2)*32.065+com(2);
                            com(3)=points(atom,3)*32.065+com(3);
                            total_mass=total_mass+32.065;
                        end
                        if molecule.Model.HeterogenAtom(atom).element=="P"
                            com(1)=points(atom,1)*30.9738+com(1);
                            com(2)=points(atom,2)*30.9738+com(2);
                            com(3)=points(atom,3)*30.9738+com(3);
                            total_mass=total_mass+30.9738;
                        end
                        if molecule.Model.HeterogenAtom(atom).element=="F"
                            com(1)=points(atom,1)*18.9984+com(1);
                            com(2)=points(atom,2)*18.9984+com(2);
                            com(3)=points(atom,3)*18.9984+com(3);
                            total_mass=total_mass+18.9984;
                        end
                    end
                    com(:)=com(:)/total_mass;
                    
                    %% find sphere of size that encompasses all heavy atoms centered about c.o.m.
                    molecule_radius=0;
                    for atom=1:num_heavy_atoms
                        %% be careful of using norm in Matlab - not 2-norm for a matrix
                        if sqrt((points(atom,1)-com(1))^2+(points(atom,2)-com(2))^2+(points(atom,3)-com(3))^2)>molecule_radius
                            molecule_radius=sqrt((points(atom,1)-com(1))^2+(points(atom,2)-com(2))^2+(points(atom,3)-com(3))^2);
                        end
                    end
                    
                    mass_density=total_mass/(4/3*pi*molecule_radius^3);  %% not point density
                    if mass_density>largest_mass_density
                        largest_mass_density=mass_density;
                        largest_mass_density_stereoisomer=fileno;
                    end
                end
            end
        end  %% num_heavy_atoms>min_heavy_atoms
        
        %% a check - this shouldn't happen
        if num_heavy_atoms<min_heavy_atoms
            largest_mass_density=.1;
        end
        
        FitnessFcn(molecule_idx)=-largest_mass_density;
        
    end
    
end

end


