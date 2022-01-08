function Parse_Results(ligand_dir,OutDockingName,idx_first)

%% input idx_first - offset in labeling in OutDockingName file
%% ligand_dir

%% example
%% ligand_dir='/home/gordoncs/Ligand_GA_MultiObj/Paxlovid_CYP_3A4/CYP_3A4_molecule_47_64';
%% idx_first=1;
%%OutDockingName='/home/gordoncs/Ligand_GA_MultiObj/Paxlovid_CYP_3A4/CYP_3A4_molecule_47_64/Summary_CYP_3A4_molecule_47_64.txt';

%% number of molecules in directory
%% should only be molecule# sub-directories in ligand_dir
number_molecules=size(dir(ligand_dir),1)-2;

%% out file
fid=fopen(OutDockingName,'w');
fprintf(fid,'%s\n',ligand_dir);
fprintf(fid,'%s\n',OutDockingName);
fprintf(fid,'molecule # offset %d\n\n',idx_first);
fprintf(fid,'%s %s %s %s %s %s %s %s %s %s %s\n\n','molecule','steiso','Score','S(PLP)','S(hbond)','S(cho)','S(metal)','DE(clash)','DE(tors)','intcor','time');

for molecule_idx=1:number_molecules
  %% determine steiso_num of a molecule_idx
  num_steiso=(size(dir(ligand_dir+"/molecule"+molecule_idx),1)-4)/3;
  
  for fileno=1:num_steiso
     log_file=ligand_dir+"/molecule"+molecule_idx+"/stereoisomer"+fileno+"/bestranking.lst";

%% parse bestranking.lst file for the docking information
%% check if an almost enpty bestranking.lst file 
    if exist(log_file,'file')>0 && dir(log_file).bytes>289
%% result
       fid_read=fopen(log_file,'r');
       for i=1:8
          text=fgetl(fid_read);
       end
       text2=strsplit(text);  %% the first in cell array is null
%% text2=strsplit(text);
       fclose(fid_read);

%% write output file in text line form, not csv
    fprintf(fid,'%d %d %s %s %s %s %s %s %s %s %s\n',molecule_idx+(idx_first-1),fileno,text2{2},text2{3},text2{3},text2{4},text2{5},text2{6},text2{7},text2{8},text2{9});
    end

%% there are 9 columns
%% Format is:
%%   Score S(PLP) S(hbond) S(cho) S(metal) DE(clash) DE(tors) intcor time 
%%   58.18 -48.99 3.61     0.00   0.00     8.19      1.03     8.60   430.347 
%% the 2nd is a text line, to be csv after loading in excel

  end  %% fileno 
end  %% molecule_idx

end  %% Parse_Results
