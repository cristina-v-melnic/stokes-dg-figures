folder_dir = "files/1011/search_1/";

listing = dir(folder_dir);
fnames = {listing.name};
n_files = numel(fnames);
sigmas = [];
names = [];

for i = 1:n_files
    if contains(fnames(i),"sigma")
      newStr = split(fnames(i),'_');
      sigmas = [sigmas, str2double(newStr(7))];
      names = [names, string(fnames(i))];
    end    
end

nr_matrices = numel(names)/4;
%cond_nrs = zeros(3, nr_matrices/3);
cond_nrs = [];

%strcat("velocity_",vel_spaces(3))
for i=1:nr_matrices
        filename = strcat(folder_dir, string(names(i)));
        A= readmatrix(filename,delimitedTextImportOptions('DataLines', [3, Inf]));
        A = cell2table(A);
        A(:, end) = [];
        
        filename_1 = strcat(folder_dir, string(names(nr_matrices+i)));
        BT = readmatrix(filename_1,delimitedTextImportOptions('DataLines', [3, Inf]));
        BT = cell2table(BT);
        BT(:, end) = [];
        
        M1 = [A, BT];
        
        filename_2 = strcat(folder_dir, string(names(2*nr_matrices+i)));
        B = readmatrix(filename_2,delimitedTextImportOptions('DataLines', [3, Inf]));
        B = cell2table(B);
        B(:, end) = [];
        
        filename_3 = strcat(folder_dir, string(names(3*nr_matrices+i)));
        C = readmatrix(filename_3,delimitedTextImportOptions('DataLines', [3, Inf]));
        C = cell2table(C);
        C(:, end) = [];
        
        M2 = [B, C];
        
        M1 = table2array(M1);
        M2 = table2array(M2);
        
        M = str2double(vertcat(M1, M2));

        cond_nrs = [cond_nrs, cond(M)];
end
%output the data(
file_out = fopen(strcat(folder_dir,"matlab_results"),'w');
%fprintf(file_out, '%6s %12s %12s %12s\r\n', 'sigma', 'TH2', 'RT1', 'BDM2');
%fprintf(file_out, '%2.e %12.2f %12.2f %12.2f\r\n',[sigmas(1:nr_matrices/3); cond_nrs(1,:); cond_nrs(2,:); cond_nrs(3,:)]);
fprintf(file_out, '%2.1e %2.2e\r\n',[sigmas(1:nr_matrices); cond_nrs]);

fclose(file_out);

%plots, uncomment if needed
%scatter(sigmas(1:nr_matrices), cond_nrs)
%set(gca,'xscale','log')
%set(gca,'yscale','log')
