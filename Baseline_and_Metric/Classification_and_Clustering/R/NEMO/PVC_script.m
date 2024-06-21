
% variables called subtype_name and expected_prefix should be defined
addpath('PVC_CODE_PATH')

A = importdata(strcat('RESULTS_PATH/', subtype_name, '/1'));
B = importdata(strcat('RESULTS_PATH/', subtype_name, '/2'));

file_names = dir(['PVC_INDICES_PATH/', subtype_name.cellstr{1}]);
num_patients = size(A.data, 2);
only_file_names = {};

for i=3:size(file_names, 1)
    only_file_names{i-2} = file_names(i).name;
end
    
for j=3:size(file_names,1)
    rand('seed',42);
    clear truth
    clear full_clustering
    file_name = file_names(j).name;
    if endsWith(file_name, 'clustering') || ismember([file_name, '_clustering'], only_file_names) || (~startsWith(file_name, expected_prefix))
        continue
    end
    
    tmp = strsplit(file_name, '_');
    num_clusters = str2num(tmp{3});
    num_clusters
    full_file_name = ['PVC_INDICES_PATH/', subtype_name.cellstr{1},'/', file_name];
    file_content = importdata(full_file_name);
    tic;
    rand_vec1 = repmat(true, num_patients, 1);
    rand_vec2 = repmat(false, num_patients, 1);
    rand_vec2(file_content) = true;
    
    common_samples_logical = rand_vec1 & rand_vec2;
    only_A_logical = rand_vec1 & ~rand_vec2;
    only_B_logical = ~rand_vec1 & rand_vec2;
    X2 = transpose(A.data(:,common_samples_logical));
    Y2 = transpose(B.data(:,common_samples_logical));
    X1 = transpose(A.data(:,only_A_logical));
    Y3 = transpose(B.data(:,only_B_logical));

    truth(1:num_patients) = randi([1,2], num_patients, 1);
    
    options = [];
    options.lamda = 0.01;
    options.latentdim = num_clusters;

    [Ux Uy P2 P1 P3 objValue C] = PVCclust(X2, Y2, X1, Y3, num_clusters, truth, options);
    
    num_common = sum(common_samples_logical);
    num_A_only = sum(only_A_logical);
    num_B_only = sum(only_B_logical);
    full_clustering(1:num_patients) = -1;
    full_clustering(common_samples_logical) = C(1:num_common);
    full_clustering(only_A_logical) = C((num_common +1 ):(num_common + num_A_only));
    full_clustering(only_B_logical) = C((num_common + num_A_only + 1):(num_common + num_A_only + num_B_only));
    
    
    time_elapsed = toc;
    dlmwrite([full_file_name, '_', num2str(time_elapsed), '_clustering'], full_clustering);

end