function load_data_bat(joints_no, shuffle_sort)
% LOAD_DATA_BAT  Create MSRC-12 train/test files from MSRC12_Skeleton.mat.
%   joints_no contains joint IDs as strings; shuffle_sort(1,:) selects training
%   subjects. Trajectories are resampled to 128 frames per joint.

load MSRC12_Skeleton;
class_no = size(SmthTrj, 1);
subject_no = size(SmthTrj, 2);
joints_num = length(joints_no);
joints = zeros(1, joints_num);
TRAJDB = cell(2, []);
TRAJSAMPLES = cell(2, []);
fixed_length = 128;
for i = 1:joints_num
    joints(1, i) = str2double(joints_no{i});
end

for i = 1:class_no
    for j = 1:subject_no
        Traj_Data = SmthTrj{i, j};
        if i == 6 && j == 18
            continue;
        end
        if ismember(j, shuffle_sort(1, :))
            for n = 1:length(Traj_Data(joints(1, 1), :)) % the number of subjects
                TRAJDB{1, end + 1} = [i, j];
                for k = 1:joints_num % the joints read
                    fprintf('loading the samples %d-%d-%d\n', i, j, k);
                    TRAJDB{2, end} = [TRAJDB{2, end} TrjResizeTime(Traj_Data{joints(1, k), n}, fixed_length)];
                    %                     TRAJDB{2,end} = [TRAJDB{2,end} Traj_Data{joints(1,k),n}];
                end
            end
        else
            for n = 1:length(Traj_Data(joints(1, 1), :)) % the number of subjects
                TRAJSAMPLES{1, end + 1} = [i, j];
                for k = 1:joints_num % the joints read
                    fprintf('loading the samples %d-%d-%d\n', i, j, k);
                    TRAJSAMPLES{2, end} = [TRAJSAMPLES{2, end} TrjResizeTime(Traj_Data{joints(1, k), n}, fixed_length)];
                    %                     TRAJSAMPLES{2,end} = [TRAJSAMPLES{2,end} Traj_Data{joints(1,k),n}];
                end
            end
        end
    end
end

save DB TRAJDB;
save SAMPLES TRAJSAMPLES;
