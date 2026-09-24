function load_IPtxt_bat(BAT_FOLDER)
% LOAD_IPTXT_BAT  Load the prepared IP hand trajectories and create train/test files.
%   The active path reads Hand_Data/Data.mat and resamples to 64 frames.
%   BAT_FOLDER is used only by the commented raw-text loading alternative.

file_ext = '.txt';
db_i = 0;
test_i = 0;
personid = [0, 2:11, 13:15, 17:22];
db_id = personid(1:10);
test_id = personid(11:end);

%% read random db data from class folder
% folder_content = dir([BAT_FOLDER,'*',file_ext]);
% ndata=size(folder_content,1);
%  TRAJDB = cell (2,8000);
%  TRAJSAMPLES = cell (2,8000);
% for k=1:ndata;
%     string= [BAT_FOLDER,folder_content(k,1).name];
%     fprintf ('Loading txt data...%s\n',string);
%     %% read the txt files
%     [trajectory,id] = readIPtxt(string);
%     if id(2) ==17 || id(2) == 18
%         continue;
%     end
%     if ismember(id(1),db_id)
%         db_i = db_i+1;
%         TRAJDB{1,db_i} = id;
%         TRAJDB{2,db_i} = trajectory;
%     else
%          test_i = test_i+1;
%          TRAJSAMPLES{1,test_i} = id;
%         TRAJSAMPLES{2,test_i} = trajectory;
%     end
% end
% TRAJDB(:,2401)=[];
% save DB TRAJDB;
% save SAMPLES TRAJSAMPLES;
% fclose('all');

%% read  db data from Data.mat given by Gy
personid = [1, 2:11, 13:15, 17:22];
db_id = personid(1:10);
test_id = personid(11:end);
TRAJDB = cell (2, 8000);
TRAJSAMPLES = cell (2, 8000);
folder = 'Hand_Data/';
load ([folder 'Data.mat']); % in Data.mat, 1. gesture id; 2. personid;
fixe_length = 64;
for n = 1:size(Data, 1)
    fprintf ('Loading the %dth data...\n', n);
    id = Data{n, 2};
    trajectory = Data{n, 1};
    if id(1) == 17 || id(1) == 18
        continue;
    end
    if ismember(id(2), db_id)
        db_i = db_i + 1;
        TRAJDB{1, db_i} = id;
        ResizeTrj1 = TrjResizeTime(trajectory(:, 1:3), fixe_length);
        ResizeTrj2 = TrjResizeTime(trajectory(:, 4:6), fixe_length);
        TRAJDB{2, db_i} = [ResizeTrj1 ResizeTrj2];
        %         TRAJDB{2,db_i} = trajectory;
    else
        test_i = test_i + 1;
        TRAJSAMPLES{1, test_i} = id;
        ResizeTrj1 = TrjResizeTime(trajectory(:, 1:3), fixe_length);
        ResizeTrj2 = TrjResizeTime(trajectory(:, 4:6), fixe_length);
        TRAJSAMPLES{2, test_i} = [ResizeTrj1 ResizeTrj2];
        %         TRAJSAMPLES{2,test_i} = trajectory;
    end
end
save DB TRAJDB;
save SAMPLES TRAJSAMPLES;
