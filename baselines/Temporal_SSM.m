function TSSM = Temporal_SSM(des, descrip_flag, slide_win, kernel, belta, c)
% TEMPORAL_SSM  Build a temporal self-similarity matrix from descriptor rows.
%   slide_win defines the temporal patch; descrip_flag selects the metric.
%   This is the legacy single-scale TSSM helper, not Multiscale SSM [23].
%   Some metric branches require external distance functions; see the review.

% Defaults restore the behaviour of the older, shorter call forms still used by
% the ASL and MSRC-12 generators, e.g. Temporal_SSM(des,flag) and
% Temporal_SSM(des,5,1,1).  belta = 0.5e-3, c = 0 are the ASL values of the
% sigmoid distance in the conference paper (they used to be hard-coded here).
if nargin < 3
    slide_win = 1;
end
if nargin < 4
    kernel = 0;
end
if nargin < 5
    belta = 0.5e-3;
end
if nargin < 6
    c = 0;
end
m = size(des, 1); % temporal length
n = size(des, 2); % dimension
d = zeros(m, m); % Similarity matrix
%% direct computing a SSM
patch = floor(slide_win / 2); % patch size
des_matrix = zeros(m - slide_win + 1, n * slide_win);
if descrip_flag == 5 || descrip_flag == 6
    for i = patch + 1:m - patch
        % consturct the patch with size = 5
        des_matrix(i - patch, 1:n * slide_win) = reshape(des(i - patch:i + patch, :), 1, slide_win * n);
    end
else
    for i = patch + 1:m - patch
        des_matrix(i - patch, 1:(n * slide_win) / 2) = reshape(des(i - patch:i + patch, 1:n / 2), 1, ...
            slide_win * (n / 2)); % consturct the patch with size = 5
        des_matrix(i - patch, (n * slide_win) / 2 + 1:end) = reshape(des(i - patch:i + patch, n / 2 + 1:end), ...
            1, slide_win * (n / 2));
    end
end

switch descrip_flag
    case {2, 3, 4}
        d = feature_dist_matching(des_matrix, des_matrix);
    case 1
        d = distance_matrix_fd(des_matrix, des_matrix);
    case {5, 6}
        if kernel == 0
            d = distance_matrix_norm2(des_matrix, des_matrix);
        elseif kernel == 1
            % the new metric
            %             belta = 1e-8;
            %             d = pdist2(des_matrix,des_matrix).^2;
            %             d = 1-exp(-belta*d);
            %             d = d - mean(mean(d));
            %             belta =sum(sum((d.^2)))/(m^2);

            % tanh function
            %             belta = 0.8e-7;%default
            %             belta = 0.9e-7;
            %             tanh_matrix = belta*(des_matrix*des_matrix');
            %             d = tanh(tanh_matrix);

            % tanh function of norm2
            distNorm = pdist2(des_matrix, des_matrix);
            %             belta = 0.5e-3; % for asl dataset
            tanh_matrix = belta * distNorm - c;
            d = tanh(tanh_matrix);

            %             % for debug
            %             if min(min(tanh_matrix))<-1 || max(max(tanh_matrix))>1
            %                 d;
            %             end
        end
    case {7, 8}
        d = distance_matrix_norm1(des_matrix, des_matrix);
    case 9
        if kernel == 0
            d = hist_cost_2(des_matrix, des_matrix);
            %             d = pdist2(des_matrix,des_matrix);
        elseif kernel == 1
            tanh_matrix = hist_cost_2(des_matrix, des_matrix);
            d = tanh(tanh_matrix);

            %             % for debug
            %             if min(min(0.8e-7*(des_matrix*des_matrix')))<-1 || max(max(0.8e-7*(des_matrix*des_matrix')))>1
            %                 d;
            %             end
        end
    otherwise
        disp('error input descrip_flag')
end
TSSM = d;
%% first computing a SSM between individual elements and then perform image filtering
% des_matrix = des;
% switch descrip_flag
%     case  {2,3,4}
%         d = feature_dist_matching(des_matrix,des_matrix);
%     case 1
%         d = distance_matrix_fd(des_matrix,des_matrix);
%     case {5,6}
%         d = distance_matrix_norm2(des_matrix,des_matrix);
%     case {7,8}
%         d = distance_matrix_norm1(des_matrix,des_matrix);
%     otherwise
%         disp('error input descrip_flag')
% end
% % H = fspecial('gaussian',slide_win,slide_win);%gaussian filter
% H = eye(slide_win);% sum filter
% TSSM = imfilter(d,H);
