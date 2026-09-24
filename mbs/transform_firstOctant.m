function Sb_trans = transform_firstOctant(Sb, inremental)
% TRANSFORM_FIRSTOCTANT  Transform points into the first octant for 2-D recognition.
%   The historical inremental argument name and cached transform are retained.
%   Known chained-comparison issues are documented in docs/CODE_REVIEW.md.

if nargin > 1 && inremental > 2
    load trasform_matrix
    Sb(:, end + 1) = 1; % homogeneous coordinates
    Sb_trans = transform_matrix * Sb'; %#ok<NODEF> % perform transformation
    Sb_trans = Sb_trans';
    Sb_trans(:, end) = [];
    Sb_trans = int32(Sb_trans);
    return;
end

transform_matrix = [1 0 0 0; ...
                    0 1 0 0; ...
                    0 0 1 0; ...
                    0 0 0 1]; %#ok<NASGU>

if size(Sb, 2) > 2
    t1 = Sb(1, 1);
    t2 = Sb(1, 2);
    t3 = Sb(1, 3);
    translation_matrix = [1 0 0 -t1; ...
                          0 1 0 -t2; ...
                          0 0 1 -t3; ...
                          0 0 0 1];
    Sb(:, end + 1) = 1; % homogeneous coordinates
    Sb_transla = translation_matrix * Sb';
    %% OXY plane rotation
    Sb_transla_xy = Sb_transla(1:2, :);
    Sb_transla_xy(:, ismember(Sb_transla_xy(1:2, :)', [0 0], 'rows')) = [];
    octant = atan2(Sb_transla_xy(2, 1:end), Sb_transla_xy(1, 1:end));
    %   octant(find(octant == 0))=[]; %#ok<FNDSB> eliminate the 0 element.

    if ~isempty(octant)
        min_octant = min(octant);
        max_octant = max(octant);

        if -pi < min_octant < -3 * pi / 4 && max_octant > 0
            trans_angle = -max_octant;
        else
            trans_angle = -min_octant;
        end

        % min_octant = -min_octant; % positive angle represent the anticlockwise rotation
        % rotate by z axis on OXY plane
        rotate_matrix = [cos(trans_angle) -sin(trans_angle) 0 0; ...
                         sin(trans_angle) cos(trans_angle) 0 0; ...
                         0 0 1 0; ...
                         0 0 0 1];
        Sb_trans = rotate_matrix * Sb_transla; % perform transformation
        transform_matrix = rotate_matrix * translation_matrix;
        %% OXZ plane rotation -- rotate Y axis on OXZ
        Sb_trans_xz = Sb_trans(1:2:3, :);
        Sb_trans_xz(:, ismember(Sb_trans_xz(1:2, :)', [0 0], 'rows')) = [];
        octant = atan2(Sb_trans_xz(2, 1:end), Sb_trans_xz(1, 1:end));
        %     octant(find(octant == 0))=[]; %#ok<FNDSB> eliminate the 0 element.
        if isempty(octant)
            Sb_trans = Sb_trans';
            Sb_trans(:, end) = [];
            Sb_trans = int32(Sb_trans);
            return;
        end

        min_octant = min(octant);
        max_octant = max(octant);
        if -pi < min_octant < -3 * pi / 4 && max_octant > 0
            trans_angle = -max_octant;
        else
            trans_angle = -min_octant;
        end

        rotate_matrix = [cos(-trans_angle) 0 sin(-trans_angle) 0; ...
                             0 1 0 0; ...
                             -sin(-trans_angle) 0 cos(-trans_angle) 0; ...
                             0 0 0 1];
        Sb_trans = rotate_matrix * Sb_trans;
        transform_matrix = rotate_matrix * transform_matrix;

        Sb_trans = Sb_trans';
        Sb_trans(:, end) = [];
        Sb_trans = int32(Sb_trans);
        return;
    else
        Sb_transla_xz = Sb_transla(1:2:3, :);
        %% OXZ plane rotation -- rotate Y axis on OXZ
        Sb_transla_xz(:, ismember(Sb_transla_xz(1:2, :)', [0 0], 'rows')) = [];
        octant = atan2(Sb_transla_xz(2, 1:end), Sb_transla_xz(1, 1:end));
        %         octant(find(octant == 0))=[]; %#ok<FNDSB> eliminate the 0 element.
        if isempty(octant)
            Sb_trans = Sb_transla';
            Sb_trans(:, end) = [];
            Sb_trans = int32(Sb_trans);
            return;
        end

        min_octant = min(octant);
        max_octant = max(octant);
        if -pi < min_octant < -3 * pi / 4 && max_octant > 0
            trans_angle = -max_octant;
        else
            trans_angle = -min_octant;
        end

        rotate_matrix = [cos(-trans_angle) 0 sin(-trans_angle) 0; ...
                         0 1 0 0; ...
                         -sin(-trans_angle) 0 cos(-trans_angle) 0; ...
                         0 0 0 1];
        Sb_trans = rotate_matrix * Sb_transla;
        transform_matrix = rotate_matrix * translation_matrix;
        Sb_trans = Sb_trans';
        Sb_trans(:, end) = [];
        Sb_trans = int32(Sb_trans);
        return;
    end
end

t1 = Sb(1, 1);
t2 = Sb(1, 2);
translation_matrix = [1 0 -t1; ...
                      0 1 -t2; ...
                      0 0 1];
Sb(:, end + 1) = 1; % homogeneous coordinates
Sb_transla = translation_matrix * Sb';
Sb_transla_xy = Sb_transla;
Sb_transla_xy(:, ismember(Sb_transla_xy(1:2, :)', [0 0], 'rows')) = [];
octant = atan2(Sb_transla_xy(2, 1:end), Sb_transla_xy(1, 1:end));
% octant(find(octant == 0))=[]; %#ok<FNDSB> eliminate the 0 element.
if isempty(octant)
    Sb_trans = Sb_transla';
    Sb_trans(:, end) = [];
    return;
end
min_octant = min(octant);
max_octant = max(octant);

if 0 >= max_octant && min_octant < 0
    trans_angle = -min_octant;
elseif 0 < max_octant && min_octant >= 0
    trans_angle = -min_octant;
else
    trans_angle = -max_octant;
end

% min_octant = -min_octant; % positive angle represent the anticlockwise rotation
rotate_matrix = [cos(trans_angle) -sin(trans_angle) 0; ...
                 sin(trans_angle) cos(trans_angle) 0; ...
                 0 0 1];

Sb_trans = rotate_matrix * Sb_transla; % perform transformation
transform_matrix = rotate_matrix * translation_matrix;
Sb_trans = Sb_trans';
Sb_trans(:, end) = [];
Sb_trans = int32(Sb_trans);
