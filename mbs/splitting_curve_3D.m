function MBS = splitting_curve_3D(Curve_xyz, Width)
% SPLITTING_CURVE_3D  Split a 3-D trajectory using its XY and XZ projections.
%   Curve_xyz is N-by-3; Width is the allowed segment thickness.
%   MBS.data contains segments and MBS.B_E_seg their endpoint indices.
%   Requires the Determine_segment MEX implementation.

k = 1;
Sb = Curve_xyz(1, :);
n = size(Curve_xyz, 1);
data = {};
parameter = [0 0 0 0];
B_E_seg = [];
MBS.data = data;
MBS.parameter = parameter;
MBS.B_E_seg = B_E_seg;

%% begin main procedure
while k < n
    k = k + 1;
    Sb(end + 1, :) = Curve_xyz(k, :);
    %% determine the OXY and OXZ plane parameters respectively
    [~, eSegment_xy] = Determine_segment(Sb(:, 1:2), Width);
    [~, eSegment_xz] = Determine_segment(Sb(:, 1:2:3), Width);
    if eSegment_xz == 0
        eSegment_xz = 1000;
    end % MAX U point = 1000 in MEX files
    if eSegment_xy == 0
        eSegment_xy = 1000;
    end % MAX U point = 1000 in MEX files
    if min(eSegment_xy, eSegment_xz) < k
        break;
    end
end
bSegment = 1;
eSegment = min(eSegment_xy, eSegment_xz);
MBS.data{end + 1} = Curve_xyz(bSegment:eSegment, :);
MBS.B_E_seg(end + 1, :) = [bSegment eSegment];
while k < n
    while k < n
        %% perry
        bSegment = bSegment + 1;
        Sb = Curve_xyz(bSegment:k, :);
        %% determine the OXY and OXZ plane parameters respectively
        [~, eSegment_xy] = Determine_segment(Sb(:, 1:2), Width);
        [~, eSegment_xz] = Determine_segment(Sb(:, 1:2:3), Width);
        if eSegment_xz == 0
            eSegment_xz = 1000;
        end % MAX U point = 1000 in MEX files
        if eSegment_xy == 0
            eSegment_xy = 1000;
        end % MAX U point = 1000 in MEX files
        if min(eSegment_xy, eSegment_xz) == size(Sb, 1)
            break;
        end
    end
    while k < n
        k = k + 1;
        Sb(end + 1, :) = Curve_xyz(k, :);
        %% determine the OXY plane parameters
        [~, eSegment_xy] = Determine_segment(Sb(:, 1:2), Width);
        [~, eSegment_xz] = Determine_segment(Sb(:, 1:2:3), Width);
        if eSegment_xz == 0
            eSegment_xz = 1000;
        end % MAX U point = 1000 in MEX files
        if eSegment_xy == 0
            eSegment_xy = 1000;
        end % MAX U point = 1000 in MEX files
        if min(eSegment_xy, eSegment_xz) < size(Sb, 1)
            break;
        end
    end
    eSegment = k - 1;
    MBS.data{end + 1} = Curve_xyz(bSegment:eSegment, :);
    MBS.B_E_seg(end + 1, :) = [bSegment eSegment];
end
