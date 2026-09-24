function MBS = splitting_curve_2D(Curve_xy, Width)
% SPLITTING_CURVE_2D  Split a 2-D trajectory into overlapping blurred segments.
%   Curve_xy is N-by-2; Width is the allowed segment thickness.
%   MBS.data contains segments and MBS.B_E_seg their endpoint indices.
%   This legacy branch needs structure/testOctant helpers not bundled here.

k = 1;
Sb = [Curve_xy(1, :)];
n = size(Curve_xy, 1);
global a;
global b;
global u;
global w;
a = 0;
b = 1;
u = 0;
w = b;
isSameOctant = true;
data = {};
parameter = [0 0 0 0];
B_E_seg = [];
MBS = structure(data, parameter, B_E_seg);
% decremental = 1;
global incremental;
incremental = 1;
%% begin main procedure
while (w - 1) / max(abs(a), abs(b)) <= Width && k < n
    k = k + 1;
    incremental = incremental + 1;
    Sb(end + 1, :) = Curve_xy(k, :); %#ok<AGROW>
    isSameOctant = testOctant(Sb);
    if isSameOctant
        Sb_trans = transform_firstOctant(Sb, 1);
        Determine_xy(Sb_trans, 1);
    else
        break;
    end
end
bSegment = 1;
eSegment = k - 1;
MBS.data{end + 1} = Curve_xy(bSegment:eSegment, :);
MBS.B_E_seg(end + 1, :) = [bSegment eSegment];
while k < n
    while (w - 1) / max(abs(a), abs(b)) > Width || (~isSameOctant)
        bSegment = bSegment + 1;
        Sb = Curve_xy(bSegment:k, :);
        isSameOctant = testOctant(Sb);
        if isSameOctant
            Sb_trans = transform_firstOctant(Sb, 1);
            Determine_xy(Sb_trans, 1);
        else
            continue;
        end
    end
    incremental = 1;
    while (w - 1) / max(abs(a), abs(b)) <= Width && k < n
        k = k + 1;
        incremental = incremental + 1;
        Sb(end + 1, :) = Curve_xy(k, :); %#ok<AGROW>
        isSameOctant = testOctant(Sb);
        if isSameOctant
            Sb_trans = transform_firstOctant(Sb, 1);
            Determine_xy(Sb_trans, 1);
        else
            break;
        end
    end
    eSegment = k - 1;
    MBS.data{end + 1} = Curve_xy(bSegment:eSegment, :);
    MBS.B_E_seg(end + 1, :) = [bSegment eSegment];
    incremental = 1;
end
