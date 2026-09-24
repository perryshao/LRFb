function Determine_xy(Sb, incremental)
% DETERMINE_XY  Update the bounding-line parameters of a 2-D blurred segment.
%   Sb contains integer point coordinates; incremental selects cached updates.
%   The bounding line is returned through globals a, b, u and w.
%   This legacy 2-D path has unresolved helpers; see docs/CODE_REVIEW.md.

global a;
global b;
global u;
global w;
% global incremental;
Sb = double(Sb);
Sb = remove_stapoint(Sb); % remove the static points
n = size(Sb, 1);
if nargin > 1 && incremental > 4
    if exist('Determine.mat', 'file')
        load Determine_xy
        Determine_xy_sub;
        return;
    end
end

if n == 1
    return;
end

%% verify whether the curve is line at beginning
L = Sb(1, :);
U = Sb(1, :); % initialize the L and U with first point of Sb
for i = 1:n - 1
    vector1 = Sb(i, :);
    vector2 = Sb(i + 1, :);
    cross_vector = vector2(1) * vector1(2) - vector2(2) * vector1(1);
    if cross_vector == 0
        %       line_flag = true;
        P = Sb(i + 1, :) - Sb(i, :);
        g = gcd(P(2), P(1));
        P = P / g;
        a = P(2);
        b = P(1);
        w = 1;
        u = 0;
        L = Sb(i + 1, :);
        if i == n - 1
            return;
        end
    else
        %         Sb=Sb(i-1:end,:);
        Sb = [Sb(1, :); Sb(i:end, :)];
        break;
    end
end
%% reestimate the length of curve
n = size(Sb, 1);
if n < 3
    P = Sb(2, :) - Sb(1, :);
    g = gcd(P(2), P(1));
    P = P / g;
    a = P(2);
    b = P(1);
    w = 1;
    u = 0;
    clear Sb;
    save Determine_xy;
    return;
else
    [Hull_dequedata] = MelkmanConvexHull(Sb(1:3, :)); % initialize the Convex Hull using 3 points
    % Hull_dequedata = int32(Hull_dequedata);
end
%% begin intialize
% isSegment = true;
% a = 0; b = 1; u = 0; w = b;
Mc = Sb(1, :); %#ok<NASGU> % first point of Sb
% L = Mc; U = Mc; % initialize the L and U with first point of Sb
i = 1; % first point of Sb
while i < n;
    i = i + 1;
    Mc = Sb(i, :);
    if i < 3
        [Hull_dequedata] = MelkmanConvexHull(Sb(1:3, :));
    else
        [Hull_dequedata] = MelkmanConvexHull(Hull_dequedata, Mc);
        if ~ismember(Hull_dequedata(1, :), Mc, 'rows')
            continue;
        end
        if ~any(ismember(Hull_dequedata, L, 'rows'))
            vector1 = Hull_dequedata(1, :) - Hull_dequedata(2, :);
            vector2 = Hull_dequedata(1, :) - L;
            cross_vector1 = vector2(1) * vector1(2) - vector2(2) * vector1(1);
            vector1 = Hull_dequedata(end, :) - Hull_dequedata(end - 1, :);
            vector2 = Hull_dequedata(end, :) - L;
            cross_vector2 = vector2(1) * vector1(2) - vector2(2) * vector1(1);
            if cross_vector1 == 0 || cross_vector2 == 0
                L = Hull_dequedata(1, :);
            end
        end
        if ~any(ismember(Hull_dequedata, U, 'rows'))
            vector1 = Hull_dequedata(1, :) - Hull_dequedata(2, :);
            vector2 = Hull_dequedata(1, :) - L;
            cross_vector1 = vector2(1) * vector1(2) - vector2(2) * vector1(1);
            vector1 = Hull_dequedata(end, :) - Hull_dequedata(end - 1, :);
            vector2 = Hull_dequedata(end, :) - L;
            cross_vector2 = vector2(1) * vector1(2) - vector2(2) * vector1(1);
            if cross_vector1 == 0 || cross_vector2 == 0
                U = Hull_dequedata(end, :);
            end
        end
        if ~any(ismember(Hull_dequedata, L, 'rows')) && ...
           ~any(ismember(Hull_dequedata, U, 'rows'))

            U = Hull_dequedata(end, :);
            L = Hull_dequedata(1, :);
        end
    end
    Hull_size = size(Hull_dequedata, 1);

    r = a * Mc(1) - b * Mc(2);
    if r == u
        U = Mc;
    end
    if r == u + w - 1
        L = Mc;
    end
    if r <= u - 1
        U = Mc;
        % Let the N the point before M in the upper convex hull
        if i < 3
            Nc = Sb(i - 1, :);
        else
            if (sign(b) < 0 && sign(a) < 0)
                Nc = Hull_dequedata(2, :); %% perry
            else
                Nc = Hull_dequedata(end - 1, :);
            end
        end
        a0 = Mc(2) - Nc(2);
        b0 = Mc(1) - Nc(1);
        a = a0 / gcd(a0, b0);
        b = b0 / gcd(a0, b0);
        u = a * Mc(1) - b * Mc(2);
        % find the first point C in the lower part of the convex hull
        % starting at L such that : slopt of [C,Cnext]>a/b
        if i >= 3
            L_index = ismember(Hull_dequedata, L, 'rows');
            L_index = find(L_index == 1);
            for j = L_index(1):-1:2
                C = Hull_dequedata(j, :);
                C_next = Hull_dequedata(j - 1, :);
                slope = C_next - C;
                slope = slope(2) / slope(1);
                if slope >= (a / b)
                    L = C;
                    break;
                end
            end
        end
        w = a * L(1) - b * L(2) - u + 1;
    elseif r >= u + w - 1
        L = Mc;
        % Let N the point before M in the lower convex hull
        if i < 3
            Nc = Sb(i - 1, :);
        else
            Nc = Hull_dequedata(1 + 1, :);
        end
        a0 = Mc(2) - Nc(2);
        b0 = Mc(1) - Nc(1);
        a = a0 / gcd(a0, b0);
        b = b0 / gcd(a0, b0);
        % find the first point C in the upper part of the convex hull
        % starting at U
        if i >= 3
            U_index = ismember(Hull_dequedata, U, 'rows');
            U_index = find(U_index == 1);

            for j = U_index(1):(Hull_size - 1)
                C = Hull_dequedata(j, :);
                C_next = Hull_dequedata(j + 1, :);
                slope = C_next - C;
                slope = slope(2) / slope(1);
                if slope < (a / b)
                    U = C;
                    break;
                end
            end
        end

        u = a * U(1) - b * U(2);
        w = a * Mc(1) - b * Mc(2) - u + 1;
    end
end
% save Determine_xy L U Hull_dequedata
