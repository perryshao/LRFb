function integral = integral_invariant(Curve_xyz, Width, Radius)
% INTEGRAL_INVARIANT  Compute local integral invariants for a trajectory.
%   Curve_xyz is N-by-3; Width sets the blurred-segment thickness.
%   Radius is measured after coordinate-range normalization.
%   Returns N-by-2 descriptors and restores repeated-point rows.

[Curve_xyz, staindex] = remove_stapoint(Curve_xyz);
curve_xyz_int = curve_grid(Curve_xyz, 1000);
Curve_xyz = normalization(Curve_xyz);
if size(curve_xyz_int, 2) > 2
    MBS = splitting_curve_3D(curve_xyz_int, Width);
else
    MBS = splitting_curve_2D(curve_xyz_int, Width);
end
n = size(curve_xyz_int, 1);
m = size(MBS.data, 2);
integral = ones(n, 2) * 0.5;
% integral1 = zeros(n,2);
%% determine the Right and Left key point through MBS
for i = 1:m
    if i == 1
        E_num = 1;
    else
        E_num = MBS.B_E_seg(i - 1, 2) + 1;
    end

    for k = E_num:MBS.B_E_seg(i, 2)
        L(k) = MBS.B_E_seg(i, 1); %#ok<AGROW>
    end

    if i == m
        B_num = n;
    else
        B_num = MBS.B_E_seg(i + 1, 1) - 1;
    end
    for k = MBS.B_E_seg(i, 1):B_num
        R(k) = MBS.B_E_seg(i, 2); %#ok<AGROW>
    end
end

%% construct the Segment curve from MBS
% index= [MBS.B_E_seg(:,1);MBS.B_E_seg(:,2);];
% index = unique(index);

% Segment_curve = Curve_xyz(index,:);
% plot3d(Segment_curve);
% Segment_curve = interpolation(Segment_curve,n*2,0);
% Curve_xyz = Segment_curve;
% n = size(Curve_xyz,1);

%% plot the original curve and initial the plot parameters

% set(gcf,'position',[50,50,1152,864]);
% view(-37.5,30);
% plot3(Curve_xyz(:,1),Curve_xyz(:,2),Curve_xyz(:,3),'.k-','LineWidth',2);grid on;hold on;

%% dynamically plot the MBS procedure
% plot_MBS;
%% compute the integral invariants
for i = 2:(n - 2)
    integral(i, :) = estimate_integral(Curve_xyz, Curve_xyz(i, :), Curve_xyz(L(i), :), Curve_xyz(R(i), :), ...
        Radius, i); % Radius of the circumcircle to [CL(i),C(i),CR(i)]
end
integral = fillstapoint(integral, staindex);

%% compare the difference between two integral invariants by different approach.
% integral1 = fillstapoint(integral1,staindex);
% plot(integral(:,1),'*-b');hold on, plot(integral1(:,1),'*-r');
% figure(2),plot(integral(:,2),'*-b');hold on, plot(integral1(:,2),'*-r');
