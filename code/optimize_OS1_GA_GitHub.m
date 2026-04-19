function result = optimize_OS1_GA_GitHub(cfg)
% OPTIMIZE_VLP16_GA_NEW
% Replaces the internal engine of the BO script with Genetic Algorithm.
% Retains robust GUI integration, plotting, and export features.

if nargin < 1
    cfg = struct();
end

% --- GUI / Plotting Setup ---
if isfield(cfg,'AxesHandle')
    ax = cfg.AxesHandle;
else
    ax = [];
end

if isfield(cfg,'EmbedPlots')
    embedPlots = cfg.EmbedPlots;
else
    embedPlots = false;
end

% --- General Settings (Same as BO) ---
clc
% close all % Don't close all if running from GUI
rng(42);

%% -------------------- Toggle knobs --------------------
SYMMETRY.ENABLE = true;
SYMMETRY.PLANE  = 'YZ';

FIX_ROLL.ENABLE = true;
FIX_ROLL.DEG    = 0;

%% -------------------- Scene / backpack path --------------------
L = 50; H = 20;
box = make_open_box(L, H);
BACKPACK_RPY_DEG = [0 0 0];
R_backpack = Rxyz_deg(BACKPACK_RPY_DEG);

MOVE_AXIS   = 'x';
Z_BACKPACK  = 1.60;
PATH_LEN    = 40;
PATH_STEP   = 1.0;
MARGIN      = 5.0;

%% -------------------- OS1-64 scan model --------------------
OS1_VERT_FOV_DEG    = 42.4;   % approx vertical FOV
OS1_VERT_BEAMS      = 64;     % channels
OS1_AZ_SAMPLES      = 180;    % rays per pose (speed knob)
OS1_SCANS           = 1;      % slices per pose (speed knob)
BEAM_DIVERGENCE_DEG = 0.18;   % small jitter
MAX_RANGE           = 120;    % practical range

% Approximate vertical channel angles (replace with real OS1 channel angles if available)
vert_angles = linspace(-OS1_VERT_FOV_DEG/2, +OS1_VERT_FOV_DEG/2, OS1_VERT_BEAMS);

%% -------------------- Mount bounds --------------------
b.x     = [-0.10 0.15];
b.y     = [-0.25 0.25];
b.z     = [-0.15 0.15];
b.pitch = [-90  90];
b.roll  = [-90  90];

%% -------------------- Objective & Settings --------------------
VOX                  = 0.5;
OVERLAP_MIN          = 0.05;
OVERLAP_MAX          = 0.10;
MIN_ANGULAR_DIV      = 30;
LAMBDA               = 2000;
BEAM_DIV_BO          = 0; % Deterministic for optimization
MASTERS_MAX          = 5; 

% --- Apply GUI Overrides ---
if isfield(cfg,'VOX'),             VOX = cfg.VOX; end
if isfield(cfg,'LAMBDA'),          LAMBDA = cfg.LAMBDA; end
if isfield(cfg,'OVERLAP_MIN'),     OVERLAP_MIN = cfg.OVERLAP_MIN; end
if isfield(cfg,'OVERLAP_MAX'),     OVERLAP_MAX = cfg.OVERLAP_MAX; end
if isfield(cfg,'MIN_ANGULAR_DIV'), MIN_ANGULAR_DIV = cfg.MIN_ANGULAR_DIV; end
if isfield(cfg,'PATH_LEN'),        PATH_LEN = cfg.PATH_LEN; end
if isfield(cfg,'PATH_STEP'),       PATH_STEP = cfg.PATH_STEP; end
if isfield(cfg,'MOVE_AXIS'),       MOVE_AXIS = cfg.MOVE_AXIS; end
if isfield(cfg,'SYMMETRY_ENABLE'), SYMMETRY.ENABLE = cfg.SYMMETRY_ENABLE; end
if isfield(cfg,'FIX_ROLL_ENABLE'), FIX_ROLL.ENABLE = cfg.FIX_ROLL_ENABLE; end

% --- Build Scene Helpers ---
meta = make_vox_meta(box, VOX);
T_B_list = build_straight_path(box, R_backpack, Z_BACKPACK, MOVE_AXIS, PATH_LEN, PATH_STEP, MARGIN);

%% -------------------- GA SETUP --------------------
% Vector Definition: 
% x = [ N_masters, (x1 y1 z1 p1 r1), ..., (xK yK zK pK rK) ]
% K = MASTERS_MAX

nvars = 1 + 5*MASTERS_MAX;

% Bounds
% First var is N_masters (1..MASTERS_MAX)
% Subsequent vars are repeating blocks of mount bounds
lb = [1,           repmat([b.x(1) b.y(1) b.z(1) b.pitch(1) b.roll(1)], 1, MASTERS_MAX)];
ub = [MASTERS_MAX, repmat([b.x(2) b.y(2) b.z(2) b.pitch(2) b.roll(2)], 1, MASTERS_MAX)];

% Integer Constraints (N_masters is integer)
IntCon = 1; 

% Helper to convert vector 'x' to 'tbl' for the objective function
vec2tbl = @(x) local_vec_to_tbl(x, MASTERS_MAX);

% Objective Function
% Note: evaluate_layout_BO returns -J (negative). GA minimizes.
% So minimizing -J maximizes J. This is consistent.
f_ga = @(x) evaluate_layout_BO( ...
    vec2tbl(x), meta, box, T_B_list, MAX_RANGE, vert_angles, ...
    OS1_AZ_SAMPLES, OS1_SCANS, BEAM_DIV_BO, ...
    OVERLAP_MIN, OVERLAP_MAX, MIN_ANGULAR_DIV, LAMBDA, ...
    SYMMETRY, FIX_ROLL);


% GA Options
gaOpts = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', 80, ...
    'MaxGenerations', 100, ... % Tuned for speed/quality balance
    'UseParallel', false, ...
    'FunctionTolerance', 1e-4);

% Hybrid Refinement (Pattern Search)
gaOpts = optimoptions(gaOpts, 'HybridFcn', @patternsearch);

% --- Plotting Integration ---
if isfield(cfg,'ConvAxes') && ~isempty(cfg.ConvAxes)
    % GUI Mode: Use custom output function to draw into specific axes
    gaOpts = optimoptions(gaOpts, 'OutputFcn', @(opt,state,flag) ga_custom_plot_gui(opt,state,flag,cfg.ConvAxes));
else
    % Script Mode: Standard plot function
    gaOpts = optimoptions(gaOpts, 'PlotFcn', {@gaplotbestf});
end

%% -------------------- RUN GA --------------------
tic;
fprintf('Starting GA Optimization...\n');
[xbest, fmin, exitflag, output] = ga(f_ga, nvars, [], [], [], [], lb, ub, [], IntCon, gaOpts);

% Convert best vector back to table for post-processing
bestTbl = vec2tbl(xbest);

%% -------------------- Rebuild Winner (Post-Processing) --------------------
% This section is identical to the trusted BO code 
[theta_list, mirrors_used] = expand_masters_to_full_set(bestTbl, SYMMETRY, FIX_ROLL);
N  = numel(theta_list);
mounts = cell(1,N);
perSensorSets = cell(1,N);
global_cov = false(meta.n,1);

for s = 1:N
    th = theta_list{s};
    M_local = theta_to_mount_local_VLP(th);
    % Re-simulate with actual divergence for final check
cov_s = simulate_vlp16_bitmap( ...
    M_local, T_B_list, box, MAX_RANGE, ...
    vert_angles, OS1_AZ_SAMPLES, OS1_SCANS, BEAM_DIVERGENCE_DEG, meta);

    perSensorSets{s} = bitmap_to_keys(meta, cov_s);
    global_cov = global_cov | cov_s;
    mounts{s} = M_local;
end

[unionKeys_final, countsMap] = union_from_per_sensor(perSensorSets);
total_union = numel(unionKeys_final);
overlap_ge2 = 0;
if ~isempty(countsMap)
    Vc = values(countsMap);
    overlap_ge2 = sum(cellfun(@(x) x>=2, Vc));
end
overlap_ratio_final = overlap_ge2 / max(1,total_union);

% Angular diversity on winner
min_ang_all = Inf;
if N >= 2
    for i = 1:N-1
        zi = mounts{i}(1:3,3);
        zi = zi/norm(zi);
        for j = i+1:N
            zj = mounts{j}(1:3,3);
            zj = zj/norm(zj);
            c  = max(-1,min(1,zi.'*zj));
            a  = acosd(c);
            if a < min_ang_all, min_ang_all = a; end
        end
    end
end

penalty_overlap = overlap_reward(overlap_ratio_final, OVERLAP_MIN, OVERLAP_MAX);
penalty_angular = angular_reward(min_ang_all, MIN_ANGULAR_DIV);
J = total_union + penalty_overlap + penalty_angular - LAMBDA * N;

fprintf('→ Winner: N=%d (sym=%d) | union=%d | overlap≈%.1f%% | minAng=%.1f° | J=%.1f\n', ...
    N, mirrors_used, total_union, 100*overlap_ratio_final, min_ang_all, J);
history = [N, total_union, J];

%% -------------------- JSON Export --------------------
out.sensor_type = 'Ouster OS1-64';
out.symmetry.enabled = SYMMETRY.ENABLE;
out.symmetry.plane   = SYMMETRY.PLANE;
out.roll.fixed       = FIX_ROLL.ENABLE;
out.roll.fixed_deg   = FIX_ROLL.DEG;
out.backpack.path_axis   = MOVE_AXIS;
out.backpack.path_length = PATH_LEN;
out.backpack.path_step   = PATH_STEP;
out.backpack.margin      = MARGIN;
out.optimization.total_coverage = total_union;
out.optimization.lambda         = LAMBDA;
out.optimization.N              = N;
out.optimization.method         = 'GA';

thetas = zeros(N,5);
for i=1:N
    M = mounts{i};
    out.mounts(i).position  = M(1:3,4).';
    out.mounts(i).euler_deg = mat2eulxyz_deg(M(1:3,1:3));
    thetas(i,:) = pose_to_theta(M);
end
out.trajectory_world = cellfun(@(T) struct('R', T(1:3,1:3), 't', T(1:3,4).'), T_B_list, 'uni', false);
json = jsonencode(out, 'PrettyPrint', true);
fid = fopen('C:................................\GA_best_mounts_os1_64_symm.json','w');
fwrite(fid,json); fclose(fid);
fprintf('✓ Exported JSON: GA_best_mounts_os1_64_symm.json\n');
toc

%% -------------------- Visualization --------------------
if embedPlots && ~isempty(ax)
    plot_multi_setup_with_path_VLP16(box, mounts, perSensorSets, meta, T_B_list, OS1_VERT_FOV_DEG, ax);
else
    plot_multi_setup_with_path_VLP16(box, mounts, perSensorSets, meta, T_B_list, OS1_VERT_FOV_DEG);

end
drawnow; shg;

%% -------------------- Per-sensor CSV --------------------
[uniq_counts, overlap_count, total_union_chk, contrib_counts] = coverage_stats(perSensorSets);
names = {'id','x','y','z','roll_deg','pitch_deg','yaw_deg','contrib','unique_only','contrib_pct'};
rows = cell(N,1);
for i = 1:N
    M = mounts{i};
    pos = M(1:3,4).';
    eul = mat2eulxyz_deg(M(1:3,1:3));
    rows{i} = [i, pos, eul, contrib_counts(i), uniq_counts(i), ...
               100*contrib_counts(i)/max(1,total_union)];
end
T = array2table(cell2mat(rows), 'VariableNames', names);
writetable(T, 'mount_summary_os1_64_GA.csv');
fprintf('✓ Saved CSV: mount_summary_os1_64_GA.csv\n');


% Return struct
result = struct('N',N,'total_unique',total_union,'objective',J,'lambda',LAMBDA, ...
    'history',history,'mount_thetas',thetas,'mounts_local',{mounts}, ...
    'trajectory',{T_B_list});

end

%% ====================== GA Helpers ======================

function tbl = local_vec_to_tbl(x, K)
    % Converts GA vector back to Table for the Objective Function
    Nm = max(1, min(K, round(x(1))));
    tbl = table; tbl.N_masters = Nm;
    ptr = 2;
    for s = 1:K
        xi = x(ptr:ptr+4);
        ptr = ptr + 5;
        tbl.(sprintf('xM_%d',s)) = xi(1);
        tbl.(sprintf('yM_%d',s)) = xi(2);
        tbl.(sprintf('zM_%d',s)) = xi(3);
        tbl.(sprintf('pM_%d',s)) = xi(4);
        tbl.(sprintf('rM_%d',s)) = xi(5);
    end
end

function [state, options, optchanged] = ga_custom_plot_gui(options, state, flag, axConv)
    % Output function to draw GA convergence into the GUI axes
    optchanged = false;
    if isempty(axConv) || ~isvalid(axConv), return; end
    
    switch flag
        case 'init'
            cla(axConv);
            hold(axConv, 'on');
            grid(axConv, 'on');
            xlabel(axConv, 'Generation');
            ylabel(axConv, 'Penalty (Neg. Objective)');
            title(axConv, 'GA Convergence');
        case 'iter'
            % Get best score history
            % 'state.Best' is usually a vector of best scores per generation
            if ~isempty(state.Best)
                gens = 0:(length(state.Best)-1);
                plot(axConv, gens, state.Best, 'b-o', 'LineWidth', 1.5);
                drawnow limitrate;
            end
        case 'done'
            % Final update
    end
end

%% ====================== Shared Objectives & Helpers ======================
% (Below are the exact same helpers as in the BO code)

function f = evaluate_layout_BO( ...
    tbl, meta, box, T_B_list, MAX_RANGE, vert_angles, ...
    VLP_AZ_SAMPLES, VLP_SCANS, BEAM_DIV_BO, ...
    OVERLAP_MIN, OVERLAP_MAX, MIN_ANGULAR_DIV, LAMBDA, ...
    SYMMETRY, FIX_ROLL)

    % Expand masters -> full sensor list
    [theta_list, ~] = expand_masters_to_full_set(tbl, SYMMETRY, FIX_ROLL);
    Nact = numel(theta_list);

    % Simulate
    cov_list = cell(Nact,1);
    M_list   = cell(Nact,1);
    for s2 = 1:Nact
        th = theta_list{s2};
        Mloc = theta_to_mount_local_VLP(th);
        cov_s = simulate_vlp16_bitmap( ...
            Mloc, T_B_list, box, MAX_RANGE, ...
            vert_angles, VLP_AZ_SAMPLES, VLP_SCANS, BEAM_DIV_BO, meta);
        cov_list{s2} = cov_s;
        M_list{s2}   = Mloc;
    end

    % Union + counts
    if isempty(cov_list)
        union_bits = false(meta.n,1);
        counts = zeros(meta.n,1,'uint8');
    else
        union_bits = cov_list{1};
        counts     = uint8(cov_list{1});
        for t = 2:numel(cov_list)
            union_bits = union_bits | cov_list{t};
            counts     = counts + uint8(cov_list{t});
        end
    end
    total_union   = sum(union_bits);
    overlap_ratio = sum(counts>=2) / max(1,total_union);

    % Angular diversity
    if Nact < 2
        min_ang = Inf;
    else
        min_ang = Inf;
        for i2=1:Nact
            zi = M_list{i2}(1:3,3)/norm(M_list{i2}(1:3,3));
            for j2=i2+1:Nact
                zj = M_list{j2}(1:3,3)/norm(M_list{j2}(1:3,3));
                c  = max(-1,min(1,dot(zi,zj)));
                a  = acosd(c);
                if a < min_ang, min_ang = a; end
            end
        end
    end

    R_overlap  = overlap_reward(overlap_ratio, OVERLAP_MIN, OVERLAP_MAX);
    R_angular  = angular_reward(min_ang, MIN_ANGULAR_DIV);
    Jbo        = double(total_union) + R_overlap + R_angular - LAMBDA * double(Nact);
    
    % GA minimizes, so return negative J
    f = -Jbo;
end

function [theta_list, mirrored] = expand_masters_to_full_set(tbl, SYM, FIXR)
    Nm = tbl.N_masters;
    masters = cell(1,Nm);
    for s = 1:Nm
        th = [tbl.(sprintf('xM_%d',s)), ...
              tbl.(sprintf('yM_%d',s)), ...
              tbl.(sprintf('zM_%d',s)), ...
              tbl.(sprintf('pM_%d',s)), ...
              tbl.(sprintf('rM_%d',s))];
        if FIXR.ENABLE, th(5) = FIXR.DEG; end
        masters{s} = th;
    end

    if ~SYM.ENABLE
        theta_list = masters;
        mirrored = 0;
        return;
    end

    theta_list = cell(1, 2*Nm);
    k = 0;
    for s = 1:Nm
        thM = masters{s};
        % Mirror across YZ plane
        thS = thM;
        thS(1) = -thM(1); % mirror X
        thS(4) = -thM(4); % flip pitch
        if FIXR.ENABLE
            thS(5) = FIXR.DEG;
        else
            thS(5) = thM(5);
        end
        k = k+1; theta_list{k} = thM;
        k = k+1; theta_list{k} = thS;
    end
    mirrored = Nm;
end

% --- Geometry / Sim / Vox helpers ---
function T_B_list = build_straight_path(box, R_backpack, z_backpack, axis, path_len, step, margin)
    switch lower(axis)
        case 'x', xs = (-path_len/2):step:(path_len/2); ys = zeros(size(xs));
        case 'y', ys = (-path_len/2):step:(path_len/2); xs = zeros(size(ys));
        otherwise, error('MOVE_AXIS must be x or y.');
    end
    xs = max(xs, box.xmin+margin); xs = min(xs, box.xmax-margin);
    ys = max(ys, box.ymin+margin); ys = min(ys, box.ymax-margin);
    T_B_list = cell(numel(xs),1);
    for i=1:numel(xs)
        pB = [xs(i), ys(i), z_backpack]';
        T_B_list{i} = [R_backpack, pB; 0 0 0 1];
    end
end

function M = theta_to_mount_local_VLP(theta)
    x=theta(1); y=theta(2); z=theta(3); pitch=theta(4); roll=theta(5);
    Rrp = rotx(roll) * roty(pitch);
    M = [Rrp, [x;y;z]; 0 0 0 1];
end

function cov = simulate_vlp16_bitmap(M_local, T_B_list, box, maxR, vert_angles_deg, az_samp, scans, div_deg, meta)
    cov = false(meta.n,1);
    nv = numel(vert_angles_deg);
    az_list = linspace(0, 360*(1-1/az_samp), az_samp);
    D_local = zeros(3, az_samp*nv);
    idx = 0;
    for a = 1:az_samp
        for v = 1:nv
            idx = idx+1;
            D_local(:,idx) = os1_angles_to_dir(az_list(a), vert_angles_deg(v));
        end
    end
    epsi = max(0.25*meta.vox, 1e-3);
    for j = 1:numel(T_B_list)
        T_B = T_B_list{j};
        M_world = T_B * M_local;
        Rw = M_world(1:3,1:3);
        tw = M_world(1:3,4);
        D_world = Rw * D_local;
        D_world = D_world ./ vecnorm(D_world);
        for s = 1:scans
            phase_cols = circshift(D_world, s-1, 2);
            if div_deg > 0
                alpha = deg2rad(div_deg)*0.5;
                noise = alpha * randn(size(phase_cols));
                phase_cols = phase_cols + noise;
                phase_cols = phase_cols ./ vecnorm(phase_cols);
            end
            [hitMask, tnear] = ray_box_vec(tw, phase_cols, box, maxR);
            if ~any(hitMask), continue; end
            P = tw + phase_cols(:,hitMask) .* tnear(hitMask);
            P(1,:) = min(max(P(1,:), box.xmin+epsi), box.xmax-epsi);
            P(2,:) = min(max(P(2,:), box.ymin+epsi), box.ymax-epsi);
            P(3,:) = min(max(P(3,:), box.zmin+epsi), box.zmax-epsi);
            lin = vox_lin_index(P, meta);
            cov(lin(lin>0)) = true;
        end
    end
end

function [hitMask, tnear] = ray_box_vec(p, D, box, maxR)
    px = p(1); py = p(2); pz = p(3);
    dx = D(1,:); dy = D(2,:); dz = D(3,:);
    eps_div = 1e-12;
    invDx = 1./max(abs(dx),eps_div) .* sign(dx);
    invDy = 1./max(abs(dy),eps_div) .* sign(dy);
    invDz = 1./max(abs(dz),eps_div) .* sign(dz);
    tx1 = (box.xmin - px) .* invDx;   tx2 = (box.xmax - px) .* invDx;
    ty1 = (box.ymin - py) .* invDy;   ty2 = (box.ymax - py) .* invDy;
    tz1 = (box.zmin - pz) .* invDz;   tz2 = (box.zmax - pz) .* invDz;
    tmin = max([min(tx1,tx2); min(ty1,ty2); min(tz1,tz2)], [], 1);
    tmax = min([max(tx1,tx2); max(ty1,ty2); max(tz1,tz2)], [], 1);
    hitMask = (tmax >= tmin) & (tmax > 0) & (tmin <= maxR);
    tnear = tmin;
    inside = (tmin < 0);
    tnear(inside) = tmax(inside);
    tnear = max(tnear, 1e-6);
    tnear(~hitMask) = 0;
end

function lin = vox_lin_index(p, m)
    i = floor((p(1,:))/m.vox) - m.i0;
    j = floor((p(2,:))/m.vox) - m.j0;
    k = floor((p(3,:))/m.vox) - m.k0;
    ok = (i>=0 & i<m.nx & j>=0 & j<m.ny & k>=0 & k<m.nz);
    lin = zeros(1,numel(i),'uint32');
    lin(ok) = uint32( double(i(ok)) + double(m.nx)*( double(j(ok)) + double(m.ny)*double(k(ok)) ) + 1 );
end

function meta = make_vox_meta(box, vox)
    meta.vox = vox;
    meta.i0 = floor((box.xmin)/vox);
    meta.j0 = floor((box.ymin)/vox);
    meta.k0 = floor((box.zmin)/vox);
    meta.nx = floor((box.xmax - box.xmin)/vox);
    meta.ny = floor((box.ymax - box.ymin)/vox);
    meta.nz = floor((box.zmax - box.zmin)/vox);
    meta.n  = meta.nx * meta.ny * meta.nz;
end

function keys = bitmap_to_keys(m, cov)
    idx = find(cov(:)).' - 1;
    if isempty(idx), keys = {}; return; end
    [i, rem] = ind2sub_vec([m.nx, m.ny, m.nz], idx);
    j = rem(:,1);
    k = rem(:,2);
    keys = cell(1,numel(i));
    for t=1:numel(i), keys{t} = sprintf('%d_%d_%d', i(t), j(t), k(t)); end
end

function [i, rest] = ind2sub_vec(sz, idx0)
    nx = sz(1); ny = sz(2);
    k = floor(idx0 / (nx*ny));
    rem2 = idx0 - k*(nx*ny);
    j = floor(rem2 / nx);      i = rem2 - j*nx;
    rest = [j(:) k(:)];
    i = i(:).';
end

function box = make_open_box(L, H)
    box.L = L; box.H = H;
    box.xmin = -L/2; box.xmax = L/2;
    box.ymin = -L/2; box.ymax = L/2;
    box.zmin = 0;     box.zmax = H;
    box.min = [box.xmin, box.ymin, box.zmin];
    box.max = [box.xmax, box.ymax, box.zmax];
end

function R = Rxyz_deg(rpy)
    r = deg2rad(rpy(1)); p = deg2rad(rpy(2)); y = deg2rad(rpy(3));
    Rx = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
    Ry = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
    Rz = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1];
    R = Rz * Ry * Rx;
end

function eul = mat2eulxyz_deg(R)
    pitch = atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2));
    roll  = atan2(R(3,2), R(3,3));
    yaw   = atan2(R(2,1), R(1,1));
    eul = rad2deg([roll, pitch, yaw]);
end

function d = os1_angles_to_dir(az_deg, el_deg)
    az = deg2rad(az_deg); el = deg2rad(el_deg);
    x = cos(el) * cos(az); y = cos(el) * sin(az); z = sin(el);
    d = [x; y; z] / norm([x y z]);
end

function [unionKeys, counts] = union_from_per_sensor(perSensorSets)
    counts = containers.Map('KeyType','char','ValueType','int32');
    ALL = {};
    for s = 1:numel(perSensorSets)
        ks = perSensorSets{s};
        ALL = [ALL, ks];
        for i = 1:numel(ks)
            k = ks{i};
            if isKey(counts,k), counts(k) = counts(k)+1; else, counts(k) = 1; end
        end
    end
    unionKeys = unique(ALL);
end

function [uniq_counts, overlap_count, total_union, contrib_counts] = coverage_stats(perSensorSets)
    S = numel(perSensorSets);
    contrib_counts = zeros(S,1);
    counts = containers.Map('KeyType','char','ValueType','int32');
    for s = 1:S
        ks = perSensorSets{s};
        contrib_counts(s) = numel(ks);
        for i = 1:numel(ks)
            k = ks{i};
            if isKey(counts,k), counts(k) = counts(k)+1; else, counts(k) = 1; end
        end
    end
    uniq_counts = zeros(S,1);
    if S>0 && ~isempty(counts)
        unique_owner = containers.Map('KeyType','char','ValueType','int32');
        K = keys(counts);
        for i=1:numel(K)
            k = K{i};
            if counts(k)==1, unique_owner(k) = 0; end
        end
        for s=1:S
            ks = perSensorSets{s};
            for i=1:numel(ks)
                k = ks{i};
                if isKey(unique_owner,k), unique_owner(k) = s; end
            end
        end
        U = values(unique_owner);
        for i=1:numel(U)
            owner = U{i};
            if owner>=1 && owner<=S, uniq_counts(owner) = uniq_counts(owner)+1; end
        end
    end
    if isempty(counts)
        overlap_count = 0; total_union = 0; return;
    end
    K = keys(counts);
    total_union = numel(K);
    overlap_count = 0;
    for i=1:total_union
        if counts(K{i}) >= 2, overlap_count = overlap_count + 1; end
    end
end

function P = keys_to_points_from_meta(K, meta)
    n = numel(K); P = zeros(n,3);
    io = double(meta.i0); jo = double(meta.j0); ko = double(meta.k0);
    v  = double(meta.vox);
    for t = 1:n
        ijk = sscanf(K{t}, '%d_%d_%d');
        P(t,1) = (double(ijk(1)) + io + 0.5) * v;
        P(t,2) = (double(ijk(2)) + jo + 0.5) * v;
        P(t,3) = (double(ijk(3)) + ko + 0.5) * v;
    end
end

function S = classify_hits(P, box, vox)
    tol = max(vox*2.0, 0.5);
    S = struct();
    if isempty(P)
        S.floor=false(0,1); S.ceil=false(0,1);
        S.xmin=false(0,1); S.xmax=false(0,1); S.ymin=false(0,1); S.ymax=false(0,1);
        S.walls=false(0,1); return;
    end
    S.floor = abs(P(:,3) - box.zmin) <= tol;
    S.ceil  = abs(P(:,3) - box.zmax) <= tol;
    S.xmin  = abs(P(:,1) - box.xmin) <= tol & ~S.floor & ~S.ceil;
    S.xmax  = abs(P(:,1) - box.xmax) <= tol & ~S.floor & ~S.ceil;
    S.ymin  = abs(P(:,2) - box.ymin) <= tol & ~S.floor & ~S.ceil;
    S.ymax  = abs(P(:,2) - box.ymax) <= tol & ~S.floor & ~S.ceil;
    S.walls = S.xmin | S.xmax | S.ymin | S.ymax;
end

function Pw = inset_wall_points(P, S, box, vox)
    eps = max(vox*0.25, 0.1);
    Pw = P;
    if ~isempty(P)
        Pw(S.xmin,1) = box.xmin + eps;
        Pw(S.xmax,1) = box.xmax - eps;
        Pw(S.ymin,2) = box.ymin + eps;
        Pw(S.ymax,2) = box.ymax - eps;
    end
end

function Pf = inset_floor_points(P, box, vox)
    if isempty(P), Pf = P; return; end
    eps = max(vox*0.25, 0.1);
    Pf = P; Pf(:,3) = max(P(:,3), box.zmin + eps);
end

function Pc = inset_ceil_points(P, box, vox)
    if isempty(P), Pc = P; return; end
    eps = max(vox*0.25, 0.1);
    Pc = P; Pc(:,3) = min(P(:,3), box.zmax - eps);
end

function plot_multi_setup_with_path_VLP16(box, mounts_local, perSensorSets, vox_meta, T_B_list, vFOVdeg, ax)
    if nargin < 7 || isempty(ax)
        fig = figure('Position',[180 80 1280 900]);
        ax  = axes('Parent', fig);
    end
    
    cla(ax);
    hold(ax,'on');
    axis(ax,'equal');
    grid(ax,'off');
    ax.XGrid = 'off'; ax.YGrid = 'off'; ax.ZGrid = 'off';
    ax.XTick = []; ax.YTick = []; ax.ZTick = [];
   xlabel(ax, 'X (m)', 'FontSize', 14, 'FontWeight', 'bold');
   ylabel(ax, 'Y (m)', 'FontSize', 14, 'FontWeight', 'bold');
   zlabel(ax, 'Z (m)', 'FontSize', 14, 'FontWeight', 'bold');
    title(ax,'Multi-OS1-64 Coverage (GA Result)');

    ax.Color = 'none';
    view(ax,45,25);

    vox = vox_meta.vox;
    S   = numel(mounts_local);
    havePath = ~isempty(T_B_list);

    % Path
    if havePath
        P = zeros(numel(T_B_list),3);
        for i=1:numel(T_B_list), P(i,:) = T_B_list{i}(1:3,4)'; end
        plot3(ax,P(:,1),P(:,2),P(:,3),'-','LineWidth',1.8,'Color',[0 0 0 0.35]);
        scatter3(ax,P(1,1),P(1,2),P(1,3),40,[0 0.5 1],'filled');
        scatter3(ax,P(end,1),P(end,2),P(end,3),40,[1 0.2 0.2],'filled');
    end

    % Sensors at start pose
    colors = lines(max(S,7));
    if havePath
        for s=1:S
            M_world = T_B_list{1} * mounts_local{s};
            draw_sensor_glyph(M_world, 0.5*colors(mod(s-1,size(colors,1))+1,:), sprintf('Sensor %d',s));
            draw_fov_cone(remove_spin_yaw(M_world), vFOVdeg, 15, colors(mod(s-1,size(colors,1))+1,:), true);
        end
    end

    % Coverage points
    cnt = containers.Map('KeyType','char','ValueType','int32');
    for s=1:S
        ks = perSensorSets{s};
        for i=1:numel(ks)
            k = ks{i};
            if isKey(cnt,k), cnt(k) = cnt(k)+1; else, cnt(k) = 1; end
        end
    end

    % Per-sensor uniques
    for s=1:S
        ks = perSensorSets{s};
        if isempty(ks), continue; end
        uniq = {};
        for i=1:numel(ks)
            k = ks{i};
            if cnt(k)==1, uniq{end+1} = k; end
        end
        if isempty(uniq), continue; end
        Puniq = keys_to_points_from_meta(uniq, vox_meta);
        Scls  = classify_hits(Puniq, box, vox);
        Pw    = inset_wall_points(Puniq, Scls, box, vox);
        Pf    = inset_floor_points(Puniq, box, vox);
        c = colors(mod(s-1,size(colors,1))+1,:);
        if any(Scls.walls)
            scatter3(ax,Pw(Scls.walls,1), Pw(Scls.walls,2), Pw(Scls.walls,3), 16, c, 'filled', 'MarkerFaceAlpha', 0.9);
        end
        if any(Scls.floor)
            scatter3(ax,Pf(Scls.floor,1), Pf(Scls.floor,2), Pf(Scls.floor,3), 8, c, 'filled', 'MarkerFaceAlpha', 0.25, 'HandleVisibility','off');
        end
        if any(Scls.ceil)
            Pc = inset_ceil_points(Puniq(Scls.ceil,:), box, vox);
            scatter3(ax,Pc(:,1), Pc(:,2), Pc(:,3), 8, c, 'filled', 'MarkerFaceAlpha', 0.25);
        end
    end

    % Overlap
    overlapK = {};
    if ~isempty(cnt)
        K = keys(cnt); vv = values(cnt);
        mask = cellfun(@(x) x>=2, vv);
        overlapK = K(mask);
    end
    if ~isempty(overlapK)
        Pover = keys_to_points_from_meta(overlapK, vox_meta);
        Sov   = classify_hits(Pover, box, vox);
        Pwo   = inset_wall_points(Pover, Sov, box, vox);
        Pfo   = inset_floor_points(Pover, box, vox);
        if any(Sov.walls)
            scatter3(ax,Pwo(Sov.walls,1), Pwo(Sov.walls,2), Pwo(Sov.walls,3), 20, [0.55 0.2 0.75], 'd', 'filled', 'MarkerFaceAlpha', 0.9, 'HandleVisibility','off');
        end
        if any(Sov.floor)
            scatter3(ax,Pfo(Sov.floor,1), Pfo(Sov.floor,2), Pfo(Sov.floor,3), 10, [0.55 0.2 0.75], 'd', 'filled', 'HandleVisibility','off');
        end
        if any(Sov.ceil)
            Pc = inset_ceil_points(Pover(Sov.ceil,:), box, vox);
            scatter3(ax,Pc(:,1), Pc(:,2), Pc(:,3), 10, [0.55 0.2 0.75], 'd', 'filled', 'HandleVisibility','off');
        end
    end
    xlim(ax,[box.xmin box.xmax]); ylim(ax,[box.ymin box.ymax]); zlim(ax,[box.zmin box.zmax]);
    if ~isa(ax,'matlab.ui.control.UIAxes')
        camlight(ax,'headlight'); lighting(ax,'gouraud');
    end

    function draw_sensor_glyph(T, color, label) %#ok<INUSD>
        r = 0.06; len = 0.12; nTheta = 36; nX = 12;
        theta = linspace(0, 2*pi, nTheta); x = linspace(-len/2, len/2, nX);
        [TH, XX] = meshgrid(theta, x);
        YY = r * cos(TH); ZZ = r * sin(TH);
        P = [XX(:), YY(:), ZZ(:), ones(numel(XX),1)]';
        PW = T * P;
        Xw = reshape(PW(1,:), size(XX)); Yw = reshape(PW(2,:), size(XX)); Zw = reshape(PW(3,:), size(XX));
        surf(ax,Xw, Yw, Zw, 'FaceColor', color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        origin = T(1:3,4); axlen = 0.35;
        quiver3(ax, origin(1),origin(2),origin(3), T(1,1)*axlen, T(2,1)*axlen, T(3,1)*axlen, 'r','LineWidth',1.8,'MaxHeadSize',0.3,'HandleVisibility','off');
        quiver3(ax, origin(1),origin(2),origin(3), T(1,2)*axlen, T(2,2)*axlen, T(3,2)*axlen, 'g','LineWidth',1.8,'MaxHeadSize',0.3,'HandleVisibility','off');
        quiver3(ax, origin(1),origin(2),origin(3), T(1,3)*axlen, T(2,3)*axlen, T(3,3)*axlen, 'b','LineWidth',1.8,'MaxHeadSize',0.3,'HandleVisibility','off');
    end

    function draw_fov_cone(T, fov_deg, range_m, color, two_sided)
        if nargin < 5, two_sided = true; end
        nTheta = 64; nS = 24;
        theta = linspace(0,2*pi,nTheta); s = linspace(0, range_m, nS);
        r = s * tan(deg2rad(fov_deg/2));
        [TH, SS] = meshgrid(theta, s);
        RR = repmat(r(:), 1, nTheta);
        X = RR.*cos(TH); Y = RR.*sin(TH); Z = SS;
        Sx = [0 0 1 0; 0 1 0 0; 1 0 0 0; 0 0 0 1];
        PW = T * (Sx * [X(:) Y(:) Z(:) ones(numel(X),1)]');
        Xw = reshape(PW(1,:),size(X)); Yw = reshape(PW(2,:),size(Y)); Zw = reshape(PW(3,:),size(Z));
        surf(ax, Xw, Yw, Zw, 'FaceColor', color, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility','off');
        if two_sided
            Zm  = -Z;
            PWm = T * (Sx * [X(:) Y(:) Zm(:) ones(numel(X),1)]');
            Xwm = reshape(PWm(1,:),size(X)); Ywm = reshape(PWm(2,:),size(Y)); Zwm = reshape(PWm(3,:),size(Zm));
            surf(ax, Xwm, Ywm, Zwm, 'FaceColor', color, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility','off');
        end
    end
end

function TnoYaw = remove_spin_yaw(T)
    R = T(1:3,1:3);
    x = R(:,1); x = x ./ norm(x);
    up = [0;0;1];
    y = cross(up, x);
    if norm(y) < 1e-8, up = [0;1;0]; y = cross(up, x); end
    y = y ./ norm(y); z = cross(x, y);
    R0 = [x y z];
    TnoYaw = T; TnoYaw(1:3,1:3) = R0;
end

function val = overlap_reward(rho, rho_min, rho_max)
    if rho < rho_min, val = -1000 * (rho_min - rho);
    elseif rho > rho_max, val = -600 * (rho - rho_max);
    else, val = 300; end
end

function val = angular_reward(min_ang, ang_min)
    if isinf(min_ang), val = 100;
    elseif min_ang < ang_min, val = -2000 * (ang_min - min_ang) / ang_min;
    else, val = 100; end
end

function th = pose_to_theta(M)
    x = M(1,4); y = M(2,4); z = M(3,4);
    e = mat2eulxyz_deg(M(1:3,1:3));
    th = [x,y,z, e(2), e(1)];
end