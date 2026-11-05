% plot_parameter_region.m
% Regions of the two theorems in (C2,C3)-space.  Theorem A regions are also
% intersected with the half-plane  C3 > -(fCC+frr) C2 + (fCC+fRR)(fCC+frr)(fRR+frr).
%
% Theorem A (phi-version):    phi_i    = f_ii / d_i,              S       = -sum(phi)
% Theorem B (theta-version):  theta_ij = (f_ii+f_jj)/(d_i+d_j),   S_theta = -sum(theta)
%
% Eigenvalues Lambda = -k^2.  The unstable wavenumbers are k = kappa1..kappa2, so
% m < sigma_1 = -kappa2^2 <= ... <= sigma_n = -kappa1^2 < M, with stable neighbours
% m = -(kappa2+1)^2 and M = -(kappa1-1)^2 (M = 0 when kappa1 = 1).
% Every sub-window [a,b] with kappa1 <= a <= b <= kappa2 is plotted as well (full window first);
% Every window, including a = b, uses the multiple-unstable-mode table (sigma_1 = sigma_n when a = b).
% Rows of msM are [m sigma_1 sigma_n M], rows of kwin are [a b].
%
% Figure 1: Theorem A, one subplot per row of msM.
% Figure 2, left: every Theorem A (steady-state) and Theorem B (Hopf) polygon on one plane.
% Figure 2, right: exact map of the zero-mode wedge, assuming the zero mode is stable.
%           It lists every window that can occur, so it does not depend on kappa1 and kappa2.

%% ---- parameters (edit these) ----
% f_ii can be any real.  f_ii < 0 gives S > 0 (always the last row of the table);
% f_ii > 0 lets S < 0 and reaches the empty, triangle, first quadrilateral and pentagon rows.
%
% With dC = drho, fcc and frhorho only change the picture through their sum; their product
% shifts it along C2.  Uncomment one set:

% Baseline: all f_ii < 0, so the Theorem B set is empty (also used by cycle_atlas.jl).
% fcc = -5.0;  frr = -1.0;  frhorho = -5.0;

% Set 1: g(0) = 0 visibly cuts the Theorem B polygons; S_theta <= -5 leaves mode 1 no window alone.
% fcc = 1.0;  frr = -4.0;  frhorho = 1.0;

% Set 2: stronger Rac decay; most Theorem B polygons survive, every table row appears for Theorem B.
% fcc = 2.0;  frr = -8.0;  frhorho = 2.0;

% Set 3: every non-empty Theorem B polygon with kappa2 <= 4 survives in full.
% fcc = 3.0;  frr = -12.0;  frhorho = 2.0;

% Set 4: same sum as Set 1, so the same picture shifted by -9 along C2.
% fcc = 4.0;  frr = -4.0;  frhorho = -2.0;

% Set 5: Rac autocatalysis with tr L(0) > 0, which violates the assumption in Figure 2 that the zero
% mode is stable.
% fcc = -2.0;  frr = 8.0;  frhorho = -3.0;

% Set 6 (with the diffusivities below): steady-state and Hopf regions both inside the wedge, and
% mode k = 1 has both a steady-state and a Hopf edge on the stable region.
% fcc = -7.5;  frr = 3.5;  frhorho = -1.0;

% Set 7: a wide range of steady-state patterns under mutual inhibition (C2 > 0, C3 < 0): each of the
% modes k = 1..8 can be the only unstable mode, and windows reach up to k = 12.
fcc = -5.0;  frr = -16.0;  frhorho = -5.0;

% Diffusivities for the baseline and Sets 1-5 and 7:
dc = 0.05;  dr = 0.5;  drho = 0.05;     % d_i > 0
% Diffusivities for Set 6:
% dc = 5.0;  dr = 0.5;  drho = 0.2;
d   = [dc dr drho];         % [dC dR drho]


phi   = [fcc/dc frr/dr frhorho/drho];     % [phiC phiR phirho]             (Theorem A)
theta = [(fcc+frr)/(dc+dr) (frhorho+frr)/(drho+dr) (fcc+frhorho)/(dc+drho)];     % [thetaCR thetaRhoR thetaCRho]  (Theorem B)

kappa1 = 1;                 % first unstable wavenumber: sigma_n = -kappa1^2  (kappa1 >= 1)
kappa2 = 4;                 % last unstable wavenumber:  sigma_1 = -kappa2^2  (kappa2 >= kappa1)
if kappa1 < 1 || kappa2 < kappa1 || any(mod([kappa1 kappa2], 1))
    error('need integers 1 <= kappa1 <= kappa2 (the unstable modes sigma_i must be nonzero).');
end


kwin = [kappa1 kappa2];                         % full window first, then every sub-window
for a = kappa1:kappa2
    for b = a:kappa2
        if ~(a == kappa1 && b == kappa2), kwin(end+1,:) = [a b]; end %#ok<SAGROW>
    end
end
msM = [-(kwin(:,2)+1).^2, -kwin(:,2).^2, -kwin(:,1).^2, -(kwin(:,1)-1).^2] + 0;   % rows of [m sigma_1 sigma_n M]; +0 turns -0 into 0

%% ---- Theorem A (phi) ----
S = -sum(phi);
f = [fcc frr frhorho];                          % f_ii in R
[VA, namesA, labA, geomA] = compute_regions(msM, kwin, @(m,sn,s1,M) region_vertices_phi(d, phi, m, sn, s1, M));
print_cases(msM, kwin, S, 'S', geomA);
[aA, bA] = halfplane(f);
fig1 = [];
if ~isempty(vertcat(VA{:}))
    plot_subplots(VA, namesA, labA, S, 'S', aA, bA);
    fig1 = gcf;
end

%% ---- Theorem B (theta) ----
S_theta = -sum(theta);
[VB, namesB, labB, geomB] = compute_regions(msM, kwin, @(m,sn,s1,M) region_vertices_theta(d, theta, m, sn, s1, M));
print_cases(msM, kwin, S_theta, 'S_theta', geomB);

%% ---- Figure 2: steady-state and Hopf windows together, and the map of cycle space ----
plot_cycle_space(f, d, phi, theta, kwin);

%% ---- sanity check: line L through every vertex with subscript 0 (comment out when done) ----
% L:  C3 = dR*phiR*C2 - dC*dR*drho*xi3.  Every region should lie on or below L.
LA = [d(2)*phi(2), -prod(d)*prod(phi)];         % [slope intercept]
for ax = findobj(fig1, 'Type', 'axes')'
    if ~isempty(findobj(ax, 'Type', 'text', 'String', 'empty set')), continue, end
    xx = xlim(ax);  xlim(ax, xx);  ylim(ax, ylim(ax));   % keep the axes fit to the regions
    plot(ax, xx, LA(1)*xx + LA(2), '--', 'Color', [0.95 0.55 0], 'LineWidth', 2, ...
         'DisplayName', 'L: C_3 = d_R\phi_R C_2 - d_Cd_Rd_\rho\xi_3');
end

%% ---- helpers ----
function [V, names, lab, geom] = compute_regions(msM, kwin, vertex_fun)
% Vertices, vertex labels, a two-line title {window + geometry; mode values} and the set
% geometry for every row of msM = [m s1 sn M].
    n = size(msM, 1);
    V = cell(n, 1);  names = cell(n, 1);  lab = cell(n, 1);  geom = cell(n, 1);
    for k = 1:n
        m = msM(k,1);  s1 = msM(k,2);  sn = msM(k,3);  M = msM(k,4);
        [V{k}, names{k}, geom{k}] = vertex_fun(m, sn, s1, M);
        lab{k} = {sprintf('k = %d..%d:  %s', kwin(k,:), geom{k}); ...
                  sprintf('m = %g, \\sigma_1 = %g, \\sigma_n = %g, M = %g', m, s1, sn, M)};
    end
end

function print_cases(msM, kwin, S, Sname, geom)
% Table of the thresholds m+sigma_n < m+M < sigma_1+sigma_n < sigma_1+M and the geometry
% reached by S, per row.
    fprintf('%s = %g\n', Sname, S);
    fprintf('  %7s | %6s %6s %6s %6s | %8s %8s %11s %8s | %s\n', 'k', 'm', 'sig_1', 'sig_n', 'M', ...
            'm+sig_n', 'm+M', 'sig_1+sig_n', 'sig_1+M', 'geometry');
    for k = 1:size(msM, 1)
        m = msM(k,1);  s1 = msM(k,2);  sn = msM(k,3);  M = msM(k,4);
        fprintf('  %7s | %6g %6g %6g %6g | %8g %8g %11g %8g | %s\n', sprintf('%d..%d', kwin(k,:)), ...
                m, s1, sn, M, m+sn, m+M, s1+sn, s1+M, geom{k});
    end
end

function [a, b] = halfplane(f)
% Half-plane  C3 > a*C2 + b  from f = [fCC fRR frr].
    a = -(f(1) + f(3));
    b = (f(1) + f(2)) * (f(1) + f(3)) * (f(2) + f(3));
end

function [xl, yl] = axis_limits(Vall)
    xl = [min(Vall(:,1)) max(Vall(:,1))];  xl = xl + 0.3*max(diff(xl), eps)*[-1 1];
    yl = [min(Vall(:,2)) max(Vall(:,2))];  yl = yl + 0.3*max(diff(yl), eps)*[-1 1];
end

function H = halfplane_shape(a, b, xl, yl)
    top = max([yl(2), a*xl + b]) + diff(yl);    % stays above the boundary line
    H = polyshape([xl(1) xl(2) xl(2) xl(1)], [a*xl(1)+b, a*xl(2)+b, top, top]);
end

function plot_subplots(V, names, lab, S, Sname, a, b)
% One subplot per row of msM for a single theorem, each with axes fit to its own region.
% Empty a, b: no half-plane.
    if isempty(vertcat(V{:}))
        warning('%s = %g: every set is empty.', Sname, S);  return
    end
    figure('Color', 'w');
    tiledlayout('flow');
    for k = 1:numel(V)
        nexttile; hold on; grid on
        xlabel('C_2'); ylabel('C_3');
        title({lab{k}{1}, sprintf('%s   (%s = %g)', lab{k}{2}, Sname, S)}, 'FontSize', 9);
        if isempty(V{k})
            text(0.5, 0.5, 'empty set', 'Units', 'normalized', 'HorizontalAlignment', 'center');
            continue
        end
        [xl, yl] = axis_limits(V{k});
        H = [];
        if ~isempty(a), H = halfplane_shape(a, b, xl, yl); end
        draw_region(V{k}, names{k}, [0.2 0.45 0.85], [0.85 0.25 0.25], H, '-', 'o');
        if ~isempty(H), plot(xl, a*xl + b, 'k-', 'LineWidth', 1.5); end
        xlim(xl); ylim(yl);
        if isempty(H), continue, end
        if k == 1
            legend({'theorem region', 'region \cap half-plane', ...
                    'C_3 = -(f_{CC}+f_{\rho\rho})C_2 + \Pi'}, 'Location', 'best');
        end
    end
end

function plot_cycle_space(f, d, phi, theta, kwin)
% Left panel: every polygon of both theorems, for the windows in kwin and those that can reach the wedge.
% Right panel: map of the (C2,C3) plane under the assumption that the zero mode is stable.
% Inside the wedge g(0) < 0, Q(0) > 0, only windows of modes in ZA (for g) and ZB (for Q) can occur,
% so this map is complete, whatever kappa1 and kappa2 are.  It shows, at every point of the wedge,
% how many eigenvalues each mode has with positive real part, and it marks which bifurcation
% happens on each edge of the region where the homogeneous equilibrium is stable.
    trL0 = sum(f);
    if trL0 >= 0
        warning('tr L(0) = %g >= 0: the zero mode is unstable, so these f_ii violate the standing assumption; skipping the map.', trL0);
        return
    end
    st = warning('off', 'MATLAB:polyshape:repairedBySimplify');
    restore = onCleanup(@() warning(st));

    e2 = f(1)*f(2) + f(1)*f(3) + f(2)*f(3);
    [a0, b0] = g_lines(f, d, 0);                % g(0) < 0:  C3 < a0*C2 + b0
    [aQ, bQ] = q_lines(f, d, 0);                % Q(0) > 0:  C3 > aQ*C2 + bQ
    tip = [e2, a0*e2 + b0];
    [ZA, ZB] = needed_modes(f, d, phi, theta, e2);
    funA = @(m,sn,s1,M) region_vertices_phi(d, phi, m, sn, s1, M);
    funB = @(m,sn,s1,M) region_vertices_theta(d, theta, m, sn, s1, M);
    [PA, WA] = window_polygons(all_windows(ZA), funA);
    [PB, WB] = window_polygons(all_windows(ZB), funB);
    [FA, WFA] = window_polygons(unique([all_windows(ZA); kwin], 'rows'), funA);
    [FB, WFB] = window_polygons(unique([all_windows(ZB); kwin], 'rows'), funB);
    colA = [0.20 0.45 0.85];  colB = [0.90 0.40 0.25];  colM = [0.55 0.25 0.70];

    % Left-panel axes: windows that can reach the wedge; if a theorem has none, its smallest polygon.
    Vfit = tip;
    for P = [PA PB], Vfit = [Vfit; P{1}.Vertices]; end %#ok<AGROW>
    if isempty(PA), Vfit = [Vfit; smallest_polygon_vertices(FA)]; end
    if isempty(PB), Vfit = [Vfit; smallest_polygon_vertices(FB)]; end

    % Clip the polygons to the wedge, then fit the axes to the part inside it.
    Vall = tip;
    for P = [PA PB], Vall = [Vall; P{1}.Vertices]; end %#ok<AGROW>
    [xl, yl] = axis_limits(Vall);
    Hzero = wedge_shape(a0, b0, aQ, bQ, xl, yl);
    [PA, WA] = clip_polygons(PA, WA, Hzero);
    [PB, WB] = clip_polygons(PB, WB, Hzero);
    Vin = vertcat(zeros(0,2), tip);
    for P = [PA PB], Vin = [Vin; P{1}.Vertices]; end %#ok<AGROW>
    if size(Vin, 1) == 1                        % no polygon meets the wedge
        x0 = e2 - max(1, abs(e2));
        Vin = [tip; x0, a0*x0 + b0; x0, aQ*x0 + bQ];
    end
    [xl, yl] = axis_limits(Vin);
    Hzero = wedge_shape(a0, b0, aQ, bQ, xl, yl);
    Hbox = polyshape([xl(1) xl(2) xl(2) xl(1)], [yl(1) yl(1) yl(2) yl(2)]);

    UA = polyshape();  for i = 1:numel(PA), UA = union(UA, PA{i}); end
    UB = polyshape();  for j = 1:numel(PB), UB = union(UB, PB{j}); end
    Hstable = subtract(subtract(Hzero, UA), UB);

    figure('Color', 'w', 'Position', [80 80 1350 600]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    nexttile; hold on; box on
    plot_all_windows(FA, WFA, FB, WFB, a0, b0, aQ, bQ, tip, colA, colB, Vfit);

    nexttile; hold on; box on
    h = gobjects(0, 1);  labs = {};
    h(end+1,1) = plot(subtract(Hbox, Hzero), 'FaceColor', [0.82 0.82 0.82], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    labs{end+1,1} = 'zero mode unstable (outside the assumption)';
    if area(Hstable) > 0
        h(end+1,1) = plot(Hstable, 'FaceColor', [0.55 0.85 0.55], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
        labs{end+1,1} = 'stable';
    end
    for i = 1:numel(PA)
        hi = plot(PA{i}, 'FaceColor', colA, 'FaceAlpha', 0.35, 'EdgeColor', colA*0.7, 'LineWidth', 0.5);
        if i == 1, h(end+1,1) = hi; labs{end+1,1} = 'one eigenvalue with Re > 0 at each mode of the window (Theorem A)'; end %#ok<AGROW>
    end
    for j = 1:numel(PB)
        Pj = subtract(PB{j}, UA);
        if area(Pj) == 0, continue, end
        hj = plot(Pj, 'FaceColor', colB, 'FaceAlpha', 0.35, 'EdgeColor', colB*0.7, 'LineWidth', 0.5);
        if ~any(strcmp(labs, 'two eigenvalues with Re > 0 at each mode of the window (Theorem B)'))
            h(end+1,1) = hj; labs{end+1,1} = 'two eigenvalues with Re > 0 at each mode of the window (Theorem B)'; %#ok<AGROW>
        end
    end
    % A point in a Theorem A and a Theorem B polygon is mixed unless the B window lies in the A window.
    for i = 1:numel(PA)
        for j = 1:numel(PB)
            if WB(j,1) >= WA(i,1) && WB(j,2) <= WA(i,2), continue, end
            Pij = intersect(PA{i}, PB{j});
            if area(Pij) == 0, continue, end
            hm = plot(Pij, 'FaceColor', colM, 'FaceAlpha', 0.5, 'EdgeColor', colM*0.7, 'LineWidth', 0.5);
            if ~any(strcmp(labs, 'mixed: one eigenvalue at some modes, two at others'))
                h(end+1,1) = hm; labs{end+1,1} = 'mixed: one eigenvalue at some modes, two at others'; %#ok<AGROW>
            end
        end
    end
    minArea = 0.004 * diff(xl) * diff(yl);      % windows too small to label are printed instead
    label_windows(PA, WA, colA*0.6, minArea);
    label_windows(PB, WB, colB*0.6, minArea);

    % Wedge boundaries and tip (double zero eigenvalue at the zero mode).
    plot(xl, a0*xl + b0, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.75, 'HandleVisibility', 'off');
    plot(xl, aQ*xl + bQ, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.75, 'HandleVisibility', 'off');
    h(end+1,1) = plot(tip(1), tip(2), 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
    labs{end+1,1} = 'wedge tip: C_2 = e_2, double zero eigenvalue';
    if xl(1) < 0 && xl(2) > 0
        plot([0 0], yl, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
        text(0, yl(2), ' C_2 = 0', 'VerticalAlignment', 'top', 'FontSize', 9);
    end

    % Edges of the stable region: it is the convex set max Q-lines < C3 < min g-lines over {0, Z}.
    kU = [0 ZA];  kL = [0 ZB];
    [aU, bU] = g_lines(f, d, kU);
    [aL, bL] = q_lines(f, d, kL);
    aAll = [aU aL];  bAll = [bU bL];
    [I, J] = find(triu(true(numel(aAll)), 1));
    xs = (bAll(J) - bAll(I)) ./ (aAll(I) - aAll(J));
    xs = unique([xl(1), xs(isfinite(xs) & xs > xl(1) & xs < xl(2)), xl(2)]);
    rngU = nan(numel(kU), 2);  rngL = nan(numel(kL), 2);
    for s = 1:numel(xs) - 1
        xm = (xs(s) + xs(s+1)) / 2;
        [U, iU] = min(aU*xm + bU);
        [Lo, iL] = max(aL*xm + bL);
        if U <= Lo, continue, end
        rngU(iU,:) = [min(rngU(iU,1), xs(s)), max(rngU(iU,2), xs(s+1))];
        rngL(iL,:) = [min(rngL(iL,1), xs(s)), max(rngL(iL,2), xs(s+1))];
    end
    edgeNames = {'zero-mode eigenvalue crosses 0: g(0) = 0', 'mode k steady (Turing) onset: g(-k^2) = 0', ...
                 'uniform Hopf: Q(0) = 0', 'mode k Hopf: Q(-k^2) = 0'};
    edgeStyle = {{'k', '-'}, {colA*0.7, '-'}, {'k', '--'}, {colB*0.8, '--'}};
    shown = false(1, 4);
    for u = find(~isnan(rngU(:,1)))'
        t = 1 + (kU(u) > 0);
        hE = plot(rngU(u,:), aU(u)*rngU(u,:) + bU(u), 'Color', edgeStyle{t}{1}, 'LineStyle', edgeStyle{t}{2}, 'LineWidth', 2.5);
        if ~shown(t), h(end+1,1) = hE; labs{end+1,1} = edgeNames{t}; shown(t) = true; end %#ok<AGROW>
        if kU(u) > 0, label_edge(rngU(u,:), aU(u), bU(u), kU(u), xl, yl); end
    end
    for u = find(~isnan(rngL(:,1)))'
        t = 3 + (kL(u) > 0);
        hE = plot(rngL(u,:), aL(u)*rngL(u,:) + bL(u), 'Color', edgeStyle{t}{1}, 'LineStyle', edgeStyle{t}{2}, 'LineWidth', 2.5);
        if ~shown(t), h(end+1,1) = hE; labs{end+1,1} = edgeNames{t}; shown(t) = true; end %#ok<AGROW>
        if kL(u) > 0, label_edge(rngL(u,:), aL(u), bL(u), kL(u), xl, yl); end
    end

    xlim(xl); ylim(yl);
    xlabel('C_2'); ylabel('C_3');
    if e2 >= 0
        note = 'e_2 \geq 0: no mode k \geq 1 can have Q < 0 or a Hopf bifurcation (Theorem E)';
    else
        note = 'e_2 < 0: the wedge lies in C_2 < e_2 < 0, so C_2 > 0 is excluded';
    end
    title({'Map of the wedge (stable zero mode)', ...
           sprintf('f = (%g, %g, %g),  tr L(0) = %g,  e_2 = %g', f, trL0, e2), note});
    legend(h, labs, 'Location', 'bestoutside');

    fprintf('Cycle-space map: e2 = %g, wedge tip (C2,C3) = (%g, %g)\n', e2, tip);
    modes = @(Z) strtrim([num2str(Z) repmat('none', 1, isempty(Z))]);
    fprintf('  modes with R < 0 (g lines needed besides k = 0):  %s\n', modes(ZA));
    fprintf('  modes with Rbar < 0 (Q lines needed besides k = 0): %s\n', modes(ZB));
    print_windows('Theorem A', PA, WA);
    print_windows('Theorem B', PB, WB);
end

function plot_all_windows(FA, WFA, FB, WFB, a0, b0, aQ, bQ, tip, colA, colB, V)
% Every Theorem A (steady-state) and Theorem B (Hopf) polygon on the current axes.  Darker fill marks
% the part inside the zero-mode wedge.  A centre marker keeps very small polygons visible.
% The axes fit the points V, so larger polygons may run off them; zoom out to see those.
    if size(V, 1) == 1
        x0 = tip(1) - max(1, abs(tip(1)));
        V = [tip; x0, a0*x0 + b0; x0, aQ*x0 + bQ];
    end
    [xl, yl] = axis_limits(V);
    Hzero = wedge_shape(a0, b0, aQ, bQ, xl, yl);
    h = gobjects(0, 1);  labs = {};
    h(end+1,1) = plot(Hzero, 'FaceColor', [0.55 0.85 0.55], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    labs{end+1,1} = 'zero-mode wedge';
    h = [h; draw_windows(FA, WFA, Hzero, colA, '-')];
    if ~isempty(FA), labs{end+1,1} = 'steady-state windows (Theorem A)'; end
    h = [h; draw_windows(FB, WFB, Hzero, colB, '--')];
    if ~isempty(FB), labs{end+1,1} = 'Hopf windows (Theorem B)'; end
    plot(xl, a0*xl + b0, 'k-',  'LineWidth', 1, 'HandleVisibility', 'off');
    plot(xl, aQ*xl + bQ, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(tip(1), tip(2), 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 9, 'HandleVisibility', 'off');
    xlim(xl); ylim(yl);
    xlabel('C_2'); ylabel('C_3');
    title({'Steady-state and Hopf windows', 'solid line: g(0) = 0,  dashed line: Q(0) = 0'});
    legend(h, labs, 'Location', 'best');
end

function V = smallest_polygon_vertices(F)
% Vertices of the smallest polygon in F (empty if F is empty).
    V = zeros(0, 2);
    if isempty(F), return, end
    [~, i] = min(cellfun(@area, F));
    V = F{i}.Vertices;
end

function hFirst = draw_windows(P, W, Hzero, col, ls)
% Fill the polygons, darken their part in the wedge, label each by its window.
    hFirst = gobjects(0, 1);
    for i = 1:numel(P)
        hi = plot(P{i}, 'FaceColor', col, 'FaceAlpha', 0.15, 'EdgeColor', col*0.7, 'LineStyle', ls);
        if i == 1, hFirst = hi; end
        Pin = intersect(P{i}, Hzero);
        if area(Pin) > 0
            plot(Pin, 'FaceColor', col, 'FaceAlpha', 0.45, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        [cx, cy] = centroid(P{i});
        plot(cx, cy, '.', 'Color', col*0.7, 'MarkerSize', 8, 'HandleVisibility', 'off');
        text(cx, cy, sprintf(' %d..%d', W(i,:)), 'Color', col*0.6, 'FontSize', 7, 'Clipping', 'on');
    end
end

function H = wedge_shape(a0, b0, aQ, bQ, xl, yl)
% The rectangle xl-by-yl intersected with g(0) < 0 and Q(0) > 0.
    H = intersect(g_halfplanes_shape(a0, b0, xl, yl), halfplane_shape(aQ, bQ, xl, yl));
end

function [ZA, ZB] = needed_modes(f, d, phi, theta, e2)
% Modes k >= 1 with R(k^2) < 0 and with Rbar(k^2) < 0.  At every other mode, the g and Q conditions
% hold throughout the wedge (Theorem D).
    delta = prod(d);
    xi2 = phi(1)*phi(2) + phi(1)*phi(3) + phi(2)*phi(3);
    pCR = d(1)+d(2);  pRr = d(2)+d(3);  pCr = d(1)+d(3);
    P = pCR*pCr*pRr;
    cbar = pCR*(e2 + f(3)^2) + pCr*f(2)^2 + pRr*(e2 + f(1)^2);
    ZA = modes_below_zero(delta, -delta*sum(phi), delta*xi2 - d(2)*e2);
    ZB = modes_below_zero(P, -P*sum(theta), cbar);
end

function Z = modes_below_zero(A, B, C)
% Integers k >= 1 with A*k^4 + B*k^2 + C < 0, for A > 0.
    disc = B^2 - 4*A*C;
    Z = zeros(1, 0);
    if disc <= 0, return, end
    qLo = (-B - sqrt(disc)) / (2*A);
    qHi = (-B + sqrt(disc)) / (2*A);
    k = 1:floor(sqrt(max(qHi, 0)));
    Z = k(k.^2 > qLo & k.^2 < qHi);
end

function W = all_windows(Z)
% Every window a..b with a, b in Z, one row [a b] per window.
    W = zeros(0, 2);
    for a = Z
        for b = Z(Z >= a)
            W(end+1,:) = [a b]; %#ok<AGROW>
        end
    end
end

function [P, W] = window_polygons(W0, vertex_fun)
% The non-empty polygon of each window in the rows of W0, and those windows.
    P = cell(1, 0);  W = zeros(0, 2);
    for i = 1:size(W0, 1)
        a = W0(i,1);  b = W0(i,2);
        V = vertex_fun(-(b+1)^2, -a^2, -b^2, -(a-1)^2 + 0);
        if isempty(V), continue, end
        [~, order] = sort(atan2(V(:,2) - mean(V(:,2)), V(:,1) - mean(V(:,1))));
        Pab = polyshape(V(order,1), V(order,2));
        if area(Pab) > 0
            P{end+1} = Pab;  W(end+1,:) = [a b]; %#ok<AGROW>
        end
    end
end

function [Pout, Wout] = clip_polygons(P, W, H)
% Intersect each polygon with H and drop those left with no area.
    Pout = cell(1, 0);  Wout = zeros(0, 2);
    for i = 1:numel(P)
        Pi = intersect(P{i}, H);
        if area(Pi) > 0
            Pout{end+1} = Pi;  Wout(end+1,:) = W(i,:); %#ok<AGROW>
        end
    end
end

function label_windows(P, W, col, minArea)
    for i = 1:numel(P)
        if area(P{i}) < minArea, continue, end
        [cx, cy] = centroid(P{i});
        text(cx, cy, sprintf('%d..%d', W(i,:)), 'Color', col, 'FontSize', 8, 'HorizontalAlignment', 'center');
    end
end

function label_edge(rng, a, b, k, xl, yl)
    xm = mean(rng);  ym = a*xm + b;
    if xm > xl(1) && xm < xl(2) && ym > yl(1) && ym < yl(2)
        text(xm, ym, sprintf(' k=%d', k), 'FontSize', 8, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
    end
end

function print_windows(name, P, W)
    if isempty(P)
        fprintf('  %s: no window meets the wedge\n', name);
        return
    end
    fprintf('  %s windows in the wedge (area):', name);
    for i = 1:numel(P), fprintf('  %d..%d (%.3g)', W(i,:), area(P{i})); end
    fprintf('\n');
end

function [a, b] = q_lines(f, d, k)
% Q(-k^2) > 0  <=>  C3 > a*C2 + b, at each wavenumber k.
    g = f(:) - d(:) * (k(:)'.^2);               % 3-by-numel(k): g_i = f_ii - d_i k^2
    a = -(g(1,:) + g(3,:));
    b = (g(1,:) + g(2,:)) .* (g(1,:) + g(3,:)) .* (g(2,:) + g(3,:));
end

function [a, b] = g_lines(f, d, k)
% g(-k^2) < 0  <=>  C3 < a*C2 + b, with a = g_R and b = -g_C*g_R*g_rho at each wavenumber k.
    g = f(:) - d(:) * (k(:)'.^2);               % 3-by-numel(k): g_i = f_ii - d_i k^2
    a = g(2,:);
    b = -prod(g, 1);
end

function H = g_halfplanes_shape(a, b, xl, yl)
% The rectangle xl-by-yl intersected with every half-plane C3 < a(j)*C2 + b(j).
    H = polyshape([xl(1) xl(2) xl(2) xl(1)], [yl(1) yl(1) yl(2) yl(2)]);
    for j = 1:numel(a)
        ends = a(j)*xl + b(j);
        low = min([yl(1), ends]) - diff(yl);    % stays below the boundary line
        H = intersect(H, polyshape([xl(1) xl(2) xl(2) xl(1)], [ends(1) ends(2) low low]));
    end
end

function hP = draw_region(V, names, colRegion, colClip, H, ls, mk)
% Fill the polygon with vertices V (light), its intersection with H (dark; skipped if H is empty),
% mark and label the vertices.  ls / mk: edge line style and vertex marker.
% Returns the handle of the region patch.
    [~, order] = sort(atan2(V(:,2) - mean(V(:,2)), V(:,1) - mean(V(:,1))));  % ccw
    P = polyshape(V(order,1), V(order,2));
    if P.NumRegions == 0                       % vertices coincide (e.g. S = m+sigma+M)
        hP = plot(V(1,1), V(1,2), mk, 'Color', colRegion, 'MarkerFaceColor', colRegion);
    else
        hP = plot(P, 'FaceColor', colRegion, 'FaceAlpha', 0.25, 'EdgeColor', colRegion, 'LineStyle', ls);
        if ~isempty(H)
            Q = intersect(P, H);
            plot(Q, 'FaceColor', colClip, 'FaceAlpha', 0.45, 'EdgeColor', colClip, 'LineStyle', ls);
        end
    end
    plot(V(:,1), V(:,2), mk, 'Color', colRegion, 'MarkerFaceColor', colRegion);
    text(V(:,1), V(:,2), names, 'Color', colRegion*0.7, 'Interpreter', 'latex', ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

function [pairs, names, geom] = vertex_table(S, m, sn, s1, M, Vsym)
% The (x,y) pairs and labels of the vertices, and the set geometry, from the
% multiple-unstable-mode table (m < sigma_1 <= sigma_n < M, thresholds
% m+sigma_n < m+M < sigma_1+sigma_n < sigma_1+M).
% Vsym is the LaTeX symbol of a vertex ('V' or '\overline{V}').
    if S <= m + sn
        pairs = [];                         tags = {};                                                  geom = 'empty';
    elseif S < m + M
        pairs = [m sn; m 0; sn 0];          tags = {'m\sigma_n', 'm0', '\sigma_n 0'};                   geom = 'triangle';
    elseif S < s1 + sn
        pairs = [m M; m sn; M 0; sn 0];     tags = {'mM', 'm\sigma_n', 'M0', '\sigma_n 0'};             geom = 'quadrilateral (I)';
    elseif S < s1 + M
        % pentagon when sigma_1 < sigma_n; V_{sigma_1 sigma_n} lies on an edge when sigma_1 = sigma_n
        pairs = [m M; m sn; s1 sn; s1 0; M 0];  tags = {'mM', 'm\sigma_n', '\sigma_1\sigma_n', '\sigma_1 0', 'M0'};  geom = 'pentagon';
    else
        pairs = [s1 M; s1 sn; m M; m sn];   tags = {'\sigma_1 M', '\sigma_1\sigma_n', 'mM', 'm\sigma_n'}; geom = 'quadrilateral (II)';
    end
    names = cellfun(@(t) ['$' Vsym '_{' t '}$'], tags, 'UniformOutput', false);
end

function [V, names, geom] = region_vertices_phi(d, phi, m, sn, s1, M)
% Theorem A: vertices in terms of phi_i = f_ii/d_i.
    S = -sum(phi);
    [pairs, names, geom] = vertex_table(S, m, sn, s1, M, 'V');
    xi2 = phi(1)*phi(2) + phi(1)*phi(3) + phi(2)*phi(3);
    xi3 = prod(phi);
    V = zeros(size(pairs,1), 2);
    for i = 1:size(pairs,1)
        x = pairs(i,1);  y = pairs(i,2);  z = S - x - y;
        C2 = d(1)*d(3) * (xi2 - x*y - (x+y)*z);
        C3 = d(1)*d(2)*d(3) * (phi(2)*C2/(d(1)*d(3)) - x*y*z - xi3);
        V(i,:) = [C2 C3];
    end
end

function [V, names, geom] = region_vertices_theta(d, th, m, sn, s1, M)
% Theorem B: vertices of the set defined by the Q inequalities, in terms of
% theta_ij = (f_ii+f_jj)/(d_i+d_j), th = [thetaCR thetaRhoR thetaCRho].
    S = -sum(th);
    [pairs, names, geom] = vertex_table(S, m, sn, s1, M, '\overline{V}');
    thCR = th(1);  thRr = th(2);  thCr = th(3);
    z2 = thCR*thCr + thCR*thRr + thCr*thRr;        % zeta_2
    z3 = thCR*thCr*thRr;                           % zeta_3
    pCR = d(1)+d(2);  pRr = d(2)+d(3);  pCr = d(1)+d(3);   % (dC+dR), (dR+drho), (dC+drho)
    V = zeros(size(pairs,1), 2);
    for i = 1:size(pairs,1)
        x = pairs(i,1);  y = pairs(i,2);  z = S - x - y;
        C2 = pCR*pRr * (z2 - x*y - (x+y)*z);
        C3 = pCR*pCr*pRr * (x*y*z + z3) - C2*pCr*thCr;
        V(i,:) = [C2 C3];
    end
end
