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
% Figure 2: Theorem A and B polygons with the half-planes g(-k^2) < 0 of Theorem B, k in kHP.
%           Also shows the zero-mode wedge g(0) < 0, Q(0) > 0, the standing assumption of theorems.tex.
% Figure 3: every Theorem A region overlaid on one plane (axes fit to Theorem A only).

%% ---- parameters (edit these) ----
% f_ii can be any real.  f_ii < 0 gives S > 0 (always the last row of the table);
% f_ii > 0 lets S < 0 and reaches the empty, triangle, first quadrilateral and pentagon rows.
%
% With dC = drho, fcc and frhorho only change the picture through their sum; their product
% shifts it along C2.  Uncomment one set:

% Baseline: all f_ii < 0, so the Theorem B set is empty (also used by cycle_atlas.jl).
% fcc = -5.0;  frr = -1.0;  frhorho = -5.0;

% Set 1: Theorem A and B polygons on one scale; g(0) = 0 visibly cuts the Theorem B polygons.
fcc = 1.0;  frr = -4.0;  frhorho = 1.0;

% Set 2: stronger Rac decay; most Theorem B polygons survive, every table row appears for Theorem B.
% fcc = 2.0;  frr = -8.0;  frhorho = 2.0;

% Set 3: every non-empty Theorem B polygon with kappa2 <= 4 survives in full.
% fcc = 3.0;  frr = -12.0;  frhorho = 2.0;

% Set 4: same sum as Set 1, so the same picture shifted by -9 along C2.
% fcc = 4.0;  frr = -4.0;  frhorho = -2.0;

% Set 5: Rac autocatalysis with tr L(0) > 0, which violates the standing assumption that the zero
% mode is stable (see theorems_general_trace.tex for this case).
% fcc = -2.0;  frr = 8.0;  frhorho = -3.0;

dc = 0.05;                  % d_i > 0
dr = 0.5;
drho= 0.05;
d   = [dc dr drho];         % [dC dR drho]


phi   = [fcc/dc frr/dr frhorho/drho];     % [phiC phiR phirho]             (Theorem A)
theta = [(fcc+frr)/(dc+dr) (frhorho+frr)/(drho+dr) (fcc+frhorho)/(dc+drho)];     % [thetaCR thetaRhoR thetaCRho]  (Theorem B)

kappa1 = 1;                 % first unstable wavenumber: sigma_n = -kappa1^2  (kappa1 >= 1)
kappa2 = 4;                 % last unstable wavenumber:  sigma_1 = -kappa2^2  (kappa2 >= kappa1)
if kappa1 < 1 || kappa2 < kappa1 || any(mod([kappa1 kappa2], 1))
    error('need integers 1 <= kappa1 <= kappa2 (the unstable modes sigma_i must be nonzero).');
end

kHP = 0:20;                 % wavenumbers of the half-planes g(-k^2) < 0 drawn in Figure 2

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

%% ---- Figure 2: Theorem A and B polygons with the half-planes of Theorem B ----
plot_polygons_and_halfplanes(VA, namesA, labA, VB, namesB, labB, f, d, kHP, S, S_theta);

%% ---- Figure 3: all Theorem A regions on one plane ----
% Vall = vertcat(VA{:}, VB{:});
Vall = vertcat(VA{:});                          % axes fit to Theorem A only
% if isempty(Vall), error('every set of both theorems is empty.'); end
if isempty(Vall)
    warning('every set of Theorem A is empty; skipping the Theorem A overlay.');
    return
end
[xl, yl] = axis_limits(Vall);
HA = halfplane_shape(aA, bA, xl, yl);

figure('Color', 'w'); hold on; grid on
n = size(msM, 1);
cols = lines(n);
hA = gobjects(n, 1);  % hB = gobjects(n, 1);
for k = 1:n
    if ~isempty(VA{k})
        hA(k) = draw_region(VA{k}, namesA{k}, cols(k,:), cols(k,:), HA, '-', 'o');
    end
    % if ~isempty(VB{k})
    %     hB(k) = draw_region(VB{k}, namesB{k}, cols(k,:), cols(k,:), [], '--', 's');
    % end
end
hlA = plot(xl, aA*xl + bA, 'k-',  'LineWidth', 1.5);
xlim(xl); ylim(yl);
xlabel('C_2'); ylabel('C_3');
% title({sprintf('All regions of both theorems,  S = %g,  S_\\theta = %g', S, S_theta), ...
%        'solid: \phi-theorem,  dashed: \theta-theorem,  dark fill = \phi-region \cap half-plane'});
title({sprintf('All regions of Theorem A,  S = %g', S), ...
       'dark fill = \phi-region \cap half-plane'});
kA = isgraphics(hA);  % kB = isgraphics(hB);
% legend([hA(kA); hB(kB); hlA], ...
%        [cellfun(@(t) ['\phi:  ' t{1}],   labA(kA), 'UniformOutput', false); ...
%         cellfun(@(t) ['\theta:  ' t{1}], labB(kB), 'UniformOutput', false); ...
%         {'half-plane boundary (\phi)'}], 'Location', 'best');
legend([hA(kA); hlA], ...
       [cellfun(@(t) ['\phi:  ' t{1}], labA(kA), 'UniformOutput', false); ...
        {'half-plane boundary (\phi)'}], 'Location', 'best');

%% ---- sanity check: line L through every vertex with subscript 0 (comment out when done) ----
% L:  C3 = dR*phiR*C2 - dC*dR*drho*xi3.  Every region should lie on or below L.
LA = [d(2)*phi(2), -prod(d)*prod(phi)];         % [slope intercept]
for ax = [findobj(fig1, 'Type', 'axes'); findobj(gcf, 'Type', 'axes')]'
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

function plot_polygons_and_halfplanes(VA, namesA, labA, VB, namesB, labB, f, d, kHP, S, S_theta)
% Theorem A polygons (solid), Theorem B polygons (dashed), the zero-mode boundaries g(0) = 0 and
% Q(0) = 0, and the smallest-slope boundary g(-k^2) = 0 among those in kHP that cross the axes.
% Grey = all half-planes g(-k^2) < 0, k in kHP.  Green = zero-mode wedge.
% Dark fill: Theorem A polygon in the wedge; Theorem B polygon in the grey set.  Dark fill needs
% tr L(0) < 0.  Only finitely many half-planes are used, so the dark Theorem B fill can only
% overestimate the Theorem B set; a warning names any other k that would cut the region shown.
    Vall = vertcat(VA{:}, VB{:});
    if isempty(Vall)
        warning('every set of both theorems is empty; skipping the half-plane figure.');
        return
    end
    [xl, yl] = axis_limits(Vall);
    [a, b] = g_lines(f, d, kHP);
    Hg = g_halfplanes_shape(a, b, xl, yl);
    trL0 = sum(f);
    [a0, b0] = g_lines(f, d, 0);                % g(0) < 0:  C3 < a0*C2 + b0
    [aQ, bQ] = halfplane(f);                    % Q(0) > 0:  C3 > aQ*C2 + bQ
    Hzero = intersect(g_halfplanes_shape(a0, b0, xl, yl), halfplane_shape(aQ, bQ, xl, yl));
    if trL0 < 0
        HA = Hzero;  HB = Hg;
    else
        HA = [];  HB = [];
        warning('tr L(0) = %g >= 0: the zero mode is unstable, so these f_ii violate the standing assumption.', trL0);
    end

    figure('Color', 'w'); hold on; grid on
    hG = plot(Hg, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    hZ = plot(Hzero, 'FaceColor', [0.2 0.7 0.3], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    n = numel(VA);
    cols = lines(n);
    hA = gobjects(n, 1);  hB = gobjects(n, 1);
    for k = 1:n
        if ~isempty(VA{k})
            hA(k) = draw_region(VA{k}, namesA{k}, cols(k,:), cols(k,:), HA, '-', 'o');
        end
        if ~isempty(VB{k})
            hB(k) = draw_region(VB{k}, namesB{k}, cols(k,:), cols(k,:), HB, '--', 's');
        end
    end
    % the zero-mode boundaries, always drawn
    h0 = plot(xl, a0*xl + b0, 'k-',  'LineWidth', 1.5);
    hQ = plot(xl, aQ*xl + bQ, 'k--', 'LineWidth', 1.5);
    % of the other g boundaries, draw only the largest- and smallest-slope ones crossing the axes
    % (k = 0 has the largest slope and is already drawn above)
    ends = [a*xl(1) + b; a*xl(2) + b];
    jCross = find(min(ends, [], 1) <= yl(2) & max(ends, [], 1) >= yl(1));
    [~, iMax] = max(a(jCross));
    [~, iMin] = min(a(jCross));
    jDraw = jCross(unique([iMax iMin]));
    jDraw = jDraw(kHP(jDraw) ~= 0);
    hL = gobjects(0);
    for j = jDraw
        hL = plot(xl, a(j)*xl + b(j), '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.75);
        % label each line where it leaves the axes
        xe = [xl(2), (yl(1) - b(j))/a(j), (yl(2) - b(j))/a(j)];
        ye = [a(j)*xl(2) + b(j), yl(1), yl(2)];
        i = find(xe >= xl(1) & xe <= xl(2) & ye >= yl(1) & ye <= yl(2), 1);
        if ~isempty(i)
            text(xe(i), ye(i), sprintf('k=%d', kHP(j)), 'Color', [0.35 0.35 0.35], 'FontSize', 8, ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        end
    end
    xlim(xl); ylim(yl);
    xlabel('C_2'); ylabel('C_3');
    title({sprintf('Theorems A and B with the half-planes g(-k^2) < 0,  k = %d..%d', kHP(1), kHP(end)), ...
           sprintf('S_\\phi = %g,  S_\\theta = %g,  tr L(0) = %g', S, S_theta, trL0), ...
           'solid: \phi-polygon (dark: in the zero-mode wedge),  dashed: \theta-polygon (dark: in the half-planes)'});
    kA = isgraphics(hA);  kB = isgraphics(hB);
    legend([hA(kA); hB(kB); hG; hZ; h0; hQ; hL(:)], ...
           [cellfun(@(t) ['\phi:  ' t{1}],   labA(kA), 'UniformOutput', false); ...
            cellfun(@(t) ['\theta:  ' t{1}], labB(kB), 'UniformOutput', false); ...
            {'g(-k^2) < 0 for all k in kHP'}; {'zero-mode wedge: g(0) < 0, Q(0) > 0'}; ...
            {'g(0) = 0'}; {'Q(0) = 0'}; repmat({'g(-k^2) = 0, smallest slope'}, numel(hL), 1)], ...
           'Location', 'best');

    % Does a half-plane outside kHP cut the region shown?
    kAll = 0:max(200, kHP(end));
    [aAll, bAll] = g_lines(f, d, kAll);
    xx = linspace(xl(1), xl(2), 400);
    [env, iMin] = min(aAll(:)*xx + bAll(:), [], 1);
    missing = setdiff(kAll(iMin(env <= yl(2))), kHP);
    if ~isempty(missing)
        warning('half-planes outside kHP also cut this figure, at k = %s; extend kHP to include them.', ...
                mat2str(missing));
    end
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
% Theorem B: vertices of the set defined by the Q inequalities (see theorems.tex), in terms of
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
