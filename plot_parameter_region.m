% plot_parameter_region.m
% Regions of the two theorems in (C2,C3)-space.  Theorem A regions are also
% intersected with the half-plane  C3 > -(fCC+frr) C2 + (fCC+fRR)(fCC+frr)(fRR+frr).
%
% Theorem A (phi-version):    phi_i    = d_i / f_ii,              S       = -sum(phi)
% % Theorem B (theta-version):  theta_ij = (f_ii+f_jj)/(d_i+d_j),   S_theta = -sum(theta)
%
% Eigenvalues Lambda = -k^2.  The unstable wavenumbers are k = kappa1..kappa2, so
% m < sigma_1 = -kappa2^2 <= ... <= sigma_n = -kappa1^2 < M, with stable neighbours
% m = -(kappa2+1)^2 and M = -(kappa1-1)^2 (M = 0 when kappa1 = 1).
% Every sub-window [a,b] with kappa1 <= a <= b <= kappa2 is plotted as well (full window first);
% Every window, including a = b, uses the multiple-unstable-mode table (sigma_1 = sigma_n when a = b).
% Rows of msM are [m sigma_1 sigma_n M], rows of kwin are [a b].
%
% Figure 1: Theorem A, one subplot per row of msM.
% % Figure 2: Theorem B, one subplot per row of msM.
% Figure 2: every Theorem A region overlaid on one plane.
%           (Theorem B is commented out; axes are fit to Theorem A only.)

%% ---- parameters (edit these) ----
% f_ii can be any nonzero real.  f_ii < 0 gives S > 0 (always the last row of the table);
% f_ii > 0 lets S < 0 and reaches the empty, triangle and first two quadrilateral rows.
fcc = -0.1;
frr = -0.1;
frhorho = -0.1;

dc = 5.0;                   % d_i > 0
dr = 1.0;
drho= 5.0;
d   = [dc dr drho];         % [dC dR drho]


phi   = [dc/fcc dr/frr drho/frhorho];     % [phiC phiR phirho]             (Theorem A)
% theta = [(fcc+frr)/(dc+dr) (frhorho+frr)/(drho+dr) (fcc+frhorho)/(dc+drho)];     % [thetaCR thetaRhoR thetaCRho]  (Theorem B)

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
f = [fcc frr frhorho];                          % f_ii in R \ {0}
if any(f == 0), error('phi_i = d_i/f_ii is undefined for f_ii = 0.'); end
[VA, namesA, labA, geomA] = compute_regions(msM, kwin, @(m,sn,s1,M) region_vertices_phi(d, phi, m, sn, s1, M));
print_cases(msM, kwin, S, 'S', geomA);
[aA, bA] = halfplane(f);
plot_subplots(VA, namesA, labA, S, 'S', aA, bA);
fig1 = gcf;

% %% ---- Theorem B (theta) ----
% S_theta = -sum(theta);
% [VB, namesB, labB, geomB] = compute_regions(msM, kwin, @(m,sn,s1,M) region_vertices_theta(d, theta, m, sn, s1, M));
% print_cases(msM, kwin, S_theta, 'S_theta', geomB);
% plot_subplots(VB, namesB, labB, S_theta, 'S_\theta', [], []);       % no half-plane

%% ---- Figure 2: all Theorem A regions on one plane ----
% Vall = vertcat(VA{:}, VB{:});
Vall = vertcat(VA{:});                          % axes fit to Theorem A only
% if isempty(Vall), error('every set of both theorems is empty.'); end
if isempty(Vall), error('every set of Theorem A is empty.'); end
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
        pairs = [m M; m sn; M 0; s1 0];     tags = {'mM', 'm\sigma_n', 'M0', '\sigma_1 0'};             geom = 'quadrilateral (I)';
    elseif S < s1 + M
        pairs = [m M; m sn; M 0; sn 0];     tags = {'mM', 'm\sigma_n', 'M0', '\sigma_n 0'};             geom = 'quadrilateral (II)';
    else
        pairs = [s1 M; s1 sn; m M; m sn];   tags = {'\sigma_1 M', '\sigma_1\sigma_n', 'mM', 'm\sigma_n'}; geom = 'quadrilateral (III)';
    end
    names = cellfun(@(t) ['$' Vsym '_{' t '}$'], tags, 'UniformOutput', false);
end

function [V, names, geom] = region_vertices_phi(d, phi, m, sn, s1, M)
% Theorem A: vertices in terms of phi_i = d_i/f_ii.
    S = -sum(phi);
    [pairs, names, geom] = vertex_table(S, m, sn, s1, M, 'V');
    xi2 = phi(1)*phi(2) + phi(1)*phi(3) + phi(2)*phi(3);
    xi3 = prod(phi);
    V = zeros(size(pairs,1), 2);
    for i = 1:size(pairs,1)
        x = pairs(i,1);  y = pairs(i,2);  z = S - x - y;
        C2 = d(1)*d(3) * (xi2 - x*y + (x+y)*z);
        C3 = d(1)*d(2)*d(3) * (phi(2)*C2/(d(1)*d(3)) - x*y*z - xi3);
        V(i,:) = [C2 C3];
    end
end

% function [V, names, geom] = region_vertices_theta(d, th, m, sn, s1, M)
% % Theorem B: vertices in terms of theta_ij = (f_ii+f_jj)/(d_i+d_j),
% % th = [thetaCR thetaRhoR thetaCRho].
%     S = -sum(th);
%     [pairs, names, geom] = vertex_table(S, m, sn, s1, M, '\overline{V}');
%     thCR = th(1);  thRr = th(2);  thCr = th(3);
%     z2 = thCR*thCr + thCR*thRr + thCr*thRr;        % zeta_2
%     z3 = thCR*thCr*thRr;                           % zeta_3
%     pCR = d(1)+d(2);  pRr = d(2)+d(3);  pCr = d(1)+d(3);   % (dC+dR), (dR+drho), (dC+drho)
%     V = zeros(size(pairs,1), 2);
%     for i = 1:size(pairs,1)
%         x = pairs(i,1);  y = pairs(i,2);  z = S - x - y;
%         C2 = pCR*pRr * (x*y + (x+y)*z - z2);
%         C3 = -pCR*pCr*pRr * (x*y*z - thCr*pCR*pRr*C2 - z3);      % as written in the theorem
%         % C3 = -pCR*pCr*pRr * (x*y*z - thCr*C2/(pCR*pRr) - z3);  % if the theta*C2 term should be divided, as in Theorem A
%         V(i,:) = [C2 C3];
%     end
% end
