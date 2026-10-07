% plot_parameter_region.m
% Region of the theorem in (C2,C3)-space, and its intersection with the
% half-plane  C3 > -(fCC+frr) C2 + (fCC+fRR)(fCC+frr)(fRR+frr).
% Convention: phi_i = d_i / f_ii.
% Figure 1: one subplot per row of msM.   Figure 2: all regions overlaid.

%% ---- parameters (edit these) ----
d   = [1 1 1];          % [dC dR drho]
phi = -[2.0 1.5 2.0];   % [phiC phiR phirho]
msM = [-3^2 -2^2 -1^2;           % rows of [m sigma M]
       -4^2 -3^2 -2^2;
       -5^2 -4^2 -3^2;
       -6^2 -5^2 -4^2];

%% ---- half-plane  C3 > a*C2 + b ----
S = -sum(phi);
f = d ./ phi;                                   % f_ii = d_i / phi_i
a = -(f(1) + f(3));
b = (f(1) + f(2)) * (f(1) + f(3)) * (f(2) + f(3));

%% ---- vertices for every row ----
n = size(msM, 1);
V = cell(n, 1);  names = cell(n, 1);  lab = cell(n, 1);
for k = 1:n
    [V{k}, names{k}] = region_vertices(d, phi, msM(k,1), msM(k,2), msM(k,3));
    lab{k} = sprintf('m = %g, \\sigma = %g, M = %g', msM(k,:));
end
Vall = vertcat(V{:});
if isempty(Vall), error('S = %g: every set is empty.', S); end
xl = [min(Vall(:,1)) max(Vall(:,1))];  xl = xl + 0.3*max(diff(xl), eps)*[-1 1];
yl = [min(Vall(:,2)) max(Vall(:,2))];  yl = yl + 0.3*max(diff(yl), eps)*[-1 1];
H  = polyshape([xl(1) xl(2) xl(2) xl(1)], [a*xl(1)+b, a*xl(2)+b, yl(2)+1e3, yl(2)+1e3]);

%% ---- figure 1: individual subplots ----
figure('Color', 'w');
tiledlayout('flow');
for k = 1:n
    nexttile; hold on; grid on
    xlabel('C_2'); ylabel('C_3');
    title(sprintf('%s   (S = %g)', lab{k}, S));
    if isempty(V{k})
        text(0.5, 0.5, 'empty set', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        continue
    end
    draw_region(V{k}, names{k}, [0.2 0.45 0.85], [0.85 0.25 0.25], H);
    plot(xl, a*xl + b, 'k-', 'LineWidth', 1.5);
    if k == 1
        legend({'theorem region', 'region \cap half-plane', ...
                'C_3 = -(f_{CC}+f_{\rho\rho})C_2 + \Pi'}, 'Location', 'best');
    end
end

%% ---- figure 2: all regions on one plane ----
figure('Color', 'w'); hold on; grid on
cols = lines(n);
h = gobjects(n, 1);
for k = 1:n
    if isempty(V{k}), continue; end
    h(k) = draw_region(V{k}, names{k}, cols(k,:), cols(k,:), H);
end
hline = plot(xl, a*xl + b, 'k-', 'LineWidth', 1.5);
xlim(xl); ylim(yl);
xlabel('C_2'); ylabel('C_3');
title(sprintf('All regions,  S = %g   (dark fill = region \\cap half-plane)', S));
keep = isgraphics(h);
legend([h(keep); hline], [lab(keep); {'C_3 = -(f_{CC}+f_{\rho\rho})C_2 + \Pi'}], 'Location', 'best');

%% ---- helpers ----
function hP = draw_region(V, names, colRegion, colClip, H)
% Fill the polygon with vertices V (light), its intersection with H (dark),
% mark and label the vertices.  Returns the handle of the region patch.
    [~, order] = sort(atan2(V(:,2) - mean(V(:,2)), V(:,1) - mean(V(:,1))));  % ccw
    P = polyshape(V(order,1), V(order,2));
    if P.NumRegions == 0                       % vertices coincide (e.g. S = m+sigma+M)
        hP = plot(V(1,1), V(1,2), 'o', 'Color', colRegion, 'MarkerFaceColor', colRegion);
    else
        hP = plot(P, 'FaceColor', colRegion, 'FaceAlpha', 0.25, 'EdgeColor', colRegion);
        Q  = intersect(P, H);
        plot(Q, 'FaceColor', colClip, 'FaceAlpha', 0.45, 'EdgeColor', colClip);
    end
    plot(V(:,1), V(:,2), 'o', 'Color', colRegion, 'MarkerFaceColor', colRegion);
    text(V(:,1), V(:,2), names, 'Color', colRegion*0.7, ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

function [V, names] = region_vertices(d, phi, m, sg, M)
% Vertices from the theorem's table.
    S = -sum(phi);
    if S <= m + sg
        pairs = [];                       names = {};
    elseif S < m + M
        pairs = [m sg; m 0; sg 0];        names = {'V_{m\sigma}','V_{m0}','V_{\sigma0}'};
    elseif S < sg + M
        pairs = [m M; m sg; M 0; sg 0];   names = {'V_{mM}','V_{m\sigma}','V_{M0}','V_{\sigma0}'};
    else
        pairs = [sg M; m M; m sg];        names = {'V_{\sigmaM}','V_{mM}','V_{m\sigma}'};
    end
    xi2 = phi(1)*phi(2) + phi(1)*phi(3) + phi(2)*phi(3);
    xi3 = prod(phi);
    V = zeros(size(pairs,1), 2);
    for i = 1:size(pairs,1)
        x = pairs(i,1);  y = pairs(i,2);  z = S - x - y;
        C2 = d(1)*d(3) * (x*y + (x+y)*z - xi2);
        C3 = d(1)*d(2)*d(3) * (phi(2)*C2/(d(1)*d(3)) - x*y*z - xi3);
        V(i,:) = [C2 C3];
    end
end
