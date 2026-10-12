function ag_inventario_s(nRuns)
% AG_INVENTARIO_S  Algoritmo genetico para la optimizacion de inventarios (5 productos).
% Cruce uniforme con reparacion, mutacion gaussiana con reintentos.
% Estudios: 1 (seleccion), 2 (iteraciones) y 3 (mutacion: Pm y sigma).
% Uso:  ag_inventario_s              -> 8 corridas por escenario
%       ag_inventario_s(2)           -> prueba rapida
%       ag_inventario_s('verificar') -> compara la forma cerrada con integral
% Individuo: 14 variables [Q1..Q4 S1..S5 K1..K5]; Q5 se despeja de la restriccion.
% Genera res_s_inventario.mat y muestra las tablas en pantalla.
if nargin >= 1 && ischar(nRuns), verificar(); return; end
if nargin < 1, nRuns = 8; end
rng(31);

cfg.n = 60;  cfg.nRuns = nRuns;  cfg.iter0 = 100;  cfg.Pm0 = 0.15;  cfg.frac0 = 0.15;
cfg.iterLista = [10 20 30 50 100 150 200];
cfg.PmLista = [0 0.02 0.1 0.25 0.5 0.9];
cfg.fracLista = [0.02 0.05 0.15 0.3 0.6];       % sigma como fraccion del rango de cada variable
metodos = {'elitismo','torneo','ruleta'};
prob = problema_inventario(cfg.frac0);

%% Estudio 1: metodo de seleccion
S1 = struct(); Rs = [];
for k = 1:3
    S1.(metodos{k}) = correr(prob, cfg.n, cfg.iter0, cfg.Pm0, metodos{k}, nRuns);
    if k == 1, Rs = S1.(metodos{k}); else, Rs(k) = S1.(metodos{k}); end
end
[~, ib] = min([Rs.media]);
mejorSel = metodos{ib};
imprimir_tabla(['Estudio 1 - inventarios (mejor seleccion: ' mejorSel ')'], metodos, Rs);

%% Estudio 2: numero de iteraciones
S2 = []; etiq = cell(1, numel(cfg.iterLista));
for k = 1:numel(cfg.iterLista)
    r = correr(prob, cfg.n, cfg.iterLista(k), cfg.Pm0, mejorSel, nRuns);
    if k == 1, S2 = r; else, S2(k) = r; end
    etiq{k} = sprintf('iter=%d', cfg.iterLista(k));
end
[mejorMedia, ib] = min([S2.media]);
k0 = find([S2.media] <= mejorMedia + S2(ib).std, 1, 'first');   % a una desviacion estandar de la mejor
iterOpt = cfg.iterLista(k0);
imprimir_tabla(['Estudio 2 - inventarios (iteraciones elegidas: ' num2str(iterOpt) ')'], etiq, S2);

%% Estudio 3: mutacion (primero Pm, luego sigma con la mejor Pm no nula)
Rpm = []; etiqPm = cell(1, numel(cfg.PmLista));
for k = 1:numel(cfg.PmLista)
    r = correr(prob, cfg.n, iterOpt, cfg.PmLista(k), mejorSel, nRuns);
    if k == 1, Rpm = r; else, Rpm(k) = r; end
    etiqPm{k} = sprintf('Pm=%.2f', cfg.PmLista(k));
end
[~, ib] = min([Rpm(2:end).media]);
pmOpt = cfg.PmLista(ib + 1);
Rsg = []; etiqSg = cell(1, numel(cfg.fracLista));
for k = 1:numel(cfg.fracLista)
    r = correr(problema_inventario(cfg.fracLista(k)), cfg.n, iterOpt, pmOpt, mejorSel, nRuns);
    if k == 1, Rsg = r; else, Rsg(k) = r; end
    etiqSg{k} = sprintf('sigma=%.2f', cfg.fracLista(k));
end
[~, ib] = min([Rsg.media]);
fracOpt = cfg.fracLista(ib);
S3 = struct(); S3.pm = Rpm; S3.sig = Rsg;
imprimir_tabla(['Estudio 3 - inventarios (Pm=' num2str(pmOpt) ', sigma=' num2str(fracOpt) ' del rango)'], [etiqPm, etiqSg], [Rpm, Rsg]);

save('res_s_inventario.mat', 'cfg', 'S1', 'S2', 'S3', 'mejorSel', 'iterOpt', 'pmOpt', 'fracOpt');
fprintf('\nResultados guardados en res_s_inventario.mat\n');
end

%% ===================== FUNCION OBJETIVO =====================
function p = params_inventario()
persistent P
if isempty(P)
    P.D  = [12000 80000 150000 50000 100000];
    P.m  = [20 220 900 120 180];
    P.mu = [300 4000 7000 2000 3500];
    P.sg = [100 1200 2500 700 1100];
    P.Smax = [200 1200 2500 800 1500];
    P.Kmax = [500 2000 5000 1500 3000];
    P.I = 8000;
    P.Qmin = 10;       % valor positivo (el intervalo es (0, inf))
    P.Qmax = P.I;      % supuesto: un pedido no supera la inversion promedio total
end
p = P;
end

function Z = fobj_inventario(Q, S, K, metodo)
% Z = sum_i D_i/(Q_i m_i) * E[(X_i - R_i)+],  R_i = K_i + S_i,  X_i ~ Normal(mu_i, sigma_i)
if nargin < 4, metodo = 'cerrada'; end
p = params_inventario();
R = K + S; Z = 0;
for i = 1:5
    if strcmp(metodo, 'integral')
        g = @(x) (x - R(i)).*exp(-(x - p.mu(i)).^2/(2*p.sg(i)^2))/(p.sg(i)*sqrt(2*pi));
        J = integral(g, R(i), Inf);
    else
        z = (R(i) - p.mu(i))/p.sg(i);
        cola = 0.5*erfc(z/sqrt(2));                 % 1 - Phi(z)
        J = (p.mu(i) - R(i))*cola + p.sg(i)*exp(-z^2/2)/sqrt(2*pi);
    end
    Z = Z + p.D(i)/(Q(i)*p.m(i))*J;
end
end

function verificar()
Q = [2000 3000 1000 1500 1000]; S = [100 600 1000 400 700]; K = [400 1500 3000 1000 2000];
fprintf('Z forma cerrada = %.6f\n', fobj_inventario(Q, S, K, 'cerrada'));
fprintf('Z con integral  = %.6f   (deben coincidir, cerca de 1473.0238)\n', fobj_inventario(Q, S, K, 'integral'));
end

function q5 = q5_inventario(x)
% Restriccion sum(Q_i/2 + S_i) = I  ->  se despeja Q5
p = params_inventario();
q5 = 2*(p.I - sum(x(5:9)) - sum(x(1:4))/2);
end

function Z = fit_inventario(x)
p = params_inventario();
q5 = q5_inventario(x);
if q5 < p.Qmin - 1e-9 || q5 > p.Qmax + 1e-9, Z = 1e6; return; end
Z = fobj_inventario([x(1:4) q5], x(5:9), x(10:14), 'cerrada');
end

function [lb, ub] = limites_inv()
p = params_inventario();
lb = [repmat(p.Qmin,1,4), zeros(1,10)];
ub = [repmat(p.Qmax,1,4), p.Smax, p.Kmax];
end

%% ===================== POBLACION, CRUCE Y MUTACION =====================
function P = pob_inventario(n)
% S y K al azar; la suma de los cinco Q queda fijada por la restriccion y se reparte al azar
p = params_inventario();
P = zeros(n, 14); k = 0;
while k < n
    S = rand(1,5).*p.Smax; K = rand(1,5).*p.Kmax;
    total = 2*(p.I - sum(S));
    w = rand(1,5); Q = total*w/sum(w);
    if all(Q >= p.Qmin) && all(Q <= p.Qmax)
        k = k + 1; P(k,:) = [Q(1:4) S K];
    end
end
end

function H = cruce_uniforme(padres, nh)
% Cruce uniforme gen a gen. Si el hijo viola la restriccion (Q5 fuera de limites),
% se repara reemplazandolo por una mezcla aritmetica de los padres (siempre factible).
p = params_inventario();
np = size(padres,1); H = zeros(nh, 14);
for k = 1:nh
    ij = randperm(np, 2); a = padres(ij(1),:); b = padres(ij(2),:);
    mascara = rand(1,14) < 0.5;
    h = a; h(mascara) = b(mascara);
    q5 = q5_inventario(h);
    if q5 < p.Qmin || q5 > p.Qmax
        lam = rand; h = lam*a + (1 - lam)*b;
    end
    H(k,:) = h;
end
end

function x = mut_inventario(x, frac)
% Mutacion gaussiana sobre una variable (sigma = frac * rango). Hasta 10 intentos para
% que Q5 siga dentro de sus limites; si ninguno sirve, el individuo no cambia.
[lb, ub] = limites_inv(); p = params_inventario();
for intento = 1:10
    j = randi(14); y = x;
    y(j) = min(max(x(j) + frac*(ub(j) - lb(j))*randn, lb(j)), ub(j));
    q5 = q5_inventario(y);
    if q5 >= p.Qmin && q5 <= p.Qmax, x = y; return; end
end
end

function prob = problema_inventario(frac)
prob.fobj = @(x) fit_inventario(x);
prob.pob = @(n) pob_inventario(n);
prob.cruce = @(P, nh) cruce_uniforme(P, nh);
prob.mut = @(x) mut_inventario(x, frac);
end

%% ===================== NUCLEO DEL ALGORITMO GENETICO =====================
function idx = seleccionar(metodo, fit, m)
fit = fit(:); n = numel(fit);
switch metodo
    case 'elitismo'
        [~, o] = sort(fit, 'ascend'); idx = o(1:m);
    case 'torneo'          % grupos de 4, gana el mejor y sale del conjunto
        libres = 1:n; idx = zeros(m,1);
        for j = 1:m
            g = libres(randperm(numel(libres), min(4, numel(libres))));
            [~, w] = min(fit(g));
            idx(j) = g(w);
            libres(libres == g(w)) = [];
        end
    case 'ruleta'          % con reposicion, probabilidad proporcional a 1/f
        c = cumsum(1./fit); c = c/c(end); idx = zeros(m,1);
        for j = 1:m, idx(j) = find(rand <= c, 1, 'first'); end
    otherwise
        error('Metodo de seleccion desconocido');
end
end

function f = evaluar(fobj, P)
f = zeros(size(P,1), 1);
for i = 1:size(P,1), f(i) = fobj(P(i,:)); end
end

function R = ag(prob, n, iter, Pm, metodoSel)
t0 = tic;
P = prob.pob(n); fit = evaluar(prob.fobj, P);
[bestFit, ib] = min(fit); best = P(ib,:);
hist = zeros(iter+1, 1); hist(1) = bestFit;
m = floor(n/2); nh = n - m;
for g = 1:iter
    idx = seleccionar(metodoSel, fit, m);
    Pp = P(idx,:); fp = fit(idx);
    if bestFit < min(fp)                 % el mejor historico siempre participa
        [~, iw] = max(fp); Pp(iw,:) = best; fp(iw) = bestFit;
    end
    H = prob.cruce(Pp, nh);
    for k = find(rand(nh,1) < Pm)'
        H(k,:) = prob.mut(H(k,:));
    end
    fh = evaluar(prob.fobj, H);
    P = [Pp; H]; fit = [fp; fh];
    [fm, ib] = min(fit);
    if fm < bestFit, bestFit = fm; best = P(ib,:); end
    hist(g+1) = bestFit;
end
R.mejor = best; R.fit = bestFit; R.hist = hist; R.t = toc(t0);
end

function res = correr(prob, n, iter, Pm, sel, nRuns)
Y = zeros(nRuns,1); T = Y; H = zeros(iter+1, nRuns); ind = cell(nRuns,1);
for r = 1:nRuns
    R = ag(prob, n, iter, Pm, sel);
    Y(r) = R.fit; T(r) = R.t; H(:,r) = R.hist; ind{r} = R.mejor;
end
[res.mejor, ib] = min(Y);
res.ind = ind{ib};
res.media = mean(Y);
res.std = std(Y);
res.t = mean(T);
res.Y = Y;
res.curva = mean(H, 2);
end

function s = fmt_ind(x)
x15 = [x(1:4) q5_inventario(x) x(5:14)];
s = ['Q=[' sprintf('%.0f ', x15(1:5)) '] S=[' sprintf('%.0f ', x15(6:10)) '] K=[' sprintf('%.0f ', x15(11:15)) ']'];
end

function imprimir_tabla(titulo, etiq, Rs)
fprintf('\n%s\n', titulo);
fprintf('| Configuracion | Tiempo (s) | Mejor Z | Media | Desv. est. | Mejor individuo |\n');
fprintf('|---|---|---|---|---|---|\n');
for k = 1:numel(Rs)
    r = Rs(k);
    fprintf('| %s | %.3f | %.4f | %.4f | %.4f | %s |\n', etiq{k}, r.t, r.mejor, r.media, r.std, fmt_ind(r.ind));
end
end