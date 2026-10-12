function ag_ackley_s(nRuns)
% AG_ACKLEY_S  Algoritmo genetico para la funcion de Ackley (d = 5 en los estudios 1 a 3).
% Cruce SBX, mutacion gaussiana con sigma como fraccion del rango.
% Estudios: 1 (seleccion), 2 (iteraciones), 3 (mutacion: Pm y sigma) y 4 (dimensiones).
% Uso:  ag_ackley_s              -> 8 corridas por escenario
%       ag_ackley_s(2)           -> prueba rapida
%       ag_ackley_s('verificar') -> comprueba la funcion objetivo
% Genera res_s_ackley.mat y muestra las tablas en pantalla.
if nargin >= 1 && ischar(nRuns), verificar(); return; end
if nargin < 1, nRuns = 8; end
rng(31);

cfg.n = 60;  cfg.nRuns = nRuns;  cfg.iter0 = 100;  cfg.Pm0 = 0.15;
cfg.d = 5;   cfg.frac0 = 0.08;   cfg.eta = 2;   % eta: indice de distribucion de SBX
cfg.iterLista = [10 20 30 50 100 150 200];
cfg.PmLista = [0 0.02 0.1 0.25 0.5 0.9];
cfg.fracLista = [0.01 0.03 0.08 0.15 0.3];      % sigma como fraccion del rango (20)
cfg.dimLista = [2 3 10 20 50 100];
ref = 3 - exp(1);
metodos = {'elitismo','torneo','ruleta'};
prob = problema_ackley(cfg.d, cfg.frac0, cfg.eta);

%% Estudio 1: metodo de seleccion
S1 = struct(); Rs = [];
for k = 1:3
    S1.(metodos{k}) = correr(prob, cfg.n, cfg.iter0, cfg.Pm0, metodos{k}, nRuns, ref);
    if k == 1, Rs = S1.(metodos{k}); else, Rs(k) = S1.(metodos{k}); end
end
[~, ib] = min([Rs.media]);
mejorSel = metodos{ib};
imprimir_tabla(['Estudio 1 - Ackley d=' num2str(cfg.d) ' (mejor seleccion: ' mejorSel ')'], metodos, Rs);

%% Estudio 2: numero de iteraciones
S2 = []; etiq = cell(1, numel(cfg.iterLista));
for k = 1:numel(cfg.iterLista)
    r = correr(prob, cfg.n, cfg.iterLista(k), cfg.Pm0, mejorSel, nRuns, ref);
    if k == 1, S2 = r; else, S2(k) = r; end
    etiq{k} = sprintf('iter=%d', cfg.iterLista(k));
end
[mejorMedia, ib] = min([S2.media]);
k0 = find([S2.media] <= mejorMedia + S2(ib).std, 1, 'first');   % a una desviacion estandar de la mejor
iterOpt = cfg.iterLista(k0);
imprimir_tabla(['Estudio 2 - Ackley (iteraciones elegidas: ' num2str(iterOpt) ')'], etiq, S2);

%% Estudio 3: mutacion (primero Pm, luego sigma con la mejor Pm no nula)
Rpm = []; etiqPm = cell(1, numel(cfg.PmLista));
for k = 1:numel(cfg.PmLista)
    r = correr(prob, cfg.n, iterOpt, cfg.PmLista(k), mejorSel, nRuns, ref);
    if k == 1, Rpm = r; else, Rpm(k) = r; end
    etiqPm{k} = sprintf('Pm=%.2f', cfg.PmLista(k));
end
[~, ib] = min([Rpm(2:end).media]);
pmOpt = cfg.PmLista(ib + 1);
Rsg = []; etiqSg = cell(1, numel(cfg.fracLista));
for k = 1:numel(cfg.fracLista)
    r = correr(problema_ackley(cfg.d, cfg.fracLista(k), cfg.eta), cfg.n, iterOpt, pmOpt, mejorSel, nRuns, ref);
    if k == 1, Rsg = r; else, Rsg(k) = r; end
    etiqSg{k} = sprintf('sigma=%.2f', cfg.fracLista(k));
end
[~, ib] = min([Rsg.media]);
fracOpt = cfg.fracLista(ib);
S3 = struct(); S3.pm = Rpm; S3.sig = Rsg;
imprimir_tabla(['Estudio 3 - Ackley (Pm=' num2str(pmOpt) ', sigma=' num2str(fracOpt) ' del rango)'], [etiqPm, etiqSg], [Rpm, Rsg]);

%% Estudio 4: numero de dimensiones
S4 = []; etiq4 = cell(1, numel(cfg.dimLista));
for k = 1:numel(cfg.dimLista)
    d = cfg.dimLista(k);
    r = correr(problema_ackley(d, fracOpt, cfg.eta), cfg.n, iterOpt, pmOpt, mejorSel, nRuns, ref);
    if k == 1, S4 = r; else, S4(k) = r; end
    etiq4{k} = sprintf('d=%d', d);
end
imprimir_tabla('Estudio 4 - Ackley vs dimensiones', etiq4, S4);

save('res_s_ackley.mat', 'cfg', 'S1', 'S2', 'S3', 'S4', 'mejorSel', 'iterOpt', 'pmOpt', 'fracOpt');
fprintf('\nResultados guardados en res_s_ackley.mat\n');
end

%% ===================== FUNCION OBJETIVO =====================
function f = fobj_ackley(x)
% Ackley con a=20, b=0.2, c=2*pi y constante final 3 en lugar de e.
a = 20; b = 0.2; c = 2*pi; d = numel(x);
t1 = -a*exp(-b*sqrt(sum(x.^2)/d));
t2 = -exp(sum(cos(c*x))/d);
f = t1 + t2 + a + 3;
end

function verificar()
fprintf('Ackley en el origen (d=5): %.6f (esperado 3-e = %.6f)\n', fobj_ackley(zeros(1,5)), 3 - exp(1));
end

%% ===================== POBLACION, CRUCE Y MUTACION =====================
function P = pob_ackley(n, d, lim)
P = lim(1) + (lim(2) - lim(1))*rand(n, d);
end

function H = cruce_sbx(padres, nh, lim, eta)
% Cruce SBX (simulated binary crossover) gen a gen; se conserva uno de los dos hijos
np = size(padres,1); L = size(padres,2); H = zeros(nh, L);
for k = 1:nh
    ij = randperm(np, 2); a = padres(ij(1),:); b = padres(ij(2),:);
    u = rand(1, L); beta = zeros(1, L);
    bajo = u <= 0.5;
    beta(bajo)  = (2*u(bajo)).^(1/(eta+1));
    beta(~bajo) = (1./(2*(1 - u(~bajo)))).^(1/(eta+1));
    if rand < 0.5
        h = 0.5*((1 + beta).*a + (1 - beta).*b);
    else
        h = 0.5*((1 - beta).*a + (1 + beta).*b);
    end
    H(k,:) = min(max(h, lim(1)), lim(2));
end
end

function x = mut_gauss(x, frac, lim)
% Mutacion gaussiana sobre un gen; sigma = frac * rango del dominio
j = randi(numel(x));
x(j) = min(max(x(j) + frac*(lim(2) - lim(1))*randn, lim(1)), lim(2));
end

function prob = problema_ackley(d, frac, eta)
lim = [-10 10];
prob.fobj = @(x) fobj_ackley(x);
prob.pob = @(n) pob_ackley(n, d, lim);
prob.cruce = @(P, nh) cruce_sbx(P, nh, lim, eta);
prob.mut = @(x) mut_gauss(x, frac, lim);
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

function res = correr(prob, n, iter, Pm, sel, nRuns, refVal)
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
res.rmse = sqrt(mean((Y - refVal).^2));
res.Y = Y;
res.curva = mean(H, 2);
end

function s = fmt_ind(x)
if numel(x) > 5
    s = ['[' sprintf('%.4f ', x(1:5)) '... (d=' num2str(numel(x)) ')]'];
else
    s = ['[' sprintf('%.4f ', x) ']'];
end
end

function imprimir_tabla(titulo, etiq, Rs)
fprintf('\n%s\n', titulo);
fprintf('| Configuracion | Tiempo (s) | Mejor f(X) | Media | Desv. est. | RMSE | Mejor individuo |\n');
fprintf('|---|---|---|---|---|---|---|\n');
for k = 1:numel(Rs)
    r = Rs(k);
    fprintf('| %s | %.3f | %.4f | %.4f | %.4f | %.4f | %s |\n', etiq{k}, r.t, r.mejor, r.media, r.std, r.rmse, fmt_ind(r.ind));
end
end