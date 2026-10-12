function ag_viajero_s(nRuns)
% AG_VIAJERO_S  Algoritmo genetico para el problema del viajero (13 ubicaciones).
% Cruce PMX, mutacion por inversion, torneo de 4 con un ganador, ruleta con reposicion.
% Estudios: 1 (seleccion), 2 (iteraciones) y 3 (probabilidad de mutacion).
% Uso:  ag_viajero_s              -> 8 corridas por escenario
%       ag_viajero_s(2)           -> prueba rapida
%       ag_viajero_s('verificar') -> comprueba la funcion objetivo y la tabla
% Genera res_s_viajero.mat y muestra las tablas en pantalla.
if nargin >= 1 && ischar(nRuns), verificar(); return; end
if nargin < 1, nRuns = 8; end
rng(31);

cfg.n = 60;  cfg.nRuns = nRuns;  cfg.iter0 = 100;  cfg.Pm0 = 0.15;
cfg.iterLista = [10 20 30 50 100 150 200];
cfg.PmLista = [0 0.02 0.1 0.25 0.5 0.9];
metodos = {'elitismo','torneo','ruleta'};
prob = problema_viajero();

%% Estudio 1: metodo de seleccion
S1 = struct(); Rs = [];
for k = 1:3
    S1.(metodos{k}) = correr(prob, cfg.n, cfg.iter0, cfg.Pm0, metodos{k}, nRuns);
    if k == 1, Rs = S1.(metodos{k}); else, Rs(k) = S1.(metodos{k}); end
end
[~, ib] = min([Rs.media]);
mejorSel = metodos{ib};
imprimir_tabla(['Estudio 1 - viajero (mejor seleccion: ' mejorSel ')'], metodos, Rs);

%% Estudio 2: numero de iteraciones
S2 = []; etiq = cell(1, numel(cfg.iterLista));
for k = 1:numel(cfg.iterLista)
    r = correr(prob, cfg.n, cfg.iterLista(k), cfg.Pm0, mejorSel, nRuns);
    if k == 1, S2 = r; else, S2(k) = r; end
    etiq{k} = sprintf('iter=%d', cfg.iterLista(k));
end
% criterio: menor numero de iteraciones cuya media queda a una desviacion estandar de la mejor
[mejorMedia, ib] = min([S2.media]);
k0 = find([S2.media] <= mejorMedia + S2(ib).std, 1, 'first');
iterOpt = cfg.iterLista(k0);
imprimir_tabla(['Estudio 2 - viajero (iteraciones elegidas: ' num2str(iterOpt) ')'], etiq, S2);

%% Estudio 3: probabilidad de mutacion
Rpm = []; etiq = cell(1, numel(cfg.PmLista));
for k = 1:numel(cfg.PmLista)
    r = correr(prob, cfg.n, iterOpt, cfg.PmLista(k), mejorSel, nRuns);
    if k == 1, Rpm = r; else, Rpm(k) = r; end
    etiq{k} = sprintf('Pm=%.2f', cfg.PmLista(k));
end
[~, ib] = min([Rpm(2:end).media]);          % se elige entre las Pm distintas de cero
pmOpt = cfg.PmLista(ib + 1);
S3 = struct(); S3.pm = Rpm;
imprimir_tabla(['Estudio 3 - viajero (Pm elegida: ' num2str(pmOpt) ')'], etiq, Rpm);

save('res_s_viajero.mat', 'cfg', 'S1', 'S2', 'S3', 'mejorSel', 'iterOpt', 'pmOpt');
fprintf('\nResultados guardados en res_s_viajero.mat\n');
end

%% ===================== FUNCION OBJETIVO =====================
function D = matriz_tiempos()
% Tabla I de la guia (minutos). 1 UCO, 2 Rionegro, 3 La Ceja, 4 El Carmen, 5 Marinilla,
% 6 Guarne, 7 El Retiro, 8 Santuario, 9 La Union, 10 San Vicente, 11 Cocorna,
% 12 El Penol, 13 Guatape. D(i,j) = tiempo de i hacia j.
D = [  0  10  35  29  16  30  40  27  56  44  67  54  81
       6   0  32  31  16  32  34  28  52  44  63  56  83
      34  32   0  29  45  58  35  57  26  73  91  83 110
      29  29  28   0  32  47  57  31  49  61  67  70  99
      18  18  46  33   0  23  52  19  68  39  54  41  69
      30  30  58  46  25   0  46  36  79  39  71  62  90
      40  37  36  59  52  54   0  64  61  80 104  89 117
      32  32  57  32  18  36  65   0  82  50  41  51  82
      55  52  23  48  67  80  56  73   0  95 116 101 130
      47  47  76  62  42  37  73  53  97   0 102  81 111
      79  79 106  80  66  87 110  59 148 101   0  97 100
      59  59  86  74  44  65  90  55 107  80  93   0  29
      89  89 116 104  74  96 119  86 137 110 103  33   0 ];
end

function T = fobj_viajero(ruta, regreso)
% Tiempo total (min) de un recorrido. Si no es una permutacion de 1..13 devuelve 10000.
if nargin < 2, regreso = false; end
persistent D
if isempty(D), D = matriz_tiempos(); end
if numel(ruta) ~= 13 || ~isequal(sort(ruta(:)'), 1:13)
    T = 10000; return
end
T = 0;
for i = 1:12, T = T + D(ruta(i), ruta(i+1)); end
if regreso, T = T + D(ruta(13), ruta(1)); end
end

function verificar()
D = matriz_tiempos();
fprintf('Ruta 1:13 sin regreso: %d (esperado 641)\n', fobj_viajero(1:13, false));
fprintf('Ruta 1:13 con regreso: %d (esperado 730)\n', fobj_viajero(1:13, true));
fprintf('Suma de toda la tabla: %d (la guia dice 9534)\n', sum(D(:)));
fprintf('Ruta invalida: %d (esperado 10000)\n', fobj_viajero([1 1 3:13], false));
end

%% ===================== POBLACION, CRUCE Y MUTACION =====================
function P = pob_viajero(n)
P = zeros(n, 13);
for i = 1:n, P(i,:) = randperm(13); end
end

function H = cruce_pmx(padres, nh)
% Cruce PMX: el hijo copia un tramo del padre 1 y completa el resto usando
% la correspondencia entre los tramos de ambos padres.
np = size(padres,1); L = size(padres,2); H = zeros(nh, L);
for k = 1:nh
    ij = randperm(np, 2); p1 = padres(ij(1),:); p2 = padres(ij(2),:);
    c = sort(randperm(L, 2)); a = c(1); b = c(2);
    h = zeros(1, L); h(a:b) = p1(a:b);
    for i = [1:a-1, b+1:L]
        v = p2(i);
        while any(h(a:b) == v)
            v = p2(find(p1 == v, 1));
        end
        h(i) = v;
    end
    H(k,:) = h;
end
end

function x = mut_inversion(x)
% Mutacion por inversion: se invierte el orden de un tramo aleatorio
ij = sort(randperm(numel(x), 2));
x(ij(1):ij(2)) = x(ij(2):-1:ij(1));
end

function prob = problema_viajero()
prob.fobj = @(x) fobj_viajero(x, false);
prob.pob = @(n) pob_viajero(n);
prob.cruce = @(P, nh) cruce_pmx(P, nh);
prob.mut = @(x) mut_inversion(x);
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

function imprimir_tabla(titulo, etiq, Rs)
fprintf('\n%s\n', titulo);
fprintf('| Configuracion | Tiempo (s) | Mejor f(X) | Media | Desv. est. | Mejor individuo |\n');
fprintf('|---|---|---|---|---|---|\n');
for k = 1:numel(Rs)
    r = Rs(k);
    fprintf('| %s | %.3f | %.4f | %.4f | %.4f | [%s] |\n', etiq{k}, r.t, r.mejor, r.media, r.std, sprintf('%d ', r.ind));
end
end