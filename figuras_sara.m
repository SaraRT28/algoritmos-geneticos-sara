% figuras_sara.m  (en la misma carpeta que res_s_viajero.mat, res_s_ackley.mat y res_s_inventario.mat)
V = load('res_s_viajero.mat'); A = load('res_s_ackley.mat'); I = load('res_s_inventario.mat');
D = {V, A, I};
nom     = {'Viajero (min)','Ackley (d=5)','Inventarios (Z)'};
metodos = {'elitismo','torneo','ruleta'};
etq     = {'Elitismo','Torneo','Ruleta'};

% ---------- 1) Selección: resultado de cada corrida y su media ----------
figure('Color','w','Position',[100 100 1000 340]);
for p = 1:3
    subplot(1,3,p); hold on;
    for k = 1:3
        y = D{p}.S1.(metodos{k}).Y;
        plot(k + 0.12*(rand(size(y)) - 0.5), y, 'o', 'Color', [0.2 0.45 0.8]);
        plot(k, mean(y), 'k_', 'MarkerSize', 22, 'LineWidth', 2);
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', etq, 'XLim', [0.5 3.5]); grid on;
    ylabel('Mejor fitness por corrida'); title(nom{p});
end
exportgraphics(gcf, 'fig_s_seleccion.png', 'Resolution', 300);

% ---------- 2) Convergencia promedio ----------
figure('Color','w','Position',[100 100 1000 340]);
for p = 1:3
    subplot(1,3,p); hold on;
    for k = 1:3
        y = D{p}.S1.(metodos{k}).curva;
        if p == 2, y = y - (3 - exp(1)) + 1e-12; end
        plot(0:numel(y)-1, y, 'LineWidth', 1.5);
    end
    if p == 2, set(gca, 'YScale', 'log'); ylabel('f - (3-e)'); else, ylabel('Mejor fitness'); end
    grid on; xlabel('Iteración'); title(nom{p});
    if p == 1, legend(etq, 'Location', 'northeast'); end
end
exportgraphics(gcf, 'fig_s_convergencia.png', 'Resolution', 300);

% ---------- 3) Iteraciones: calidad (arriba) y tiempo (abajo) ----------
figure('Color','w','Position',[100 100 1000 520]);
for p = 1:3
    R = D{p}.S2; it = D{p}.cfg.iterLista;
    subplot(2,3,p);   errorbar(it, [R.media], [R.std], '-o', 'LineWidth', 1.2); grid on;
    ylabel('Mejor fitness'); title(nom{p});
    subplot(2,3,p+3); plot(it, [R.t], '-s', 'LineWidth', 1.2); grid on;
    xlabel('Iteraciones'); ylabel('Tiempo (s)');
end
exportgraphics(gcf, 'fig_s_iteraciones.png', 'Resolution', 300);

% ---------- 4) Probabilidad de mutación ----------
figure('Color','w','Position',[100 100 1000 340]);
for p = 1:3
    R = D{p}.S3.pm;
    subplot(1,3,p); errorbar(D{p}.cfg.PmLista, [R.media], [R.std], '-o', 'LineWidth', 1.2); grid on;
    xlabel('P_m'); ylabel('Mejor fitness'); title(nom{p});
end
exportgraphics(gcf, 'fig_s_pm.png', 'Resolution', 300);

% ---------- 5) Magnitud de la mutación (fracción del rango) ----------
figure('Color','w','Position',[100 100 700 340]);
R = A.S3.sig; subplot(1,2,1);
errorbar(A.cfg.fracLista, [R.media], [R.std], '-o', 'LineWidth', 1.2); grid on;
xlabel('\sigma / rango'); ylabel('Mejor fitness'); title('Ackley (d=5)');
R = I.S3.sig; subplot(1,2,2);
errorbar(I.cfg.fracLista, [R.media], [R.std], '-o', 'LineWidth', 1.2); grid on;
xlabel('\sigma / rango'); ylabel('Mejor fitness'); title('Inventarios');
exportgraphics(gcf, 'fig_s_sigma.png', 'Resolution', 300);

% ---------- 6) Dimensiones ----------
figure('Color','w','Position',[100 100 700 340]);
subplot(1,2,1); semilogx(A.cfg.dimLista, [A.S4.rmse], '-o', 'LineWidth', 1.2); grid on;
xlabel('Dimensiones d'); ylabel('RMSE');
subplot(1,2,2); semilogx(A.cfg.dimLista, [A.S4.t], '-s', 'LineWidth', 1.2); grid on;
xlabel('Dimensiones d'); ylabel('Tiempo (s)');
exportgraphics(gcf, 'fig_s_dimensiones.png', 'Resolution', 300);

% ---------- 7) Función de Ackley: curva en d=1 y mapa de contorno en d=2 ----------
a = 20; b = 0.2; c = 2*pi;
x = linspace(-10, 10, 2000); f1 = -a*exp(-b*abs(x)) - exp(cos(c*x)) + a + 3;
[X, Y] = meshgrid(linspace(-10, 10, 300));
Z = -a*exp(-b*sqrt((X.^2 + Y.^2)/2)) - exp((cos(c*X) + cos(c*Y))/2) + a + 3;
figure('Color','w','Position',[100 100 900 380]);
subplot(1,2,1); plot(x, f1, 'LineWidth', 1.2); grid on;
xlabel('x'); ylabel('f(x)'); title('(a) Ackley, d = 1');
subplot(1,2,2); contourf(X, Y, Z, 30); colorbar;
xlabel('x_1'); ylabel('x_2'); title('(b) Ackley, d = 2');
exportgraphics(gcf, 'fig_s_ackley.png', 'Resolution', 300);

% ---------- 8) Diagrama de flujo ----------
figure('Color','w','Position',[100 100 520 640]);
axes('Position', [0 0 1 1]); axis off; hold on; xlim([0 1]); ylim([0 1]);
paso = {'1. Población inicial aleatoria','2. Evaluar el fitness','3. Seleccionar n/2 padres', ...
        '4. Cruzar parejas de padres','5. Mutar hijos con prob. P_m','6. Evaluar hijos y actualizar el mejor', ...
        '¿Se completaron las iteraciones?','7. Resultado: mejor individuo'};
y = linspace(0.93, 0.07, 8);
for k = 1:8
    rectangle('Position', [0.12 y(k)-0.035 0.62 0.07], 'Curvature', 0.25, 'FaceColor', [0.93 0.97 0.9]);
    text(0.43, y(k), paso{k}, 'HorizontalAlignment', 'center', 'FontSize', 9);
    if k < 8, quiver(0.43, y(k)-0.035, 0, y(k+1)-y(k)+0.07, 0, 'k', 'MaxHeadSize', 1.2); end
end
plot([0.74 0.88 0.88 0.74], [y(7) y(7) y(3) y(3)], 'k-');
quiver(0.82, y(3), -0.08, 0, 0, 'k', 'MaxHeadSize', 2);
text(0.9, mean([y(7) y(3)]), 'No', 'Rotation', 90, 'HorizontalAlignment', 'center');
text(0.45, y(7)-0.057, 'Sí');
xlim([0 1]); ylim([0 1]);
exportgraphics(gcf, 'fig_s_flujo.png', 'Resolution', 300);

% ---------- individuos completos de Ackley (para el anexo) ----------
for k = 1:numel(A.cfg.dimLista)
    fprintf('d = %d, f = %.4f\n', A.cfg.dimLista(k), A.S4(k).mejor);
    fprintf('%.4f ', A.S4(k).ind); fprintf('\n\n');
end