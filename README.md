# Algoritmos genéticos: viajero, Ackley e inventarios

Trabajo de Sara Ramirez Toro, Bioingeniería, Unidad 3 (Optimización metaheurística), Universidad Católica de Oriente.

Implementación en MATLAB de un algoritmo genético aplicado a tres problemas de minimización: el problema del viajero (13 ubicaciones), la función de Ackley y un problema de inventarios con restricción de presupuesto.

## Archivos

- `ag_viajero_s.m`: problema del viajero (cruce PMX, mutación por inversión).
- `ag_ackley_s.m`: función de Ackley (cruce SBX, mutación gaussiana).
- `ag_inventario_s.m`: inventarios (cruce uniforme con reparación, mutación gaussiana).
- `figuras_sara.m`: genera las figuras del informe a partir de los resultados.
- `salida_s_*.txt`: tablas impresas por cada corrida real.

## Cómo ejecutarlo

Requiere MATLAB (la función `exportgraphics` necesita R2020a o posterior). Con los archivos en una misma carpeta:

    ag_viajero_s('verificar')    % comprueba la función objetivo
    ag_viajero_s                 % corrida real (8 corridas por escenario)
    ag_ackley_s
    ag_inventario_s
    figuras_sara                 % figuras del informe

Cada script guarda sus resultados en `res_s_*.mat` e imprime las tablas en pantalla. Con `ag_..._s(2)` se hace una prueba rápida de 2 corridas.

## Configuración

Población de 60 individuos, 8 corridas por escenario y semilla `rng(31)`.
