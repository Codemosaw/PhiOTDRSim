function tren_triangular = TrenTriangular(d)
    if mod(d, 2) == 0
        error('El número de puntos debe ser impar.');
    end
    
    % Calcular el número de puntos para cada triangular
    puntos_triangular = (d + 1) / 5;
    
    % Generar cada triangular y concatenarlas
    tren_triangular = [];
    for i = 1:5
        tren_triangular = [tren_triangular, triangularSignal(puntos_triangular)];
    end
end

