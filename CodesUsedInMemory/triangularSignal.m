function y = triangularSignal(p)
   
    if mod(p, 2) == 0
        error('El número de puntos debe ser impar.');
    end
    
    % Calcular el número de puntos a cada lado del pico de la triangular
    n_points_side = (p - 1) / 2;
    
    % Generar puntos para la parte ascendente de la triangular
    x_ascend = linspace(-1, 0, n_points_side + 1);
    y_ascend = x_ascend + 1;
    
    % Generar puntos para la parte descendente de la triangular
    x_descend = linspace(0, 1, n_points_side + 1);
    y_descend = 1 - x_descend;
    
    % Combinar las partes ascendente y descendente para formar la triangular completa
    x = [x_ascend(1:end-1) x_descend];
    y = [y_ascend(1:end-1) y_descend];
    
    % Normalizar para que la integral bajo la curva sea igual a 1
    y = y / trapz(x, y);
end
