function bloch(N, M, t, H, J_x, J_y, J_z, psi0, rho_A)
% 
%   
    quantum_base;

    % Calculate values
    theta_steps = 15;
    phi_steps = 30;
    q = zeros(theta_steps, phi_steps);
    for i = 1:theta_steps
        for j = 1:phi_steps

            theta = pi * (i - 1) / (theta_steps - 1);
            phi = 2 * pi * (j - 1) / (phi_steps - 1);

            q(i, j) = abs(Q(theta, phi));

        end
    end



    % Bloch space
    theta_matrix = linspace(0, pi, theta_steps);
    phi_matrix = linspace(0, 2*pi, phi_steps);
    [thetamesh, phimesh] = meshgrid(theta_matrix, phi_matrix);



    % Calculate eigenvalues of J_y for P(J_y)
    [V, D] = eig(J_y);
    P_y = zeros(N + 1, 1);

    for i = 1:N + 1
        % ith eigenvector
        v = V(:,i);

        P_y(i) = v' * rho_A * v;
    end

    % Bloch sphere

    radius = N / 2;
    x = radius * sin(thetamesh) .* cos(phimesh);
    y = radius * sin(thetamesh) .* sin(phimesh);
    z = radius * cos(thetamesh);

    surf(x, y, z, q')
    shading 'interp'
    view([0,1,0])
    axis on
    axis tight
    axis equal
    grid off;
    xlabel('Jz')
    zlabel('Jy')
    set(gca, 'XTickLabel', '');
    set(gca, 'YTickLabel', '');
    set(gca, 'ZTickLabel', '');

    %subplot(2, 1, 2)
    % Histogram

    %bar(-N/2:N/2, P_y, 'k');
    %xl = xlabel('J_y');
    %yl = ylabel('P(J_y)');
    %set([xl yl], 'interpreter', 'tex');

end

