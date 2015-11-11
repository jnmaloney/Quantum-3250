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



    figure
    title('Husimi Q-Function of a coherent spin state')
    colormap(spring);

    contourf(phimesh, -thetamesh, q')
    ylabel('\phi phi', 'FontSize', 14);
    xlabel('\theta theta', 'FontSize', 14);
    set(gca,'YTick',[0 pi/2 pi])
    set(gca,'YTickLabel',{'0','\pi/2','\pi'})
    set(gca,'XTick',[0  pi/2 pi  3*pi/2 2*pi])
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})


    figure

    %map = linspace(0.3, 1.0, 100);
    %map = [0 * map; 0.45 * map; 0.74 * map]';
    %colormap(map);
    colormap(spring);
    %title(sprintf('Husimi-Q function and J_y projection for N=%l atoms, M=%i photons and gt=%i', N, M, t));



    subplot(2, 1, 1)
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

    subplot(2, 1, 2)
    % Histogram

    bar(-N/2:N/2, P_y, 'k');
    xl = xlabel('J_y');
    yl = ylabel('P(J_y)');
    set([xl yl], 'interpreter', 'tex');



    % subplot(2,2,2)
    %     
    % contourf(phimesh, -thetamesh, q')
    % colormap(spring);
    % %imagesc(q);
    % ylabel('\phi phi', 'FontSize', 14);
    % xlabel('\theta theta', 'FontSize', 14);
    % set(gca,'YTick',[0 pi/2 pi])
    % set(gca,'YTickLabel',{'0','\pi/2','\pi'})
    % set(gca,'XTick',[0  pi/2 pi  3*pi/2 2*pi])
    % set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})



    % subplot(2,2,4)
    % 
    % %histogram(varJy);
    % 
    % gt_steps = 100;
    % F_A = zeros(gt_steps, 1);
    % F_AB = zeros(gt_steps, 1);
    % gt_space = linspace(0, 0.15, gt_steps);
    % for i = 1:gt_steps
    %     gt = gt_space(i);
    %     psi = U(H_m, gt) * alpha(0, 0);
    %     %F_A = 4 * var(G_A, psi);
    %     F_A(i) = 4 * var(J_y, psi);
    %     F_AB(i) = 4 * (var(J_y^2, psi) - var(J_y, psi)^2);
    % end
    % 
    % plot(F_A, 'r');
    % hold on
    % plot(F_AB, 'b');



end

