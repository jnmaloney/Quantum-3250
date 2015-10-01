function example_discord()
% discord Quantum Discord
%
%
%
%

% ref: PhysRevLett.88.017901.pdf
% Quantum Discord: A Measure of the Quantumness of Correlations


    N = 2;
    quantum_base;

    
    
    z_range = linspace(0.1, 0.9, 42);
    theta_range = linspace(-pi, pi, 42);
    [theta_mesh, z_mesh] = meshgrid(theta_range, z_range);
    q = zeros(42, 42);

    for z = 1:42
        for theta = 1:42
    
            % Create rho_A
            %rho_A = 0.5 * [1 0; 0 1];
            %rho_A = [1 0; 0 0];
            psi_0 = [1; 1] / sqrt(2);
            rho_A = psi_0 * psi_0';
            %rho_A = [0.5 0.1; 0.1 0.5];


              
            % Create rho_AB
            rho_AB = 0.5 * [1 0; 0 1] +...
                     z_range(z)/2 * [0 1; 1 0];
                 
                 
           % Calculate discord
           q(z, theta) = discord(rho_A, rho_AB, theta_range(theta));
           
           
        end
    end
    
    figure;
    colormap(spring);
    surf(theta_mesh, z_mesh, q)


end

