    Re_list = logspace(2, log10(2000), 20);
    Re_list = Re_list(1);
    delta_list = logspace(-6, 0, 20);
    delta_list = delta_list(1);
    B = eye(9);
    num_sims = 10; 
    norm_u_results = zeros(length(Re_list), length(delta_list));
    norm_d_results = zeros(length(Re_list), length(delta_list));
    n = [1;0;0;0;0;0;0;0;-1];
    %     gamma_results = zeros(length(Re_list), length(delta_list));
    %     test_u_delta = zeros(length(Re_list), length(delta_list));
    %     test_u_d_gamma = zeros(length(Re_list), length(delta_list));
    %     reachable_delta = zeros(length(Re_list), length(delta_list));
    %     lambda_min_results = zeros(length(Re_list), length(delta_list));
    Lx = 1.75*pi;
    Lz = 1.2*pi;
    global alpha
    alpha = (2*pi)/Lx;
    global Beta
    Beta = pi/2;
    global Gamma
    Gamma = 2*pi/Lz;
    global KBG
    KBG = sqrt(Beta^2+Gamma^2);
    global KAG
    KAG = sqrt(alpha^2+Gamma^2);
    global KABG
    KABG = sqrt(alpha^2+Beta^2+Gamma^2);
    deltaf = zeros(length(Re_list), length(delta_list));
    delta_max = zeros(1, length(Re_list));
    [RHS_J_mean_shear, nonlinear, u] = nonliner();
    delete(gcp('nocreate'));
    parpool(4);
    parfor ind_Re = 1:length(Re_list)
        Re = Re_list(ind_Re);
        [local_deltaf,max_norm_u_d]=LMI_Re(Re, delta_list,RHS_J_mean_shear, nonlinear, u, ind_Re, norm_u_results, norm_d_results, num_sims, B);
        deltaf(ind_Re, :) = local_deltaf;
        delta_max(ind_Re) = max(local_deltaf);
%         test_u_delta(ind_Re, :) = local_test(ind_Re, :);       
    end
    
    %     if any(delta_max == 0)
    %         ind_Re = find(delta_max);
    %     end
    
%     log_Re_list = log10(Re_list);
%     log_delta_max = log10(delta_max);
%     
%     % line fitting to the data
%     coeffs = polyfit(log_Re_list, log_delta_max, 1);
%     sigma = coeffs(1);
%     A_log = coeffs(2);
%     
%     figure;
%     loglog(Re_list, delta_max, 'o'); hold on;
%     fit_line = 10.^polyval(coeffs, log_Re_list);
%     loglog(Re_list, fit_line, '-');
%     xlabel('Re');
%     ylabel('\delta_{p}');
%     title(['Scaling exponent \sigma = ', num2str(sigma)]);
%     
%     disp(['Scaling exponent (sigma): ', num2str(sigma)]);
%     disp(['Intercept (A): ', num2str(A_log)]);

%     num_sims = 1; 
%     norm_u_results = zeros(length(Re_list), length(delta_list));
%     norm_d_results = zeros(length(Re_list), length(delta_list));


%     for ind_Re = 1:length(Re_list)
%         Re = Re_list(ind_Re);
%         for ind_delta = 1:length(delta_list)
%             delta = delta_list(ind_delta);
%             gamma = gamma_results(ind_Re, ind_delta);
%             if isnan(gamma)
%                 test_u_delta(ind_Re, ind_delta) = NaN;
%                 test_u_d_gamma(ind_Re, ind_delta) = NaN;
%                 continue;  % Skip if the LMI was infeasible
%             end
% 
%             norm_u_simulations = zeros(num_sims, 1);
%             norm_d_simulations = zeros(num_sims, 1);
%             T = 10000;
%             for sim = 1:num_sims
%                 d_forcing_func = generate_d_forcing(delta, gamma, T);
%                 [norm_u_simulations(sim), u] = compute_norm_u(Re,RHS_R_viscous ,B,RHS_J_mean_shear, d_forcing_func, T);
%                 norm_d_simulations(sim) = norm(d_forcing_func(0));  
%             end
% 
%             norm_u_results(ind_Re, ind_delta) = max(norm_u_simulations);
%             norm_d_results(ind_Re, ind_delta) = max(norm_d_simulations);
%             disp(['Re: ', num2str(Re), ', delta: ', num2str(delta), ', Norm u: ', num2str(norm_u_results(ind_Re, ind_delta)), ', Norm d: ', num2str(norm_d_results(ind_Re, ind_delta)), ' , Optimal Gamma: ', num2str(gamma)]);
%             
%             if norm_u_results(ind_Re, ind_delta) < delta
%                 test_u_delta(ind_Re, ind_delta) = 1;
%             else 
%                 test_u_delta(ind_Re, ind_delta) = 0;
%             end
%         end
%     end




function [local_deltaf,max_norm_u_d]=LMI_Re(Re,delta_list,RHS_J_mean_shear, nonlinear, u, ind_Re, norm_u_results, norm_d_results, num_sims, B)
    global Beta
    global Gamma
    global alpha
    global KBG
    global KABG
    global KAG
    local_deltaf = zeros(1, length(delta_list));
    
    linear_term = [(Beta^2)/Re;
        ((4*Beta^2)/3 + Gamma^2)/Re;
        (Beta^2+Gamma^2)/Re;
        (3*alpha^2+4*Beta^2)/(3*Re);
        (alpha^2+Beta^2)/Re;
        (3*alpha^2+4*Beta^2+3*Gamma^2)/(3*Re);
        (alpha^2+Beta^2+Gamma^2)/Re;
        (alpha^2+Beta^2+Gamma^2)/Re;
        (9*Beta^2)/Re];
    RHS_R_viscous = diag(linear_term);
    A = RHS_J_mean_shear - RHS_R_viscous;

    for ind_delta = 1:length(delta_list)
        yalmip('clear');
        delta = delta_list(ind_delta);
        delta2 = delta^2;

        [F_square] = F__square(nonlinear, u);

        %             s = sdpvar(length(F_square), 1);
        %             s_bound = zeros(size(A));
        %             for m_ind = 1:length(F_square)
        %                  s_bound = s_bound + s(m_ind) * delta2 * double(F_square{m_ind});
        %             end
        %
        %             if length(s) <= length(A)
        %                 diag_s = diag(s);
        %             else
        %                 diag_s = diag(s(1:length(A))) + eye(size(A)) * s(length(A) + 1);
        %             end
        
        [result_yalmip, P, Gamma2] = LMI(A, B, F_square, delta2);
        
        if result_yalmip.problem == 0 && all(all(value(Gamma2) > 0))
                gamma_optimal = sqrt(value(Gamma2));
                gamma_results(ind_Re, ind_delta) = gamma_optimal;
                reachable_delta(ind_Re,ind_delta)= delta;
        else
                gamma_results(ind_Re, ind_delta) = NaN; 
                reachable_delta(ind_Re,ind_delta)=NaN;
        end
            
        disp(['Re: ', num2str(Re), ', delta: ', num2str(delta), ', Optimal gamma: ', num2str(gamma_results(ind_Re, ind_delta))]);
       
        gamma = gamma_results(ind_Re, ind_delta);
        
        if isnan(gamma)
                test_u_delta(ind_Re, ind_delta) = NaN;
                %test_u_d_gamma(ind_Re, ind_delta) = NaN;
                continue;  % Skip if the LMI was infeasible
        end
        
       max_norm_u_d_tmp = test_u_delta_func(gamma, Re, delta, num_sims, norm_u_results, norm_d_results, A, B, P);
       max_norm_u_d{ind_Re, ind_delta} = max_norm_u_d_tmp;


%         if result_yalmip.problem == 0 && all(all(value(Gamma2) > 0))
%                 gamma_optimal = sqrt(value(Gamma2));
%                 gamma_results(ind_Re, ind_delta) = gamma_optimal;
%                 reachable_delta(ind_Re,ind_delta)= delta;
%         else
%                 gamma_results(ind_Re, ind_delta) = NaN; 
%                 reachable_delta(ind_Re,ind_delta)=NaN;
%         end
%             
%         disp(['Re: ', num2str(Re), ', delta: ', num2str(delta), ', Optimal gamma: ', num2str(gamma_results(ind_Re, ind_delta))]);
%         
%         gamma = gamma_results(ind_Re, ind_delta);
%         if isnan(gamma)
%                 test_u_delta(ind_Re, ind_delta) = NaN;
%                 %test_u_d_gamma(ind_Re, ind_delta) = NaN;
%                 continue;  % Skip if the LMI was infeasible
%         end
% 
%         norm_u_simulations = zeros(num_sims, 1);
%         norm_d_simulations = zeros(num_sims, 1);
%         T = 10000;
%         for sim = 1:num_sims
%                 d_forcing_func = generate_d_forcing(delta, gamma, T);
%                 [norm_u_simulations(sim), u] = compute_norm_u(Re,RHS_R_viscous ,B,RHS_J_mean_shear, d_forcing_func, T);
%                 norm_d_simulations(sim) = norm(d_forcing_func(0));  
%         end
% 
%         norm_u_results(ind_Re, ind_delta) = max(norm_u_simulations);
%         norm_d_results(ind_Re, ind_delta) = max(norm_d_simulations);
%         disp(['Re: ', num2str(Re), ', delta: ', num2str(delta), ', Norm u: ', num2str(norm_u_results(ind_Re, ind_delta)), ', Norm d: ', num2str(norm_d_results(ind_Re, ind_delta)), ' , Optimal Gamma: ', num2str(gamma)]);
%             
%         if norm_u_results(ind_Re, ind_delta) < delta
%                 test_u_delta(ind_Re, ind_delta) = 1;
%         else 
%                 test_u_delta(ind_Re, ind_delta) = 0;
%         end
        
    end

end

function [max_norm_u_d] = test_u_delta_func(gamma, Re, delta, num_sims, norm_u_results, norm_d_results, A, B, P)
        norm_u_simulations = zeros(num_sims, 1);
        norm_d_simulations = zeros(num_sims, 1);
        T = 10000;
        omgea_list = logspace(-3, 3, 20);
        for ind_omega = 1:length(omgea_list)
            omega = omgea_list(ind_omega);
            d_forcing_func = generate_d_forcing(delta, gamma, T, omega);
            [~, ~, ~, ~, max_SVD, ~] = compute_norm_u(Re, A, B, d_forcing_func, T, omega);
            max_SVD_list(ind_omega) = max_SVD;
            
            for sim = 1:num_sims
                    d_forcing_func = generate_d_forcing(delta, gamma, T, omega);
                    
                    [norm_u_simulations(sim), u, ~, ~,~, t] = compute_norm_u(Re, A, B, d_forcing_func, T, omega);
                    
                    dt = diff(t);
                    dt = dt(1);
                    diag(u*value(P)*u');
                    
                    norm_d_simulations(sim) = sqrt(sum(sum(d_forcing_func(t').^2))*dt);
                    norm_u_d(sim) = norm_u_simulations(sim)/norm_d_simulations(sim);
                    
%                     figure(3)
%                     plot(t, u(:,1)); hold on
%                     plot(t, u(:,2)); hold on
            end
            %norm_u_d(find(norm_u_d== inf)) = -inf;
            max_norm_u_d(ind_omega) = max(norm_u_d);
        end

%          norm_u_results(ind_Re, ind_delta) = max(norm_u_simulations);
%          norm_d_results(ind_Re, ind_delta) = max(norm_d_simulations);
%          disp(['Re: ', num2str(Re), ', delta: ', num2str(delta), ', Norm u: ', num2str(norm_u_results(ind_Re, ind_delta)), ', Norm d: ', num2str(norm_d_results(ind_Re, ind_delta)), ' , Optimal Gamma: ', num2str(gamma)]);
%                  
%          if norm_u_results(ind_Re, ind_delta) < delta
%                      test_u_delta(ind_Re, ind_delta) = 1;
%          else 
%                      test_u_delta(ind_Re, ind_delta) = 0;
%          end

%        figure(1)
%        plot(log10(omgea_list), log10(max_SVD_list)); hold on
%        plot(log10(omgea_list), log10(max_norm_u_d));

end

function [RHS_J_mean_shear,nonlinear, u] = nonliner()
    global Beta
    global Gamma
    global alpha
    global KBG
    global KABG
    global KAG
    u = sym('u', [9, 1]);
    term1 = -sqrt(3/2)*Beta*Gamma*u(6)*u(8)/KABG+sqrt(3/2)*Beta*Gamma*u(2)*u(3)/KBG;
    term2 = (5/3)*sqrt(2/3)*Gamma^2*u(4)*u(6)/KAG-Gamma^2*u(5)*u(7)/(sqrt(6)*KAG) ...
        -alpha*Beta*Gamma*u(5)*u(8)/(sqrt(6)*KAG*KABG)-sqrt(3/2)*Beta*Gamma*u(1)*u(3)/KBG-sqrt(3/2)*Beta*Gamma*u(3)*u(9)/KBG;
    term3 = 2*alpha*Beta*Gamma*(u(4)*u(7)+u(5)*u(6))/(sqrt(6)*KAG*KBG)+(Beta^2*(3*alpha^2+Gamma^2)-3*Gamma^2*(alpha^2+Gamma^2))*u(4)*u(8)/(sqrt(6)*KAG*KBG*KABG);
    term4 = -alpha*u(1)*u(5)/sqrt(6)-10*alpha^2*u(2)*u(6)/(3*sqrt(6)*KAG)  ...
        -sqrt(3/2)*alpha*Beta*Gamma*u(3)*u(7)/KAG*KBG-sqrt(3/2)*alpha^2*Beta^2*u(3)*u(8)/KAG*KBG*KABG-alpha*u(5)*u(9)/sqrt(6);
    term5 =  alpha*u(1)*u(4)/sqrt(6)+alpha^2*u(2)*u(7)/(sqrt(6)*KAG)-alpha*Beta*Gamma*u(2)*u(8)/(sqrt(6)*KAG*KABG)+alpha*u(4)*u(9)/sqrt(6)+2*alpha*Beta*Gamma*u(3)*u(6)/(sqrt(6)*KAG*KBG);
    term6 =  alpha*u(1)*u(7)/sqrt(6)+sqrt(3/2)*Beta*Gamma*u(1)*u(8)/KABG  ...
        +10*(alpha^2-Gamma^2)*u(2)*u(4)/(KAG*3*sqrt(6))-2*sqrt(2/3)*u(3)*u(5)*alpha*Beta*Gamma/(KAG*KBG)+alpha*u(7)*u(9)/sqrt(6)+sqrt(3/2)*Beta*Gamma*u(8)*u(9)/KABG;
    term7 = -alpha*(u(1)*u(6)+u(6)*u(9))/sqrt(6)+(Gamma^2-alpha^2)*u(2)*u(5)/(sqrt(6)*KAG)+alpha*Beta*Gamma*u(3)*u(4)/(sqrt(6)*KAG*KBG);
    term8 = 2*alpha*Beta*Gamma*u(2)*u(5)/(sqrt(6)*KAG*KABG)+Gamma^2*(3*alpha^2-Beta^2+3*Gamma^2)*u(3)*u(4)/(sqrt(6)*KAG*KBG*KABG);
    term9 = sqrt(3/2)*Beta*Gamma*u(2)*u(3)/KBG-sqrt(3/2)*Beta*Gamma*u(6)*u(8)/KABG;
    nonlinear= [term1;
        term2;
        term3;
        term4;
        term5;
        term6;
        term7;
        term8;
        term9];
    
    
    
    nonlinear_gradient = sym(zeros(length(u), length(u)));
    for i = 1:length(nonlinear)
        nonlinear_gradient(i, :) = gradient(nonlinear(i), u);
    end
    %nonlinear_gradinet = gradient(nonlinear, u);
    a_bar = [1;0;0;0;0;0;0;0;0];
    nonlinear_gradient_sub = double(subs(nonlinear_gradient, u, a_bar));
    RHS_J_mean_shear = nonlinear_gradient_sub;
end



function [F_square] = F__square(nonlinear, u)

    F_square = cell(1, length(nonlinear));
    for n_ind = 1:length(nonlinear)
        %initialize F{n_ind}
        F{n_ind} = sym(zeros(length(u), length(u)));
        for x_ind = 1:length(u)
            for y_ind = 1:length(u)
                F{n_ind}(x_ind, y_ind) = 1/2 * diff(diff(nonlinear(n_ind), u(x_ind)), u(y_ind));
            end
        end
    
        %convert F{n_ind} to double if it's symbolic
        [V, D] = eig(double(F{n_ind}));
        %store F_square as double
        F_square{n_ind} = double(V * D^2 /(V));
    end

end



function[result_yalmip, P, Gamma2] = LMI(A, B, F_square, delta2)
            yalmip('clear');
            sdp_options = sdpsettings('solver', 'mosek');
            options.bisection.absgaptol = 1e-16;
            [s, diag_s, s_bound] = s_diag_s(F_square, A, delta2);

            Gamma2 = sdpvar(1,1);
            P = sdpvar(size(A, 1), size(A, 2));
            G11 = A' * P + P * A + s_bound;
            G12 = P;
            G13 = P * B;
            G22 = -diag_s;
            G33 = -eye(9);
            G = [G11, G12, G13;
                 G12', G22, zeros(9, 9);
                 G13', zeros(9, 9), G33];
            V_ineq = P - Gamma2 * eye(9);
            F = [V_ineq >= 0, G <= 0, s >= 0];
            
            % Solve the optimization problem
            result_yalmip = optimize(F, -Gamma2, sdp_options);
end



function[s, diag_s, s_bound] = s_diag_s(F_square, A, delta2)
    s = sdpvar(length(F_square), 1);
    s_bound = zeros(size(A));
    for m_ind = 1:length(F_square)
        s_bound = s_bound + s(m_ind) * delta2 * double(F_square{m_ind});
    end
    
    if length(s) <= length(A)
        diag_s = diag(s);
    else
        diag_s = diag(s(1:length(A))) + eye(size(A)) * s(length(A) + 1);
    end
end




function d_forcing_func = generate_d_forcing(delta, gamma, T, omega)
        phi = rand(1,1)*2*pi;
        d_forcing_func = @(t) (delta*gamma)/sqrt(T)*[sin(omega*t)*cos(phi);sin(omega*t)*sin(phi);sin(omega*t)*cos(phi); sin(omega*t)*sin(phi);sin(omega*t)*cos(phi); sin(omega*t)*sin(phi); sin(omega*t)*cos(phi); sin(omega*t)*sin(phi);sin(omega*t)*cos(phi)]; 
end


function [norm_u, u, G, sys, max_SVD, t] = compute_norm_u(Re, A, B, d_forcing_func, T, omega)
        
        u0 = zeros(size(A, 1), 1);
        tspan = linspace(0, T, 10000);  
   
        I = eye(size(A));
        C = eye(9);  
        D = zeros(9, 9);  
    
        G = C * inv(1i * omega * I - A) * B;
        max_SVD = max(svd(G));
        sys = ss(A, B, C, D);
        %ode_system = @(t, u) A * u + B * d_forcing_func(t); 
        %u0 = zeros(9,1);
    
        dt = diff(tspan);
        dt = dt(1);
    
        %[t, u] = ode45(@(t,u) ode_system_123(t, u,A,B, d_forcing_func), tspan, u0);
        [t, u] = rk4_solver(@(t,u) ode_system_123(t, u,A,B, d_forcing_func), tspan, u0);

        
        %compute the time-varying average norm_u
        %integrand = @(t) interp1(tspan, sum(u.^2, 2), t, 'linear', 'extrap');
        %norm_u = sqrt(integral(integrand, 0, tspan(end), 'ArrayValued', true) / tspan(end));
        
        %old norm integrate over time
        % dt=diff(tspan);
        % dt=dt(1);
        % u_square=sum(u.^2,2);
        % norm_u=sqrt(sum(u_square(2:end).*dt)/max(tspan));
    
        %1 norm maximal over time
        u_square=sqrt(sum(u.^2,2));
        norm_u = sqrt(sum(u_square.^2)*dt);
        
        %norm_u=max(u_square);    
 end


 function [du_dt] = ode_system_123(t, u, A, B, d_forcing_func)
global Beta
global Gamma
global alpha
global KBG
global KABG
global KAG

% force=randn(9,1);
% forcing=u_upper_bound*force/norm(force);

term1 = -sqrt(3/2)*Beta*Gamma*u(6)*u(8)/KABG+sqrt(3/2)*Beta*Gamma*u(2)*u(3)/KBG;
term2 = (5/3)*sqrt(2/3)*Gamma^2*u(4)*u(6)/KAG-Gamma^2*u(5)*u(7)/(sqrt(6)*KAG) ...
    -alpha*Beta*Gamma*u(5)*u(8)/(sqrt(6)*KAG*KABG)-sqrt(3/2)*Beta*Gamma*u(1)*u(3)/KBG-sqrt(3/2)*Beta*Gamma*u(3)*u(9)/KBG;
term3 = 2*alpha*Beta*Gamma*(u(4)*u(7)+u(5)*u(6))/(sqrt(6)*KAG*KBG)+(Beta^2*(3*alpha^2+Gamma^2)-3*Gamma^2*(alpha^2+Gamma^2))*u(4)*u(8)/(sqrt(6)*KAG*KBG*KABG);
term4 = -alpha*u(1)*u(5)/sqrt(6)-10*alpha^2*u(2)*u(6)/(3*sqrt(6)*KAG)  ...
    -sqrt(3/2)*alpha*Beta*Gamma*u(3)*u(7)/KAG*KBG-sqrt(3/2)*alpha^2*Beta^2*u(3)*u(8)/KAG*KBG*KABG-alpha*u(5)*u(9)/sqrt(6);
term5 =  alpha*u(1)*u(4)/sqrt(6)+alpha^2*u(2)*u(7)/(sqrt(6)*KAG)-alpha*Beta*Gamma*u(2)*u(8)/(sqrt(6)*KAG*KABG)+alpha*u(4)*u(9)/sqrt(6)+2*alpha*Beta*Gamma*u(3)*u(6)/(sqrt(6)*KAG*KBG);
term6 =  alpha*u(1)*u(7)/sqrt(6)+sqrt(3/2)*Beta*Gamma*u(1)*u(8)/KABG  ...
    +10*(alpha^2-Gamma^2)*u(2)*u(4)/(KAG*3*sqrt(6))-2*sqrt(2/3)*u(3)*u(5)*alpha*Beta*Gamma/(KAG*KBG)+alpha*u(7)*u(9)/sqrt(6)+sqrt(3/2)*Beta*Gamma*u(8)*u(9)/KABG;
term7 = -alpha*(u(1)*u(6)+u(6)*u(9))/sqrt(6)+(Gamma^2-alpha^2)*u(2)*u(5)/(sqrt(6)*KAG)+alpha*Beta*Gamma*u(3)*u(4)/(sqrt(6)*KAG*KBG);
term8 = 2*alpha*Beta*Gamma*u(2)*u(5)/(sqrt(6)*KAG*KABG)+Gamma^2*(3*alpha^2-Beta^2+3*Gamma^2)*u(3)*u(4)/(sqrt(6)*KAG*KBG*KABG);
term9 = sqrt(3/2)*Beta*Gamma*u(2)*u(3)/KBG-sqrt(3/2)*Beta*Gamma*u(6)*u(8)/KABG;
nonlinear= [term1;
    term2;
    term3;
    term4;
    term5;
    term6;
    term7;
    term8;
    term9];

du_dt = A * u + B * d_forcing_func(t)+ nonlinear;
end


function [t, u] = rk4_solver(ode_func, tspan, u0)
    % Number of time steps
    N = length(tspan);
    % Preallocate u matrix (number of time steps by 9 states)
    u = zeros(N, length(u0));
    u(1, :) = u0;  % Set initial condition

    % Time stepping
    for i = 1:N-1
        h = tspan(i+1) - tspan(i);  % Calculate step size
        % Perform one RK4 step and store result in the next row
        u(i+1, :) = rk4_step(u(i, :), h, ode_func);
    end

    t = tspan;  % Output the time vector
end

function u_next = rk4_step(u_current, h, ode_func)
    % Calculate the RK4 coefficients
    k1 = h * ode_func(0, u_current')';
    k2 = h * ode_func(0, (u_current' + 0.5 * k1'))';
    k3 = h * ode_func(0, (u_current' + 0.5 * k2'))';
    k4 = h * ode_func(0, (u_current' + k3'))';

    % Update the state using RK4 formula
    u_next = u_current + (k1 + 2*k2 + 2*k3 + k4) / 6;
end
