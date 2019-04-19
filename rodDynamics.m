function [g_next,xi_next,eta_next,mu_next,lambda_next] = rodDynamics(g,xi,eta,mu,lambda,F)
    %the variational integrator with implicit stepping
    %appears like there might be a subtle point in stepping g one step
    %ahead of the others
    
    
    %expecting columns of a and b to be the values
    [six,N] = size(eta);
    %size(b);
    
    D = 1e-2;
    A = pi/4*D^2;
    I = pi/64*D^4;
    J = 2*I;
    E = 100e3;
    G = E/3;
    K = diag([E*I,E*I,G*J,G*A,G*A,E*A]);
    rho = 1000;
    M = rho*diag([I,I,J,A,A,A]);
    vis = 300;
    V = vis*diag([3*I,3*I,J,A,A,3*A]);
    
    %a cable force
    %F = -1;
    N_act = length(F);
    r = @(i) [D/4*cos(i*2*pi/N_act);D/4*sin(i*2*pi/N_act);0];
    %r = [D/4;0;0];  
    
    g0 = eye(4);
    eta0 = [0;0;0;0;0;0];
    WL = [0;0;0;0;0;0];
    for k=1:N_act
        WL = WL+[-skew([0;0;F(k)])*r(k);0;0;F(k)];
    end
    xi_ref = [0;0;0;0;0;1];
    xi_tip = K\WL+xi_ref; %free end
    
    L = 10e-2;
    ds = L/(N-1);
    dt = 0.1*ds/sqrt(E/rho);
%     if isnan(dt) 
%         dt = 0.1*ds/sqrt(E/rho);
%     end
    
    %max_wave = max(max(xi_dot))
    
%     for i=1:N
%         (M*eta(:,i))\(K*(xi(:,i)-xi_ref))
%     end
    
    %dt = ds/max_wave
%     (xi-xi_prev)/dt
%     xi_dot = zeros(6,N);
%     for i=1:N-1
%         xi_dot(:,i) = (eta(:,i+1)-eta(:,i))/ds + ad(xi(:,i))*eta(:,i);
%     end
%     xi_dot
%     max_wave = max(max(xi_dot))
%     dt = ds/max_wave
%     if max_wave == 0 
%         dt = 0.1*ds/sqrt(E/rho);
%     end
    
    dTauInv = @(x) eye(6)-ad(x)/2;
    tau = @(x) expm(se(x));
    itau = @(x) unse(logm(x));
    tauSO = @(x) expm(skew(x));
    itauSO = @(x) unskew(logm(x));
    
    %grav = [-9.81;0;0];
    grav = [0;0;-9.81];
    %grav = [0;0;0];
    
    
%     %compute g_next from eta
%     g_next = zeros(12,N);
%     g_next(:,1) = [1,0,0,0,1,0,0,0,1,0,0,0];
%     for i=2:N %can skip base since it is fixed
%         G = [reshape(g(1:9,i),3,3),g(10:12,i);0,0,0,1];
%         G_n = G*tau(dt*eta(:,i));
%         g_next(:,i) = [reshape(G_n(1:3,1:3),9,1);G_n(1:3,4)];
%     end
%     
%     %get xi_next from g
%     xi_next = zeros(6,N);
%     xi_next(:,end) = xi_tip;
%     for i=1:N-1
%         G_n = [reshape(g(1:9,i+1),3,3),g(10:12,i+1);0,0,0,1];
%         G = [reshape(g(1:9,i),3,3),g(10:12,i);0,0,0,1];
%         xi_next(:,i) = unse(logm(G\G_n)/ds);
%     end
    
    %step xi
    xi_next = zeros(6,N);
    xi_next(:,end) = xi_tip;
    for i=1:N-1
        %omega_next = itauSO(tauSO(-dt*eta(1:3,i))*tauSO(ds*xi(1:3,i))*tauSO(dt*eta(1:3,i+1)))/ds;
        %xi_next(:,i) = [omega_next;xi_ref(4:6)];
        xi_next(:,i) = itau(tau(-dt*eta(:,i))*tau(ds*xi(:,i))*tau(dt*eta(:,i+1)))/ds;
    end
        
    %step g
    g_next = zeros(12,N);
    g_next(:,1) = [1,0,0,0,1,0,0,0,1,0,0,0];
    G = eye(4);
    for i=2:N
        G = G*tau(ds*xi_next(:,i-1));
        g_next(:,i) = [reshape(G(1:3,1:3),9,1);G(1:3,4)];
    end
    
    %conjugate varibles (mommentum and stress)
    %slightly more efficient to pass around, but for testing not important
    %mu = zeros(6,N);
    lambda_next = zeros(6,N);
    for j=1:N
        %mu(:,j) = dTauInv(dt*eta(:,j))'*M*eta(:,j);
        lambda_next(:,j) = dTauInv(ds*xi_next(:,j))'*K*(xi_next(:,j)-xi_ref);
    end

    %step momentum and velocity
    eta_next = zeros(6,N);
    eta_next(:,1) = eta0;
    %options = optimset('TolFun',1e-10,'Display','off');
    %options = optimoptions('fsolve','Algorithm','trust-region','Display','off','SpecifyObjectiveGradient',true,'TolFun',1e-14);
    W_bar_curr = [0;0;0;0;0;0];
    W_bar_back = [0;0;0;0;0;0];
    W_bar_prev = [0;0;0;0;0;0];
    W_bar_back = [0;0;0;reshape(g_next(1:9,1),3,3)'*grav]*rho*A;
    W_bar_back = W_bar_back + V*(xi_next(:,1)-xi(:,1))/dt;
    mu_next = zeros(6,N);
    for j=2:N     
        
        xi_der = (xi_next(:,j-1)-xi(:,j-1))/dt;
        %xi_der = (eta(:,j)-eta(:,j-1))/ds+ad(xi(:,j-1))*eta(:,j-1);
        W_bar = -V*xi_der;
        R = reshape(g(1:9,j-1),3,3);
        W_bar = W_bar + [0;0;0;R'*grav]*rho*A;
        W_bar = W_bar;
        %cable force
        nu = xi(4:6,j-1);
        omega = xi(1:3,j-1);
        nu_der = (xi(4:6,j)-nu)/ds;
        omega_der = (xi(1:3,j)-omega)/ds;
        for k=1:N_act
            pa_der = R*(nu-skew(r(k))*omega);
            pa_dder = R*(skew(omega)*(nu+skew(omega)*r(k))+nu_der+skew(omega_der)*r(k));
            ta_prime = -skew(pa_der)^2*pa_dder/norm(pa_der)^3;
            W_cable = -F(k)*[skew(r(k))*R'*ta_prime;R'*ta_prime];
            W_bar = W_bar+W_cable;
        end
        %W_bar = W_bar + V*xi_dot(:,j); %viscosity
        mu_next(:,j) = Ad(tau(dt*eta(:,j)))'*mu(:,j)+dt/ds*(lambda_next(:,j)-Ad(tau(ds*xi_next(:,j-1)))'*lambda_next(:,j-1))+dt*W_bar;
        %mu_next = dt/dt_prev*Ad(tau(dt_prev*eta(:,j)))'*mu(:,j)+dt/ds*(lambda(:,j)-Ad(tau(ds*xi(:,j-1)))'*lambda(:,j-1))+dt*W_bar;
        %mu_next = Ad(tau(dt*eta(:,j)))'*mu(:,j)+dt/ds*(lambda(:,j)-Ad(tau(ds*xi(:,j-1)))'*lambda(:,j-1))+dt*(W_bar_curr);%+W_bar_back+W_bar_prev);
        %mu_next = Ad(tau(dt*eta(:,j)))'*mu(:,j)+dt/ds*(lambda(:,j)-Ad(tau(ds*xi(:,j-1)))'*lambda(:,j-1));%+W_bar_back+W_bar_prev);
        %eta_next(:,j) = fsolve(@(eta) mu_next-dTauInv(dt*eta)'*M*eta,eta(:,j),options);
        eta_next(:,j) = newtonSolve(@(x) eta_equation(dt,mu_next(:,j),M,x), eta(:,j), 1e-14);
        %eta_next(:,j) = newtonSolve(@(x) eta_equation(dt,mu_next+dt*(V*((x-eta_next(:,j-1))/ds+ad(xi_next(:,j))*x)),M,x), eta(:,j), 1e-14);
        %eta_next(:,j) = fsolve(@(x) eta_equation(dt,mu_next,M,x),eta(:,j),options);
        
        W_bar_back = W_bar_curr;
    end
    %approximate energy
    T = 0;
    U = 0;
    for j=1:N
        T = T+1/2*eta_next(:,j)'*M*eta_next(:,j);
        U = U+1/2*xi_next(:,j)'*K*xi_next(:,j);
    end
    Energy = T+U;
    
%     %the last one is slightly shifted
%     W_bar = [0;0;0;0;0;0]; %temp for now
%     mu_next = Ad(tau(dt*eta(:,j)))'*mu(:,j)+2*dt/ds*(lambda(:,j)-Ad(tau(ds*xi(:,j-1)))'*lambda(:,j-1))+dt*W_bar;
%     eta_next(:,j) = fsolve(@(eta) mu_next-dTauInv(dt*eta)'*M*eta,eta(:,j));

    %possibly compute wave speed
%     xi_dot = zeros(6,N);
%     for i=1:N
%         %xi_dot(:,i) = (eta_next(:,i+1)-eta_next(:,i))/ds + ad(xi_next(:,i))*eta_next(:,i);
%         xi_dot(:,i) = (xi_next(:,i)-xi(:,i))/dt;
%     end
    
    %approx wave speeds
%     wave_speed = 0;
%     for i=2:N
%         xi_prime = (xi_next(:,i)-xi_next(:,i-1))/ds;
%         eta_dot = (eta_next(:,i)-eta(:,i))/dt;
%         wave_speeds = sqrt(abs(eta_dot./xi_prime));
%         wave_speed = max([wave_speed;wave_speeds]);
%     end
%     wave_speed;
    %rough_speed = sqrt(abs(((eta_next(2,2)-eta(2,2))/dt)/((xi_next(2,2)-xi_next(2,1))/ds)))
%     rough_speed = sqrt(abs(((eta_next(2,end)-eta(2,end))/dt)/((xi_next(2,end)-xi_next(2,end-1))/ds)));
%     dt_next = ds/(15*rough_speed);
%     dt_low = ds/(10*sqrt(E/rho));
%     dt_high = 1/30;
%     if dt_next < dt_low
%         dt_next = dt_low;
%     elseif dt_next > dt_high
%         dt_next = dt_high;
%     elseif isnan(dt_next)
%         dt_next = dt_low;
%     end
%     dt_next
end

function [F,J] = eta_equation(dt,mu,M,eta)
    %computes the values of mu-dTauInv(dt*eta)'*M*eta and the jacobian for
    %potentially faster convergence
    F = mu-(eye(6)-ad(dt*eta)/2)'*M*eta;
    if nargout > 1
        M1 = M(1,1);
        M2 = M(2,2);
        M3 = M(3,3);
        M4 = M(4,4);
        M5 = M(5,5);
        M6 = M(6,6);
        n1 = eta(1);
        n2 = eta(2);
        n3 = eta(3);
        n4 = eta(4);
        n5 = eta(5);
        n6 = eta(6);
        J = [-M1, (M2*dt*n3)/2 - (M3*dt*n3)/2, (M2*dt*n2)/2 - (M3*dt*n2)/2, 0, (M5*dt*n6)/2 - (M6*dt*n6)/2, (M5*dt*n5)/2 - (M6*dt*n5)/2; (M3*dt*n3)/2 - (M1*dt*n3)/2, -M2, (M3*dt*n1)/2 - (M1*dt*n1)/2, (M6*dt*n6)/2 - (M4*dt*n6)/2, 0, (M6*dt*n4)/2 - (M4*dt*n4)/2; (M1*dt*n2)/2 - (M2*dt*n2)/2, (M1*dt*n1)/2 - (M2*dt*n1)/2, -M3, (M4*dt*n5)/2 - (M5*dt*n5)/2, (M4*dt*n4)/2 - (M5*dt*n4)/2, 0; 0, -(M6*dt*n6)/2, (M5*dt*n5)/2, -M4, (M5*dt*n3)/2, -(M6*dt*n2)/2; (M6*dt*n6)/2, 0, -(M4*dt*n4)/2, -(M4*dt*n3)/2, -M5, (M6*dt*n1)/2; -(M5*dt*n5)/2, (M4*dt*n4)/2, 0, (M4*dt*n2)/2, -(M5*dt*n1)/2, -M6];

    end
end