function [g,xi,eta,mu,lambda] = initRod(N)
    %just create the initial matrices for testing
    D = 1e-2;
    A = pi/4*D^2;
    I = pi/64*D^4;
    J = 2*I;
    E = 100e3;
    G = E/3;
    rho = 1000;
    K = diag([E*I,E*I,G*J,G*A,G*A,E*A]);
    xi_ref = [0;0;0;0;0;1];
    L = 10e-2;
    ds = L/(N-1);
    eta = zeros(6,N);
    mu = zeros(6,N);
    lambda = zeros(6,N);
    %a(:,1) = [0;0;0;0;0.1;0]; %a little velocity at the base to start
    %b = repmat([0;0;0;0;0;ds*0.99],1,N-1);
    xi = repmat([0;0;0;0;0;1],1,N);
    %xi = repmat([0;0;0;0;0;1],1,N);
    g = zeros(12,N);
    G = eye(4);
    g(:,1) = [reshape(G(1:3,1:3),9,1);G(1:3,4)];
    for i=2:N
        G = G*expm(se(ds*xi(:,i-1)));
        g(:,i) = [reshape(G(1:3,1:3),9,1);G(1:3,4)];
    end
    
    for i=1:N
        lambda(:,i) = (eye(6)-ad(xi(:,i))/2)'*K*(xi(:,i)-xi_ref);
    end
    
end