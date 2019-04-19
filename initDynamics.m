function [a,b,g] = initDynamics(N)
    %just create the initial matrices for testing
    L = 10e-2;
    ds = L/(N-1);
    a = zeros(6,N);
    %a(:,1) = [0;0;0;0;0.1;0]; %a little velocity at the base to start
    %b = repmat([0;0;0;0;0;ds*0.99],1,N-1);
    b = repmat([0;1/N;0;0;0;ds],1,N);
    g = zeros(12,N);
    G = eye(4);
    g(:,1) = [reshape(G(1:3,1:3),9,1);G(1:3,4)];
    for i=2:N
        G = G*expm(se(b(:,i-1)));
        g(:,i) = [reshape(G(1:3,1:3),9,1);G(1:3,4)];
    end
end