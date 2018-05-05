function examples(t)
%%% summary: emgrtest (run emgr demos)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2018)
%$
    rand('seed',1009);
    randn('seed',1009);

    switch(lower(t))

        case 'hnm' % Hyperbolic Network Model

            sys.M = 1;								% Number of inputs
            sys.N = 64;								% Number of states
            sys.Q = 1;								% Number of outputs
            sys.h = 0.01;							% Time step
            sys.T = 1.0;							% Time horizon

            A = sqrt(sys.N) * gallery('tridiag',sys.N,1,-2,1);			% System matrix
            B = ones(sys.N,sys.M);						% Input matrix
            C = B'; 								% Output matrix

            sys.f = @(x,u,p,t) A*tanh(p.*x) + B*u;				% Vector field
            sys.g = @(x,u,p,t) C*x;						% Output functional
            sys.p = ones(sys.N,1) * [0.5,1.0];					% Training parameter range
            sys.q = 0.5 + 0.5 * rand(sys.N,5);					% Test parameter

            curios(sys,'combined-reduction','minimality-based',{'active','linpar'});

        case 'isp' % Inverse Sylvester Procedure

            sys.M = 1;								% Number of inputs
            sys.N = 64;								% Number of states
            sys.Q = 1;								% Number of outputs
            sys.h = 0.01;							% Time step
            sys.T = 1.0;							% Time horizon

            a = 1e-1;
            b = 1e+1;
            WX = -diag( a*((b/a).^rand(sys.N,1)) );				% Balanced cross gramian
            B = randn(sys.N,sys.M);						% Balanced input and output matrix
            A = sylvester(WX,WX,B*B') - sqrt(eps)*speye(sys.N);			% Balanced system matrix
            Q = orth(randn(sys.N,sys.N));					% Unbalancing transformation
            A = Q'*A*Q;								% Unbalanced system matrix
            B = Q'*B;								% Unbalanced input matrix
            C = B';								% Unbalanced output matrix

            sys.f = @(x,u,p,t) A*x + B*u;					% Vector field
            sys.g = @(x,u,p,t) C*x;						% Output functional

            curios(sys,'state-reduction','nonlinear-balanced-truncation');

        case 'fss' % Flexible Space Structures

            K = 32;								% Number of subsystems

            sys.M = 1;								% Number of inputs
            sys.N = 2*K;							% Number of states
            sys.Q = 1;								% Number of outputs
            sys.h = 0.01;							% Time step
            sys.T = 1.0;							% Time horizon

            xi = rand(1,K) * 0.001;						% Sample damping ratio
            omega = rand(1,K) * 10.0;						% Sample natural frequencies
            A_k = cellfun(@(p) sparse([-2.0*p(1)*p(2),-p(2);p(2),0]), ...
                                num2cell([xi;omega],1),'UniformOutput',0);	% Subsystem blocks
            A = blkdiag(A_k{:});						% System matrix
            B = kron(rand(K,sys.M),[1;0]);					% Input matrix
            C = 10.0 * rand(sys.Q,2*K);						% Output matrix

            sys.f = @(x,u,p,t) A*x + B*u;					% Vector field
            sys.g = @(x,u,p,t) C*x;						% Output functional

            curios(sys,'state-reduction','nonlinear-direct-truncation');

        case 'nrc' % Nonlinear Resistor-Capacitor cascade

            sys.M = 1;								% Number of inputs
            sys.N = 64;								% Number of states
            sys.Q = 1;								% Number of outputs
            sys.h = 0.01;							% Time step
            sys.T = 1.0;							% Time horizon

            g = @(x) exp(x) + x - 1.0;						% Diode nonlinearity
            A0 = sparse(1,1,1,sys.N,sys.N);
            A1 = spdiags(ones(sys.N-1,1),-1,sys.N,sys.N) - speye(sys.N);
            A1(1,1) = 0;
            A2 = spdiags([ones(sys.N-1,1);0],0,sys.N,sys.N) - spdiags(ones(sys.N,1),1,sys.N,sys.N);
            B = sparse(1,1,1,sys.N,1);

            sys.f = @(x,u,p,t) -g(A0*x) + g(A1*x) - g(A2*x) + B*u;		% Vector field
            sys.g = @(x,u,p,t) x(1);						% Output functional
            sys.v = @(t) ones(sys.M,1)*(t<=0.5*sys.T);				% Test input

            curios(sys,'state-reduction','nonlinear-direct-truncation');

        case 'lte' % Linear Transport Equation

            sys.M = 1;								% Number of inputs
            sys.N = 256;							% Number of states
            sys.Q = 1;								% Number of outputs
            sys.h = 1.0./sys.N;							% Time step
            sys.T = 1.0;							% Time horizon

            A = sys.N * spdiags(ones(sys.N,1)*[1,-1],[-1,0],sys.N,sys.N);	% System matrix
            B = sparse(1,1,sys.N,sys.N,1);					% Input matrix
            C = sparse(1,sys.N,1.0,1,sys.N);					% Output matrix

            sys.f = @(x,u,p,t) p*A*x + B*u;					% Vector field
            sys.F = @(x,u,p,t) p*A'*x + C'*u;					% Adjoint vector field
            sys.g = @(x,u,p,t) C*x;						% Output functional
            sys.ut = Inf;
            sys.v = @(t) exp(((t-0.1).^2)./(-0.001));				% Input function
            sys.p = 1.4;							% Transport velocity
            sys.q = sys.p;

            curios(sys,'state-reduction','linear-direct-truncation');

        case 'fbc' % Five-Body Choreography

            n = 5;								% Number of bodies

            sys.M = 0;								% Number of inputs
            sys.N = 4*n;							% Number of states
            sys.Q = 2*n;							% Number of outputs
            sys.h = 0.01;							% Time step
            sys.T = 1.0;							% Time horizon

            sys.f = @(x,u,p,t) [x((2*n)+1:end);acc(x(1:2*n),u,p)];		% Vector field
            sys.g = @(x,u,p,t)  x(1:2*n);					% output functional
            sys.xs = [1.449;  0.0;    0.400; -0.345; -1.125;...			% Initial condition
                      0.448; -1.125; -0.448;  0.400;  0.345;...
                      0.0;   -0.922; -1.335;  0.810; -0.919;...
                     -0.349;  0.919; -0.349;  1.335;  0.810];
            sys.p = ones(n,1);							% Parameters

            curios(sys,'state-reduction','observability-truncation',{'position'});

        case 'qso' % Quasi-Stable Orbits inside black holes

            sys.M = 0;								% Number of inputs
            sys.N = 4;								% Number of states
            sys.Q = 3;								% Number of outputs
            sys.h = 0.005;							% Time step
            sys.T = 5.0;							% Time horizon

            sys.f = @(x,u,p,t) orbit(x,u,p,t);					% Vector field
            sys.g = @(x,u,p,t) bl2c(x,u,p,t);					% Output functional

            fprintf('\nParameters: E, L, Q, a, e, epsilon, mu\n\n');

            % Fermion
            sys.xs = [0.4;pi/2;0;0];						% Initial state
            sys.p = [0.568;1.13;0.13;0.9982;0.05;0;1] * [0.9,1.1];		% Parameter range

            curios(sys,'sensitivity-analysis','input-output-based');

            % Photon
            EE = 10.5;
            sys.xs = [0.2;pi/2;0;0];						% Initial state
            sys.p = [EE;1.38*EE;0.03*EE*EE;0.9982;0.05;0;0]*[0.9,1.1];		% Parameter range

            curios(sys,'sensitivity-analysis','input-output-based',{'hold'});
    end
end

function y = acc(x,u,p,t) % Acceleration vector field component

    N = numel(x)/2;
    A = reshape(x,[2,N]);
    y = zeros(2,N);

    for n = 1:N
        B = bsxfun(@minus,A,A(:,n));
        Z = p'./(sqrt(1e-6+(sum(B.^2))).^3);
        B = bsxfun(@times,B,Z);
        y(:,n) = sum(B,2);
    end

    y = y(:);
end

function x = orbit(x,u,p,t) % Generalized orbit vector-field

    E = p(1); % E
    L = p(2); % L
    Q = p(3); % Q
    a = p(4); % a
    e = p(5); % e
    ep = p(6); % epsilon
    mu = p(7); % mu

    D  = x(1)^2 - 2*x(1) + a^2 + e^2;
    S  = x(1)^2 + a^2*cos(x(2))^2;
    P  = E*(x(1)^2 + a^2) + e*ep*x(1) - a*L;
    Vt = Q - cos(x(2))^2*(a^2*(mu^2 - E^2) + L^2*sin(x(2))^(-2) );
    Vr = P^2 - D*(mu^2*x(1)^2 + (L - a*E)^2 + Q);

    x = abs([ sqrt(Vr) ; ...
              sqrt(Vt) ; ...
              L * sin(x(2))^(-2) + a*(P/D-E) ; ...
              a * (L - a * E * sin(x(2))^2) + P/D * (x(1)^2 + a^2) ]./S);
end

function y = bl2c(x,u,p,t) % Boyer-Lindquist to Cartesian coordinate conversion

    y = [ sqrt(x(1)^2 + p(4)^2) * sin(x(2)) * cos(x(3)); ...
          sqrt(x(1)^2 + p(4)^2) * sin(x(2)) * sin(x(3)); ...
          x(1) * cos(x(2)) ];
end
