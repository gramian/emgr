function estDemo(t)
%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% version: 5.9 (2021-01-21)
%%% authors: Christian Himpe (0000-0003-2194-6754)
%%% license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%%% summary: estDemo - run emgr examples via est

    switch lower(t)

        case 'hnm', hnm(); % Hyperbolic Network Model

        case 'isp', isp(); % Inverse Sylvester Procedure

        case 'fss', fss(); % Flexible Space Structures

        case 'nrc', nrc(); % Nonlinear Resistor-Capacitor Cascade

        case 'rqo', rqo(); % Random Diagonal System with Quadratic Output

        case 'lte', lte(); % Linear Transport Equation

        case 'aps', aps(); % All Pass System

        case 'fbc', fbc(); % Five-Body Choreography

        case 'qso', qso(); % Quasi-Stable Orbits Inside Black Holes

        otherwise,  error('Unknown example code!');
    end%switch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HYPERBOLIC NETWORK MODEL

function hnm()

    name = 'Hyperbolic Network Model';

    disp(['Example: ',name]);

    sys.M = 1;									% Number of inputs
    sys.N = 64;								% Number of states
    sys.Q = 1;									% Number of outputs

    A = sqrt(sys.N) * gallery('tridiag',sys.N,1,-2,1);				% System matrix
    B = tanh(1:sys.N)';							% Input matrix
    C = B'; 									% Output matrix

    sys.f = @(x,u,p,t) A * tanh(p.*x) + B * u;					% Vector field
    sys.g = @(x,u,p,t) C * x;							% Output functional
    sys.dt = 0.01;								% Time step
    sys.Tf = 1.0;								% Time horizon
    sys.p = ones(sys.N,1) * [0.5,1.0];						% Training parameter range

    task.type = 'combined_reduction';
    task.method = 'minimality';
    task.variant = 'dominant_subspaces';

    config.test = true;
    config.num_test_param = 10;
    config.pcentering = 'linear';
    config.extra_input = 'yes';
    config.skip_x = 3;
    config.skip_p = 3;

    [R,S] = est(sys,task,config);

    MORscore = S

    figure('Name',name,'NumberTitle','off');
    h = surf(R{1},R{2},R{5});
    xlabel('State Dimension');
    ylabel('Parameter Dimension');
    zlabel('Relative Error');
    xlim([R{1}(1),R{1}(end)]);
    ylim([R{2}(1),R{2}(end)]);
    zl = zlim();
    zlim([zl(1),1]);
    set(gca,'ZScale','log');
    set(h,'CData',log10(get(h,'CData')));
    set(gca,'CLim',log10(get(gca,'ZLim')));
    view(135,15);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INVERSE SYLVESTER PROCEDURE

function isp()

    name = 'Inverse Sylvester Procedure';

    disp(['Example: ',name]);

    rand('seed',1009);
    randn('seed',1009);

    sys.M = 1;									% Number of inputs
    sys.N = 64;								% Number of states
    sys.Q = 1;									% Number of outputs

    a = 1e-1;									% Minimum cross gramian singular value
    b = 1e+1;									% Maximum cross gramian singular value
    WX = -diag( a*((b/a).^rand(sys.N,1)) );					% Balanced cross gramian
    B = randn(sys.N,sys.M);							% Balanced input and output matrix
    A = sylvester(WX,WX,B*B') - sqrt(eps)*speye(sys.N);			% Balanced system matrix
    Q = orth(randn(sys.N,sys.N));						% Unbalancing transformation
    A = Q'*A*Q;								% Unbalanced system matrix
    B = Q'*B;									% Unbalanced input matrix
    C = B';									% Unbalanced output matrix

    sys.f = @(x,u,p,t) A*x + B*u;						% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional
    sys.F = @(x,u,p,t) A'*x + C'*u;						% Adjoint vector field
    sys.dt = 0.01;								% Time step
    sys.Tf = 1.0;								% Time horizon

    task.type = 'model_reduction';
    task.method = 'dominant_subspaces';
    task.variant = 'minimality';

    config.test = true;
    config.linearity = 'linear';

    [R,S] = est(sys,task,config);

    MORscore = S

    errorplot(name,R{1},R{3},R{4},R{5},R{6});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLEXIBLE SPACE STRUCTURES

function fss()

    name = 'Flexible Space Structure';

    disp(['Example: ',name]);

    K = 32;									% Number of subsystems

    sys.M = 1;									% Number of inputs
    sys.N = 2*K;								% Number of states
    sys.Q = 1;									% Number of outputs

    xi = rand(1,K) * 0.001;							% Sample damping ratio
    omega = rand(1,K) * 10.0;							% Sample natural frequencies
    A_k = cellfun(@(p) sparse([-2.0*p(1)*p(2),-p(2);p(2),0]), ...		% Subsystem blocks
                       num2cell([xi;omega],1),'UniformOutput',false);
    A = blkdiag(A_k{:});							% System matrix
    B = kron(rand(K,sys.M),[1;0]);						% Input matrix
    C = 10.0 * rand(sys.Q,2*K);						% Output matrix

    sys.f = @(x,u,p,t) A*x + B*u;						% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional
    sys.F = @(x,u,p,t) A'*x + C'*u;						% Adjoint vector field
    sys.dt = 0.01;								% Time step
    sys.Tf = 1.0;								% Time horizon

    task.type = 'model_reduction';
    task.method = 'dominant_subspaces';
    task.variant = 'minimality';

    config.test = true;
    config.linearity = 'linear';

    [R,S] = est(sys,task,config);

    MORscore = S

    errorplot(name,R{1},R{3},R{4},R{5},R{6});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NONLINEAR RESISTOR CAPACITOR CASCADE

function nrc()

    name = 'Nonlinear Resistor-Capacitor Cascade';

    disp(['Example: ',name]);

    sys.M = 1;									% Number of inputs
    sys.N = 64;								% Number of states
    sys.Q = 1;									% Number of outputs

    g = @(x) exp(x) + x - 1.0;							% Diode nonlinearity
    A0 = sparse(1,1,1,sys.N,sys.N);
    A1 = spdiags(ones(sys.N-1,1),-1,sys.N,sys.N) - speye(sys.N);
    A1(1,1) = 0;
    A2 = spdiags([ones(sys.N-1,1);0],0,sys.N,sys.N) - spdiags(ones(sys.N,1),1,sys.N,sys.N);
    B = sparse(1,1,1,sys.N,1);

    sys.f = @(x,u,p,t) -g(A0*x) + g(A1*x) - g(A2*x) + B*u;			% Vector field
    sys.g = @(x,u,p,t) x(1);							% Output functional
    sys.dt = 0.01;								% Time step
    sys.Tf = 1.0;								% Time horizon

    %sys.v = @(t) ones(sys.M,1)*(t<=0.5*sys.Tf);				% Test input

    task.type = 'model_reduction';
    task.method = 'dominant_subspaces';
    task.variant = 'minimality';

    config.test = true;

    [R,S] = est(sys,task,config);

    MORscore = S

    errorplot(name,R{1},R{3},R{4},R{5},R{6});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RANDOM DIAGONAL SYSTEM WITH QUADRATIC OUTPUT

function rqo()

    name = 'Random Diagonal System with Quadratic Output';

    disp(['Example: ',name]);

    rand('seed',1009);

    sys.M = 1;
    sys.N = 64;
    sys.Q = 1;

    A = spdiags(-rand(sys.N,1),0,sys.N,sys.N);
    B = rand(sys.N,1);

    sys.f = @(x,u,p,t) A*x + B*u;
    sys.g = @(x,u,p,t) norm(x);
    sys.dt = 0.01;
    sys.Tf = 1.0;

    task.type = 'model_reduction';
    task.method = 'dominant_subspaces';
    task.variant = 'minimality';

    config.test = true;
    config.kernel = 'gauss';

    [R,S] = est(sys,task,config);

    MORscore = S

    errorplot(name,R{1},R{3},R{4},R{5},R{6});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LINEAR TRANSPORT EQUATION

function lte()

    name = 'Linear Transport Equation';

    disp(['Example: ',name]);

    sys.M = 1;									% Number of inputs
    sys.N = 256;								% Number of states
    sys.Q = 1;									% Number of outputs

    A = spdiags(sys.N*ones(sys.N,1)*[1,-1],[-1,0],sys.N,sys.N);			% System matrix
    B = sparse(1,1,sys.N,sys.N,1);						% Input matrix
    C = sparse(1,sys.N,1.0,1,sys.N);						% Output matrix

    sys.f = @(x,u,p,t) p*A*x + B*u;						% Vector field
    sys.F = @(x,u,p,t) p*A'*x + C'*u;						% Adjoint vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional
    sys.dt = 0.5./sys.N;							% Time step
    sys.Tf = 1.5;								% Time horizon
    sys.p = 1;									% Transport velocity

    task.type = 'model_reduction';
    task.method = 'dominant_subspaces';
    task.variant = 'minimality';

    config.test = true;
    config.linearity = 'linear';
    config.skip_x = 4;

    [R,S] = est(sys,task,config);

    MORscore = S

    errorplot(name,R{1},R{3},R{4},R{5},R{6});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALL PASS SYSTEM

function aps()

    name = 'All Pass System';

    disp(['Example: ',name]);

    sys.M = 1;									% Number of inputs
    sys.N = 64;								% Number of states
    sys.Q = 1;									% Number of outputs

    A = gallery('tridiag',sys.N,-1,0,1); A(1,1) = -0.5;			% System matrix
    B = sparse(1,1,1,sys.N,1);							% Input matrix
    C = -B';									% Output matrix

    sys.f = @(x,u,p,t) A*x + B*u;						% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional
    sys.dt = 0.01;								% Time step
    sys.Tf = 1.0;								% Time horizon

    task.type = 'model_reduction';
    task.method = 'dmd_galerkin';
    task.variant = 'controllability';

    config.test = true;

    [R,S] = est(sys,task,config);

    MORscore = S

    errorplot(name,R{1},R{3},R{4},R{5},R{6});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIVE BODY CHOREOGRAPHY

function fbc()

    name = 'Five-Body Choreography';

    disp(['Example: ',name]);

    n = 5;									% Number of bodies

    sys.M = 0;									% Number of inputs
    sys.N = 4*n;								% Number of states
    sys.Q = 2*n;								% Number of outputs
    sys.f = @(x,u,p,t) [x((2*n)+1:end);acc(x(1:2*n),u,p)];			% Vector field
    sys.g = @(x,u,p,t)  x(1:2*n);						% Output functional
    sys.dt = 0.01;								% Time step
    sys.Tf = 1.0;								% Time horizon
    sys.p = ones(n,1);								% Parameters
    sys.x0 = [1.449;  0.0;    0.400; -0.345; -1.125;...			% Initial condition
              0.448; -1.125; -0.448;  0.400;  0.345;...
              0.0;   -0.922; -1.335;  0.810; -0.919;...
             -0.349;  0.919; -0.349;  1.335;  0.810];

    task.type = 'model_reduction';
    task.method = 'poor_man';
    task.variant = 'observability';

    config.centering = 'mean';
    config.kernel = @(x,y) x(size(x,1)/2+2:end,:) * pinv(y(:,size(y,2)/2+1:end-1)'); % velocity dmd kernel
    config.test = true;

    [R,S] = est(sys,task,config);

    MORscore = S

    errorplot(name,R{1},R{3},R{4},R{5},R{6});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUASI STABLE ORBITS

function qso()

    name = 'Quasi-Stable Orbits Inside Black Hole';

    disp(['Example: ',name]);

    sys.M = 1;									% Number of inputs
    sys.N = 4;									% Number of states
    sys.Q = 3;									% Number of outputs
    sys.f = @orbit;								% Vector field
    sys.g = @bl2c;								% Output functional
    sys.dt = 0.005;								% Time step
    sys.Tf = 5.0;								% Time horizon

    % Fermion
    sys.xs = [0.4;0.5*pi;0;0];							% Initial state
    sys.p = [[0.568;1.13;0.13;0.9982;0.05;1.0] * [0.9,1.1]; [-0.01,0.01]];	% Parameter range

    task.type = 'parameter_sensitivity';
    task.method = 'minimality';

    config = struct();

    RF = est(sys,task,config);

    % Photon
    EE = 10.5;
    sys.xs = [0.2;0.5*pi;0;0];							% Initial state
    sys.p = [[EE;1.38*EE;0.03*EE*EE;0.9982;0.05] * [0.9,1.1]; [-0.01,0.01;-0.01,0.01]];	% Parameter range

    RP = est(sys,task,config);

    figure('Name',name,'NumberTitle','off');
    set(gca,'XLim',[0,8],'XTickLabel',{' ','E','L','Q','a','e','\mu','\epsilon',' '},'YScale','log','YGrid','on','NextPlot','add');
    bar([RF,RP]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS

function errorplot(name,dims,l0,l1,l2,l8)

    figure('Name',name,'NumberTitle','off');
    set(gca,'XLim',[dims(1),dims(end)],'YLim',[1e-16,1.1],'YScale','log','Ytick',[1e-16,1e-8,1e-0],'YGrid','on','NextPlot','add');
    plot(dims,l1+eps,'r','LineWidth',3);
    plot(dims,l2+eps,'g','LineWidth',3);
    plot(dims,l8+eps,'b','LineWidth',3);
    plot(dims,l0+eps,'k--','LineWidth',3);
    legend('L_1 Error','L_2 Error','L_\infty Error','L_0 Error','location','NorthEast');
    xlabel('State Dimension');
    ylabel('Relative Error');
    pbaspect([2,1,1]);
end

function a = acc(x,u,p,t) % N-body acceleration vector field component

    N = numel(x)/2;
    A = reshape(x,[2,N]);
    y = zeros(2,N);

    for n = 1:N
        B = bsxfun(@minus,A,A(:,n));
        Z = p'./(sqrt(1e-6+(sum(B.^2))).^3);
        B = bsxfun(@times,B,Z);
        y(:,n) = sum(B,2);
    end%for

    a = y(:);
end

function x = orbit(x,u,p,t) % Generalized orbit vector-field

    E = p(1);  % E
    L = p(2);  % L
    Q = p(3);  % Q
    a = p(4);  % a
    e = p(5);  % e
    mu = p(6); % mu
    ep = p(7); % epsilon

    D  = x(1)^2 - 2.0*x(1) + a^2 + e^2;
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

