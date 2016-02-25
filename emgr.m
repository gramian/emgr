function W = emgr(f,g,s,t,w,pr,nf,ut,us,xs,um,xm)
% emgr - Empirical Gramian Framework ( Version: 3.9 )
% Copyright (c) 2013-2016 Christian Himpe ( gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% SYNTAX:
%    W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm]);
%
% SUMMARY:
%    emgr - EMpirical GRamian framemwork,
%    computation of empirical gramians for model reduction,
%    system identification and uncertainty quantification.
%    Enables gramian-based nonlinear model order reduction.
%    Compatible with OCTAVE and MATLAB.
%
% ARGUMENTS:
%   (func handle)  f - system function handle; signature: xdot = f(x,u,p)
%   (func handle)  g - output function handle; signature:    y = g(x,u,p)
%        (vector)  s - system dimensions [inputs,states,outputs]
%        (vector)  t - time discretization [step,stop]
%          (char)  w - gramian type:
%            * 'c' : empirical controllability gramian (WC)
%            * 'o' : empirical observability gramian (WO)
%            * 'x' : empirical cross gramian (WX)
%            * 'y' : empirical linear cross gramian (WY)
%            * 's' : empirical sensitivity gramian (WS)
%            * 'i' : empirical identifiability gramian (WI)
%            * 'j' : empirical joint gramian (WJ)
% (matrix,vector,scalar) [pr = 0] - parameters, each column is one set
%        (vector,scalar) [nf = 0] - options, 12 components:
%            + zero(0),init(1),steady(2),mean(3),median(4),midr(5),rms(6) center
%            + linear(0), log(1), geom(2), single(3), sparse(4) input scales
%            + linear(0), log(1), geom(2), single(3), sparse(4) state scales
%            + unit(0), reciproce(1), dyadic(2), single(3) input rotations
%            + unit(0), reciproce(1), dyadic(2), single(3) state rotations
%            + single(0), double(1), scaled(2) run
%            + regular(0), non-symmetric(1) cross gramian; only: WX, WJ
%            + plain(0), robust(1) parameters; only: WC, WY
%            + active(0), passive(1) parameter; only: WI, WJ
%            + none(0), linear(1), logarithmic(2) parameter centering
%            + default(0), exclusive options:
%                  * use rms-centering(1); only: WS
%                  * use Schur-complement(1); only: WI
%                  * use detailed Schur-complement(1); only: WJ
%            + assume(0), enforce(1) gramian symmetry
%  (matrix,vector,scalar) [ut = 1] - input; default: delta impulse
%         (vector,scalar) [us = 0] - steady-state input
%         (vector,scalar) [xs = 0] - steady-state and initial state x0
%  (matrix,vector,scalar) [um = 1] - input scales
%  (matrix,vector,scalar) [xm = 1] - initial-state scales
%
% RETURNS:
%            (matrix)  W - Gramian Matrix (only: WC, WO, WX, WY)
%              (cell)  W - {State-,Parameter-} Gramian (only: WS, WI, WJ)
%
% CITATION:
%    C. Himpe (2016). emgr - Empirical Gramian Framework (Version 3.9)
%    [Software]. Available from http://gramian.de . doi:10.5281/zenodo.46523 .
%
% SEE ALSO:
%    gram
%
% KEYWORDS:
%    model reduction, empirical gramian, cross gramian, mor
%
% Further information: <http://gramian.de>
%*

    % Custom ODE Solver
    global ODE;
    if(isa(ODE,'function_handle')==0), ODE = @rk2; end;

    % Version Info
    if( (nargin==1) && strcmp(f,'version') ), W = 3.9; return; end;

    % Default Arguments
    if( (nargin<6)  || isempty(pr) ), pr = 0.0; end;
    if( (nargin<7)  || isempty(nf) ), nf = 0;   end;
    if( (nargin<8)  || isempty(ut) ), ut = 1.0; end;
    if( (nargin<9)  || isempty(us) ), us = 0.0; end;
    if( (nargin<10) || isempty(xs) ), xs = 0.0; end;
    if( (nargin<11) || isempty(um) ), um = 1.0; end;
    if( (nargin<12) || isempty(xm) ), xm = 1.0; end;

    % System Dimensions
    J = s(1);               % number of inputs
    N = s(2);               % number of states
    O = s(3);               % number of outputs
    M = 0;                  % internal variable used by WS, WI, WJ
    if(numel(s)==4)
        M = s(4);
    end;

    h = t(1);               % width of time step
    T = floor(t(2)/h) + 1;  % number of time steps plus initial value

    w = lower(w);           % ensure lower case gramian type

    P = size(pr,1);         % number of parameters
    Q = size(pr,2);         % number of parameter sets

    % Linear Chirp Input
    if( isnumeric(ut) && numel(ut)==1 && ut==Inf )
        ut = @(t) 0.5*cos(pi*(t+10*t.*t))+0.5;
    end;

    % Discretize Procedural Input
    if(isa(ut,'function_handle'))
        uf = ut;
        ut = zeros(J,T);
        for l=1:T
            ut(:,l) = uf(l*h);
        end;
    end;

    % Lazy Arguments
    if( isnumeric(g) && g==1 ), g = @(x,u,p) x; O = N; end;

    if(numel(nf)<12), nf(12)    = 0;  end;
    if(numel(ut)==1), ut(1:J,1) = (1.0/h)*ut; end;
    if(numel(us)==1), us(1:J,1) = us; end;
    if(numel(xs)==1), xs(1:N,1) = xs; end;
    if(numel(um)==1), um(1:J,1) = um; end;
    if(numel(xm)==1), xm(1:N+(w=='y')*(J-N),1) = xm; end;

    if(size(ut,2)==1), ut(:,2:T) = 0.0; end;
    if(size(us,2)==1), us = repmat(us,[1,T]); end;
    if(size(um,2)==1), um = scales(um,nf(2),nf(4)); end;
    if(size(xm,2)==1), xm = scales(xm,nf(3),nf(5)); end;

%% PARAMETRIC SETUP

    if( (nf(8) && w~='o') || w=='s' || w=='i' || w=='j' )

        if(Q==1), error('ERROR! emgr: min and max parameter required!'); end;

        pmin = min(pr,[],2);
        pmax = max(pr,[],2);

        if( nf(8) || w=='s' ) % assemble (controllability) parameter scales
            pn = size(um,2);
        else                  % assemble (observability) parameter scales
            pn = size(xm,2);
        end;

        pl = (1.0:floor(pn/2))./floor(pn/2);
        pu = (1.0:ceil(pn/2))./ceil(pn/2);

        switch(nf(10)) % parameter centering

            case 1, % linear
                pr = mean(pr,2);
                pm = [(pmin - pr)*pl , (pmax - pr)*pu];

            case 2, % logarithmic
                lpmin = log(pmin);
                lpmax = log(pmax);
                lpavg = 0.5*(lpmax - lpmin);
                pr = pmin.*exp(lpavg);
                pm = [bsxfun(@times,exp((lpavg - lpmin)*pl),pmin), ...
                      bsxfun(@times,exp((lpmax - lpavg)*pu),pr) ];
                pm = bsxfun(@minus,pm,pr);

            otherwise, % none
                pr = pmin;
                pm = (pmax - pmin)*((1:pn)./pn);
        end;

        Q = 1;
    end;

%% STATE-SPACE SETUP

    if( w=='c' || w=='o' || w=='x' || w=='y' )

        C = size(um,2); % number of input scales
        D = size(xm,2); % number of state scales

        switch(nf(1)) % residual types

            case 1, % initial state
                res = @(d) d(:,1);

            case 2, % steady state
                res = @(d) d(:,end);

            case 3, % mean state
                res = @(d) mean(d,2);

            case 4, % median state
                res = @(d) median(d,2);

            case 5, % midrange
                res = @(d) 0.5*(min(d,2)+max(d,2));

            case 6, % rms
                res = @(d) sqrt(sum(d.*d,2));

            otherwise, % zero state
                res = @(d) zeros(size(d,1),1);
        end;

        switch(nf(6)) % scaled runs

            case 1, % preconditioned run
                nf(6) = 0;
                WT = emgr(f,g,s,t,w,pr,nf,ut,us,xs,um,xm);
                TX = sqrt(diag(WT));
                TX = TX(1:(N-(M>0 && w~='c')*P));
                tx = 1.0./TX;
                F = f; f = @(x,u,p) TX.*F(tx.*x,u,p);
                G = g; g = @(x,u,p)     G(tx.*x,u,p);

            case 2, % steady state (input) scaled run
                TU = us(:,1);
                TX = xs;
                TX = TX(1:(N-(M>0 && w~='c')*P));
                TU(TU==0) = 1.0; tu = 1.0./TU;
                TX(TX==0) = 1.0; tx = 1.0./TX;
                F = f; f = @(x,u,p) TX.*F(tx.*x,tu.*u,p);
                G = g; g = @(x,u,p)     G(tx.*x,tu.*u,p);
        end;

        if(nf(8)) % robust parameter
            J = J + P;
            ut = [ut;ones(P,T)];
            us = [us;repmat(pr,[1,T])];
            um = [um;pm];
            if(w=='y'), xm = [xm;pm]; end;
            F = f; f = @(x,u,p) F(x,u(1:J-P),u(J-P+1:J));
            G = g; g = @(x,u,p) G(x,u(1:J-P),u(J-P+1:J));
        end;

        m = N - P*(M>0 && w=='x'); % non-zero rows if joint gramian
        W = zeros(m,N);            % preallocate gramian
    end;

%% GRAMIAN COMPUTATION

    switch(w) % empirical gramian types

        case 'c', % controllability gramian
            for q=1:Q
                pp = pr(:,q);
                for c=1:C
                    for j=1:J % parfor
                        if(um(j,c)==0), continue; end;
                        if(M>0)
                            up = pr + sparse(M,1,um(j,c),P,1);
                            x = ODE(f,1,t,xs,us,up);
                        else
                            uu = us + bsxfun(@times,ut,um(j,c)*(1:J==j)');
                            x = ODE(f,1,t,xs,uu,pp);
                        end;
                        x = bsxfun(@minus,x,res(x));
                        x = x * (1.0./um(j,c));
                        W = W + (x*x'); % offload
                    end;
                end;
            end;
            W = W * (h/(C*Q));

        case 'o', % observability gramian
            o = zeros(O*T,N);
            for q=1:Q
                pp = pr(:,q);
                for d=1:D
                    for n=1:N % parfor
                        if(xm(n,d)==0), continue; end;
                        xx = xs + xm(n,d)*(1:N==n)';
                        if(M>0 && n>M && nf(9))
                            y = ODE(f,g,t,xx(1:M),us+ut,xx(M+1:end));
                        elseif(M>0)
                            y = ODE(f,g,t,xx(1:M),us,xx(M+1:end));
                        else
                            y = ODE(f,g,t,xx,us,pp);
                        end;
                        y = bsxfun(@minus,y,res(y));
                        y = y * (1.0/xm(n,d));
                        o(:,n) = y(:);
                    end;
                    W = W + (o'*o); % offload
                end;
            end;
            W = W * (h/(D*Q));

        case 'x', % cross gramian
            if(J~=O && nf(7)==0), error('ERROR! emgr: non-square system!'); end;
            o = zeros(O,T,N);
            for q=1:Q
                pp = pr(:,q);
                for d=1:D
                    for n=1:N % parfor
                        if(xm(n,d)==0), continue; end;
                        xx = xs + xm(n,d)*(1:N==n)';
                        if(M>0 && n>M && nf(9))
                            y = ODE(f,g,t,xx(1:M),us+ut,xx(M+1:end));
                        elseif(M>0)
                            y = ODE(f,g,t,xx(1:M),us,xx(M+1:end));
                        else
                            y = ODE(f,g,t,xx,us,pp);
                        end;
                        y = bsxfun(@minus,y,res(y));
                        y = y * (1.0/xm(n,d));
                        o(:,:,n) = y;
                    end;
                    o = permute(o,[2,3,1]); % generalized transposition
                    if(nf(7))
                        o(:,:,1) = sum(o,3);
                    end;
                    for c=1:C
                        for j=1:J % parfor
                            if(um(j,c)==0), continue; end;
                            uu = us + bsxfun(@times,ut,um(j,c)*(1:J==j)');
                            if(M>0)
                                x = ODE(f,1,t,xs(1:M),uu,xs(M+1:end));
                            else
                                x = ODE(f,1,t,xs,uu,pp);
                            end;
                            x = bsxfun(@minus,x,res(x));
                            x = x * (1.0./um(j,c));
                            if(nf(7)) % non-symmetric cross gramian
                                W = W + (x*o(:,:,1)); % offload
                            else      % regular cross gramian
                                W = W + (x*o(:,:,j)); % offload
                            end;
                        end;
                    end;
                    o = reshape(o,O,T,N); % reset
                end;
            end;
            W = W * (h/(C*D*Q));

        case 'y', % linear cross gramian
            if(J~=O && nf(8)==0), error('ERROR! emgr: non-square system!'); end;
            for q=1:Q
                pp = pr(:,q);
                for c=1:C
                    for j=1:J % parfor
                        if(um(j,c)==0 || xm(j,c)==0), continue; end;
                        uu = us + bsxfun(@times,ut,um(j,c)*(1:J==j)');
                        x = ODE(f,1,t,xs,uu,pp);
                        x = bsxfun(@minus,x,res(x));
                        x = x * (1.0./um(j,c));
                        uu = us + bsxfun(@times,ut,xm(j,c)*(1:J==j)');
                        z = ODE(g,1,t,xs,uu,pp);
                        z = bsxfun(@minus,z,res(z));
                        z = z * (1.0./xm(j,c));
                        W = W + (x*z'); % offload
                    end;
                end;
            end;
            W = W * (h/(C*Q));

        case 's', % sensitivity gramian
            W = cell(1,2);
            ps = sparse(P,1);
            nf(8) = 0;
            W{1} = emgr(f,g,[J,N,O],t,'c',ps,nf,ut,us,xs,um,xm);
            W{2} = zeros(P,1);
            for p=1:P
                V = emgr(f,g,[1,N,O,p],t,'c',pr,nf,ut,us,xs,pm(p,:),xm);
                W{1} = W{1} + V;        % approximate controllability gramian
                W{2}(p) = trace(V);
            end;
            if(nf(11))
                W{2} = W{2} - mean(W{2});
            end;
            W{2} = spdiags(W{2},0,P,P); % sensitivity gramian

        case 'i', % identifiability gramian
            W = cell(1,2);
            ps = sparse(P,1);
            V = emgr(f,g,[J,N+P,O,N],t,'o',ps,nf,ut,us,[xs;pr],um,[xm;pm]);
            W{1} = V(1:N,1:N);         % observability gramian
            W{2} = V(N+1:N+P,N+1:N+P); % identifiability gramian
            if(nf(11))
                W{2} = W{2} - V(N+1:N+P,1:N)*ainv(W{1})*V(1:N,N+1:N+P);
            end;

        case 'j', % joint gramian
            W = cell(1,2);
            ps = sparse(P,1);
            V = emgr(f,g,[J,N+P,O,N],t,'x',ps,nf,ut,us,[xs;pr],um,[xm;pm]);
            W{1} = V(1:N,1:N); % cross gramian
            %W{2} = zeros(P,P); % cross-identifiability gramian
            if(nf(11))
                W{2} = -0.5*V(1:N,N+1:N+P)'*pinv(W{1}+W{1}')*V(1:N,N+1:N+P);
            else
                W{2} = -0.5*V(1:N,N+1:N+P)'*ainv(W{1}+W{1}')*V(1:N,N+1:N+P);
            end;

        otherwise,
            error('ERROR! emgr: unknown gramian type!');
    end;

    if(nf(12) && (w=='c' || w=='o' || w=='y' || w=='x') ) % enforce symmetry
        W(1:m,1:m) = 0.5*(W(1:m,1:m) + W(1:m,1:m)');
    end;
end

%% ======== SCALES SELECTOR ========
function s = scales(s,d,e)

    switch(d)

        case 0, % linear
            s = s*[0.25,0.50,0.75,1.0];

        case 1, % logarithmic
            s = s*[0.001,0.01,0.1,1.0];

        case 2, % geometric
            s = s*[0.125,0.25,0.5,1.0];

        case 4, % sparse
            s = s*[0.38,0.71,0.92,1.0];

        otherwise, % single
            %s = s;
    end;

    switch(e)

        case 1, % reciproce
            s = [1.0./s,s];

        case 2, % dyadic
            s = s*s';

        case 3, % single
            %s = s;

        otherwise, % unit
            s = [-s,s];
    end;
end

%% ======== FAST APPROXIMATE INVERSION ========
function x = ainv(m)

    d = diag(m);
    d(d~=0) = 1.0./d(d~=0);
    n = numel(d);
    x = bsxfun(@times,m,-d);
    x = bsxfun(@times,x,d');
    x(1:n+1:end) = d;
end

%% ======== DEFAULT ODE INTEGRATOR ========
function x = rk2(f,g,t,z,u,p)

    if(isnumeric(g) && g==1), g = @(x,u,p) x; end;

    h = t(1);
    L = floor(t(2)/h) + 1;

    x(:,1) = g(z,u(:,end),p);
    x(end,L) = 0; % preallocate trajectory

    for l=2:L % 2nd order Ralston's Runge-Kutta Method
        k1 = h*f(z,u(:,l-1),p);
        k2 = h*f(z + 0.666666666666667*k1,u(:,l-1),p);
        z = z + 0.25*k1 + 0.75*k2;
        x(:,l) = g(z,u(:,l-1),p);
    end;
end

