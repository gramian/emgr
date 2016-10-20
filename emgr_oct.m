function W = emgr_oct(f,g,s,t,w,pr=0,nf=0,ut=1,us=0,xs=0,um=1,xm=1,pm=[])
### project: emgr - Empirical Gramian Framework ( http://gramian.de )
### version: 5.0-oct ( 2016-10-20 )
### authors: Christian Himpe ( 0000-0003-2194-6754 )
### license: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
#
## SYNTAX:
#   W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[pm]);
#
## SUMMARY:
#   empirical gramian matrix computation for model reduction,
#   decentralized control, sensitivity analysis, parameter identification,
#   uncertainty quantification and combined state and parameter reduction
#   of large-scale input-output systems. Enables data-driven analysis of
#   input-output coherence and gramian-based nonlinear model order reduction.
#   Compatible with OCTAVE and MATLAB.
#
## ARGUMENTS:
#   (handle) f : system function handle; signature: xdot = f(x,u,p,t)
#   (handle) g : output function handle; signature:    y = g(x,u,p,t)
#   (vector) s : system dimensions [inputs,states,outputs]
#   (vector) t : time discretization [step,stop]
#   (string) w : character encoding gramian type:
#      * 'c' empirical controllability gramian (WC)
#      * 'o' empirical observability gramian (WO)
#      * 'x' empirical cross gramian (WX aka WCO and XCG)
#      * 'y' empirical linear cross gramian (WY)
#      * 's' empirical sensitivity gramian (WS)
#      * 'i' empirical identifiability gramian (WI)
#      * 'j' empirical joint gramian (WJ)
#   (matrix) pr = 0 : parameters, each column is one set
#   (vector) nf = 0 : option flags, ten components:
#      * centering: no(0), init(1), steady(2), mean(3), rms(4), midrange(5)
#      * input scales: single(0), linear(1), geometric(2), log(3), sparse(4)
#      * state scales: single(0), linear(1), geometric(2), log(3), sparse(4)
#      * input rotations: unit(0), single(1)
#      * state rotations: unit(0), single(1)
#      * scale type: no(0), Jacobi preconditioner(1), steady-state(2)
#      * cross gramian type (only: WX,WY,WJ): regular(0), non-symmetric(1) (WZ)
#      * extra input (only: WO,WX,WS,WI,WJ): no(0), param(1), state(2), both(3)
#      * parameter centering (only: WS,WI,WJ): no(0), linear(1), logarithmic(2)
#      * Schur-complement (only: WI,WJ): detailed(0), approximate(1)
#   (handle) ut = 1 : input function handle; default: delta impulse
#   (vector) us = 0 : steady-state input
#   (vector) xs = 0 : steady-state and initial state x0
#   (matrix) um = 1 : input scales
#   (matrix) xm = 1 : initial-state scales
#   (matrix) pm = [] : parameter scales (reserved)
#
## RETURNS:
#   (matrix) W : Gramian Matrix (only: WC,WO,WX,WY)
#     (cell) W : [State-, Parameter-] Gramian (only: WS,WI,WJ)
#
## CITATION:
#   C. Himpe (2016). emgr - Empirical Gramian Framework (Version 5.0)
#   [Software]. Available from http://gramian.de . doi: 10.5281/zenodo.162135 .
#
## SEE ALSO:
#   gram
#
## KEYWORDS:
#   model reduction, empirical gramians, cross gramian matrix, MOR
#
# Further information: <http://gramian.de>
#$
    global DOT; # Inner Product Handle
    global ODE; # Integrator Handle
    global DWX; # Distributed Cross Gramian (Column Width, Index)

    if(isa(DOT,'function_handle')==0), DOT = @mtimes; end;
    if(isa(ODE,'function_handle')==0), ODE = @ssp2; end;

    # Version Info
    if(strcmp(f,'version'))
        W = 5.0;
        return;
    end;

## GENERAL SETUP

    # System Dimensions
    M = s(1);                   # Number of inputs
    N = s(2);                   # Number of states
    Q = s(3);                   # Number of outputs
    if(numel(s)==4)             # Number of augmented parameter-states / inputs
        A = s(4);
    else
        A = 0;
    end;
    P = size(pr,1);             # Dimension of parameter
    K = size(pr,2);             # Number of parameters
    h = t(1);           	# Width of time step
    L = floor(t(2)/h) + 1;      # Number of time steps plus initial value

    w = lower(w);               # Ensure lower case gramian type

    # Lazy Arguments
    if(isnumeric(g) && g==1)
        g = @(x,u,p,t) x;
        Q = N;
    end;

    if(numel(nf)<10)
        nf(10) = 0;
    end;

    if(isnumeric(ut) && numel(ut)==1)
        if(ut==Inf) # Linear Chirp Input
            mh = ones(M,1);
            sh = (1.0/L - 0.1/h) / L;
            ut = @(t) 0.5 + mh*0.5*cos(2.0*pi*(((0.1/h)+0.5*sh*t).*t));
        else # Delta Impulse Input
            mh = ones(M,1)./h;
            ut = @(t) mh*(t<=h);
        end;
    end;

    if(numel(us)==1), us = us*ones(M,1); end;
    if(numel(xs)==1), xs = xs*ones(N,1); end;
    if(numel(um)==1), um = um*ones(M,1); end;
    if(numel(xm)==1), xm = xm*ones(N-(w=='y')*(N-M),1); end;

    if(size(um,2)==1), um = scales(um,nf(2),nf(4)); end;
    if(size(xm,2)==1), xm = scales(xm,nf(3),nf(5)); end;

## GRAMIAN SETUP

    if( w=='c' || w=='o' || w=='x' || w=='y' )

        C = size(um,2); # Number of input scales
        D = size(xm,2); # Number of state scales

        if(isempty(pm))
            pm = sparse(1,max(C,D));
        end;

        switch(nf(6))

            case 1 # Preconditioned run
                nf(6) = 0;
                WT = emgr_oct(f,g,s,t,w,pr,nf,ut,us,xs,um,xm,pm);
                TX = sqrt(WT(1:N+1:end))';
                tx = 1.0./TX;
                F = f; f = @(x,u,p,t) TX.*F(tx.*x,u,p,t);
                G = g; g = @(x,u,p,t)     G(tx.*x,u,p,t);

            case 2 # Steady-state scaled run
                TU = us(:,1); TU(TU==0) = 1.0;
                tu = 1.0./TU;
                TX = xs; TX(TX==0) = 1.0;
                tx = 1.0./TX;
                F = f; f = @(x,u,p,t) TX.*F(tx.*x,tu.*u,p,t);
                G = g; g = @(x,u,p,t)     G(tx.*x,tu.*u,p,t);
        end;

        if(nf(8)==1 || nf(8)==3) # Extra input for parameter perturbations
            up = @(t) us + ut(t);
        else
            up = @(t) us;
        end

        if(nf(8)==2 || nf(8)==3) # Extra input for state perturbations
            ux = @(t) us + ut(t);
        else
            ux = @(t) us;
        end
    end;

## GRAMIAN COMPUTATION

    switch(w) # Empirical gramian types

        case 'c' # Controllability gramian
            W = zeros(N,N);
            o = zeros(A,1);
            for k = 1:K
                for c = 1:C
                    for m = find(um(:,c))' # parfor
                        uu = @(t) us + ut(t) .* sparse(m,1,um(m,c),M,1);
                        x = ODE(f,@(x,u,p,t) x,t,xs,uu,pr(:,k));
                        x -= avg(x,nf(1));
                        x *= 1.0/um(m,c);
                        W += DOT(x,x');
                    end;
                    for m = find(pm(:,c))' # parfor
                        pp = pr(:,k) + sparse(m,1,pm(m,c),P,1);
                        x = ODE(f,@(x,u,p,t) x,t,xs,up,pp);
                        x -= avg(x,nf(1));
                        x *= 1.0/pm(m,c);
                        o(m) += sum(sum(x.*x));
                    end;
                end;
            end;
            W *= h/(C*K);
            W = 0.5 * (W + W');
            if(A>0), W = {W,o}; end;

        case 'o' # Observability gramian
            W = zeros(N+A,N+A);
            o = zeros(Q*L,N+A);
            for k = 1:K
                for d = 1:D
                    for n = find(xm(:,d))' # parfor
                        xx = xs + sparse(n,1,xm(n,d),N,1);
                        y = ODE(f,g,t,xx,ux,pr(:,k));
                        y -= avg(y,nf(1));
                        y *= 1.0/xm(n,d);
                        o(:,n) = y(:);
                    end;
                    for n = find(pm(:,d))' # parfor
                        pp = pr(:,k) + sparse(n,1,pm(n,d),P,1);
                        y = ODE(f,g,t,xs,up,pp);
                        y -= avg(y,nf(1));
                        y *= 1.0/pm(n,d);
                        o(:,N+n) = y(:);
                    end;
                    W += DOT(o',o);
                end;
            end;
            W *= h/(D*K);
            W = 0.5 * (W + W');

        case 'x' # Cross gramian
            if(M~=Q && nf(7)==0), error('emgr: non-square system!'); end;

            if(isempty(DWX)) # Full cross gramian
                n0 = 0;
                a0 = -N;
                W = zeros(N,N+A);
                o = zeros(L,N+A,Q);

            else # Distributed cross gramian
                if(any(ceil(DWX)~=floor(DWX)) || any(DWX<=0))
                    error('emgr: non-positive or non-integer values in DWX!');
                end;
                i0 = (DWX(2)-1)*DWX(1) + 1;
                if(i0<=N) # State-space columns setup
                    i1 = min(i0+DWX(1)-1,N);
                    xm([1:i0-1,i1+1:end],:) = 0;
                    pm = zeros(1,D);
                    n0 = i0 - 1;
                else # Parameter-space columns setup
                    i0 = i0 - (ceil(N/DWX(1))*DWX(1) - N);
                    i1 = min(i0+DWX(1)-1,N+A);
                    xm = zeros(1,D);
                    pm([1:i0-N-1,i1-N+1:end],:) = 0;
                    a0 = i0 - N - 1;
                end;
                W = zeros(N,i1-i0+1);
                o = zeros(L,i1-i0+1,Q);

                if(i0>i1), W = 0; return; end;
            end;

            for k = 1:K
                for d = 1:D
                    for n = find(xm(:,d))' # parfor
                        xx = xs + sparse(n,1,xm(n,d),N,1);
                        y = ODE(f,g,t,xx,ux,pr(:,k));
                        y -= avg(y,nf(1));
                        y *= 1.0/xm(n,d);
                        o(:,n-n0,:) = y';
                    end;
                    for n = find(pm(:,d))' # parfor
                        pp = pr(:,k) + sparse(n,1,pm(n,d),P,1);
                        y = ODE(f,g,t,xs,up,pp);
                        y -= avg(y,nf(1));
                        y *= 1.0/pm(n,d);
                        o(:,n-a0,:) = y';
                    end;
                    if(nf(7)) # Non-symmetric cross gramian: cache average
                        o(:,:,1) = sum(o,3);
                    end;
                    for c = 1:C
                        for m = find(um(:,c))' # parfor
                            uu = @(t) us + ut(t) .* sparse(m,1,um(m,c),M,1);
                            x = ODE(f,@(x,u,p,t) x,t,xs,uu,pr(:,k));
                            x -= avg(x,nf(1));
                            x *= 1.0/um(m,c);
                            if(nf(7)) # Non-symmetric cross gramian
                                W += DOT(x,o(:,:,1));
                            else      # Regular cross gramian
                                W += DOT(x,o(:,:,m));
                            end;
                        end;
                    end;
                end;
            end;
            W *= h/(C*D*K);

        case 'y' # Linear cross gramian
            if(M~=Q && nf(7)==0), error('emgr: non-square system!'); end;
            W = zeros(N,N);
            o = zeros(N,L,M);
            for k = 1:K
                for c = 1:C
                    for m = find(um(:,c))' # parfor
                        uu = @(t) us + ut(t) .* sparse(m,1,um(m,c),M,1);
                        x = ODE(f,@(x,u,p,t) x,t,xs,uu,pr(:,k));
                        x -= avg(x,nf(1));
                        o(:,:,m) = x * (1.0/um(m,c));
                    end;
                    if(nf(7)) # Non-symmetric cross gramian: cache average
                        o(:,:,1) = sum(o,3);
                    end;
                    for m = find(xm(:,c))' # parfor
                        uu = @(t) us + ut(t) .* sparse(m,1,xm(m,c),Q,1);
                        z = ODE(g,@(x,u,p,t) x,t,xs,uu,pr(:,k));
                        z -= avg(z,nf(1));
                        z *= 1.0/xm(m,c);
                        if(nf(7)) # Non-symmetric cross gramian
                            W += DOT(o(:,:,1),z');
                        else      # Regular cross gramian
                            W += DOT(o(:,:,m),z');
                        end;
                    end;
                end;
            end;
            W *= h/(C*K);

        case 's' # Sensitivity gramian
            [pr,pm] = pscales(pr,size(um,2),nf(9));
            W = emgr_oct(f,g,[M,N,Q,P],t,'c',pr,nf,ut,us,xs,um,xm,pm);
            # W{1} # Controllability gramian
            # W{2} # Sensitivity gramian diagonal

        case 'i' # Identifiability gramian
            [pr,pm] = pscales(pr,size(xm,2),nf(9));
            V = emgr_oct(f,g,[M,N,Q,P],t,'o',pr,nf,ut,us,xs,um,xm,pm);
            W{1} = V(1:N,1:N);  # Observability gramian
            WM = V(1:N,N+1:N+P);
            if(nf(10))          # Identifiability gramian
                W{2} = V(N+1:N+P,N+1:N+P);
            else
                W{2} = V(N+1:N+P,N+1:N+P) - (WM'*ainv(W{1})*WM);
            end;

        case 'j' # Joint gramian
            [pr,pm] = pscales(pr,size(xm,2),nf(9));
            V = emgr_oct(f,g,[M,N,Q,P],t,'x',pr,nf,ut,us,xs,um,xm,pm);
            if(isempty(DWX)==0) # Distributed cross gramian
                W = V;
                return;
            end;
            W{1} = V(1:N,1:N);  # Cross gramian
            WM = V(1:N,N+1:N+P);
            if(nf(10))          # Cross-identifiability gramian
                W{2} = -0.5 * (WM'*WM);
            else
                W{2} = -0.5 * (WM'*ainv(W{1}+W{1}')*WM);
            end;

        otherwise
            error('emgr: unknown gramian type!');
    end;
end

## ======== INPUT AND STATE SCALES ========
function s = scales(s,d,e)
### summary: scales (Input and initial state perturbation scales)
#$

    switch(d)

        case 1 # Linear
            s *= [0.25,0.50,0.75,1.0];

        case 2 # Geometric
            s *= [0.125,0.25,0.5,1.0];

        case 3 # Logarithmic
            s *= [0.001,0.01,0.1,1.0];

        case 4 # Sparse
            s *= [0.38,0.71,0.92,1.0];
    end;

    if(e==0)
        s = [-s,s];
    end;
end

## ======== PARAMETER SCALES ========
function [pr,pm] = pscales(p,n,e)
### summary: pscales (Parameter perturbation scales)
#$

    if(size(p,2)==1), error('emgr: min + max parameter required!'); end;

    pmin = min(p,[],2);
    pmax = max(p,[],2);

    switch(e) # Parameter centering

        case 1 # Linear
            pr = 0.5 * (pmax + pmin);
            pm = (pmax - pmin) * linspace(0,1.0,n);
            pm = pm + pmin - pr;

        case 2 # Logarithmic
            lmin = log(pmin);
            lmax = log(pmax);
            pr = real(exp(0.5 * (lmax + lmin)));
            pm = (lmax - lmin) * linspace(0,1.0,n);
            pm = pm + lmin;
            pm = real(exp(pm)) + (pmin - pr);

        otherwise # None
            pr = pmin;
            pm = (pmax - pmin)*linspace(0,1.0,n);
    end;
end

## ======== TRAJECTORY AVERAGE ========
function m = avg(x,d)
### summary: avg (State and output trajectory centering)
#$
    switch(d)

        case 1 # Initial state / output
            m = x(:,1);

        case 2 # Steady state / output
            m = x(:,end);

        case 3 # Mean state / output
            m = mean(x,2);

        case 4 # Root-mean-square state / output
            m = sqrt(sum(x.*x,2));

        case 5 # Midrange state / output
            m = 0.5*(max(x,[],2)-min(x,[],2));

        otherwise # None
            m = zeros(size(x,1),1);
    end;
end

## ======== FAST APPROXIMATE INVERSION ========
function x = ainv(m)
### summary: ainv (Quadratic complexity approximate inverse matrix)
#$
    d = diag(m);
    d(d~=0) = 1.0./d(d~=0);
    n = numel(d);
    x = m .* (-d);
    x = x .* (d');
    x(1:n+1:end) = d;
end

## ======== DEFAULT ODE INTEGRATOR ========
function y = ssp2(f,g,t,x0,u,p)
### summary: ssp2 (Low-Storage Stability Preserving Runge-Kutta SSP32)
#$
    global STAGES;

    if(isscalar(STAGES)==0), STAGES = 3; end;

    h = t(1);
    K = floor(t(2)/h) + 1;

    x = x0;
    y = g(x,u(0),p,0);
    y(end,K) = 0; # Preallocate trajectory

    hs = h/(STAGES-1);
    xk = x;

    for k = 1:(K-1)

        tk = k*h;
        uk = u(tk);
        for s=1:(STAGES-1)

            xk += hs*f(xk,uk,p,tk);
            tk += hs;
        end
        xk = ((STAGES-1)*xk + x + h*f(xk,uk,p,tk))./STAGES;
        x = xk;
        y(:,k+1) = g(x,uk,p,tk);
    end;
end

