function W = emgr_oct(f,g,s,t,w,pr=0,nf=0,ut=1,us=0,xs=0,um=1,xm=1,pm=[])
# emgr - Empirical Gramian Framework ( Version: 4.0.oct )
# Copyright (c) 2013-2016 Christian Himpe ( gramian.de )
# released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
#
# SYNTAX:
#    W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[pm]);
#
# SUMMARY:
#    emgr - EMpirical GRamian framework,
#    computation of empirical gramians for model reduction,
#    system identification and uncertainty quantification.
#    Enables gramian-based nonlinear model order reduction.
#    Compatible with OCTAVE and MATLAB.
#
# ARGUMENTS:
#   (func handle)  f - system function handle; signature: xdot = f(x,u,p)
#   (func handle)  g - output function handle; signature:    y = g(x,u,p)
#        (vector)  s - system dimensions [inputs,states,outputs]
#        (vector)  t - time discretization [step,stop]
#          (char)  w - gramian type:
#            * 'c' : empirical controllability gramian (WC)
#            * 'o' : empirical observability gramian (WO)
#            * 'x' : empirical cross gramian (WX aka WCO or XCG)
#            * 'y' : empirical linear cross gramian (WY)
#            * 's' : empirical sensitivity gramian (WS)
#            * 'i' : empirical identifiability gramian (WI)
#            * 'j' : empirical joint gramian (WJ)
# (matrix,vector,scalar) [pr = 0] - parameters, each column is one set
#        (vector,scalar) [nf = 0] - options, ten components:
#            + zero(0), init(1), steady(2), mean(3), rms(4) centering
#            + single(0), linear(1), geom(2), log(3), sparse(4) input scales
#            + single(0), linear(1), geom(2), log(3), sparse(4) state scales
#            + unit(0), single(1) input rotations
#            + unit(0), single(1) state rotations
#            + single(0), preconditioned double(1), steady-state scaled(2) run
#            + regular(0), non-symmetric(1) cross gramian; only: WX, WY, WJ
#            + active(0), passive(1) parameter; only: WS, WI, WJ
#            + none(0), linear(1), geom(2) parameter centering; only: WS, WI, WJ
#            + approximate(0), detailed(1) Schur-complement; only: WI, WJ
#  (matrix,vector,scalar) [ut = 1] - input; default: delta impulse
#         (vector,scalar) [us = 0] - steady-state input
#         (vector,scalar) [xs = 0] - steady-state and initial state x0
#  (matrix,vector,scalar) [um = 1] - input scales
#  (matrix,vector,scalar) [xm = 1] - initial-state scales
#                (matrix) [pm = 0] - parameter scales (reserved)
#
# RETURNS:
#            (matrix)  W - Gramian Matrix (only: WC, WO, WX, WY)
#              (cell)  W - {State-,Parameter-} Gramian (only: WS, WI, WJ)
#
# CITATION:
#    C. Himpe (2016). emgr - Empirical Gramian Framework (Version 4.0)
#    [Software]. Available from http://gramian.de . doi:10.5281/zenodo.56502 .
#
# SEE ALSO:
#    gram
#
# KEYWORDS:
#    model reduction, empirical gramian, cross gramian, mor
#
# Further information: <http://gramian.de>
#*
    global DOT; # Inner Product Handle
    global ODE; # Integrator Handle
    global DWX; # Distributed Cross Gramian (Width,Index)

    if(isa(DOT,'function_handle')==0), DOT = @mtimes; end;
    if(isa(ODE,'function_handle')==0), ODE = @rk2; end;

    # Version Info
    if( (nargin==1) && strcmp(f,'version') ), W = 4.0; return; end;

    # System Dimensions
    J = s(1);              # number of inputs
    N = s(2);              # number of states
    O = s(3);              # number of outputs
    M = 0;                 # number of parameter-inputs or parameter-states
    if(numel(s)>3), M = s(4); end;
    P = size(pr,1);        # dimension of parameter
    Q = size(pr,2);        # number of parameters
    h = t(1);              # width of time step
    L = floor(t(2)/h) + 1; # number of time steps plus initial value
    w = lower(w);          # ensure lower case gramian type

    # Decaying Linear Chirp Input 
    if( isnumeric(ut) && numel(ut)==1 && ut==Inf )
        ut = @(t) 0.5 + 0.5*cos(2.0*pi*(((0.1/h)+0.5*((1.0/L-0.1/h)/L)*t).*t));
    end;

    # Discretize Procedural Input
    if(isa(ut,'function_handle'))
        uf = ut;
        ut = zeros(J,L);
        for l=1:L
            ut(:,l) = uf(l*h);
        end;
    end;

    # Lazy Arguments
    if( isnumeric(g) && g==1 ), g = @(x,u,p) x; O = N; end;

    if(numel(nf)<10), nf(10) = 0; end;
    if(numel(ut)==1), ut = (ut/h)*ones(J,1); end;
    if(numel(us)==1), us = us*ones(J,1); end;
    if(numel(xs)==1), xs = xs*ones(N,1); end;
    if(numel(um)==1), um = um*ones(J,1); end;
    if(numel(xm)==1), xm = xm*ones(N-(w=='y')*(N-J),1); end;

    if(size(ut,2)==1), ut = [ut,zeros(J,L-1)]; end;
    if(size(us,2)==1), us = repmat(us,[1,L]); end;
    if(size(um,2)==1), um = scales(um,nf(2),nf(4)); end;
    if(size(xm,2)==1), xm = scales(xm,nf(3),nf(5)); end;

    # State-Space Setup
    if( w=='c' || w=='o' || w=='x' || w=='y' )

        C = size(um,2); # number of input scales
        D = size(xm,2); # number of state scales

        if(isempty(pm)), pm = sparse(1,max(C,D)); end;

        switch(nf(6))

            case 1, # preconditioned run
                nf(6) = 0;
                WT = emgr_oct(f,g,s,t,w,pr,nf,ut,us,xs,um,xm,pm);
                TX = sqrt(WT(1:size(WT,1)+1:N*size(WT,1)));
                tx = 1.0./TX;
                F = f; f = @(x,u,p) TX.*F(tx.*x,u,p);
                G = g; g = @(x,u,p)     G(tx.*x,u,p);

            case 2, # steady state (input) scaled run
                TU = us(:,1); TU(TU==0) = 1.0;
                tu = 1.0./TU;
                TX = xs; TX(TX==0) = 1.0;
                tx = 1.0./TX;
                F = f; f = @(x,u,p) TX.*F(tx.*x,tu.*u,p);
                G = g; g = @(x,u,p)     G(tx.*x,tu.*u,p); 
        end;

        ua = us;  # active parameters
        if(nf(8)) # passive parameters
            ua = ut + us; 
        end;
    end;

## GRAMIAN COMPUTATION

    switch(w) # empirical gramian types

        case 'c', # controllability gramian
            W = zeros(N,N);
            o = zeros(N,ceil(M/N));
            for q=1:Q
                for c=1:C 
                    for j = find(um(:,c))' # parfor
                        uu = us;
                        uu(j,:) += ut(j,:) * um(j,c);
                        x = ODE(f,1,t,xs,uu,pr(:,q));
                        x -= avg(x,nf(1));
                        x *= 1.0./um(j,c);
                        W += DOT(x,x');
                    end;
                    for m = find(pm(:,c))' # parfor
                        pp = pr;
                        pp(m,1) += pm(m,c);
                        x = ODE(f,1,t,xs,ua,pp);
                        x -= avg(x,nf(1));
                        x *= 1.0./pm(m,c);
                        W += DOT(x,x');
                        o(m) += sum(sum(x.*x));
                    end;
                end;
            end;
            W *= h/(C*Q);
            W = 0.5*(W+W');
            W = [W,o];

        case 'o', # observability gramian
            W = zeros(N+M,N+M);
            o = zeros(O*L,N+M);
            for q=1:Q
                for d=1:D
                    for n = find(xm(:,d))' # parfor
                        xx = xs;
                        xx(n,1) += xm(n,d);
                        y = ODE(f,g,t,xx,us,pr(:,q));
                        y -= avg(y,nf(1));
                        y *= 1.0/xm(n,d);
                        o(:,n) = y(:);
                    end;
                    for m = find(pm(:,d))' # parfor
                        pp = pr;
                        pp(m,1) += pm(m,d);
                        y = ODE(f,g,t,xs,ua,pp);
                        y -= avg(y,nf(1));
                        y *= 1.0/pm(m,d);
                        o(:,N+m) = y(:);
                    end;
                    W += DOT(o',o);
                end;
            end;
            W *= h/(D*Q);
            W = 0.5*(W+W');

        case 'x', # cross gramian
            if(J~=O && nf(7)==0), error('ERROR! emgr: non-square system!'); end;

            if(isempty(DWX)) # full cross gramian
                n0 = 0;
                m0 = -N;
                W = zeros(N,N+M);
                o = zeros(L,N+M,O,D);

            else # distributed cross gramian column selection
                i0 = round((DWX(2)-1)*DWX(1)) + 1;
                if(i0<=N)
                    i1 = min(i0+DWX(1)-1,N);
                    xm([1:i0-1,i1+1:end],:) = 0;
                    pm = zeros(1,D);
                    n0 = i0 - 1;
                else
                    i0 = i0 - (ceil(N/DWX(1))*DWX(1) - N);
                    i1 = min(i0+DWX(1)-1,N+M);
                    xm = zeros(1,D);
                    pm([1:i0-N-1,i1-N+1:end],:) = 0;
                    m0 = i0 - N - 1;
                end;
                W = zeros(N,i1-i0+1);
                o = zeros(L,i1-i0+1,O,D);

                if(i0>i1), W = 0; return; end;
            end;

            for q=1:Q
                for d=1:D
                    for n = find(xm(:,d))' # parfor
                        xx = xs;
                        xx(n,1) += xm(n,d);
                        y = ODE(f,g,t,xx,us,pr(:,q));
                        y -= avg(y,nf(1));
                        y *= 1.0/xm(n,d);
                        o(:,n-n0,:,d) = y';
                    end;
                    for m = find(pm(:,d))' # parfor
                        pp = pr;
                        pp(m,1) += pm(m,d);
                        y = ODE(f,g,t,xs,ua,pp);
                        y -= avg(y,nf(1));
                        y *= 1.0/pm(m,d);
                        o(:,m-m0,:,d) = y';
                    end;
                end;
                if(nf(7)) # non-symmetric cross gramian: cache average
                    o(:,:,1,:) = sum(o,3);
                end;
                for c=1:C
                    for j = find(um(:,c))' # parfor
                        uu = us;
                        uu(j,:) += ut(j,:) * um(j,c);
                        x = ODE(f,1,t,xs,uu,pr(:,q));
                        x -= avg(x,nf(1));
                        x *= 1.0/um(j,c);
                        for d=1:D
                            if(nf(7)) # non-symmetric cross gramian
                                W += DOT(x,o(:,:,1,d));
                            else      # regular cross gramian
                                W += DOT(x,o(:,:,j,d));
                            end;
                        end;
                    end;
                end;
            end;
            W *= h/(C*D*Q);

        case 'y', # linear cross gramian
            if(J~=O && nf(7)==0), error('ERROR! emgr: non-square system!'); end;
            W = zeros(N,N);
            o = zeros(N,L,J);
            for q=1:Q
                for c=1:C
                    for j = find(um(:,c))' # parfor
                        uu = us;
                        uu(j,:) += ut(j,:) * um(j,c);
                        x = ODE(f,1,t,xs,uu,pr(:,q));
                        x -= avg(x,nf(1));
                        o(:,:,j) = x * (1.0./um(j,c));
                    end;
                    if(nf(7)) # non-symmetric cross gramian: cache average
                        o(:,:,1) = sum(o,3);
                    end;
                    for j = find(xm(:,c))' # parfor
                        uu = us;
                        uu(j,:) += ut(j,:) * xm(j,c);
                        z = ODE(g,1,t,xs,uu,pr(:,q));
                        z -= avg(z,nf(1));
                        z *= 1.0./xm(j,c);
                        if(nf(7)) # non-symmetric cross gramian
                            W += DOT(o(:,:,1),z');
                        else      # regular cross gramian
                            W += DOT(o(:,:,j),z');
                        end;
                    end;
                end;
            end;
            W *= h/(C*Q);

        case 's', # sensitivity gramian
            [pr,pm] = pscales(pr,size(um,2),nf(9));
            V = emgr_oct(f,g,[J,N,O,P],t,'c',pr,nf,ut,us,xs,um,xm,pm);
            W{1} = V(1:N,1:N);     # robust controllability gramian
            W{2} = V(N*N+1:N*N+P); # sensitivity gramian (diagonal)

        case 'i', # identifiability gramian
            [pr,pm] = pscales(pr,size(xm,2),nf(9));
            V = emgr_oct(f,g,[J,N,O,P],t,'o',pr,nf,ut,us,xs,um,xm,pm);
            W{1} = V(1:N,1:N);         # observability gramian
            W{2} = V(N+1:N+P,N+1:N+P); # identifiability gramian
            if(nf(10))
                W{2} -= V(N+1:N+P,1:N)*pinv(W{1})*V(1:N,N+1:N+P);
            end;

        case 'j', # joint gramian
            [pr,pm] = pscales(pr,size(xm,2),nf(9));
            V = emgr_oct(f,g,[J,N,O,P],t,'x',pr,nf,ut,us,xs,um,xm,pm);
            if(~isempty(DWX)), W = V; return; end;
            W{1} = V(1:N,1:N); # cross gramian
            if(nf(10))         # cross-identifiability gramian
                W{2} = -0.5*V(1:N,N+1:N+P)'*pinv(W{1}+W{1}')*V(1:N,N+1:N+P);
            else
                W{2} = -0.5*V(1:N,N+1:N+P)'*ainv(W{1}+W{1}')*V(1:N,N+1:N+P);
            end;

        otherwise,
            error('ERROR! emgr: unknown gramian type!');
    end;
end

## ======== PARAMETER SCALES ========
function [pr,pm] = pscales(p,n,e)

    if(size(p,2)==1), error('ERROR! emgr: min + max parameter required!'); end;

    pmin = min(p,[],2);
    pmax = max(p,[],2);

    switch(e) # parameter centering

        case 1, # linear
            pr = 0.5*(pmax + pmin);
            pm = (pmin - pr) + (pmax - pmin).*linspace(0,1.0,n);

        case 2, # logarithmic
            lmin = log(pmin);
            lmax = log(pmax);
            lmid = 0.5*(lmax + lmin);
            pr = real(exp(lmid));
            pm = real(exp(lmin + (lmax - lmin).*linspace(0,1.0,n))) - pr;

        otherwise, # none
            pr = pmin;
            pm = (pmax - pmin).*linspace(0,1.0,n);
    end;
end

## ======== INPUT AND STATE SCALES ========
function s = scales(s,d,e)

    switch(d)

        case 1, # linear
            s = s*[0.25,0.50,0.75,1.0];

        case 2, # geometric
            s = s*[0.125,0.25,0.5,1.0];

        case 3, # logarithmic
            s = s*[0.001,0.01,0.1,1.0];

        case 4, # sparse
            s = s*[0.38,0.71,0.92,1.0];
    end;

    if(e==0)
        s = [-s,s];
    end;
end

## ======== TRAJECTORY AVERAGE ========
function m = avg(x,d)

    switch(d)

        case 1, # initial state (output)
            m = x(:,1);

        case 2, # steady state (output)
            m = x(:,end);

        case 3, # mean state (output)
            m = mean(x,2);

        case 4, # root mean square state (output)
            m = sqrt(sum(x.*x,2));

        otherwise, # zero
            m = zeros(size(x,1),1);
    end;
end

## ======== FAST APPROXIMATE INVERSION ========
function x = ainv(m)

    d = diag(m);
    d(d~=0) = 1.0./d(d~=0);
    n = numel(d);
    x = bsxfun(@times,m,-d);
    x = bsxfun(@times,x,d');
    x(1:n+1:end) = d;
end

## ======== DEFAULT ODE INTEGRATOR ========
function x = rk2(f,g,t,z,u,p)

    if(isnumeric(g) && g==1), g = @(x,u,p) x; end;

    h = t(1);
    L = floor(t(2)/h) + 1;

    x(:,1) = g(z,u(:,end),p);
    x(end,L) = 0; # preallocate trajectory

    for l=2:L # 2nd order Heun's Runge-Kutta Method
        k1 = f(z,u(:,l-1),p);
        k2 = f(z + h*k1,u(:,l-1),p);
        z += (0.5*h)*(k1 + k2);
        x(:,l) = g(z,u(:,l-1),p);
    end;
end

