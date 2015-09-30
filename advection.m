function advection(o)
% finite difference discretized pde transport equation reduction
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        disp('emgr framework is required. Download at http://gramian.de/emgr.m');
        return;
    else
        global ODE;
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

    %% Setup
    J = 0;
    N = 100;
    O = N;
    R = 20;
    T = [0.0,0.01,1.0];
    L = (T(3)-T(1))/T(2);
    U = zeros(1,L);
    H = 0.1;
    X = 0:H:9.9; X = exp(-(X-1).^2)';

    P = 5;
    A = spdiags((1.0/H)*[ones(N,1),-ones(N,1)],[-1,0],N,N);

    LIN = @(x,u,p) p*A*x;
    ADJ = @(x,u,p) p*A'*x + u;
    OUT = @(x,u,p) x;

    %% Main
    Y = ODE(LIN,OUT,T,X,U,P); % Full Order

    tic;
    WO = emgr(ADJ,1,[N,N,O],T,'c',P);
    [UU,D,VV] = svd(WO); UU = UU(:,1:R); VV = VV(:,1:R)';
    a = VV*A*UU;
    c = UU;
    x = VV*X;
    lin = @(x,u,p) p*a*x;
    out = @(x,u,p) c*x;
    OFFLINE = toc

    y = ODE(lin,out,T,x,U,P);

    %% Output
    if(nargin==0), return; end
    cmax = max(y(:));
    figure();
    imagesc(sparse(y)); caxis([0 cmax]);
    colormap(antijet);
    set(gca,'YTick',0,'xtick',[]); ylabel('X');
    pbaspect([2,1,1]);
    if(o==1), print('-dpng',[mfilename(),'.png']); end;
end

%% ======== Colormap ========
function m = antijet(n)
% antijet colormap
% by Christian Himpe 2014-2015
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(nargin<1 || isempty(n)), n = 256; end;
    L = linspace(0,1,n);

    R = -0.5*sin( L*(1.37*pi)+0.13*pi )+0.5;
    G = -0.4*cos( L*(1.5*pi) )+0.4;
    B = 0.3*sin( L*(2.11*pi) )+0.3;

    m = [R;G;B]';
end
