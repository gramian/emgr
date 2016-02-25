function decentral(o)
% decentral (-ized control)
% by Christian Himpe, 2013-2016 (http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

%% SETUP
    J = 4;
    N = 16;
    O = J;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = [ones(J,1),zeros(J,L-1)];
    X = ones(N,1);

    rand('seed',1009);
    A = rand(N,N);
    A(1:N+1:end) = -0.55*N;
    B = rand(N,J);
    C = rand(O,N);

    LIN = @(x,u,p) A*x + B*u;
    OUT = @(x,u,p) C*x;

%% OFFLINE
    tic;
    WX = cell(J,O);
    for V=1:J
        for W=1:O
            lin = @(x,u,p) A*x + B(:,V)*u;
            out = @(x,u,p) C(W,:)*x;
            WX{V,W} = emgr(lin,out,[1,N,1],T,'x');
        end
    end

    PM = zeros(J,O); % Participation Matrix
    PM = cellfun(@(w) trace(w),WX);
    PM = PM./sum(sum(PM));

    OFFLINE = toc

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    imagesc(PM);
    colormap(antijet); colorbar;
    set(gca,'XTick',1:1:J,'YTick',1:1:O);
    if(nargin>0 && o==1), print('-dpng',[mfilename(),'.png']); end;
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
