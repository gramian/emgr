function decentral(o)
%%% summary: decentral (decentralized control)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2013--2017)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

%% SYSTEM SETUP
    M = 4;			% number of inputs
    N = 16;			% number of states
    Q = M;			% number of outputs
    h = 0.01;			% time step size
    T = 1.0;			% time horizon

    rand('seed',1009);
    A = rand(N,N);		% \
    A(1:N+1:end) = -0.55*N;	%  system matrix
    B = rand(N,M);		% input matrix
    C = rand(Q,N);		% output matrix

%% SUB-SYSTEM GRAMIAN COMPUTATION

    dp = @(x,y) diag(sum(x.*y',2)); % diagonal-only pseudo-kernel

    tic;
    PM = zeros(M,Q); % Participation Matrix

    for V=1:M
        for W=1:Q
            lin = @(x,u,p,t) A*x + B(:,V)*u; 					% SISO vector field
            out = @(x,u,p,t) C(W,:)*x;						% SISO output functional
            io = emgr(lin,out,[1,N,1],[h,T],'x',0,0,1,0,0,1,1,dp);		% SISO cross gramian
            PM(V,W) = trace(io)^2 / trace(io.^2);				% interaction measure
        end
    end

    PM = PM./sum(PM(:));	% Normalize participation matrix
    OFFLINE = toc

%% PLOT PARTICIPATION MATRIX
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    imagesc(PM); colorbar;
    set(gca,'XTick',1:1:M,'YTick',1:1:Q);
    if(nargin>0 && o==1), print('-dpng',[mfilename(),'.png']); end;
end

