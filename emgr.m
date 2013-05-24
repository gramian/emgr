function W = emgr(f,g,q,p,t,w,nf,ut,us,xs,um,xm,yd)
% emgr - Empirical Gramian Framework (Version 1.1)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%
% SYNOPSIS:
%	W = emgr(f,g,q,p,t,w,[nf],[ut],[us],[xs],[um],[xm],[yd]);
%
% ABOUT:
%	emgr - Empirical Gramian Framemwork, computating empirical gramians
%	for model reduction and system identification.
%	Compatible with Octave and Matlab.
%
% INPUTS:
%		  (func)  f - system function handle, signature: xdot = f(x,u,p)
%		  (func)  g - output function handle, signature: y = g(x,u,p)
%		(vector)  q - system dimensions [inputs,states,outputs]
%		(vector)  p - parameters, 0 if none
%		(vector)  t - time [start,step,stop]
%		  (char)  w - gramian type:
%			* 'C' or 'c' : empirical controllability gramian (WC)
%			* 'O' or 'o' : empirical observability gramian (WO)
%			* 'X' or 'x' : empirical cross gramian (WX / WCO)
%			* 'S' or 's' : empirical sensitivity gramian (WS)
%			* 'I' or 'i' : empirical identifiability gramian (WI)
%			* 'J' or 'j' : empirical joint gramian (WJ)
%		(vector) [nf = 0] - options, 10 components
%			+ residualize against zero(0), average(1), last(2), steady-state(3), pca(4)
%			+ unit-normal(0), pca(1) directions
%			+ linear(0), log(1), rombergseq(2), single(3) input scale spacing
%			+ linear(0), log(1), rombergseq(2), single(3) init-state scale spacing
%			+ unit(0), [factorial(1)], single(3) rotations of input directions
%			+ unit(0), [factorial(1)], single(3) rotations of init-state directions
%			+ single(0), double(0) run
%			+ disable(0), enable(1)
%				* robust parameters (WC, WS only)
%				* data-driven pca (WO, WI only)
%				* less scales (WX, WJ only)
%			+ disable(0), enable(1) data-driven gramians
%			+ euler-method-1(0), adams-bashforth-2(1)
%  (matrix,vector,scalar) [ut = 1] - input; default: delta impulse
%         (vector,scalar) [us = 0] - steady-state input
%         (vector,scalar) [xs = 0] - steady-state
%  (matrix,vector,scalar) [um = 1] - input scales
%  (matrix,vector,scalar) [xm = 1] - init-state scales
%                (matrix) [yd = 0] - observed data
%
% OUTPUT:
%	        (matrix)  W - Gramian matrix (WC, WO, WX only)
%	          (cell)  W - {State,Parameter} Gramian Matrices (WS, WI, WJ only)
%
% TODO:
% 	factorial transformations
%
% For further information see http://gramian.de
%*

global x;		%Make x available to ode1x
global o;		%Make y available to ode1y

J = q(1);		%Number of Inputs
N = q(2);		%Number of States
O = q(3);		%Number of Outputs
p = p(:);		%Ensure parameters are in column vector
P = numel(p);		%Number of parameters
h = t(2);		%Time step
T = (t(3)-t(1))/h;	%Number of time steps
w = lower(w);		%Force lower case gramian type

if (nargin<7) ||(isempty(nf)) nf = 0; end;		%If options empty, set to zero
if (nargin<8) ||(isempty(ut)) ut = 1; end;		%If input empty, set to one
if (nargin<9) ||(isempty(us)) us = 0; end;		%If steady-state input empty, set to zero
if (nargin<10)||(isempty(xs)) xs = 0; end;		%If steady-state empty, set to zero
if (nargin<11)||(isempty(um)) um = 1; end; 	%If input scales empty, set to one
if (nargin<12)||(isempty(xm)) xm = 1; end; 	%If state scales empty, set to one
if (nargin<13)||(isempty(yd)) yd = cell(2,1); end;

if (numel(nf)<10) nf(10) = 0; end;		%If options scalar, extend to vector
if (numel(ut)==1) ut = ones(J,1)*ut; end;		%If input scalar, extend to vector
if (numel(us)==1) us = ones(J,1)*us; end;		%If steady-state input scalar, extend to vector
if (numel(xs)==1) xs = ones(N,1)*xs; end;		%If steady-state scalar, extend to vector
if (numel(um)==1) um = ones(J,1)*um; end;		%If input scales scalar, extend to vector
if (numel(xm)==1) xm = ones(N,1)*xm; end;		%If state scales scalar, extend to vector

if(w=='c' || w=='o' || w=='x')
	if(w=='c'&&nf(8)~=0)
		J = J+P;
		if(size(um,1)==J) um = [um;ones(P,1)]; end;
		F = f; f = @(x,u,p) F(x,u(1:J),u(J+1:J+P));
		G = g; g = @(x,u,p) G(x,u(1:J),u(J+1:J+P));
	end;

	if(size(um,2)==1) um = scales(um,nf(3),nf(5),w=='x'&&nf(8)==0); end;	%Generate Scales
	if(size(xm,2)==1) xm = scales(xm,nf(4),nf(6),w=='x'&&nf(8)==0); end;	%Generate Scales
	C = size(um,2);								%Number input scales
	D = size(xm,2);								%Number state scales

	if(size(us,2)==1) us = us*ones(1,T); end;					%Expand input steady state for each step
	if(size(ut,2)==1) k  = (1/h); ut = [ut,zeros(J,T-1)]; else k = 1; end;	%TODO ut sparse when octave compat

	if(nf(1)==4)&&(w=='o')&&(nf(8)~=0) nf(1) = 5; end;
	if(nf(2)==1)&&(w~='o') dx = pod(ut);                     else dx = 0; end;	%Set input directions
	if(nf(2)==1)&&(w~='c') dy = pod(ode1x(f,N,h,T,xs,us,p)); else dy = 0; end;	%Set state directions

	o = cell(N,1);	%Preallocate array of outputs
	x = zeros(N,T);	%Preallocate array pf states
	W = zeros(N,N);	%Preallocate gramian

	switch(nf(1))
		case 0
			X = f(zeros(N,1),zeros(J,1),p);
			Y = g(zeros(N,1),zeros(J,1),p);
		case 3
			X = f(xs,us(:,1),p);
			Y = g(xs,us(:,1),p);
		case 5
			X = 0;
			Y = pod(yd);
		otherwise
			X = 0;
			Y = 0;
	end;

	if(nf(7)==1)	%Double Run
		nf(7) = 0;
		A = emgr(f,g,q,p,t,w,nf,ut,us,xs,um,xm,yd);
		A = sqrt(diag(A));
		B = diag(1./A);
		A = diag(A);

		F = f;
		G = g;
		f = @(x,u,p) A*F(B*x,u,p);
		g = @(x,u,p)   G(B*x,u,p);
	end;

	sn = size(yd,1)==1;
end;

switch(w)												%Switch by gramian type
	case 'c'												%Type: controllability gramian
		for c=1:C										%For all input scales
			for j=1:J									%For all input components
				uu = us + bsxfun(@times,ut,dirs(j,J,dx)*(um(j,c)*k));			%Set up input
				if(nf(9)==0) odex(f,h,T,xs,uu,p,nf(10)); else x = yd{1,c}; end;		%Simulate (nonlinear) system
				x = bsxfun(@minus,x,steady(nf(1),x,X))*(1.0/um(j,c));			%Subtract scaled steady state
				W = W + x*x';								%Vectorized sum of dyadic products x(:,t)'*x(:,t) for all t
			end;
		end;
		W = W*(h/C);										%Symmetrize and normalize by number of scales
	case 'o'												%Type: observability gramian
		for d=1:D										%For all scales
			for n=1:N									%For all state components
				xx = xs + dirs(n,N,dy)*xm(n,d);						%Set up initial value
				if(nf(9)==0) odey(f,g,h,T,xx,us,p,n,nf(10)); else o{n} = yd{2-sn,d}; end;	%Simulate (nonlinear) system
				o{n} = bsxfun(@minus,o{n},steady(nf(1),o{n},Y))*(1.0/xm(n,d));		%Subtract scaled steady state
			end;
			for n=1:N									%For each row
				for m=1:N								%For each column
					W(n,m) = W(n,m) + o{n}(:)'*o{m}(:);				%Vectorized dot product of o{n}(:,t)*o{m}(:,t)' for all t
				end;
			end;
		end;
		W = W*(h/D);										%Symmetrize and normalize by number of scales
	case 'x'												%Type: cross gramian
		if(J~=O) error('ERROR: non-square system!'); end;						%error if non square
		for d=1:D										%For all state scales
			for n=1:N 									%For all state components
				xx = xs + dirs(n,N,dy)*xm(n,d);						%Set up initial value
				if(nf(9)==0) odey(f,g,h,T,xx,us,p,n,nf(10)); else o{n} = yd{2,d}; end;	%Simulate (nonlinear) system
				o{n} = bsxfun(@minus,o{n},steady(nf(1),o{n},Y))*(1.0/xm(n,d));		%Subtract scaled steady state
			end;
			for c=1:C									%For all input scales
				for j=1:J								%For all input components
					uu = us + bsxfun(@times,ut,dirs(j,J,dx)*(um(j,c)*k));		%Set up input
					if(nf(9)==0) odex(f,h,T,xs,uu,p,nf(10)); else x = yd{1,c}; end;	%Simulate (nonlinear) system
					x = bsxfun(@minus,x,steady(nf(1),x,X))*(1.0/um(j,c));		%Subtract scaled steady state
					for m=1:N							%For each column
						W(:,m) = W(:,m) + x*o{m}(j,:)';				%Sum product of control and observe components
					end;
				end;
			end;
		end;
		W = W*(h/(C*D));								%Symmetrize and normalize by number of scales
	case 's'										%Type: sensitivity gramian
		W = cell(2,1);								%Allocate return type
		W{1} = emgr(f,g,[J N O],zeros(P,1),t,'c',nf,ut,us,xs,um,xm); 		%Compute parameterless controllability gramian
		W{2} = eye(P);								%Allocate return type
		F = @(x,u,p) f(x,zeros(J,1),p*u);						%Augment system function with constant parameter input
		G = @(x,u,p) g(x,zeros(J,1),p*u);						%Adapt Output function
		for q=1:P								%For each parameter
			V = emgr(F,G,[1 N O],(1:P==q)*p(q),t,'c',nf,0,p(q),xs,1,xm);	%Compute parameter controllability Gramian
			W{2}(q,q) = trace(V);						%Trace of current parameter controllability gramian
			W{1} = W{1} + V;							%Accumulate controllability gramian
		end;
	case 'i'										%Type: identifiability gramian
		if(size(xm,1)==N) xm = [xm;ones(P,1)]; end;				%Augment state scales
		W = cell(2,1);								%Allocate return type
		F = @(x,u,p) [f(x(1:N),u,x(N+1:N+P));zeros(P,1)];				%Augment system function with constant parameter states
		G = @(x,u,p)  g(x(1:N),u,x(N+1:N+P));					%Adapt Output function
		V = emgr(F,G,[J N+P O],0,t,'o',nf,ut,us,[xs;p],um,xm); 			%Compute Observability Gramian of augmented system
		W{1} = V(1:N,1:N);							%Extract Observability gramian
		V = chol(V+speye(N+P));							%Use shortcut to compute schur complement of augmented with ensured positive definiteness
		V = V(N+1:N+P,N+1:N+P);							%Extract parameter observability gramian using cholesky factors
		W{2} = V'*V - speye(P);	 						%Compute identifiability gramian
	case 'j'										%Type: joint gramian
		if(size(um,1)==J) um = [um;ones(P,1)]; end;				%Augment input scales
		if(size(xm,1)==N) xm = [xm;ones(P,1)]; end;				%Augment state scales
		up = [p+(p==0),zeros(P,size(ut,2)-1)];
		W = cell(2,1);								%Allocate return type
		F = @(x,u,p) [f(x(1:N),u(1:J),x(N+1:N+P));u(J+1:J+P)];			%Augment system function with constant parameter states
		G = @(x,u,p) [g(x(1:N),u(1:J),x(N+1:N+P));x(J+1:J+P)];			%Adapt Output function
		V = emgr(F,G,[J+P N+P O+P],0,t,'x',nf,[ut;up],[us;zeros(P,1)],[xs;p],um,xm);%Compute Cross Gramian of double augmented system
		S = norm(V,1);								%Bound largest eigenvalue
		W{1} = V(1:N,1:N);							%Extract Cross gramian
		V = chol(V+S*speye(N+P));							%Use shortcut to compute schur complement of augmented wit ensured positive definiteness
		V = V(N+1:N+P,N+1:N+P);							%Extract parameter cross gramian using cholesky factors
		W{2} = V'*V - S*speye(P);							%Compute cross identifiability gramian
	otherwise
		error('ERROR: unknown gramian type!');
end;

if(w=='c' || w=='o' || w=='x') W = 0.5*(W+W'); end; 	%Enforce Symmetry

end	%end emgr

%********

function d = dirs(n,N,e)

	switch(e)
		case 0    %Unit-Normal
			d = (1:N==n)';
		otherwise %PCA
			d = e;
	end;

end	%end dirs

%********

function s = scales(s,d,e,f)

	switch(d)
		case 0 %Linear
			if(f) s = s*[1/2,1];   else s = s*[1/4,1/2,3/4,1]; end;
		case 1 %Logarithmic
			if(f) s = s*[1/100,1]; else s = s*[1/1000,1/100,1/10,1]; end;
		case 2 %Romberg Sequence
			if(f) s = s*[1/4,1];   else s = s*[1/8,1/4,1/2,1]; end;
		case 3 %Single
			%s = s;
	end;

	switch(e)
		case 0 %Unit
			s = [-s,s];
		case 1 %Factorial
			s = s*(2*(dec2bin(0:2^q-1)-'0')'-1)./sqrt(2^q);	%not beautiful but working
		case 2 %Ideas: All, Dyadic

		case 3 %Single
			%s = s;
	end;

end	%end scales

%********

function y = steady(v,d,e)

	switch(v)
		case 0 %Zero
			y = e;
		case 1 %Average
			y = mean(d,2);
		case 2 %Last
			y = d(end,:);
		case 3 %Steady State
			y = e;
		case 4 %PCA
			y = pod(d);
		case 5 %Data PCA
			y = e;
	end;

end	%end steady

%********

function p = pod(y)

	p = svd(cov(bsxfun(@minus,y,mean(y,2))'),'econ');

end	%end pod

%********

function odex(f,h,T,z,u,p,q)

global x;

switch(q)
	case 0 %Eulers Method
		for t=1:T
			z = z + h*f(z,u(:,t),p);
			x(:,t) = z;
		end;
	case 1 %Adams-Bashforth Method
		m = h*f(z + 0.5*h*f(z,u(:,1),p),u(:,1),p);
		z = z + m;
		x(:,1) = z;
		m = 0.5*m;

		for t=2:T
			k = 0.5*h*f(z,u(:,t),p);
			z = z + 3.0*k-m;
			x(:,t) = z;
			m = k;
		end;
	end;

end	%end odex

%********

function odey(f,g,h,T,z,u,p,n,q)

global o;

switch(q)
	case 0 %Eulers Method
		for t=1:T
			z = z + h*f(z,u(:,t),p);
			o{n}(:,t) = g(z,u(:,t),p);
		end;
	case 1 %Adams-Bashforth Method
		m = h*f(z + 0.5*h*f(z,u(:,1),p),u(:,1),p);
		z = z + m;
		o{n}(:,1) = g(z,u(:,1),p);
		m = 0.5*m;

		for t=2:T
			k = 0.5*h*f(z,u(:,t),p);
			z = z + 3.0*k-m;
			o{n}(:,t) = g(z,u(:,t),p);
			m = k;
		end;
	end;

end	%end odey
