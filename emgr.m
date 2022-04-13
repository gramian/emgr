function W = emgr(f,g,s,t,w,pr,nf,ut,us,xs,um,xm,dp)
%% emgr - EMpirical GRamian Framework
%
%  project: emgr ( https://gramian.de )
%  version: 5.99 (2022-04-13)
%  authors: C. Himpe (0000-0003-2194-6754)
%  license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%  summary: Empirical system Gramians for (nonlinear) input-output systems.
%
% DESCRIPTION:
%
%  Empirical gramian matrix and empirical covariance matrix computation
%  for model reduction, decentralized control, nonlinearity quantification,
%  sensitivity analysis, parameter identification, uncertainty quantification &
%  combined state and parameter reduction of large-scale input-output systems.
%  Data-driven analysis of input-output coherence and system-gramian-based
%  nonlinear model order reduction. Compatible with OCTAVE and MATLAB.
%
% BRIEF:
%
%  Unsupervised learning of I/O system properties for data-driven control.
%
% ALGORITHM:
%
%  C. Himpe (2018). emgr - The Empirical Gramian Framework. Algorithms 11(7):91
%  <https://doi.org/10.3390/a11070091 doi:10.3390/a11070091>
%
% USAGE:
%
%  W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[dp])
%
% MANDATORY ARGUMENTS:
%
%   f {handle} vector field: x' = f(x,u,p,t)
%   g {handle} output functional: y = g(x,u,p,t)
%   s {vector} system dimensions: [inputs, states, outputs]
%   t {vector} time discretization: [time-step, time-horizon]
%   w  {char}  empirical gramian type:
%    * 'c' empirical controllability gramian (Wc)
%    * 'o' empirical observability gramian (Wo)
%    * 'x' empirical cross gramian (Wx aka Wco)
%    * 'y' empirical linear cross gramian (Wy)
%    * 's' empirical sensitivity gramian (Ws)
%    * 'i' empirical identifiability gramian (Wi)
%    * 'j' empirical joint gramian (Wj)
%
% OPTIONAL ARGUMENTS:
%
%  pr {matrix|0} parameter vector(s), each column is one parameter sample
%  nf {vector|0} option flags, thirteen component vector, default all zero:
%    * centering: none(0), steady(1), last(2), mean(3), rms(4), midrange(5)
%    * input scales: single(0), linear(1), geometric(2), log(3), sparse(4)
%    * state scales: single(0), linear(1), geometric(2), log(3), sparse(4)
%    * input rotations: unit(0), single(1)
%    * state rotations: unit(0), single(1)
%    * normalization (only: Wc, Wo, Wx, Wy): none(0), steady(1), Jacobi(2)
%    * state gramian variant:
%      * controllability gramian type (only: Wc, Ws): regular(0), output(1)
%      * observability gramian type (only: Wo, Wi): regular(0), averaged(1)
%      * cross gramian type (only: Wx, Wy, Wj): regular(0), non-symmetric(1)
%    * extra input (only: Wo, Wx, Ws, Wi, Wj): no(0), yes(1)
%    * parameter centering (only: Ws, Wi, Wj): none(0), lin(1), log(2), nom(3)
%    * parameter gramian variant:
%      * averaging type (only: Ws): input-state(0), input-output(1)
%      * Schur-complement (only: Wi, Wj): approx(0), coarse(1), exact(2)
%    * cross gramian partition size (only: Wx, Wj): full(0), partitioned(<N)
%    * cross gramian partition index (only: Wx, Wj): partition(>0)
%    * weighting: none(0), linear(1), squared(2), state(3), scale(4), rsqrt(5)
%  ut {handle|'i'} input function: u_t = ut(t) or character:
%    * 'i' delta impulse input
%    * 's' step input / load vector / source term
%    * 'h' havercosine decaying exponential chirp input
%    * 'a' sinc (cardinal sine) input
%    * 'r' pseudo-random binary input
%  us {vector|0} steady-state input (1 or #inputs rows)
%  xs {vector|0} steady-state and nominal initial state x_0 (1 or #states rows)
%  um {matrix|1} input scales (1 or #inputs rows)
%  xm {matrix|1} initial-state scales (1 or #states rows)
%  dp {handle|@mtimes} inner product or kernel: xy = dp(x,y)
%
% RETURNS:
%
%  W {matrix} State-space system Gramian Matrix (for: Wc, Wo, Wx, Wy)
%  W  {cell}  {State, Parameter}-space system Gramian (for: Ws, Wi, Wj)
%
% CITE AS:
%
%  C. Himpe (2022). emgr - EMpirical GRamian Framework (Version 5.99)
%  [Software]. Available from https://gramian.de . doi:10.5281/zenodo.6457616
%
% KEYWORDS:
%
%  model reduction, system gramians, empirical gramians, cross gramian, MOR
%
% SEE ALSO: gram (Control System Toolbox)
%
% COPYRIGHT: Christian Himpe
%
% For more information, see: <https://gramian.de>

  % Set Integrator Handle (i.e. for custom solvers)
  global ODE;
  if not(isa(ODE,'function_handle')), ODE = @ssp2; end%if

  % Version Info (and export default local integratorvia ODE)
  if isequal(f,'version'), W = 5.99; return; else, fState = f; end%if

  % Default Arguments
  if (nargin <  6) || isempty(pr), pr = 0.0; end%if
  if (nargin <  7) || isempty(nf), nf = 0;   end%if
  if (nargin <  8) || isempty(ut), ut = 'i'; end%if
  if (nargin <  9) || isempty(us), us = 0.0; end%if
  if (nargin < 10) || isempty(xs), xs = 0.0; end%if
  if (nargin < 11) || isempty(um), um = 1.0; end%if
  if (nargin < 12) || isempty(xm), xm = 1.0; end%if
  if (nargin < 13) || isempty(dp), dp = @mtimes; end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP

  % System Dimensions
  nInputs = s(1);   % Number of system inputs / controls
  nStates = s(2);   % Number of system states / degrees of freedom (DoF)
  nOutputs = s(3);  % Number of system outputs / quantities of interest (QoI)

  % Parameter Dimensions
  nParams = size(pr,1);        % Number of parameters / parameter dimension
  nParamSamples = size(pr,2);  % Number of parameter samples

  % Time Discretization
  tStep = t(1);                        % Time-step width
  tFinal = t(2);                       % Time horizon
  nSteps = floor(tFinal / tStep) + 1;  % Number of time-steps

  % Gramian Type
  gramianType = lower(w);

  % Flag Vector
  flags = [nf(:)',zeros(1,max(13,13-numel(nf)))];

  % Built-in Input Functions
  if not(isa(ut,'function_handle'))

    a0 = (pi / (2.0 * tStep)) * tFinal / log(4.0 * (tStep / tFinal));
    b0 = (4.0 * (tStep / tFinal)) ^ (1.0 / tFinal);

    switch lower(ut)
      case 'i',  fExcite = @(t) (t <= tStep) / tStep;				% Impulse Input
      case 's',  fExcite = @(t) 1.0;						% Step Input
      case 'h',  fExcite = @(t) 0.5 * cos(a0 * (b0 ^ t - 1.0)) + 0.5;		% Havercosine Chirp Input
      case 'a',  fExcite = @(t) sin(t / tStep) / ((t / tStep) + (t == 0));	% Sinc Input
      case 'r',  fExcite = @(t) randi([0,1],1,1);				% Pseudo-Random Binary Input
      otherwise, error(' emgr: Unknown input ut!');
    end%switch
  else
    fExcite = ut;
  end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONFIGURATION

  % Output Function
  if (isnumeric(g) && (1 == g)) || ...
     (strcmp(gramianType,'c') && not(flags(7))) || strcmp(gramianType,'y')
    fOutput = @id;
    fAdjoint = g;
  else
    fOutput = g;
  end%if

  % Trajectory Weighting
  tInstances = [0.5*tStep,tStep:tStep:tFinal];
  switch flags(13)
    case 1,    fWeight = @(traj) traj .* sqrt(tInstances);			% Linear time-weighting
    case 2,    fWeight = @(traj) traj .* (tInstances ./ sqrt(2.0));		% Quadratic time-weighting
    case 3,    fWeight = @(traj) traj ./ max(sqrt(eps),vecnorm(traj,2,1));	% State-weighting
    case 4,    fWeight = @(traj) traj ./ max(sqrt(eps),vecnorm(traj,Inf,2));	% Scale-weighting
    case 5,    fWeight = @(traj) traj ./ (pi*tInstances).^0.25;		% Reciprocal square-root time-weighting
    otherwise, fWeight = @(traj) traj;
  end%switch

  % Trajectory Centering
  switch flags(1)
    case 1,    fCenter = @(traj,xs) traj - xs;					% Steady state / output
    case 2,    fCenter = @(traj,xs) traj - traj(:,end);			% Final state / output
    case 3,    fCenter = @(traj,xs) traj - mean(traj,2);			% Temporal mean of state / output
    case 4,    fCenter = @(traj,xs) traj - sqrt(mean(traj .* traj,2));	% Temporal root-mean-square of state / output
    case 5,    fCenter = @(traj,xs) traj - 0.5*(max(traj,[],2)+min(traj,[],2));% Temporal mid-range of state / output
    otherwise, fCenter = @(traj,xs) traj;
  end%switch

  % Steady State
  vSteadyInput = repmat(us,iif(isscalar(us),nInputs,1),1);
  vSteadyState = repmat(xs,iif(isscalar(xs),nStates,1),1);

  % Gramian Normalization
  if ismember(flags(6),[1,2]) && ismember(gramianType,{'c','o','x','y'})

    if 2 == flags(6)  % Jacobi-type preconditioner
        NF = nf;
        NF(6) = 0;
        if isequal(w,'c'), NF(7) = 0; end%if
        PR = mean(pr,2);
        DP = @(x,y) sum(x(1:nStates,:) .* y(:,1:nStates)',2);  % Diagonal-only pseudo-kernel
        TX = sqrt(abs(emgr(f,g,s,t,w,PR,NF,ut,us,xs,um,xm,DP)));
    else              % Steady-state preconditioner
        TX = vSteadyState;
    end%if

    TX(abs(TX)<sqrt(eps)) = 1.0;

    vSteadyState = vSteadyState ./ TX;
    fState = @(x,u,p,t) f(TX .* x,u,p,t) ./ TX;
    fAdjoint = @(x,u,p,t) g(TX .* x,u,p,t) ./ TX;
    fOutput = @(x,u,p,t) g(TX .* x,u,p,t);
  end%if

  % Output Averaging
  nPages = iif(flags(7),1,nOutputs);

  % Extra Input (for control explicit observability)
  fSteady = iif(flags(8),@(t) vSteadyInput + fExcite(t),@(t) vSteadyInput);

  % Perturbation Scales
  vInputMax = repmat(um,iif(isscalar(um),nInputs,1),1);
  vStateMax = repmat(xm,iif(isscalar(xm),nStates,1),1);
  vOutputMax = repmat(xm,iif(isscalar(xm),nOutputs,1),1);

  mInputScales = vInputMax * iif(1 == size(um,2),scales(flags(2),flags(4)),1);
  mStateScales = vStateMax * iif(1 == size(xm,2),scales(flags(3),flags(5)),1);
  mOutputScales = vOutputMax * iif(1 == size(xm,2),scales(flags(2),flags(4)),1);

  nInputScales = size(mInputScales,2);
  nStateScales = size(mStateScales,2);
  nOutputScales = size(mOutputScales,2);

  nTotalStates = size(mStateScales,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAMIAN COMPUTATION

  W = 0.0;

  switch gramianType

    % Common Layout:
    %   For each {parameter, scale, input/state/parameter component}:
    %     Perturb, simulate, (weight, center), normalize, accumulate
    %   Output and adjoint trajectories are cached to prevent recomputation
    %   Parameter gramians 's', 'i', 'j' call state gramians 'c', 'o', 'x'

%% Empirical Controllability Gramian

    case 'c'

      for k = 1:nParamSamples
        vParam = pr(:,k);
        vSteadyOutput = fOutput(vSteadyState,vSteadyInput,vParam,0);
        for c = 1:nInputScales
          for m = 1:nInputs  % (parallelizable with `parfor`)
            sPerturb = mInputScales(m,c);
            if not(0 == sPerturb)
              vUnit = sparse(m,1,sPerturb,nInputs,1);
              fInput = @(t) vSteadyInput + vUnit * fExcite(t);
              mTraj = ODE(fState,fOutput,t,vSteadyState,fInput,vParam);
              mTraj = fWeight(fCenter(mTraj,vSteadyOutput)) ./ sPerturb;
              W = W + dp(mTraj,mTraj');
            end%if
          end%for
        end%for
      end%for
      W = W * (tStep / (nInputScales * nParamSamples));

%% Empirical Observability Gramian

    case 'o'

      obsCache = zeros(nPages*nSteps,nTotalStates);
      for k = 1:nParamSamples
        vParam = pr(:,k);
        for d = 1:nStateScales
          for n = 1:nTotalStates  % (parallelizable with `parfor`)
            sPerturb = mStateScales(n,d);
            if not(0 == sPerturb)
              vUnit = sparse(n,1,sPerturb,nTotalStates,1);
              vInit = vSteadyState + vUnit(1:nStates);
              vParamInit = vParam;
              if nTotalStates > nStates
                vParamInit = vParamInit + vUnit(nStates+1:end);
              end%if
              vSteadyOutput = fOutput(vSteadyState,vSteadyInput,vParamInit,0);
              mTraj = ODE(fState,fOutput,t,vInit,fSteady,vParamInit);
              mTraj = fWeight(fCenter(mTraj,vSteadyOutput)) ./ sPerturb;
              if flags(7)
                obsCache(:,n) = sum(mTraj,1)';
              else
                obsCache(:,n) = reshape(mTraj',[],1);
              end%if
            end%if
          end%for
          W = W + dp(obsCache',obsCache);
        end%for
      end%for
      W = W * (tStep / (nStateScales * nParamSamples));

%% Empirical Cross Gramian

    case 'x'

      assert((nInputs == nOutputs) || flags(7),' emgr: non-square system!');

      colFirst = 1;            % Start partition column index
      colLast = nTotalStates;  % Final partition column index

      % Partitioned Cross Gramian
      if flags(11) > 0
        parSize = round(flags(11));   % Partition size
        parIndex = round(flags(12));  % Partition index
        colFirst = colFirst + (parIndex - 1) * parSize;
        colLast = min(colFirst + (parSize - 1),nStates);
        if colFirst > nStates
          colFirst = colFirst - (ceil(nStates / parSize) * parSize - nStates);
          colLast = min(colFirst + parSize - 1,nTotalStates);
        end%if

        if (parIndex < 1) || (colFirst > colLast) || (colFirst < 0)
          return;
        end%if
      end%if

      obsCache = zeros(nSteps*nPages,colLast-colFirst+1);
      for k = 1:nParamSamples
        vParam = pr(:,k);
        for d = 1:nStateScales
          for n = 1:(colLast-colFirst+1)  % (parallelizable with `parfor`)
            sPerturb = mStateScales(colFirst+n-1,d);
            if not(0 == sPerturb)
              vUnit = sparse(colFirst+n-1,1,sPerturb,nTotalStates,1);
              vInit = vSteadyState + vUnit(1:nStates);
              vParamInit = vParam;
              if nTotalStates > nStates
                vParamInit = vParamInit + vUnit(nStates+1:end);
              end%if
              vSteadyOutput = fOutput(vSteadyState,vSteadyInput,vParamInit,0);
              mTraj = ODE(fState,fOutput,t,vInit,fSteady,vParamInit);
              mTraj = fWeight(fCenter(mTraj,vSteadyOutput)) ./ sPerturb;
              if flags(7)
                obsCache(:,n) = sum(mTraj,1)';
              else
                obsCache(:,n) = reshape(mTraj',[],1);
              end%if
            end%if
          end%for
          for c = 1:nInputScales  % (parallelizable with `parfor`)
            for m = 1:nInputs
              sPerturb = mInputScales(m,c);
              if not(0 == sPerturb)
                vUnit = sparse(m,1,sPerturb,nInputs,1);
                fInput = @(t) vSteadyInput + vUnit * fExcite(t);
                mTraj = ODE(fState,@id,t,vSteadyState,fInput,vParam);
                mTraj = fWeight(fCenter(mTraj,vSteadyInput)) ./ sPerturb;
                nBlock = iif(flags(7),0,(m - 1) * nSteps);
                W = W + dp(mTraj,obsCache(nBlock+1:nBlock+nSteps,:));
              end%if
            end%for
          end%for
        end%for
      end%for
      W = W * (tStep / (nInputScales * nStateScales * nParamSamples));

%% Empirical Linear Cross Gramian

    case 'y'

      assert((nInputs == nOutputs) || flags(7),' emgr: non-square system!');
      assert(nInputScales == nOutputScales,' emgr: scale count mismatch!');

      adjCache = zeros(nSteps,nStates,nPages);
      for k = 1:nParamSamples
        vParam = pr(:,k);
        for c = 1:nInputScales
          for q = 1:nOutputs  % (parallelizable with `parfor`)
            sPerturb = mOutputScales(q,c);
            if not(0 == sPerturb)
              vUnit = sparse(q,1,sPerturb,nOutputs,1);
              fInput = @(t) vSteadyInput + vUnit * fExcite(t);
              mTraj = ODE(fAdjoint,@id,t,vSteadyState,fInput,vParam);
              mTraj = fWeight(fCenter(mTraj,vSteadyInput)) ./ sPerturb;
              adjCache(:,:,q) = mTraj';
            end%if
          end%for
          if flags(7)
            adjCache(:,:,1) = sum(adjCache,3);
          end%if
          for m = 1:nInputs  % (parallelizable with `parfor`)
            sPerturb = mInputScales(m,c);
            if not(0 == sPerturb)
              vUnit = sparse(m,1,sPerturb,nInputs,1);
              fInput = @(t) vSteadyInput + vUnit * fExcite(t);
              mTraj = ODE(fState,@id,t,vSteadyState,fInput,vParam);
              mTraj = fWeight(fCenter(mTraj,vSteadyInput)) ./ sPerturb;
              W = W + dp(mTraj,adjCache(:,:,iif(flags(7),1,m)));
            end%if
          end%for
        end%for
      end%for
      W = W * (tStep / (nInputScales * nParamSamples));

%% Empirical Sensitivity Gramian

    case 's'

      % Controllability Gramian
      [pr,mParamScales] = paramScales(pr,flags(9),nInputScales);
      WC = emgr(f,g,s,t,'c',pr,flags,ut,us,xs,um,xm,dp);

      if not(flags(10))  % Input-state sensitivity gramian
        DP = @(x,y) sum(sum(x .* y'));       % Trace pseudo-kernel
      else               % Input-output sensitivity gramian
        DP = @(x,y) y;                       % Custom pseudo-kernel
        flags(7) = 1;
        Y = emgr(f,g,s,t,'o',pr,flags,ut,us,xs,um,xm,DP);
        flags(7) = 0;
        DP = @(x,y) abs(sum(y(:) .* Y(:)));  % Custom pseudo-kernel
      end%if

      % (Diagonal) Sensitivity Gramian
      WS = zeros(nParams,1);
      for p = 1:nParams  % (parallelizable with `parfor`)
        paramSamples = repmat(pr,[1,size(mParamScales,2)]);
        paramSamples(p,:) = paramSamples(p,:) + mParamScales(p,:);
        WS(p) = emgr(f,g,s,t,'c',paramSamples,flags,ut,us,xs,um,xm,DP);
      end%for

      W = {WC,WS};

%% Empirical Identifiability Gramian

    case 'i'

      % Augmented Observability Gramian
      [pr,mParamScales] = paramScales(pr,flags(9),nStateScales);
      V = emgr(f,g,s,t,'o',pr,flags,ut,us,xs,um,[mStateScales;mParamScales],dp);

      % Return Augmented Observability Gramian
      if flags(11), W = V; return; end%if

      WO = V(1:nStates, 1:nStates);          % Observability gramian
      WM = V(1:nStates, nStates+1:end);      % Mixed block
      WI = V(nStates+1:end, nStates+1:end);  % Parameter gramian

      % Identifiability Gramian via Schur Complement
      switch flags(10)
        case 0,    WI = WI - (WM' * ainv(WO) * WM);
        case 2,    WI = WI - (WM' * pinv(WO) * WM);
      end%switch

      W = {WO,WI};

%% Empirical Joint Gramian

    case 'j'

      % Joint Gramian
      [pr,mParamScales] = paramScales(pr,flags(9),nStateScales);
      V = emgr(f,g,s,t,'x',pr,flags,ut,us,xs,um,[mStateScales;mParamScales],dp);

      % Return Joint Gramian (Partition)
      if flags(11), W = V; return; end%if

      WX = V(1:nStates, 1:nStates);      % Cross gramian
      WM = V(1:nStates, nStates+1:end);  % Mixed block

      % Cross-Identifiability Gramian via Schur Complement
      switch flags(10)
        case 1,    WI = 0.5 * (WM' * WM);
        case 2,    WI = 0.5 * (WM' * pinv(WX + WX') * WM);
        otherwise, WI = 0.5 * (WM' * ainv(WX + WX') * WM);
      end%switch

      W = {WX,WI};

    otherwise

      error(' emgr: unknown empirical gramian type!');
  end%switch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: iif

function r = iif(pre,con,alt)
% summary: inline if

  if pre, r = con; else, r = alt; end%if
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: id

function x = id(x,u,p,t)
% summary: (Output) identity functional

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: stateScales

function mScales = scales(flScales,flRot)
% summary: Input and initial state perturbation scales

  switch flScales
    case 1,    mScales = [0.25, 0.50, 0.75, 1.0];  % Linear
    case 2,    mScales = [0.125, 0.25, 0.5, 1.0];  % Geometric
    case 3,    mScales = [0.001, 0.01, 0.1, 1.0];  % Logarithmic
    case 4,    mScales = [0.01, 0.50, 0.99, 1.0];  % Sparse
    otherwise, mScales = 1.0;                      % One
  end%switch

  if 0 == flRot, mScales = [-mScales,mScales]; end%if
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: paramScales

function [vParamSteady,mParamScales] = paramScales(p,flScales,nParamScales)
% summary: Parameter perturbation scales

  [vParamMin,vParamMax] = bounds(p,2);

  switch flScales
    case 1     % Linear centering and scaling
      assert(size(p,2) >= 2,' emgr: min and max parameter required!');
      vParamSteady = 0.5 * (vParamMax + vParamMin);
      vScales = linspace(0.0,1.0,nParamScales);
    case 2     % Logarithmic centering and scaling
      assert(size(p,2) >= 2,' emgr: min and max parameter required!');
      vParamSteady = sqrt(vParamMax .* vParamMin);
      vParamMin = log(vParamMin);
      vParamMax = log(vParamMax);
      vScales = linspace(0.0,1.0,nParamScales);
    case 3     % Nominal centering and scaling
      assert(size(p,2) == 3,' emgr: min, nom, max parameter required!');
      vParamSteady = p(:,2);
      vParamMin = p(:,1);
      vParamMax = p(:,3);
      vScales = linspace(0.0,1.0,nParamScales);
    otherwise  % No centering and linear scaling
      assert(size(p,2) >= 2,' emgr: min and max parameter required!');
      vParamSteady = vParamMin;
      vParamMin = ones(size(p,1),1)./nParamScales;
      vScales = linspace(1.0/nParamScales,1.0,nParamScales);
  end%switch

  mParamScales = (vParamMax - vParamMin) * vScales + vParamMin;
  if 2 == flScales, mParamScales = exp(mParamScales); end%if
  mParamScales = mParamScales - vParamSteady;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: ainv

function x = ainv(m)
% summary: Quadratic complexity approximate inverse matrix

  % Based on truncated Neumann series: X = D^-1 - (D^-1 (M - D) D^-1)
  D = diag(m);
  k = find(abs(D) > sqrt(eps));
  D(k) = 1.0 ./ D(k);
  x = (m .* (-D)) .* D';
  x(1:numel(D) + 1:end) = D;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: ssp2

function y = ssp2(f,g,t,x0,u,p)
% summary: Low-Storage Strong-Stability-Preserving Second-Order Runge-Kutta

  global STAGES;  % Configurable number of stages for enhanced stability
  if not(isscalar(STAGES)), nStages = 3; else, nStages = STAGES; end%if

  tStep = t(1);
  nSteps = floor(t(2) / tStep) + 1;
  y = g(x0,u(0),p,0);
  y(:,nSteps) = 0.0;  % Pre-allocate trajectory

  xk1 = x0;
  for k = 2:nSteps
    xk2 = xk1;
    tCurr = (k - 1.5) * tStep;
    uCurr = u(tCurr);
    for s = 2:nStages
      xk1 = xk1 + (tStep / (nStages - 1)) * f(xk1,uCurr,p,tCurr);
    end%for
    xk1 = (xk1 * (nStages - 1) + xk2 + tStep * f(xk1,uCurr,p,tCurr)) / nStages;
    y(:,k) = g(xk1,uCurr,p,tCurr);
  end%for
end

