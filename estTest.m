function estTest(t)
%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% version: 5.9 (2021-01-21)
%%% authors: Christian Himpe (0000-0003-2194-6754)
%%% license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%%% summary: estTest - test est functionality

    rand('seed',1009);								% Seed uniform random number generator
    randn('seed',1009');

    M = 4;
    N = 16;

    A = full(-gallery('tridiag',N));
    A(1,1) = 1;
    B = double(1:N==1)'*linspace(1/M,1,M);
    C = B';

    sys = struct('M',M, ...
                 'N',N, ...
                 'Q',M, ...
                 'f',@(x,u,p,t) A*x + B*u + p, ...
                 'g',@(x,u,p,t) C*x, ...
                 'F',@(x,u,p,t) (x'*A + u'*C)', ...
                 'p',zeros(N,1), ...
                 'dt',0.01, ...
                 'Tf',1.0);

    switch lower(t)

        case 'matrix_equation', matrix_equation(sys);

        case 'singular_values', singular_values(sys);

        case 'model_reduction', model_reduction(sys);

        case 'parameter_reduction', parameter_reduction(sys);

        case 'combined_reduction', combined_reduction(sys);

        case 'decentralized_control', decentralized_control(sys);

        case 'state_sensitivity', state_sensitivity(sys);

        case 'parameter_sensitivity', parameter_sensitivity(sys);

        case 'parameter_identifiability', parameter_identifiability(sys);

        case 'uncertainty_quantification', uncertainty_quantification(sys);

        case 'nonlinearity_quantification', nonlinearity_quantification(sys);

        case 'gramian_index', gramian_index(sys);

        case 'system_index', system_index(sys);

        case 'system_norm', system_norm(sys);

        case 'tau_function', tau_function(sys);

        otherwise error('Unknown test!');

    end%switch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST MATRIX EQUATION

function matrix_equation(sys)

    task = struct('type','matrix_equation');

    method = {'lyapunov', ...
              'sylvester'};

    config = struct('linearity','linear');

    disp([char(10),'Task: ',task.type,char(10)]);

    er = zeros(numel(method),1);

    for k = 1:numel(method)

        name{k} = [method{k}];

        disp(['Testing: ',name{k}]);

        r = est(sys,setfield(task,'method',method{k}),config);

        s = est(setfield(sys,'dt',sys.dt * 0.1),setfield(task,'method',method{k}),config);

        t = est(setfield(sys,'dt',sys.dt * 0.01),setfield(task,'method',method{k}),config);

        err{k} = [norm(r - t,'Fro'),norm(s - t,'Fro')];

        EOC(k,1) = log10(norm(r - t,'Fro') / norm(s - t,'Fro'));
    end%for

    EOC

    figure('Name',task.type,'NumberTitle','off');
    set(gca,'YScale','log','YGrid','on','XLim',[1,2],'NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c,'LineWidth',3),err);
    legend(name(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST SINGULAR VALUES

function singular_values(sys)

    task = struct('type','singular_values');

    method = {'controllability', ...
              'observability', ...
              'minimality'};

    linearity = {'linear','nonlinear'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    sv = cell(numel(method),numel(linearity));
    ms = zeros(numel(method),numel(linearity));

    for k = 1:numel(method)
        for l = 1:numel(linearity)

            name{k,l} = [method{k},' (',linearity{l},')'];

            disp(['Testing: ',name{k,l}]);

            [r,s] = est(sys,setfield(task,'method',method{k}), ...
                            setfield(config,'linearity',linearity{l}));
            sv{k,l} = r;
            ms(k,l) = s;
        end%for
    end%for

    MORscore = ms

    figure('Name',task.type,'NumberTitle','off');
    set(gca,'YScale','log','YGrid','on','XLim',[1,sys.N],'NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c./max(c),'LineWidth',3),sv);
    legend(name(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST MODEL REDUCTION

function model_reduction(sys)

    task = struct('type','model_reduction');

    method = {'poor_man', ...
              'dominant_subspaces', ...
              'approx_balancing', ...
              'balanced_pod', ...
              'balanced_truncation' ...
              'dmd_galerkin'};

    variant = {'observability', ...
               'minimality'};

    linearity = {'linear','nonlinear'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    sv = cell(numel(method),numel(variant),numel(linearity));
    ms = zeros(numel(method),numel(variant),numel(linearity));

    for k = 1:numel(method)
        for l = 1:numel(variant)
            for m = 1:numel(linearity)

                name{k,l,m} = [method{k},' (',variant{l},',',linearity{m},')'];

                disp(['Testing: ',name{k,l,m}]);

                [r,s] = est(sys,setfield(setfield(task,'method',method{k}),'variant',variant{l}), ...
                                setfield(setfield(config,'linearity',linearity{m}),'test',true));

                mr{k,l,m} = r;
                ms(k,l,m) = s;
            end%for
        end%for
    end%for 

    MORscore_PM = squeeze(ms(1,:,:))
    MORscore_DS = squeeze(ms(2,:,:))
    MORscore_AB = squeeze(ms(3,:,:))
    MORscore_BT = squeeze(ms(4,:,:))
    MORscore_DG = squeeze(ms(5,:,:))

    figure('Name',task.type,'NumberTitle','off');
    yl = {'L_1','L_2','L_\infty','L_0'};
    for k = 1:numel(method)
        for l = 1:4
            subplot(4,numel(method),(l-1)*numel(method)+k);
            set(gca,'YScale','log','Ytick',[1e-16,1e-8,1e-0],'YGrid','on','XLim',[1,sys.N-1],'YLim',[1e-16,1.1],'NextPlot','add','ColorOrder',lines12);
            cellfun(@(c) plot(c{l+2},'LineWidth',3),mr(k,:,:));
            if isequal(k,1), ylabel(yl{l}); end%if
            if isequal(l,1), title(method{k},'Interpreter','none'); end%if
        end%for
    end%for
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST PARAMETER REDUCTION

function parameter_reduction(sys)

    task = struct('type','parameter_reduction');

    method = {'observability', ...
              'minimality'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    sv = cell(1,numel(method));
    ms = zeros(1,numel(method));

    for k = 1:numel(method)

        disp(['Testing: ',method{k}]);

        [r,s] = est(setfield(sys,'p',repmat([0,1],sys.N,1)), ...
                    setfield(task,'method',method{k}), ...
                    setfield(setfield(config,'test',true),'num_test_param',5));

        pr{k} = r;
        ms(k) = s;
    end%for

    MORscore = ms

    figure('Name',task.type,'NumberTitle','off');
    yl = {'L_1','L_2','L_\infty','L_0'};
    for l = 1:4
        subplot(4,1,l);
        set(gca,'YScale','log','Ytick',[1e-16,1e-8,1],'YGrid','on','XLim',[1,sys.N-1],'YLim',[1e-16,1.1],'NextPlot','add','ColorOrder',lines12);
        for k = 1:numel(method)
            plot(pr{k}{l+2},'LineWidth',3);
        end%for
        ylabel(yl{l});
    end%for
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST COMBINED REDUCTION

function combined_reduction(sys)

    task = struct('type','combined_reduction');

    variant = {'poor_man', ...
               'dominant_subspaces', ...
               'approx_balancing', ...
               'balanced_pod', ...
               'balanced_truncation'};

    method = {'observability', ...
              'minimality'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)
        for l = 1:numel(variant)

            disp(['Testing: ',method{k},' (',variant{l},')']);

            [r,s] = est(setfield(sys,'p',repmat([0,1],sys.N,1)), ...
                        setfield(setfield(task,'method',method{k}),'variant',variant{l}), ...
                        setfield(setfield(config,'test',true),'num_test_param',5));

            cr{k,l} = r;
            ms(k,l) = s;
        end%for
    end%for

    MORscore = ms'

    figure('Name',task.type,'NumberTitle','off');
    for m = 1:numel(cr)
        subplot(numel(variant),numel(method),m);
        h = surf(cr{m}{1},cr{m}{2},cr{m}{5});
        xlabel('x');
        ylabel('p');
        set(gca,'ZScale','log','CLim',log10(get(gca,'ZLim')));
        set(h,'CData',log10(get(h,'CData')));
        view(135,15);
    end%for
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST DECENTRALIZED CONTROL

function decentralized_control(sys)

    task = struct('type','decentralized_control');

    method = {'relative_gain_array', ...
              'io_pairing', ...
              'hardy_inf', ...
              'participation_matrix', ...
              'hankel_interaction', ...
              'rms_hsv', ...
              'hardy_2', ...
              'io_coherence'};

    linearity = {'linear','nonlinear'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)
        for l = 1:numel(linearity)

            disp(['Testing: ',method{k},' (',linearity{l},')']);

            dc{l,k} = est(sys,setfield(task,'method',method{k}), ...
                              setfield(config,'linearity',linearity{l}));
        end%for
    end%for

    figure('Name',task.type,'NumberTitle','off');
    for m = 1:numel(dc)
        subplot(numel(method),numel(linearity),m);
        imagesc(dc{m});
        set(gca,'Ytick',[],'Xtick',[],'Ztick',[]);
    end%for
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST STATE SENSITIVITY

function state_sensitivity(sys)

    task = struct('type','state_sensitivity');

    method = {'controllability', ...
              'observability', ...
              'minimality'};

    linearity = {'linear','nonlinear'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)
        for l = 1:numel(linearity)

            disp(['Testing: ',method{k},' (',linearity{l},')']);

            ds{l,k} = est(sys,setfield(task,'method',method{k}), ...
                              setfield(config,'linearity',linearity{l}));
        end%for
    end%for

    figure('Name',task.type,'NumberTitle','off');
    for k = 1:numel(linearity)
        subplot(1,numel(linearity),k);
        set(gca,'XLim',[1,sys.N],'NextPlot','add','ColorOrder',lines12);
        cellfun(@(c) plot(c,'LineWidth',3), ds(k,:));
    end%for
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST PARAMETER SENSITIVITY

function parameter_sensitivity(sys)

    task = struct('type','parameter_sensitivity');

    method = {'controllability', ...
              'observability', ...
              'minimality'};

    linearity = {'linear','nonlinear'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)

        disp(['Testing: ',method{k}]);

        ps{k} = est(setfield(sys,'p',repmat([0,1],sys.N,1)), ...
                    setfield(task,'method',method{k}), ...
                    config);
    end%for

    figure('Name',task.type,'NumberTitle','off');
    set(gca,'YScale','log','XLim',[1,sys.N],'NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c,'LineWidth',3),ps);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST PARAMETER IDENTIFIABILITY

function parameter_identifiability(sys)

    task = struct('type','parameter_identifiability');

    method = {'observability', ...
              'minimality'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)

        disp(['Testing: ',method{k}]);

        pj{k} = est(setfield(sys,'p',repmat([0,1],sys.N,1)), ...
                    setfield(task,'method',method{k}), ...
                    config);
    end%for

    figure('Name',task.type,'NumberTitle','off');
    set(gca,'YScale','log','XLim',[1,sys.N],'NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c./max(c),'LineWidth',3),pj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST UNCERTAINTY QUANTIFICATION

function uncertainty_quantification(sys)

    task = struct('type','uncertainty_quantification');

    method = {'controllability', ...
              'observability'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)

        disp(['Testing: ',method{k}]);

        uq{k} = est(sys,setfield(task,'method',method{k}),config);
    end%for

    figure('Name',task.type,'NumberTitle','off');
    subplot(1,2,1);
    set(gca,'XLim',[1,sys.N],'YScale','log','NextPlot','add','ColorOrder',lines12);
    plot(uq{1}./max(uq{1}),'LineWidth',3);
    subplot(1,2,2);
    set(gca,'XLim',[1,sys.Q],'YScale','log','NextPlot','add','ColorOrder',lines12);
    plot(uq{2}./max(uq{2}),'LineWidth',3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST NONLINEARITY QUANTIFICATION

function nonlinearity_quantification(sys)

    task = struct('type','nonlinearity_quantification');

    method = {'controllability', ...
              'observability', ...
              'minimality', ...
              'correlation'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)

        disp(['Testing: ',method{k}]);

        nq{k} = est(setfield(sys,'f',@(x,u,p,t) sys.f(atan(x),u,p,t)), ...
                    setfield(task,'method',method{k}), ...
                    config);
    end%for

    figure('Name',task.type,'NumberTitle','off');
    set(gca,'XLim',[1,10],'YScale','log','NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c,'LineWidth',3),nq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST GRAMIAN INDICES

function gramian_index(sys)

    task = struct('type','gramian_index');

    method = {'sigma_min', ...
              'harmonic_mean', ...
              'geometric_mean', ...
              'energy_fraction', ...
              'operator_norm', ...
              'sigma_max', ...
              'log_det', ...
              'storage_efficiency', ...
              'unobservability_index', ...
              'performance_index'};

    variant = {'controllability', ...
               'observability', ...
               'minimality'};

    config = struct('linearity','linear');

    sys.proj = est(sys, ...
                   struct('type','model_reduction', ...
                          'method','dominant_subspaces', ...
                          'variant','minimality'), ...
                   config);

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)
        for l = 1:numel(variant)

            disp(['Testing: ',method{k},'(',variant{l},')']);

            gi{k,l} = est(sys,setfield(setfield(task,'method',method{k}),'variant',variant{l}),config);
        end%for
    end%for

    figure('Name',task.type,'NumberTitle','off');
    for l = 1:numel(variant)

        subplot(1,numel(variant),l);
        set(gca,'XLim',[1,sys.N],'YScale','log','NextPlot','add','ColorOrder',lines12);
        cellfun(@(c) plot(c./max(c),'LineWidth',3),gi(:,l));
        title(variant{l},'Interpreter','None');
        legend(method,'Location','SouthOutside','Interpreter','None');
    end%for
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST SYSTEM INDICES

function system_index(sys)

    task = struct('type','system_index');

    method = {'cauchy_index', ...
...
              'system_entropy', ...
              'gramian_distance', ...
...
              'system_symmetry', ...
              'io_coherence', ...
              'system_gain', ...
              'geometric_mean_hsv', ...
              'network_sensitivity', ...
              'rv_coefficient'};

    config = struct('linearity','linear');

    sys.proj = est(sys, ...
                   struct('type','model_reduction', ...
                          'method','dominant_subspaces', ...
                          'variant','minimality'), ...
                   config);

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)

        disp(['Testing: ',method{k}]);

        si{k} = est(sys,setfield(task,'method',method{k}),config);
    end%for

    figure('Name',task.type,'NumberTitle','off');

    subplot(1,3,1);
    set(gca,'XLim',[1,sys.N],'NextPlot','add','ColorOrder',lines12);
    plot(si{1},'LineWidth',3);
    title('Discrete Indicators','Interpreter','None');
    legend(method(1),'Location','SouthOutside','Interpreter','None');

    subplot(1,3,2);
    set(gca,'XLim',[1,sys.N],'NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c./max(c),'LineWidth',3),si(2:3));
    title('Linear Indicators','Interpreter','None');
    legend(method(2:3),'Location','SouthOutside','Interpreter','None');

    subplot(1,3,3);
    set(gca,'XLim',[1,sys.N],'YScale','log','NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c./max(c),'LineWidth',3),si(4:end));
    title('Logarithmic Indicators','Interpreter','None');
    legend(method(4:end),'Location','SouthOutside','Interpreter','None');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST SYSTEM NORM

function system_norm(sys)

    task = struct('type','system_norm');

    method = {'hardy_2_norm', ...
              'hardy_inf_norm', ...
              'hilbert_schmidt_hankel_norm', ...
              'hankel_norm', ...
             };

    config = struct('linearity','linear');

    sys.proj = est(sys, ...
                   struct('type','model_reduction', ...
                          'method','dominant_subspaces', ...
                          'variant','minimality'), ...
                   config);

    disp([char(10),'Task: ',task.type,char(10)]);

    for k = 1:numel(method)

        disp(['Testing: ',method{k}]);

        sn{k} = est(sys,setfield(task,'method',method{k}),config);
    end%for

    figure('Name',task.type,'NumberTitle','off');
    set(gca,'XLim',[1,sys.N],'YScale','log','NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c./max(c),'LineWidth',3),sn);
    title('Norms','Interpreter','None');
    legend(method,'Location','SouthOutside','Interpreter','None');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST TAU FUNCTION

function tau_function(sys)

    task = struct('type','tau_function');

    linearity = {'linear','nonlinear'};

    config = struct;

    disp([char(10),'Task: ',task.type,char(10)]);

    for l = 1:numel(linearity)

        name{l} = ['tau (',linearity{l},')'];

        disp(['Testing: ',name{l}]);

        au{l} = est(sys,task,setfield(config,'linearity',linearity{l}));
    end%for

    figure('Name',task.type,'NumberTitle','off');
    set(gca,'YGrid','on','NextPlot','add','ColorOrder',lines12);
    cellfun(@(c) plot(c,'LineWidth',3),au);
    legend(name(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTION

function r = lines12(n)

    if nargin<1 || isempty(n), n = 12; end%if

    x = [166,206,227; ...
         31,120,180;  ...
         178,223,138; ...
         51,160,44;   ...
         251,154,153; ...
         227,26,28;   ...
         253,191,111; ...
         255,127,0;   ...
         202,178,214; ...
         106,61,154; ...
         255,255,153; ...
         177,89,40]./256;

    r = x(rem(0:(n-1),size(x,1))+1,:);
end

