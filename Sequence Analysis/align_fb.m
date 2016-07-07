function [mod_inds, mod_type, lvl_accum, P, ks] = align_fb(model_prediction, model_stds, lvls, dts, dboff, probs)

    if nargin < 5
        dboff = 100;
    end
    
    if nargin < 6
        % so first we need observation (emission) probabilities
        % assuming constant stddev for each level
        % obsstd = 2.0; % pA
        stayprob = 0.001;
        fwdprob = 1;
        backprob = 0.0001;
        skipprob = 0.01;

        % other transition probs
        noiseprob = 0.001;
        deepprob = 0.001;
        deepdeep = 0.01;
    else
        obsstd = probs(1);
        stayprob = probs(2);
        fwdprob = probs(3);
        backprob = probs(4);
        skipprob = probs(5);
        noiseprob = probs(6);
        deepprob = probs(7);
        deepdeep = probs(8);
    end

    if nargin < 3
        dts = 0*lvls + 0.1;
    end
    
    % rows of our matrix is model_prediction
    % columns are all of our levels
    % then we also have one for 'spurious noise' states
    % and another one for 'deep blockage' states
    M = numel(model_prediction);
    N = numel(lvls);
    
    [Nmat,Mmat] = meshgrid(lvls,model_prediction);
    stds = repmat(model_stds',numel(lvls),1)'; % orientation of this matrix has been checked and is correct
    
    %E_model = normpdf(Mmat,Nmat,obsstd);
    %E_model = normpdf(Mmat,Nmat,min(2,stds));% + 0.001*lorentzianpdf(Mmat,Nmat,stds/2).^10;
    E_model = normpdf(Mmat,Nmat,0.5+stds);
    %E_model = lorentzianpdf(Mmat,Nmat,stds/2).^10;
    % use exponential constant for prob of being random noise
    E_noise = exppdf(repmat(dts',[M,1]),0.0001);%1+0*E_model;
    %E_noise = exppdf(repmat(dts',[M,1]),0.001).*E_model;
    %E_noise = exppdf(repmat(dts',[M,1]),median(dts)/log(2)/200);
    % and add dboff pA for deep blockages
    %E_deep = normpdf(Mmat,Nmat+dboff,obsstd);
    %E_deep = ( normpdf(Mmat,Nmat+dboff,min(2,stds))) + 1e-5*exppdf(Nmat,5) )/(1+1e-5);
    E_deep = lorentzianpdf(Mmat,Nmat+dboff,stds/2).^8;%.*E_model;
    
    % put them all into one emission matrix
    E = [E_model;E_noise;E_deep];
    E = E ./ repmat(sum(E),[size(E,1) 1]);
    % do log-likelihoods
    E = log(E);
    
    % and build a transition matrix
    % ROW IS START STATE, COLUMN IS END STATE
    dn = 5;
    
    dprobs = [backprob*skipprob.^(dn-1:-1:0) stayprob fwdprob*skipprob.^(0:dn-1)];
    %dprobs = [(backprob*skipprob).^(5*((dn-1):-1:0)) stayprob (fwdprob*skipprob).^(5*(0:dn-1))];
    %dprobs = [backprob*backprob.^(10*(dn-1):-10:0) stayprob fwdprob*skipprob.^(0:10:10*(dn-1))];
    T0 = zeros(M);
    for i=1:numel(dprobs)
        di = i-dn-1;
        T0 = T0 + diag(dprobs(i)*ones(M-abs(di),1),di);
    end
    
    TI = eye(M);
    
    % now build up a big 3x3 block transition matrix
    % [ 0->0,  0->NOISE,   0->DEEP;
    %   N->0,      N->N,   N->D;
    %   D->0,      D->N,   D->D]
    
    T = [       T0,         noiseprob*TI,           deepprob*T0;
                T0,       noiseprob^2*TI,           deepprob*TI;
                T0,         noiseprob*TI,           deepdeep*TI];
    
    % normalize each row to sum to 1
    Tn = sum(T,2);
    T = T ./ repmat(Tn,[1,3*M]);
    
    T = log(T+1e-300);
    
    % now run the viterbi iteration stuff
    % how do we set initial probabilities? dunno, exponentially decaying
    % at the start or something
    P = 0*E;
    P(:,1) = E(:,1) - (2*(1:3*M)');
    P(:,1) = P(:,1) - log(sum(exp(P(:,1))));
    Pnorms = zeros(N,1);
    K = 0*E;
    
    for i=2:N
        % rows are previous state
        % columns are current state (with emissions)
        p = P(:,i-1);
        e = E(:,i);
        % outer sum matrix
        A = repmat(e',size(p))+repmat(p,size(e'));
        % add (multiply) by transition matrix
        A = A + T;
        [p_new,pi] = max(A);
        pn = log(sum(exp(p_new)));
        Pnorms(i) = pn;
        p_new = p_new - pn;
        P(:,i) = p_new;
        K(:,i) = pi;
    end
    
    % run backtrace
    ks = zeros(N,1);
    [~,k] = max(P(:,end));
    for i=N:-1:1
        ks(i) = k;
        k = K(k,i);
    end
    
%     % run forward trace
%     ks = zeros(N,1);
%     [~,k] = max(P(:,1));
%     for i=1:N
%         ks(i) = k;
%         k = K(k,i);
%     end
    
    % so now ks is the state index of each seen level
    % kinds is the model index (which is state index modulo number of mod
    % levels)
    kinds = mod(ks-1,M)+1;
    mod_inds = kinds;
    % and model_type is what type of level it is (1=normal, 2=noise,
    % 3=deep)
    mod_type = 1 + floor((ks-1)/M);
    
    % create adjusted levels for deep blockages
%     lvls_adj = lvls;
%     lvls_adj(mod_type==3) = lvls_adj(mod_type==3)+dboff;

%     lvl_accum = accumarray(kinds(mod_type~=2),lvls_adj(mod_type~=2),[M 1],@mean,nan);
    lvl_accum = accumarray(kinds(mod_type==1),lvls(mod_type==1),[M 1],@mean,nan);
    
    if nargout >= 1
        return
    end

    modlvls = model_prediction(kinds);
    
    figure(5)
    clf(5)
    
    subplot(3,1,1);
    plot(model_prediction,'o-','LineWidth',2)
    hold on
    plot(lvl_accum,'o-','LineWidth',2);
    legend('Model','Data')
    ylabel('Current (pA)')
    xlabel('Model level')
    set(gca,'FontSize',20)
    title('Best fit of data to model')
    
    subplot(3,1,2);
    modlvls(mod_type==3) = modlvls(mod_type==3) - dboff;
    plot(modlvls)
    for i=1:N
        text(i,lvls(i),num2str(kinds(i)),'FontSize',14);
    end
    hold on
    xx = 1:numel(modlvls);
    plot(xx(mod_type==2),lvls(mod_type==2),'rx','MarkerSize',10)
    plot(lvls,'-');
    ylabel('Current (pA)')
    xlabel('Measured level')
    set(gca,'FontSize',20)
    title('Matching each measured level to a model state')
    
    subplot(3,1,3);
    imagesc(1.03.^(P/2)) % scales so image shows up well
    hold on
    plot(1:N,ks,'r','LineWidth',3);
    title('Probability Matrix')
    ylabel('State')
    xlabel('Level Number')
    set(gca,'FontSize',20)
    %axis equal
    
end