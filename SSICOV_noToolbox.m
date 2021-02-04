function [fn,zeta,phi,varargout] = SSICOV_noToolbox(y,dt,varargin)
%
% -------------------------------------------------------------------------
% [fn,zeta,phi,varargout] = SSICOV_noToolbox(y,dt,varargin) identifies the modal
% parameters of the M-DOF system whose response histories are located in
% the matrix y, sampled with a time step dt.
% -------------------------------------------------------------------------
% Input:
% y: time series of ambient vibrations: matrix of size [MxN]
% dt : scalar: Time step
% Varargin: contains additional optaional parameters:
%	'Ts': scalar : time lag for covariance calculation
%	'methodCOV': scalar: method for COV estimate ( 1 or 2)
%	'Nmin': scalar: minimal number of model order
%	'Nmax': scalar: maximal number of model order
%	'eps_freq': scalar: frequency accuracy
%	'eps_zeta': scalar: % damping accuracy
%	'eps_MAC': scalar: % MAC accuracy
%	'eps_cluster': scalar: % maximal distance inside each cluster
% -------------------------------------------------------------------------
% Output:
% fn: eigen frequencies identified
% zeta:  modal damping ratio identified
% phi:mode shape identified
% varargout: structure data useful for stabilization diagram
% -------------------------------------------------------------------------
%  Syntax:
% [fn,zeta,phi] = SSICOV_noToolbox(y,dt,'Ts',30) specifies that the time lag
% has to be 30 seconds.
%
% [fn,zeta,phi] = SSICOV_noToolbox(y,dt,'Ts',30,'Nmin',5,'Nmax',40) specifies that the
% time lag has to be 30 seconds, with a system order ranging from 5 to 40.
%
% [fn,zeta,phi] = SSICOV_noToolbox(y,dt,'eps_cluster',0.05) specifies that the
% max distance inside each cluster is 0.05 hz.
%
% [fn,zeta,phi] = SSICOV_noToolbox(y,dt,'eps_freq',1e-2,'eps_MAC'.1e-2) changes the
% default accuracy for the stability checking procedure
%
% -------------------------------------------------------------------------
% Organization of the function:
% 6 steps:
% 1 - Claculation of cross-correlation function
% 2 - Construction of the block Toeplitz matrix and SVD of it
% 3 - Modal identification procedure
% 4 - Stability checking procedure
% 5 - Selection of stable poles only
% 6 - Cluster Algorithm
% -------------------------------------------------------------------------
% References:
% Magalhaes, F., Cunha, A., & Caetano, E. (2009).
% Online automatic identification of the modal parameters of a long span arch
% bridge. Mechanical Systems and Signal Processing, 23(2), 316-329.
%
% Magalhães, F., Cunha, Á., & Caetano, E. (2008).
% Dynamic monitoring of a long span arch bridge. Engineering Structures,
% 30(11), 3034-3044.
% -------------------------------------------------------------------------
% Author: E Cheynet, UiS/UiB - Norway
% Last modified: 06/12/2019
% -------------------------------------------------------------------------
%
% see also plotStabDiag.m

%%
% options: default values
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Ts',500*dt);
p.addOptional('methodCOV',1);
p.addOptional('Nmin',2);
p.addOptional('Nmax',30);
p.addOptional('eps_freq',1e-2);
p.addOptional('eps_zeta',4e-2);
p.addOptional('eps_MAC',5e-3);
p.addOptional('eps_cluster',0.2);
p.parse(varargin{:});

% Number of outputs must be >=3 and <=4.
nargoutchk(3,4)
% size of the input y
[Nyy,N]= size(y);

% shorthen the variables name
eps_freq = p.Results.eps_freq ;
eps_zeta = p.Results.eps_zeta ;
eps_MAC = p.Results.eps_MAC ;
eps_cluster = p.Results.eps_cluster ;
Nmin = p.Results.Nmin ;
Nmax = p.Results.Nmax ;

%  Natural Excitation Technique (NeXT)
[IRF,~] = NExT(y,dt,p.Results.Ts,p.Results.methodCOV);
% Block Hankel computations
[U,S,~] = blockToeplitz(IRF);
if isnan(U)
    fn = nan;
    zeta = nan;
    phi = nan;
    if nargout==4
        varargout = {nan};
    end
    return
end
% Stability check
kk=1;
for ii=Nmax:-1:Nmin % decreasing order of poles
    if kk==1
        [fn0,zeta0,phi0] = modalID(U,S,ii,Nyy,dt);
    else
        [fn1,zeta1,phi1] = modalID(U,S,ii,Nyy,dt);
        [a,b,c,d,e] = stabilityCheck(fn0,zeta0,phi0,fn1,zeta1,phi1);
        fn2{kk-1}=a;
        zeta2{kk-1}=b;
        phi2{kk-1}=c;
        MAC{kk-1}=d;
        stablity_status{kk-1}=e;
        fn0=fn1;
        zeta0=zeta1;
        phi0=phi1;
    end
    kk=kk+1;
end

% sort for increasing order of poles
stablity_status=fliplr(stablity_status);
fn2=fliplr(fn2);
zeta2=fliplr(zeta2);
phi2=fliplr(phi2);
MAC=fliplr(MAC);

% get only stable poles
[fnS,zetaS,phiS,MACS] = getStablePoles(fn2,zeta2,phi2,MAC,stablity_status);

if isempty(fnS)
    warning('No stable poles found');
    fn = nan;
    zeta = nan;
    phi = nan;
    if nargout==4
        varargout = {nan};
    end
    return
end

% Hierarchical cluster
[fn3,zeta3,phi3] = myClusterFun(fnS,zetaS,phiS);
if isnumeric(fn3)
    warning('Hierarchical cluster failed to find any cluster');
    fn = nan;
    zeta = nan;
    phi = nan;
    if nargout==4
        varargout = {nan};
    end
    return
end
% average the clusters to get the frequency and mode shapes


% Up to Nmax parameters are identified
fn = zeros(1,Nmax);
zeta = zeros(1,Nmax);
phi = zeros(Nmax,Nyy);
for ii=1:numel(fn3)
    fn(ii)=nanmean(fn3{ii});
    zeta(ii)=nanmean(zeta3{ii});
    phi(ii,:)=nanmean(phi3{ii},2);
end

phi(fn==0,:)=[];
zeta(fn==0)=[];
fn(fn==0)=[];

% sort the eigen frequencies
[fn,indSort]=sort(fn);
zeta = zeta(indSort);
phi = phi(indSort,:);

% varargout for stabilization diagram
if nargout==4
    paraPlot.status=stablity_status;
    paraPlot.Nmin = Nmin;
    paraPlot.Nmax = Nmax;
    paraPlot.fn = fn2;
    varargout = {paraPlot};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [U,S,V] = blockToeplitz(h)
        %
        % [U,S,V] = blockToeplitz(h) calculate the shifted block Toeplitz matrix T1 and
        % the  result from the SVD of T1
        %
        % Input:
        % h: 3D-matrix of cross-correlation functions
        %
        % Outputs
        % U : result from SVD of H0
        % S : result from SVD of H0
        % V : result from SVD of H0
        %%
        if or(size(h,1)~=size(h,2),ndims(h)~=3)
            error('the IRF must be a 3D matrix with dimensions <M x M x N> ')
        end
        % get block Toeplitz matrix
        N1 = round(size(h,3)/2)-1;
        M = size(h,2);
        clear H0
        for oo=1:N1
            for ll=1:N1
                T1((oo-1)*M+1:oo*M,(ll-1)*M+1:ll*M) = h(:,:,N1+oo-ll+1);
            end
        end
        if or(any(isinf(T1(:))),any(isnan(T1(:))))
            warning('Input to SVD must not contain NaN or Inf. ')
            U=nan;
            S=nan;
            V=nan;
            return
        else
            try
                [U,S,V] = svd(T1);
            catch exception
                warning(' SVD of the block-Toeplitz failed ');
                U=nan;
                S=nan;
                V=nan;
                return
            end
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [IRF,t] = NExT(x,dt,Ts,method)
        %
        % [IRF,t] = NExT(x,dt,Ts,method) implements the Natural Excitation Technique to
        % retrieve the Impulse Response Function (IRF) from the cross-correlation
        % of the measured output y.
        %
        % Input:
        % x: time series of ambient vibrations: vector of size [1xN]
        % dt : Time step
        % Ts: Duration of subsegments (T<dt*(numel(y)-1))
        % method = 1 use the fft without zero padding.
        % method = 2 calls the function xcov with zero padding.
        %
        % Output
        % IRF: impulse response function
        % t: time vector asociated with the IRF
        %
        %%
        if nargin<4, method = 2; end % the fastest method is the default method
        if ~ismatrix(x), error('Error: x must be a vector or a matrix'),end
        
        if size(x,1)>size(x,2)
            x=x';
            [Nxx,~]=size(x);
        else
            [Nxx,~]=size(x);
        end
        
        % get the maximal segment length fixed by T
        M = round(Ts/dt);
        switch method
            case 1
                
                IRF = zeros(Nxx,Nxx,M);
                for oo=1:Nxx
                    for jj=1:Nxx
                        y1 = fft(x(oo,:));
                        y2 = fft(x(jj,:));
                        h0 = ifft(y1.*conj(y2));
                        IRF(oo,jj,:) = h0(1:M);
                    end
                end
                % get time vector t associated to the IRF
                t = (0:1:M-1)*dt;
                if Nxx==1,IRF = squeeze(IRF)';end
                
            case 2
                IRF = zeros(Nxx,Nxx,M+1);
                for oo=1:Nxx
                    for jj=1:Nxx
                        [dummy,lag]=xcov(x(oo,:),x(jj,:),M,'unbiased');
                        IRF(oo,jj,:) = dummy(end-round(numel(dummy)/2)+1:end);
                    end
                end
                if Nxx==1, IRF = squeeze(IRF)'; end
                % get time vector t associated to the IRF
                t = dt.*lag(end-round(numel(lag)/2)+1:end);
        end
        % normalize the IRF
        if Nxx==1,            IRF = IRF./IRF(1);        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [fn,zeta,phi] = modalID(U,S,Nmodes,Nyy,dt)
        % [fn,zeta,phi] = modalID(H1,U,S,V,N,M) identify the modal propeties of the
        % system.
        %
        % Input:
        % U: matrix of size [N1 x N1] obtained from te function blockToeplitz
        % S: matrix of size [N1 x N1] obtained from te function blockToeplitz
        % Nmodes: Number of modes (or poles)  scalar [1x1]
        % Nyy: Number of nodes (or sensors) along the line-like structure   scalar [1x1]
        % de: time step: scalar [1x1]
        %
        % Outputs
        % fn : Identified eigen frequencies
        % zeta : Identified damping ratios
        % phi : IDentified mode shapes
        
        
        if Nmodes>=size(S,1)
            warning(['Nmodes is larger than the numer of row of S. Nmodes is reduced to ',num2str(size(S,1))]);
            % extended observability matrix
            Nmodes = size(S,1);
        end
        
        O = U(:,1:Nmodes)*sqrt(S(1:Nmodes,1:Nmodes));
        % Get A and its eigen decomposition
        
        IndO = min(Nyy,size(O,1));
        C = O(1:IndO,:);
        jb = round(size(O,1)./IndO);
        A = pinv(O(1:IndO*(jb-1),:))*O(end-IndO*(jb-1)+1:end,:);
        [Vi,Di] = eig(A);
        
        mu = log(diag(Di))./dt; % poles
        fn = abs(mu(2:2:end))./(2*pi);% eigen-frequencies
        zeta = -real(mu(2:2:end))./abs(mu(2:2:end)); % modal amping ratio
        phi = real(C(1:IndO,:)*Vi); % mode shapes
        phi = phi(:,2:2:end);
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [fn,zeta,phi,MAC,stablity_status] = stabilityCheck(fn0,zeta0,phi0,fn1,zeta1,phi1)
        % [fn,zeta,phi,MAC,stablity_status] = stabilityCheck(fn0,zeta0,phi0,fn1,zeta1,phi1)
        % calculate the stability status of each mode obtained for
        % two adjacent poles (i,j).
        %
        % Input:
        % fn0: eigen frequencies calculated for pole i: vetor of N-modes [1 x N]
        % zeta0: modal damping ratio for pole i: vetor of N-modes [1 x N]
        % phi0: mode shape for pole i: vetor of N-modes [Nyy x N]
        % fn1: eigen frequencies calculated for pole j: vetor of N-modes [1 x N+1]
        % zeta1: modal damping ratio for pole j: vetor of N-modes [1 x N+1]
        % phi1: mode shape for pole j: vetor of N-modes [Nyy x N+1]
        %
        % Output:
        % fn: eigen frequencies calculated for pole j
        % zeta:  modal damping ratio for pole i
        % phi:mode shape for pole i
        % MAC: Mode Accuracy
        % stablity_status: stabilitystatus
        %%
        
        % Preallocation
        stablity_status = [];
        fn = [];
        zeta = [];
        phi = [];
        MAC=[];
        % frequency stability
        N0 = numel(fn0);
        N1 = numel(fn1);
        for rr=1:N0
            for jj=1:N1
                stab_fn = errCheck(fn0(rr),fn1(jj),eps_freq);
                stab_zeta = errCheck(zeta0(rr),zeta1(jj),eps_zeta);
                [stab_phi,dummyMAC] = getMAC(phi0(:,rr),phi1(:,jj),eps_MAC);
                % get stability status
                if stab_fn==0,
                    stabStatus = 0; % new pole
                elseif stab_fn == 1 & stab_phi == 1 & stab_zeta == 1,
                    stabStatus = 1; % stable pole
                elseif stab_fn == 1 & stab_zeta ==0 & stab_phi == 1,
                    stabStatus = 2; % pole with stable frequency and vector
                elseif stab_fn == 1 & stab_zeta == 1 & stab_phi ==0,
                    stabStatus = 3; % pole with stable frequency and damping
                elseif stab_fn == 1 & stab_zeta ==0 & stab_phi ==0,
                    stabStatus = 4; % pole with stable frequency
                else
                    error('Error: stablity_status is undefined')
                end
                fn = [fn,fn1(jj)];
                zeta = [zeta,zeta1(jj)];
                phi = [phi,phi1(:,jj)];
                MAC = [MAC,dummyMAC];
                stablity_status = [stablity_status,stabStatus];
            end
        end
        
        [fn,ind] = sort(fn);
        zeta = zeta(ind);
        phi = phi(:,ind);
        MAC = MAC(ind);
        stablity_status = stablity_status(ind);
        
        function y = errCheck(x0,x1,eps)
            if or(numel(x0)>1,numel(x1)>1),
                error('x0 and x1 must be a scalar');
            end
            if abs(1-x0./x1)<eps % if frequency for mode i+1 is almost unchanged
                y =1;
            else
                y = 0;
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [fnS,zetaS,phiS,MACS] = getStablePoles(fn,zeta,phi,MAC,stablity_status)
        fnS = [];zetaS = [];phiS=[];MACS = [];
        for oo=1:numel(fn)
            for jj=1:numel(stablity_status{oo})
                if stablity_status{oo}(jj)==1
                    fnS = [fnS,fn{oo}(jj)];
                    zetaS = [zetaS,zeta{oo}(jj)];
                    phiS = [phiS,phi{oo}(:,jj)];
                    MACS = [MACS,MAC{oo}(jj)];
                end
            end
        end
        
        % remove negative damping
        fnS(zetaS<=0)=[];
        phiS(:,zetaS<=0)=[];
        MACS(zetaS<=0)=[];
        zetaS(zetaS<=0)=[];
        
        % Normalized mode shape
        for oo=1:size(phiS,2)
            phiS(:,oo)= phiS(:,oo)./max(abs(phiS(:,oo)));
            if diff(phiS(1:2,oo))<0
                phiS(:,oo)=-phiS(:,oo);
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [fn,zeta,phi] = myClusterFun(fn0,zeta0,phi0)
        
        [~,Nsamples] = size(phi0);
        pos = zeros(Nsamples,Nsamples);
        for i1=1:Nsamples
            for i2=1:Nsamples
                [~,MAC0] = getMAC(phi0(:,i1),phi0(:,i2),eps_MAC); % here, eps_MAC is not important.
                pos(i1,i2) = abs((fn0(i1)-fn0(i2))./fn0(i2)) +1-MAC0; % compute MAC number between the selected mode shapes
            end
            
        end
        
        if numel(pos)==1,
            warning('At least one distance (two observations) are required');
            fn = nan;
            zeta = nan;
            phi = nan;
            return
        else

        Tree = PHA_Clustering(pos);
        [~, myClus0, Number] = Cluster2(Tree,'Limit',eps_cluster);
            
            
            
            Ncluster = numel(myClus0);
            
            ss=1;
            fn = {}; zeta = {}; phi = {};
            for rr=1:Ncluster
                myClus = myClus0{rr};
                if numel(myClus)>5
                    dummyZeta = zeta0(myClus);
                    dummyFn = fn0(myClus);
                    dummyPhi = phi0(:,myClus);
                    valMin = max(0,(quantile(dummyZeta,0.25) - abs(quantile(dummyZeta,0.75)-quantile(dummyZeta,0.25))*1.5));
                    valMax =quantile(dummyZeta,0.75) + abs(quantile(dummyZeta,0.75)-quantile(dummyZeta,0.25))*1.5;
                    dummyFn(or(dummyZeta>valMax,dummyZeta<valMin)) = [];
                    dummyPhi(:,or(dummyZeta>valMax,dummyZeta<valMin)) = [];
                    dummyZeta(or(dummyZeta>valMax,dummyZeta<valMin)) = [];
                    fn{ss} = dummyFn;
                    zeta{ss} = dummyZeta;
                    phi{ss} = dummyPhi;
                    ss=ss+1;
                end
            end
            if isempty(fn)
                fn = nan;
                zeta = nan;
                phi = nan;
                return
            end
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [y,dummyMAC] = getMAC(x0,x1,eps)
        Num = abs(x0(:)'*x1(:)).^2;
        D1= x0(:)'*x0(:);
        D2= x1(:)'*x1(:);
        dummyMAC = Num/(D1.*D2);
        if dummyMAC >(1-eps)
            y = 1;
        else
            y = 0;
        end
    end
% Analysis of tree
    function [Roots, Clusters, Number] = Cluster2(Tree,Parameter,Value)
        
        % Empty roots and nodes vectors
        Roots = [];
        Nodes = [];
        
        % Vector length
        N = max(max(Tree(:,1:2)))/2+1;
        
        % Clustering parameter
        switch lower(Parameter)
            
            % Number of clusters
            case 'number'
                Number = Value;
                Limit = Tree(end-Number+1+1,3);
                if Limit == 0
                    return
                end
                
                % Dissimilarity limit
            case 'limit'
                Limit = Value;
                Number = N-find(Tree(:,3)>=Limit,1,'first')+1;
                
        end
        
        % Clusters
        Clusters = cell(Number,1);
        
        % Cluster number
        Cluster = 0;
        
        % Exploration of the node
        ExplorationDown(Tree(end,1),Tree(end,3));
        ExplorationDown(Tree(end,2),Tree(end,3));
        
        % Exploration of the subnodes of a node
        function ExplorationDown(Node,Dissimilarity)
            
            if ismember(Node,Nodes) || ismember(Node,Roots)
                return
            end
            
            % Adding of the current node in nodes list
            Nodes = [Nodes Node];
            
            if Node <= N
                
                % Root
                Roots = [Roots Node];
                
                % Root whose distance is higher than limit
                [n,~]=find(Tree(:,1:2)==Node);
                if Tree(n,3) >= Limit
                    Cluster = Cluster+1;
                end
                
                % Cluster
                Clusters{Cluster} = [Clusters{Cluster} Node];
                
            else
                
                % Nodes
                Node = Node-N;
                
                % Subnodes
                N1 = Tree(Node,1);
                N2 = Tree(Node,2);
                
                % Cluster index increment
                if Tree(Node,3) < Limit && Dissimilarity >= Limit
                    Cluster = Cluster+1;
                end
                
                % Dissimilarity of current node
                Dissimilarity = Tree(Node,3);
                
                if N1 <= N && N2 <= N
                    
                    % Roots subnodes
                    if N1 < N2
                        ExplorationDown(N1,Dissimilarity);
                        ExplorationDown(N2,Dissimilarity);
                    else
                        ExplorationDown(N2,Dissimilarity);
                        ExplorationDown(N1,Dissimilarity);
                    end
                    
                else
                    
                    % Exploration of the subnodes
                    ExplorationDown(N1,Dissimilarity);
                    ExplorationDown(N2,Dissimilarity);
                    
                end
                
            end
            
        end
        
    end


    function [Z, totalPotential, parents] = PHA_Clustering(dMatrix, S)
        % ---Purpose---
        % Performs hierarchical clustering using the PHA method
        % The function will produce a hierarchical cluster tree (Z) from the input distance matrix
        % The output Z is similar to the output by the Matlab function 'linkage'
        %
        % ---INPUT---
        % dMatrix: distance matrix (numPts X numPts) defining distances between objects
        % S: (optional) Scale factor for determining parameter delta. The default value is S=10.
        %     If two points are closer than delta, they don't have attractive force.
        %
        % ---OUTPUT---
        % Z: hierarchical cluster tree  which is represented as a matrix with size (numPts-1 X 3)
        % totalPotential: total potential values
        % parents: the parent index of each data
        %
        % ---HOW TO USE---
        % Z = PHA_Clustering(dMatrix);
        % T = cluster(Z,'maxclust',k);
        %
        % ---Author---
        % Yonggang Lu (ylu@lzu.edu.cn)
        %
        % ---Reference---
        % Yonggang Lu, Yi Wan. (2013). ¡°PHA: A Fast Potential-based Hierarchical Agglomerative
        % Clustering Method, Pattern Recognition, Vol. 46(5), pp. 1227-1239.
        %
        [numPts,numPts2] = size(dMatrix); % numPts is the number of points
        if (numPts ~= numPts2)
            error(' PotentialHierachyWDistMatrix: distance matrix should be a square matrix! ');
        end
        if (nargin < 2)
            S = 10;
        end
        % compute the dalta automatically
        minDist = zeros(numPts, 1);
        for i = 1:numPts
            mask = (dMatrix(i,:)~=0);
            minDist(i) = min(dMatrix(i, mask));
        end
        delta = mean(minDist)/S;
        totalPotential = zeros(1, numPts);
        for i = 1:numPts  % for each point
            distToAll = dMatrix(i, :);
            selIdxes  = find(distToAll >= delta);
            totalP = sum(1./dMatrix(i, selIdxes));
            totalP = totalP + (numPts-length(selIdxes)-1)*(1/delta); % for points within delta, potential = 1/delta
            totalPotential(i)= - totalP;
        end
        [sortedP, sortedIdx] = sort(totalPotential);
        parents = [1:numPts]; % stores the parent information
        distToParent = zeros(1, numPts); % stores the distance to the parent point
        for pi = 2:numPts
            centerIdx = sortedIdx(pi);
            visitedPtsIdx = sortedIdx(1:pi-1);
            distToVisited = dMatrix(centerIdx, visitedPtsIdx);
            
            [minDist, minIdx] = min(distToVisited);
            parents(centerIdx) = visitedPtsIdx(minIdx);
            distToParent(centerIdx) = minDist;
        end
        % Z returns a (numPts-1) by 3 Matrix (same as linkage)
        Z = zeros(numPts-1, 3);
        [sortedDist, sortedIdx2] = sort(distToParent);
        linkIdx = [1:numPts]; % remember the index after each layer of merging
        for i = 1:numPts-1
            mergeIdx = sortedIdx2(i+1);
            Z(i, 1) = linkIdx(parents(mergeIdx));
            Z(i, 2) = linkIdx(mergeIdx);
            Z(i, 3) = sortedDist(i+1);
            linkIdx(linkIdx==Z(i, 2)) = numPts+i;
            linkIdx(linkIdx==Z(i, 1)) = numPts+i;
        end
    end

end
