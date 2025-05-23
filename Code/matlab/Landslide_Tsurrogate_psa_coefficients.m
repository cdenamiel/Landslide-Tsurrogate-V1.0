function Landslide_Tsurrogate_psa_coefficients

disp(' ')
disp('-------------------------')
disp('Landslide-Tsurrogate v1.0')
disp('-------------------------')
disp(' ')
disp('contact: Clea Denamiel')
disp('email: clea.denamiel@live.fr')
disp(' ')
disp('For more information: ')
disp(' - User Manual: ')
disp(' - Publication: ')
disp(' ')
disp('--------------------------------------------------------------------')
disp('STEP 5: generation of the deterministic coefficients of the gPCE PSA')
disp('--------------------------------------------------------------------')
disp(' ')
disp(' ')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('Before performing this step:')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('STEP 1: Edit and run Landslide_Tsurrogate_user_input to generate the file: ../results/output_users.mat')
disp('STEP 2: Run Landslide_Tsurrogate_input_parameters to generate the file: ../results/output_param.mat')
disp('STEP 3: Run the deterministic simulations outside Landslide-Tsurrogate v1.0 and generate: ../data/input_simus.mat')
disp('        input_simus.mat: contains zeta_max_surf[nsim,nx,ny] (maximum elevation), velo_max_surf[nsim,nx,ny] (maximum speed)')
disp('        ---------------  and time_max_surf[nsim,nx,ny] (time of arrival) with nsim the total number of simulations')
disp('                         corresponding to the maximum total order and [nx,ny] the spatial dimensions of the domain used to')
disp('                         perform the deterministic simulations')
disp('        Prepare the locations where to create the surrogates and generate: ../data/surrogate_models_locations.mat')
disp('        surrogate_models_locations.mat: contains index_coast[nl] the indices where to build the surrogate models')
disp('        ------------------------------  with nl the number of surrogate models to build')
disp('STEP 4: Run Landslide_Tsurrogate_format_input to generate the file: ../results/output_model.mat')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp(' ')

%--------------
% Load the data
%--------------

my_file='../results/output_users.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

my_file='../results/output_param.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

my_file='../results/output_model.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

%---------------------------------------------------
% Create the PSA coefficient of the surrogate models
%---------------------------------------------------

coeff = surrogate_model_gauss_patterson_psa_coeff(param,model,a,b,maxdeg+1,option);
save ../results/output_coeff.mat coeff
disp(' ')
disp('The results have been saved in: ../results/output_coeff.mat')
disp(' ')
    
end

function coeff = surrogate_model_gauss_patterson_psa_coeff(param,model,a,b,maxdeg,option)
%% Pseudo Spectral Approximation with Gauss-Patterson Sparse grids 
% based on:
% (1) Florian Heiss, Viktor Winschel (2008): "Likelihood approximation by 
% numerical integration on sparse grids".Journal of Econometrics,
% Vol. 144, pp. 62-80.
%--------------------------------------------------------------------------
% Input:
% ------
%
%       a,b:     limits of the uniform distribution such as U([a,b])
%                dimension: (nmodes)
%
%    option:     0 = normal or 1 = delayed Gauss-Patterson sequencies 
%
%
%     param:     structure containing the paramters 
%                dimension: (maxdeg+1)
%--------------------------------------------
%         Z:     chosen samples of the stochastic variables 
%                dimensions: (nmodes) x (nwj) => dimension of the set of samples
%     index:     index corresponding to the unique model run 
%                dimension: (nwj)
%
%     model:     structure containing the model results 
%                dimension: (nwj)
%--------------------------------------------
%    Zeta_Z:     extracted elevation corresponding to Z at stations
%                dimensions: (nx)
%    Velo_Z:     extracted speed corresponding to Z at stations
%                dimensions: (nx)
%    Time_Z:     extracted time of arrival corresponding to Z at stations
%                dimensions: (nx)
%--------------------------------------------------------------------------
% Output:
%--------
%
% coeff:         structure containing the coefficients of the PCE
%                dimension: (maxdeg+1)
%--------------------------------------------
%     alpha:     unique set of multi-indexes for dregree maxdeg
%                dimension: (nl) x (nmodes)
%      zeta:     polynomial coefficient for zeta at nx stations
%                dimension: (nx) x (nl) 
%      velo:     polynomial coefficient for velo at nx stations
%                dimension: (nx) x (nl) 
%      time:     polynomial coefficient for time at nx stations
%                dimension: (nx) x (nl) 
%--------------------------------------------------------------------------
%
%
% Initialization of the matrices
%
nmodes      = length(a);
nx          = length(model(1).Zeta_Z);
%
% Legendre Polynomials up to maxdeg
%
Le    = cell(maxdeg,1);
Le{1} = 1;       % L_1 = 1
Le{2} = [1 0];   % L_2 = x
for n = 3:maxdeg
   Le{n} = ((2*n-3)/(n-1))*[Le{n-1} 0] - ((n-2)/(n-1))*[0 0 Le{n-2}];
end
%
% PCE computation
%
coeff = struc([]);
Z     = param(maxdeg).Z;
index = param(maxdeg).index;
for alpha_norm1 = 0:maxdeg-1
    %
    disp('-------------------------------')
    disp(['||alpha||_1 = ' num2str(alpha_norm1)])
    disp('-------------------------------')  
    % Reset index to 1 
    start = 1;   
    % multi-indexes alphal with sum equal to alpha_norm1
    alpha       = get_seq(nmodes,nmodes+alpha_norm1);
    for l=1:size(alpha,1)
        disp(['Completion: ',num2str(100*l/size(alpha,1)),'%'])
        [Znew,Wnew] = gauss_patterson_nested_rule_nd(alpha(l,:),a,b,option);
        Znew        = single(Znew);
        % Find the corresponding simulations
        Zeta_Z = nan(nx,length(Wnew));
        Velo_Z = nan(nx,length(Wnew));
        Time_Z = nan(nx,length(Wnew));
        for ss = 1:length(Wnew)
            [~,~,ind_new] = intersect(Znew(ss,:),Z,'rows');
            Zeta_Z(:,ss) = log(model(index(ind_new)).Zeta_Z);
            Velo_Z(:,ss) = log(model(index(ind_new)).Velo_Z);
            Time_Z(:,ss) = log(model(index(ind_new)).Time_Zeta_Z);
        end
        % multi-indexes alpha_sum varying between [0 0 ... 0] and alphal 
        alpha_sum = permn(0:max(alpha(l,:)),nmodes);
        inds = find(sum(alpha_sum<=repmat(alpha(l,:),size(alpha_sum,1),1),2)==nmodes);                     
        alpha_sum = alpha_sum(inds,:); %#ok<*FNDSB>
        zeta_temp = zeros(nx,size(alpha_sum,1));
        velo_temp = zeros(nx,size(alpha_sum,1));
        time_temp = zeros(nx,size(alpha_sum,1));
        for s = 1:size(alpha_sum,1)
            % Calculate the quadrature coefficients
            qrule_coeff = Wnew;
            for n = 1:nmodes
                qrule_coeff   = qrule_coeff.*polyval(Le{alpha_sum(s,n)+1},(2.*Znew(:,n)-a(n)-b(n))./(b(n)-a(n)))...
                    ./norm2_legendre_polynomials(alpha_sum(s,n),a(n),b(n));
            end      
            % Quadrature rule - Coefficients of the PCE
            for i = 1:nx
                zeta_temp(i,s) = Zeta_Z(i,:)*qrule_coeff;
                velo_temp(i,s) = Velo_Z(i,:)*qrule_coeff;
                time_temp(i,s) = Time_Z(i,:)*qrule_coeff;
            end
         end
         stop=start+size(alpha_sum,1)-1;
         alpha_final(start:stop,:) = alpha_sum;
         zeta_final(:,start:stop)= zeta_temp;
         velo_final(:,start:stop)= velo_temp;
         time_final(:,start:stop)= time_temp;
         start=stop+1;
    end
    coeff(alpha_norm1+1).alpha = alpha_final;
    coeff(alpha_norm1+1).zeta  = zeta_final;
    coeff(alpha_norm1+1).velo  = velo_final;
    coeff(alpha_norm1+1).time  = time_final;
    clear alpha_final zeta_final time_final velo_final
end
return
end

function [M, I] = permn(V, N, K)
% PERMN - permutations with repetition
%   Using two input variables V and N, M = PERMN(V,N) returns all
%   permutations of N elements taken from the vector V, with repetitions.
%   V can be any type of array (numbers, cells etc.) and M will be of the
%   same type as V.  If V is empty or N is 0, M will be empty.  M has the
%   size numel(V).^N-by-N. 
%
%   When only a subset of these permutations is needed, you can call PERMN
%   with 3 input variables:  M = PERMNSUB(V,N,K). 
%   M will return only the K-ths permutations.  The output is the same as M
%   = PERMN(V,N), followed by M = M(K,:), but it avoids memory issues by
%   generating all possible combinations first.  This is particulary useful
%   when you only need one, or a small subset of all permutations at a
%   given time. If V or K is empty, or N is zerp, M will be empty. M has
%   the size numel(K)-by-N. 
%
%   [M, I] = PERMN(...) also returns an index matrix I so that M = V(I).
%
%   Examples:
%     M = permn([1 2 3],2) % returns the 9-by-2 matrix:
%              1     1
%              1     2
%              1     3
%              2     1
%              2     2
%              2     3
%              3     1
%              3     2
%              3     3
%
%     M = permn([99 7],4) % returns the 16-by-4 matrix:
%              99     99    99    99
%              99     99    99     7
%              99     99     7    99
%              99     99     7     7
%              ...
%               7      7     7    99
%               7      7     7     7
%
%     M = permn({'hello!' 1:3},2) % returns the 4-by-2 cell array
%             'hello!'        'Ahello!'
%             'hello!'        [1x3 double]
%             [1x3 double]    'hello!'
%             [1x3 double]    [1x3 double]
%
%     V = 11:15, N = 3, K = [2 124 21 99]
%     M = permn(V, N, K) % returns the 4-by-3 matrix:
%     %        1  1  2
%     %        5  5  4
%     %        1  5  1
%     %        4  5  4
%     % which are the 2nd, 124th, 21st and 99th combinations
%     % Check with PERMN
%     M2 = permn(V,N) ; isequal(M2(K,:),M)
%     % Note that M2 is a 125-by-3 matrix
%
%     % PERMN can be used generate a binary table
%     B = permn([0 1],5)  
%
%   NB Matrix sizes increases exponentially at rate (n^N)*N.
%
%   See also PERMS, NCHOOSEK
%            ALLCOMB, PERMPOS on the File Exchange

% tested in Matlab R13, R14, 2010b, 2014a
% version 6.0 (may 2015)
% (c) Jos van der Geest
% Matlab File Exchange Author ID: 10584
% email: samelinoa@gmail.com

% History
% 1.1 updated help text
% 2.0 new faster algorithm
% 3.0 (aug 2006) implemented very fast algorithm
% 3.1 (may 2007) Improved algorithm Roger Stafford pointed out that for some values, the floor
%   operation on floating points, according to the IEEE 754 standard, could return
%   erroneous values. His excellent solution was to add (1/2) to the values
%   of A.
% 3.2 (may 2007) changed help and error messages slightly
% 4.0 (may 2008) again a faster implementation, based on ALLCOMB, suggested on the
%   newsgroup comp.soft-sys.matlab on May 7th 2008 by "Helper". It was
%   pointed out that COMBN(V,N) equals ALLCOMB(V,V,V...) (V repeated N
%   times), ALLCMOB being faster. Actually version 4 is an improvement
%   over version 1 ...
% 4.1 (jan 2010) removed call to FLIPLR, using refered indexing N:-1:1
%   (is faster, suggestion of Jan Simon, jan 2010), removed REPMAT, and
%   let NDGRID handle this
% 4.2 (apr 2011) corrrectly return a column vector for N = 1 (error pointed
%    out by Wilson).
% 4.3 (apr 2013) make a reference to COMBNSUB
% 5.0 (may 2015) NAME CHANGED (COMBN -> PERMN) and updated description,
%   following comment by Stephen Obeldick that this function is misnamed
%   as it produces permutations with repetitions rather then combinations.
% 5.1 (may 2015) always calculate M via indices
% 6.0 (may 2015) merged the functionaly of permnsub (aka combnsub) and this
%   function
narginchk(2,3) ;
if fix(N) ~= N || N < 0 || numel(N) ~= 1 ;
    error('permn:negativeN','Second argument should be a positive integer') ;
end
nV = numel(V) ;
if nargin==2, % PERMN(V,N) - return all permutations
    if nV==0 || N == 0,
        M = zeros(nV,N) ;
        I = zeros(nV,N) ;
    elseif N == 1,
        % return column vectors
        M = V(:) ;
        I = (1:nV).' ;
    else
        % this is faster than the math trick used for the call with three
        % arguments.
        [Y{N:-1:1}] = ndgrid(1:nV) ;
        I = reshape(cat(N+1,Y{:}),[],N) ;
        % I = local_allcomb(1:nV, N) ;
        M = V(I) ;
    end
else % PERMN(V,N,K) - return a subset of all permutations
    nK = numel(K) ;
    if nV == 0 || N == 0 || nK == 0
        M = zeros(numel(K), N) ;
        I = zeros(numel(K), N) ;
    elseif nK < 1 || any(K<1) || any(K ~= fix(K))
        error('permn:InvalidIndex','Third argument should contain positive integers.') ;
    else
        V = reshape(V,1,[]) ; % v1.1 make input a row vector
        nV = numel(V) ;
        Npos = nV^N ;
        if any(K > Npos)
            warning('permn:IndexOverflow', ...
                'Values of K exceeding the total number of combinations are saturated.')
            K = min(K, Npos) ;
        end             
        % The engine is based on version 3.2 of COMBN  with the correction
        % suggested by Roger Stafford. This approaches uses a single matrix
        % multiplication.
        B = nV.^(1-N:0) ;
        I = ((K(:)-.5) * B) ; % matrix multiplication
        I = rem(floor(I),nV) + 1 ;
        M = V(I) ;
    end
end
end

function PsiN2 = norm2_legendre_polynomials(order,a,b)
Le    = cell(order+1,1);
Le{1} = 1;       % L_1 = 1
Le{2} = [1 0];   % L_2 = x
for n = 3:order+1
   Le{n} = ((2*n-3)/(n-1))*[Le{n-1} 0] - ((n-2)/(n-1))*[0 0 Le{n-2}];
end
fun = @(x) polyval(Le{order+1},(2*x-a-b)/(b-a)).*polyval(Le{order+1},(2*x-a-b)/(b-a));
PsiN2 = integral(@(x)fun(x),a,b);
return
end

function [Z,W] = gauss_patterson_nested_rule_nd(alphai,a,b,option)
%% Gauss-Patterson nested rule
%-------------------------------------------------------------------------- 
% Input:
%-------
%
% alphai:     multi-index of polynomial degrees 
%             dimension: (nmodes) = number of stochastic variables 
%
% a,b:        limits of the uniform distribution such as U([a,b])
%             dimension: (nmodes) 
%
% option:     0 = normal or 1 = delayed Gauss-Patterson sequencies 
%
%-------------------------------------------------------------------------- 
% Output:
%--------
%
% Z:        quadrature points of the multivariate polynomial basis
%           dimension: (nmodes) x (nwj) => dimension of the set of samples
%
% W:        associated weights of the multivariate polynomial basis  
%           dimension: (nwj) => dimension of the set of samples
%
% -------------------------------------------------------------------------
%
% Initialization of the matrices:
[~,nmodes] = size(alphai);
z = cell(nmodes,1);
w = cell(nmodes,1);
% order of the delayed Gauss-Patterson rule up to polynomial degree: max_degree = 100
max_degree = 100;
alpha_nest_delayed = ones(max_degree+1,1);
for ll=2:max_degree+1
    p=5; 
    o=3; 
    while(p<2*(ll-1)+1)
        p=(2*p)+1;
        o=(2*o)+1;
    end 
    alpha_nest_delayed(ll)=o;
end
%
for n=1:nmodes
% For a given alpha_nest depending on the degree alpha of a given stochastic variable  
% find the points and the weights of the univariate Gauss-Patterson
% quadrature.
% 
% Calculation of the level of accuracy for quadrature
level = alphai(n)+1;
if(option==0)
% ACCURATE METHOD:
    alpha_nest = (2.^level) -1 ;    
else
% DELAYED METHOD:
    alpha_nest = alpha_nest_delayed(level);
end 
% Gauss-Patterson sequences
    [z{n},w{n}] = gauss_patterson_nested_rule_1d(alpha_nest,a(n),b(n));
end
%
% For all the stochastic variables, calculate the tensor product of the
% Pseudospectral Approximation. In other words build the tensor grids of  
% points (Z) and weights (W)
%
Z=z{1}; W=w{1};
for n=2:nmodes
    Z=[kron(Z,ones(length(z{n}),1)) kron(ones(size(Z,1),1),z{n})];
    W=kron(W,w{n});
end
return
end

function alpha = get_seq(nmodes,alpha_norm)
%% Get Sequence of multi-indexes
%-------------------------------------------------------------------------- 
% Input:
%-------  
% 
% nmodes:     number of stochastic variables
%
% alpha_norm: row sum of elements each of the rows has to have
%
%-------------------------------------------------------------------------- 
% Output:
%-------  
% 
% alpha:     matrix with nmodes columns 
%            Each row represents one vector with all elements >=1
%            and the sum of elements == norm
%
%--------------------------------------------------------------------------   
% 
seq = zeros(1,nmodes);
a=alpha_norm-nmodes;
seq(1)=a;
alpha = seq;
c=1;
while seq(nmodes)<a
    if (c==nmodes) 
        for i=(c-1):-1:1
            c=i;
            if seq(i)~=0
                break
            end
        end
    end
    seq(c) = seq(c)-1;
    c=c+1;
    seq(c) = a - sum(seq(1:(c-1)));
    if (c<nmodes) 
        seq((c+1):nmodes)=zeros(1,nmodes-c);
    end
    alpha = [alpha;seq]; %#ok<*AGROW>
end 
alpha = alpha; %#ok<*ASGSL>
end

function [zj,wj] = gauss_patterson_nested_rule_1d(alphajj,ajj,bjj)
%% Gauss-Patterson quadrature based on:
%
% (1)Prem Kythe, Michael Schaeferkotter,
%    Handbook of Computational Methods for Integration, Chapman and Hall, 
%    2004,ISBN: 1-58488-428-2, LC: QA299.3.K98.
%
% (2)Thomas Patterson,
%    The Optimal Addition of Points to Quadrature Formulae,Mathematics of 
%    Computation,Volume 22, Number 104, October 1968, pages 847-856.
%
%  Parameters:
%
%    Input, integer N, the order of the rule.
%    ORDER must be 1, 3, 7, 15, 31, 63, 127, 255 or 511.
%
%    Output, real X(N), the abscissas of the rule.
%
%    Output, real W(N), the weights of the rule.
%    The weights are positive, symmetric and should sum to 2.
%
%-------------------------------------------------------------------------- 
% Input:
%-------
%
% alphajj: univariate polynomial degree, alphajj <= 511                
%
%-------------------------------------------------------------------------- 
% Output:
%-------
%
% zjj:     nodes of the univariate rule
% wjj:     weights of the univariate rule
%--------------------------------------------------------------------------
% Discussion:
%------------
%
%    The integration interval is [ -1, 1 ]. The weight function is w(x) = 1.0.
%    The integral to approximate: Integral ( -1 <= X <= 1 ) F(X) dX
%    The quadrature rule: Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
%    The zeroth rule, of order 1, is the standard Gauss-Legendre rule.
%    The first rule, of order 3, is the standard Gauss-Legendre rule.
%    The second rule, of order 7, includes the abscissas of the previous rule.
%    Each subsequent rule is nested in a similar way.  
%    Rules are available of orders 1, 3, 7, 15, 31, 63, 127, 255 or 511.
%    The data for N = 511 was supplied by Dirk Laurie.
%
% -------------------------------------------------------------------------
%
zj = zeros (alphajj,1);
wj = zeros (alphajj,1);
switch alphajj
  case 1
    zj(1) = 0.0;
    wj(1) = 2.0;
  case 3
    zj(1) = -0.77459666924148337704;
    zj(2) =  0.0;
    zj(3) =  0.774596669241483037704;
    wj(1) = 0.555555555555555555556;
    wj(2) = 0.888888888888888888889;
    wj(3) = 0.555555555555555555556;
  case 7
    zj(1) = -0.96049126870802028342;
    zj(2) = -0.77459666924148337704;
    zj(3) = -0.43424374934680255800;
    zj(4) =  0.0;
    zj(5) =  0.43424374934680255800;
    zj(6) =  0.77459666924148337704;
    zj(7) =  0.96049126870802028342;
    wj(1) = 0.104656226026467265194;
    wj(2) = 0.268488089868333440729;
    wj(3) = 0.401397414775962222905;
    wj(4) = 0.450916538658474142345;
    wj(5) = 0.401397414775962222905;
    wj(6) = 0.268488089868333440729;
    wj(7) = 0.104656226026467265194;
  case 15
    zj( 1) = -0.99383196321275502221;
    zj( 2) = -0.96049126870802028342;
    zj( 3) = -0.88845923287225699889;
    zj( 4) = -0.77459666924148337704;
    zj( 5) = -0.62110294673722640294;
    zj( 6) = -0.43424374934680255800;
    zj( 7) = -0.22338668642896688163;
    zj( 8) =  0.0;
    zj( 9) =  0.22338668642896688163;
    zj(10) =  0.43424374934680255800;
    zj(11) =  0.62110294673722640294;
    zj(12) =  0.77459666924148337704;
    zj(13) =  0.88845923287225699889;
    zj(14) =  0.96049126870802028342;
    zj(15) =  0.99383196321275502221;
    wj( 1) = 0.0170017196299402603390;
    wj( 2) = 0.0516032829970797396969;
    wj( 3) = 0.0929271953151245376859;
    wj( 4) = 0.134415255243784220360;
    wj( 5) = 0.171511909136391380787;
    wj( 6) = 0.200628529376989021034;
    wj( 7) = 0.219156858401587496404;
    wj( 8) = 0.225510499798206687386;
    wj( 9) = 0.219156858401587496404;
    wj(10) = 0.200628529376989021034;
    wj(11) = 0.171511909136391380787;
    wj(12) = 0.134415255243784220360;
    wj(13) = 0.0929271953151245376859;
    wj(14) = 0.0516032829970797396969;
    wj(15) = 0.0170017196299402603390;
case 31
    zj( 1) = -0.99909812496766759766;
    zj( 2) = -0.99383196321275502221;
    zj( 3) = -0.98153114955374010687;
    zj( 4) = -0.96049126870802028342;
    zj( 5) = -0.92965485742974005667;
    zj( 6) = -0.88845923287225699889;
    zj( 7) = -0.83672593816886873550;
    zj( 8) = -0.77459666924148337704;
    zj( 9) = -0.70249620649152707861;
    zj(10) = -0.62110294673722640294;
    zj(11) = -0.53131974364437562397;
    zj(12) = -0.43424374934680255800;
    zj(13) = -0.33113539325797683309;
    zj(14) = -0.22338668642896688163;
    zj(15) = -0.11248894313318662575;
    zj(16) =  0.0;
    zj(17) =  0.11248894313318662575;
    zj(18) =  0.22338668642896688163;
    zj(19) =  0.33113539325797683309;
    zj(20) =  0.43424374934680255800;
    zj(21) =  0.53131974364437562397;
    zj(22) =  0.62110294673722640294;
    zj(23) =  0.70249620649152707861;
    zj(24) =  0.77459666924148337704;
    zj(25) =  0.83672593816886873550;
    zj(26) =  0.88845923287225699889;
    zj(27) =  0.92965485742974005667;
    zj(28) =  0.96049126870802028342;
    zj(29) =  0.98153114955374010687;
    zj(30) =  0.99383196321275502221;
    zj(31) =  0.99909812496766759766;
    wj( 1) = 0.00254478079156187441540;
    wj( 2) = 0.00843456573932110624631;
    wj( 3) = 0.0164460498543878109338;
    wj( 4) = 0.0258075980961766535646;
    wj( 5) = 0.0359571033071293220968;
    wj( 6) = 0.0464628932617579865414;
    wj( 7) = 0.0569795094941233574122;
    wj( 8) = 0.0672077542959907035404;
    wj( 9) = 0.0768796204990035310427;
    wj(10) = 0.0857559200499903511542;
    wj(11) = 0.0936271099812644736167;
    wj(12) = 0.100314278611795578771;
    wj(13) = 0.105669893580234809744;
    wj(14) = 0.109578421055924638237;
    wj(15) = 0.111956873020953456880;
    wj(16) = 0.112755256720768691607;
    wj(17) = 0.111956873020953456880;
    wj(18) = 0.109578421055924638237;
    wj(19) = 0.105669893580234809744;
    wj(20) = 0.100314278611795578771;
    wj(21) = 0.0936271099812644736167;
    wj(22) = 0.0857559200499903511542;
    wj(23) = 0.0768796204990035310427;
    wj(24) = 0.0672077542959907035404;
    wj(25) = 0.0569795094941233574122;
    wj(26) = 0.0464628932617579865414;
    wj(27) = 0.0359571033071293220968;
    wj(28) = 0.0258075980961766535646;
    wj(29) = 0.0164460498543878109338;
    wj(30) = 0.00843456573932110624631;
    wj(31) = 0.00254478079156187441540;  
case 63
    zj( 1) = -0.99987288812035761194;
    zj( 2) = -0.99909812496766759766;
    zj( 3) = -0.99720625937222195908;
    zj( 4) = -0.99383196321275502221;
    zj( 5) = -0.98868475754742947994;
    zj( 6) = -0.98153114955374010687;
    zj( 7) = -0.97218287474858179658;
    zj( 8) = -0.96049126870802028342;
    zj( 9) = -0.94634285837340290515;
    zj(10) = -0.92965485742974005667;
    zj(11) = -0.91037115695700429250;
    zj(12) = -0.88845923287225699889;
    zj(13) = -0.86390793819369047715;
    zj(14) = -0.83672593816886873550;
    zj(15) = -0.80694053195021761186;
    zj(16) = -0.77459666924148337704;
    zj(17) = -0.73975604435269475868;
    zj(18) = -0.70249620649152707861;
    zj(19) = -0.66290966002478059546;
    zj(20) = -0.62110294673722640294;
    zj(21) = -0.57719571005204581484;
    zj(22) = -0.53131974364437562397;
    zj(23) = -0.48361802694584102756;
    zj(24) = -0.43424374934680255800;
    zj(25) = -0.38335932419873034692;
    zj(26) = -0.33113539325797683309;
    zj(27) = -0.27774982202182431507;
    zj(28) = -0.22338668642896688163;
    zj(29) = -0.16823525155220746498;
    zj(30) = -0.11248894313318662575;
    zj(31) = -0.056344313046592789972;
    zj(32) =  0.0;
    zj(33) =  0.056344313046592789972;
    zj(34) =  0.11248894313318662575;
    zj(35) =  0.16823525155220746498;
    zj(36) =  0.22338668642896688163;
    zj(37) =  0.27774982202182431507;
    zj(38) =  0.33113539325797683309;
    zj(39) =  0.38335932419873034692;
    zj(40) =  0.43424374934680255800;
    zj(41) =  0.48361802694584102756;
    zj(42) =  0.53131974364437562397;
    zj(43) =  0.57719571005204581484;
    zj(44) =  0.62110294673722640294;
    zj(45) =  0.66290966002478059546;
    zj(46) =  0.70249620649152707861;
    zj(47) =  0.73975604435269475868;
    zj(48) =  0.77459666924148337704;
    zj(49) =  0.80694053195021761186;
    zj(50) =  0.83672593816886873550;
    zj(51) =  0.86390793819369047715;
    zj(52) =  0.88845923287225699889;
    zj(53) =  0.91037115695700429250;
    zj(54) =  0.92965485742974005667;
    zj(55) =  0.94634285837340290515;
    zj(56) =  0.96049126870802028342;
    zj(57) =  0.97218287474858179658;
    zj(58) =  0.98153114955374010687;
    zj(59) =  0.98868475754742947994;
    zj(60) =  0.99383196321275502221;
    zj(61) =  0.99720625937222195908;
    zj(62) =  0.99909812496766759766;
    zj(63) =  0.99987288812035761194;
    wj( 1) = 0.000363221481845530659694;
    wj( 2) = 0.00126515655623006801137;
    wj( 3) = 0.00257904979468568827243;
    wj( 4) = 0.00421763044155885483908;
    wj( 5) = 0.00611550682211724633968;
    wj( 6) = 0.00822300795723592966926;
    wj( 7) = 0.0104982469096213218983;
    wj( 8) = 0.0129038001003512656260;
    wj( 9) = 0.0154067504665594978021;
    wj(10) = 0.0179785515681282703329;
    wj(11) = 0.0205942339159127111492;
    wj(12) = 0.0232314466399102694433;
    wj(13) = 0.0258696793272147469108;
    wj(14) = 0.0284897547458335486125;
    wj(15) = 0.0310735511116879648799;
    wj(16) = 0.0336038771482077305417;
    wj(17) = 0.0360644327807825726401;
    wj(18) = 0.0384398102494555320386;
    wj(19) = 0.0407155101169443189339;
    wj(20) = 0.0428779600250077344929;
    wj(21) = 0.0449145316536321974143;
    wj(22) = 0.0468135549906280124026;
    wj(23) = 0.0485643304066731987159;
    wj(24) = 0.0501571393058995374137;
    wj(25) = 0.0515832539520484587768;
    wj(26) = 0.0528349467901165198621;
    wj(27) = 0.0539054993352660639269;
    wj(28) = 0.0547892105279628650322;
    wj(29) = 0.0554814043565593639878;
    wj(30) = 0.0559784365104763194076;
    wj(31) = 0.0562776998312543012726;
    wj(32) = 0.0563776283603847173877;
    wj(33) = 0.0562776998312543012726;
    wj(34) = 0.0559784365104763194076;
    wj(35) = 0.0554814043565593639878;
    wj(36) = 0.0547892105279628650322;
    wj(37) = 0.0539054993352660639269;
    wj(38) = 0.0528349467901165198621;
    wj(39) = 0.0515832539520484587768;
    wj(40) = 0.0501571393058995374137;
    wj(41) = 0.0485643304066731987159;
    wj(42) = 0.0468135549906280124026;
    wj(43) = 0.0449145316536321974143;
    wj(44) = 0.0428779600250077344929;
    wj(45) = 0.0407155101169443189339;
    wj(46) = 0.0384398102494555320386;
    wj(47) = 0.0360644327807825726401;
    wj(48) = 0.0336038771482077305417;
    wj(49) = 0.0310735511116879648799;
    wj(50) = 0.0284897547458335486125;
    wj(51) = 0.0258696793272147469108;
    wj(52) = 0.0232314466399102694433;
    wj(53) = 0.0205942339159127111492;
    wj(54) = 0.0179785515681282703329;
    wj(55) = 0.0154067504665594978021;
    wj(56) = 0.0129038001003512656260;
    wj(57) = 0.0104982469096213218983;
    wj(58) = 0.00822300795723592966926;
    wj(59) = 0.00611550682211724633968;
    wj(60) = 0.00421763044155885483908;
    wj(61) = 0.00257904979468568827243;
    wj(62) = 0.00126515655623006801137;
    wj(63) = 0.000363221481845530659694;
  case 127
    zj(  1) = -0.99998243035489159858;
    zj(  2) = -0.99987288812035761194;
    zj(  3) = -0.99959879967191068325;
    zj(  4) = -0.99909812496766759766;
    zj(  5) = -0.99831663531840739253;
    zj(  6) = -0.99720625937222195908;
    zj(  7) = -0.99572410469840718851;
    zj(  8) = -0.99383196321275502221;
    zj(  9) = -0.99149572117810613240;
    zj( 10) = -0.98868475754742947994;
    zj( 11) = -0.98537149959852037111;
    zj( 12) = -0.98153114955374010687;
    zj( 13) = -0.97714151463970571416;
    zj( 14) = -0.97218287474858179658;
    zj( 15) = -0.96663785155841656709;
    zj( 16) = -0.96049126870802028342;
    zj( 17) = -0.95373000642576113641;
    zj( 18) = -0.94634285837340290515;
    zj( 19) = -0.93832039777959288365;
    zj( 20) = -0.92965485742974005667;
    zj( 21) = -0.92034002547001242073;
    zj( 22) = -0.91037115695700429250;
    zj( 23) = -0.89974489977694003664;
    zj( 24) = -0.88845923287225699889;
    zj( 25) = -0.87651341448470526974;
    zj( 26) = -0.86390793819369047715;
    zj( 27) = -0.85064449476835027976;
    zj( 28) = -0.83672593816886873550;
    zj( 29) = -0.82215625436498040737;
    zj( 30) = -0.80694053195021761186;
    zj( 31) = -0.79108493379984836143;
    zj( 32) = -0.77459666924148337704;
    zj( 33) = -0.75748396638051363793;
    zj( 34) = -0.73975604435269475868;
    zj( 35) = -0.72142308537009891548;
    zj( 36) = -0.70249620649152707861;
    zj( 37) = -0.68298743109107922809;
    zj( 38) = -0.66290966002478059546;
    zj( 39) = -0.64227664250975951377;
    zj( 40) = -0.62110294673722640294;
    zj( 41) = -0.59940393024224289297;
    zj( 42) = -0.57719571005204581484;
    zj( 43) = -0.55449513263193254887;
    zj( 44) = -0.53131974364437562397;
    zj( 45) = -0.50768775753371660215;
    zj( 46) = -0.48361802694584102756;
    zj( 47) = -0.45913001198983233287;
    zj( 48) = -0.43424374934680255800;
    zj( 49) = -0.40897982122988867241;
    zj( 50) = -0.38335932419873034692;
    zj( 51) = -0.35740383783153215238;
    zj( 52) = -0.33113539325797683309;
    zj( 53) = -0.30457644155671404334;
    zj( 54) = -0.27774982202182431507;
    zj( 55) = -0.25067873030348317661;
    zj( 56) = -0.22338668642896688163;
    zj( 57) = -0.19589750271110015392;
    zj( 58) = -0.16823525155220746498;
    zj( 59) = -0.14042423315256017459;
    zj( 60) = -0.11248894313318662575;
    zj( 61) = -0.084454040083710883710;
    zj( 62) = -0.056344313046592789972;
    zj( 63) = -0.028184648949745694339;
    zj( 64) =  0.0;
    zj( 65) =  0.028184648949745694339;
    zj( 66) =  0.056344313046592789972;
    zj( 67) =  0.084454040083710883710;
    zj( 68) =  0.11248894313318662575;
    zj( 69) =  0.14042423315256017459;
    zj( 70) =  0.16823525155220746498;
    zj( 71) =  0.19589750271110015392;
    zj( 72) =  0.22338668642896688163;
    zj( 73) =  0.25067873030348317661;
    zj( 74) =  0.27774982202182431507;
    zj( 75) =  0.30457644155671404334;
    zj( 76) =  0.33113539325797683309;
    zj( 77) =  0.35740383783153215238;
    zj( 78) =  0.38335932419873034692;
    zj( 79) =  0.40897982122988867241;
    zj( 80) =  0.43424374934680255800;
    zj( 81) =  0.45913001198983233287;
    zj( 82) =  0.48361802694584102756;
    zj( 83) =  0.50768775753371660215;
    zj( 84) =  0.53131974364437562397;
    zj( 85) =  0.55449513263193254887;
    zj( 86) =  0.57719571005204581484;
    zj( 87) =  0.59940393024224289297;
    zj( 88) =  0.62110294673722640294;
    zj( 89) =  0.64227664250975951377;
    zj( 90) =  0.66290966002478059546;
    zj( 91) =  0.68298743109107922809;
    zj( 92) =  0.70249620649152707861;
    zj( 93) =  0.72142308537009891548;
    zj( 94) =  0.73975604435269475868;
    zj( 95) =  0.75748396638051363793;
    zj( 96) =  0.77459666924148337704;
    zj( 97) =  0.79108493379984836143;
    zj( 98) =  0.80694053195021761186;
    zj( 99) =  0.82215625436498040737;
    zj(100) =  0.83672593816886873550;
    zj(101) =  0.85064449476835027976;
    zj(102) =  0.86390793819369047715;
    zj(103) =  0.87651341448470526974;
    zj(104) =  0.88845923287225699889;
    zj(105) =  0.89974489977694003664;
    zj(106) =  0.91037115695700429250;
    zj(107) =  0.92034002547001242073;
    zj(108) =  0.92965485742974005667;
    zj(109) =  0.93832039777959288365;
    zj(110) =  0.94634285837340290515;
    zj(111) =  0.95373000642576113641;
    zj(112) =  0.96049126870802028342;
    zj(113) =  0.96663785155841656709;
    zj(114) =  0.97218287474858179658;
    zj(115) =  0.97714151463970571416;
    zj(116) =  0.98153114955374010687;
    zj(117) =  0.98537149959852037111;
    zj(118) =  0.98868475754742947994;
    zj(119) =  0.99149572117810613240;
    zj(120) =  0.99383196321275502221;
    zj(121) =  0.99572410469840718851;
    zj(122) =  0.99720625937222195908;
    zj(123) =  0.99831663531840739253;
    zj(124) =  0.99909812496766759766;
    zj(125) =  0.99959879967191068325;
    zj(126) =  0.99987288812035761194;
    zj(127) =  0.99998243035489159858;
    wj(  1) = 0.0000505360952078625176247;
    wj(  2) = 0.000180739564445388357820;
    wj(  3) = 0.000377746646326984660274;
    wj(  4) = 0.000632607319362633544219;
    wj(  5) = 0.000938369848542381500794;
    wj(  6) = 0.00128952408261041739210;
    wj(  7) = 0.00168114286542146990631;
    wj(  8) = 0.00210881524572663287933;
    wj(  9) = 0.00256876494379402037313;
    wj( 10) = 0.00305775341017553113613;
    wj( 11) = 0.00357289278351729964938;
    wj( 12) = 0.00411150397865469304717;
    wj( 13) = 0.00467105037211432174741;
    wj( 14) = 0.00524912345480885912513;
    wj( 15) = 0.00584344987583563950756;
    wj( 16) = 0.00645190005017573692280;
    wj( 17) = 0.00707248999543355546805;
    wj( 18) = 0.00770337523327974184817;
    wj( 19) = 0.00834283875396815770558;
    wj( 20) = 0.00898927578406413572328;
    wj( 21) = 0.00964117772970253669530;
    wj( 22) = 0.0102971169579563555237;
    wj( 23) = 0.0109557333878379016480;
    wj( 24) = 0.0116157233199551347270;
    wj( 25) = 0.0122758305600827700870;
    wj( 26) = 0.0129348396636073734547;
    wj( 27) = 0.0135915710097655467896;
    wj( 28) = 0.0142448773729167743063;
    wj( 29) = 0.0148936416648151820348;
    wj( 30) = 0.0155367755558439824399;
    wj( 31) = 0.0161732187295777199419;
    wj( 32) = 0.0168019385741038652709;
    wj( 33) = 0.0174219301594641737472;
    wj( 34) = 0.0180322163903912863201;
    wj( 35) = 0.0186318482561387901863;
    wj( 36) = 0.0192199051247277660193;
    wj( 37) = 0.0197954950480974994880;
    wj( 38) = 0.0203577550584721594669;
    wj( 39) = 0.0209058514458120238522;
    wj( 40) = 0.0214389800125038672465;
    wj( 41) = 0.0219563663053178249393;
    wj( 42) = 0.0224572658268160987071;
    wj( 43) = 0.0229409642293877487608;
    wj( 44) = 0.0234067774953140062013;
    wj( 45) = 0.0238540521060385400804;
    wj( 46) = 0.0242821652033365993580;
    wj( 47) = 0.0246905247444876769091;
    wj( 48) = 0.0250785696529497687068;
    wj( 49) = 0.0254457699654647658126;
    wj( 50) = 0.0257916269760242293884;
    wj( 51) = 0.0261156733767060976805;
    wj( 52) = 0.0264174733950582599310;
    wj( 53) = 0.0266966229274503599062;
    wj( 54) = 0.0269527496676330319634;
    wj( 55) = 0.0271855132296247918192;
    wj( 56) = 0.0273946052639814325161;
    wj( 57) = 0.0275797495664818730349;
    wj( 58) = 0.0277407021782796819939;
    wj( 59) = 0.0278772514766137016085;
    wj( 60) = 0.0279892182552381597038;
    wj( 61) = 0.0280764557938172466068;
    wj( 62) = 0.0281388499156271506363;
    wj( 63) = 0.0281763190330166021307;
    wj( 64) = 0.0281888141801923586938;
    wj( 65) = 0.0281763190330166021307;
    wj( 66) = 0.0281388499156271506363;
    wj( 67) = 0.0280764557938172466068;
    wj( 68) = 0.0279892182552381597038;
    wj( 69) = 0.0278772514766137016085;
    wj( 70) = 0.0277407021782796819939;
    wj( 71) = 0.0275797495664818730349;
    wj( 72) = 0.0273946052639814325161;
    wj( 73) = 0.0271855132296247918192;
    wj( 74) = 0.0269527496676330319634;
    wj( 75) = 0.0266966229274503599062;
    wj( 76) = 0.0264174733950582599310;
    wj( 77) = 0.0261156733767060976805;
    wj( 78) = 0.0257916269760242293884;
    wj( 79) = 0.0254457699654647658126;
    wj( 80) = 0.0250785696529497687068;
    wj( 81) = 0.0246905247444876769091;
    wj( 82) = 0.0242821652033365993580;
    wj( 83) = 0.0238540521060385400804;
    wj( 84) = 0.0234067774953140062013;
    wj( 85) = 0.0229409642293877487608;
    wj( 86) = 0.0224572658268160987071;
    wj( 87) = 0.0219563663053178249393;
    wj( 88) = 0.0214389800125038672465;
    wj( 89) = 0.0209058514458120238522;
    wj( 90) = 0.0203577550584721594669;
    wj( 91) = 0.0197954950480974994880;
    wj( 92) = 0.0192199051247277660193;
    wj( 93) = 0.0186318482561387901863;
    wj( 94) = 0.0180322163903912863201;
    wj( 95) = 0.0174219301594641737472;
    wj( 96) = 0.0168019385741038652709;
    wj( 97) = 0.0161732187295777199419;
    wj( 98) = 0.0155367755558439824399;
    wj( 99) = 0.0148936416648151820348;
    wj(100) = 0.0142448773729167743063;
    wj(101) = 0.0135915710097655467896;
    wj(102) = 0.0129348396636073734547;
    wj(103) = 0.0122758305600827700870;
    wj(104) = 0.0116157233199551347270;
    wj(105) = 0.0109557333878379016480;
    wj(106) = 0.0102971169579563555237;
    wj(107) = 0.00964117772970253669530;
    wj(108) = 0.00898927578406413572328;
    wj(109) = 0.00834283875396815770558;
    wj(110) = 0.00770337523327974184817;
    wj(111) = 0.00707248999543355546805;
    wj(112) = 0.00645190005017573692280;
    wj(113) = 0.00584344987583563950756;
    wj(114) = 0.00524912345480885912513;
    wj(115) = 0.00467105037211432174741;
    wj(116) = 0.00411150397865469304717;
    wj(117) = 0.00357289278351729964938;
    wj(118) = 0.00305775341017553113613;
    wj(119) = 0.00256876494379402037313;
    wj(120) = 0.00210881524572663287933;
    wj(121) = 0.00168114286542146990631;
    wj(122) = 0.00128952408261041739210;
    wj(123) = 0.000938369848542381500794;
    wj(124) = 0.000632607319362633544219;
    wj(125) = 0.000377746646326984660274;
    wj(126) = 0.000180739564445388357820;
    wj(127) = 0.0000505360952078625176247;
  case 255 
    zj(  1) = -0.99999759637974846462E+00;
    zj(  2) = -0.99998243035489159858E+00;
    zj(  3) = -0.99994399620705437576E+00;
    zj(  4) = -0.99987288812035761194E+00;
    zj(  5) = -0.99976049092443204733E+00;
    zj(  6) = -0.99959879967191068325E+00;
    zj(  7) = -0.99938033802502358193E+00;
    zj(  8) = -0.99909812496766759766E+00;
    zj(  9) = -0.99874561446809511470E+00;
    zj( 10) = -0.99831663531840739253E+00;
    zj( 11) = -0.99780535449595727456E+00;
    zj( 12) = -0.99720625937222195908E+00;
    zj( 13) = -0.99651414591489027385E+00;
    zj( 14) = -0.99572410469840718851E+00;
    zj( 15) = -0.99483150280062100052E+00;
    zj( 16) = -0.99383196321275502221E+00;
    zj( 17) = -0.99272134428278861533E+00;
    zj( 18) = -0.99149572117810613240E+00;
    zj( 19) = -0.99015137040077015918E+00;
    zj( 20) = -0.98868475754742947994E+00;
    zj( 21) = -0.98709252795403406719E+00;
    zj( 22) = -0.98537149959852037111E+00;
    zj( 23) = -0.98351865757863272876E+00;
    zj( 24) = -0.98153114955374010687E+00;
    zj( 25) = -0.97940628167086268381E+00;
    zj( 26) = -0.97714151463970571416E+00;
    zj( 27) = -0.97473445975240266776E+00;
    zj( 28) = -0.97218287474858179658E+00;
    zj( 29) = -0.96948465950245923177E+00;
    zj( 30) = -0.96663785155841656709E+00;
    zj( 31) = -0.96364062156981213252E+00;
    zj( 32) = -0.96049126870802028342E+00;
    zj( 33) = -0.95718821610986096274E+00;
    zj( 34) = -0.95373000642576113641E+00;
    zj( 35) = -0.95011529752129487656E+00;
    zj( 36) = -0.94634285837340290515E+00;
    zj( 37) = -0.94241156519108305981E+00;
    zj( 38) = -0.93832039777959288365E+00;
    zj( 39) = -0.93406843615772578800E+00;
    zj( 40) = -0.92965485742974005667E+00;
    zj( 41) = -0.92507893290707565236E+00;
    zj( 42) = -0.92034002547001242073E+00;
    zj( 43) = -0.91543758715576504064E+00;
    zj( 44) = -0.91037115695700429250E+00;
    zj( 45) = -0.90514035881326159519E+00;
    zj( 46) = -0.89974489977694003664E+00;
    zj( 47) = -0.89418456833555902286E+00;
    zj( 48) = -0.88845923287225699889E+00;
    zj( 49) = -0.88256884024734190684E+00;
    zj( 50) = -0.87651341448470526974E+00;
    zj( 51) = -0.87029305554811390585E+00;
    zj( 52) = -0.86390793819369047715E+00;
    zj( 53) = -0.85735831088623215653E+00;
    zj( 54) = -0.85064449476835027976E+00;
    zj( 55) = -0.84376688267270860104E+00;
    zj( 56) = -0.83672593816886873550E+00;
    zj( 57) = -0.82952219463740140018E+00;
    zj( 58) = -0.82215625436498040737E+00;
    zj( 59) = -0.81462878765513741344E+00;
    zj( 60) = -0.80694053195021761186E+00;
    zj( 61) = -0.79909229096084140180E+00;
    zj( 62) = -0.79108493379984836143E+00;
    zj( 63) = -0.78291939411828301639E+00;
    zj( 64) = -0.77459666924148337704E+00;
    zj( 65) = -0.76611781930376009072E+00;
    zj( 66) = -0.75748396638051363793E+00;
    zj( 67) = -0.74869629361693660282E+00;
    zj( 68) = -0.73975604435269475868E+00;
    zj( 69) = -0.73066452124218126133E+00;
    zj( 70) = -0.72142308537009891548E+00;
    zj( 71) = -0.71203315536225203459E+00;
    zj( 72) = -0.70249620649152707861E+00;
    zj( 73) = -0.69281376977911470289E+00;
    zj( 74) = -0.68298743109107922809E+00;
    zj( 75) = -0.67301883023041847920E+00;
    zj( 76) = -0.66290966002478059546E+00;
    zj( 77) = -0.65266166541001749610E+00;
    zj( 78) = -0.64227664250975951377E+00;
    zj( 79) = -0.63175643771119423041E+00;
    zj( 80) = -0.62110294673722640294E+00;
    zj( 81) = -0.61031811371518640016E+00;
    zj( 82) = -0.59940393024224289297E+00;
    zj( 83) = -0.58836243444766254143E+00;
    zj( 84) = -0.57719571005204581484E+00;
    zj( 85) = -0.56590588542365442262E+00;
    zj( 86) = -0.55449513263193254887E+00;
    zj( 87) = -0.54296566649831149049E+00;
    zj( 88) = -0.53131974364437562397E+00;
    zj( 89) = -0.51955966153745702199E+00;
    zj( 90) = -0.50768775753371660215E+00;
    zj( 91) = -0.49570640791876146017E+00;
    zj( 92) = -0.48361802694584102756E+00;
    zj( 93) = -0.47142506587165887693E+00;
    zj( 94) = -0.45913001198983233287E+00;
    zj( 95) = -0.44673538766202847374E+00;
    zj( 96) = -0.43424374934680255800E+00;
    zj( 97) = -0.42165768662616330006E+00;
    zj( 98) = -0.40897982122988867241E+00;
    zj( 99) = -0.39621280605761593918E+00;
    zj(100) = -0.38335932419873034692E+00;
    zj(101) = -0.37042208795007823014E+00;
    zj(102) = -0.35740383783153215238E+00;
    zj(103) = -0.34430734159943802278E+00;
    zj(104) = -0.33113539325797683309E+00;
    zj(105) = -0.31789081206847668318E+00;
    zj(106) = -0.30457644155671404334E+00;
    zj(107) = -0.29119514851824668196E+00;
    zj(108) = -0.27774982202182431507E+00;
    zj(109) = -0.26424337241092676194E+00;
    zj(110) = -0.25067873030348317661E+00;
    zj(111) = -0.23705884558982972721E+00;
    zj(112) = -0.22338668642896688163E+00;
    zj(113) = -0.20966523824318119477E+00;
    zj(114) = -0.19589750271110015392E+00;
    zj(115) = -0.18208649675925219825E+00;
    zj(116) = -0.16823525155220746498E+00;
    zj(117) = -0.15434681148137810869E+00;
    zj(118) = -0.14042423315256017459E+00;
    zj(119) = -0.12647058437230196685E+00;
    zj(120) = -0.11248894313318662575E+00;
    zj(121) = -0.098482396598119202090E+00;
    zj(122) = -0.084454040083710883710E+00;
    zj(123) = -0.070406976042855179063E+00;
    zj(124) = -0.056344313046592789972E+00;
    zj(125) = -0.042269164765363603212E+00;
    zj(126) = -0.028184648949745694339E+00;
    zj(127) = -0.014093886410782462614E+00;
    zj(128) =  0.0E+00;
    zj(129) =  0.014093886410782462614E+00;
    zj(130) =  0.028184648949745694339E+00;
    zj(131) =  0.042269164765363603212E+00;
    zj(132) =  0.056344313046592789972E+00;
    zj(133) =  0.070406976042855179063E+00;
    zj(134) =  0.084454040083710883710E+00;
    zj(135) =  0.098482396598119202090E+00;
    zj(136) =  0.11248894313318662575E+00;
    zj(137) =  0.12647058437230196685E+00;
    zj(138) =  0.14042423315256017459E+00;
    zj(139) =  0.15434681148137810869E+00;
    zj(140) =  0.16823525155220746498E+00;
    zj(141) =  0.18208649675925219825E+00;
    zj(142) =  0.19589750271110015392E+00;
    zj(143) =  0.20966523824318119477E+00;
    zj(144) =  0.22338668642896688163E+00;
    zj(145) =  0.23705884558982972721E+00;
    zj(146) =  0.25067873030348317661E+00;
    zj(147) =  0.26424337241092676194E+00;
    zj(148) =  0.27774982202182431507E+00;
    zj(149) =  0.29119514851824668196E+00;
    zj(150) =  0.30457644155671404334E+00;
    zj(151) =  0.31789081206847668318E+00;
    zj(152) =  0.33113539325797683309E+00;
    zj(153) =  0.34430734159943802278E+00;
    zj(154) =  0.35740383783153215238E+00;
    zj(155) =  0.37042208795007823014E+00;
    zj(156) =  0.38335932419873034692E+00;
    zj(157) =  0.39621280605761593918E+00;
    zj(158) =  0.40897982122988867241E+00;
    zj(159) =  0.42165768662616330006E+00;
    zj(160) =  0.43424374934680255800E+00;
    zj(161) =  0.44673538766202847374E+00;
    zj(162) =  0.45913001198983233287E+00;
    zj(163) =  0.47142506587165887693E+00;
    zj(164) =  0.48361802694584102756E+00;
    zj(165) =  0.49570640791876146017E+00;
    zj(166) =  0.50768775753371660215E+00;
    zj(167) =  0.51955966153745702199E+00;
    zj(168) =  0.53131974364437562397E+00;
    zj(169) =  0.54296566649831149049E+00;
    zj(170) =  0.55449513263193254887E+00;
    zj(171) =  0.56590588542365442262E+00;
    zj(172) =  0.57719571005204581484E+00;
    zj(173) =  0.58836243444766254143E+00;
    zj(174) =  0.59940393024224289297E+00;
    zj(175) =  0.61031811371518640016E+00;
    zj(176) =  0.62110294673722640294E+00;
    zj(177) =  0.63175643771119423041E+00;
    zj(178) =  0.64227664250975951377E+00;
    zj(179) =  0.65266166541001749610E+00;
    zj(180) =  0.66290966002478059546E+00;
    zj(181) =  0.67301883023041847920E+00;
    zj(182) =  0.68298743109107922809E+00;
    zj(183) =  0.69281376977911470289E+00;
    zj(184) =  0.70249620649152707861E+00;
    zj(185) =  0.71203315536225203459E+00;
    zj(186) =  0.72142308537009891548E+00;
    zj(187) =  0.73066452124218126133E+00;
    zj(188) =  0.73975604435269475868E+00;
    zj(189) =  0.74869629361693660282E+00;
    zj(190) =  0.75748396638051363793E+00;
    zj(191) =  0.76611781930376009072E+00;
    zj(192) =  0.77459666924148337704E+00;
    zj(193) =  0.78291939411828301639E+00;
    zj(194) =  0.79108493379984836143E+00;
    zj(195) =  0.79909229096084140180E+00;
    zj(196) =  0.80694053195021761186E+00;
    zj(197) =  0.81462878765513741344E+00;
    zj(198) =  0.82215625436498040737E+00;
    zj(199) =  0.82952219463740140018E+00;
    zj(200) =  0.83672593816886873550E+00;
    zj(201) =  0.84376688267270860104E+00;
    zj(202) =  0.85064449476835027976E+00;
    zj(203) =  0.85735831088623215653E+00;
    zj(204) =  0.86390793819369047715E+00;
    zj(205) =  0.87029305554811390585E+00;
    zj(206) =  0.87651341448470526974E+00;
    zj(207) =  0.88256884024734190684E+00;
    zj(208) =  0.88845923287225699889E+00;
    zj(209) =  0.89418456833555902286E+00;
    zj(210) =  0.89974489977694003664E+00;
    zj(211) =  0.90514035881326159519E+00;
    zj(212) =  0.91037115695700429250E+00;
    zj(213) =  0.91543758715576504064E+00;
    zj(214) =  0.92034002547001242073E+00;
    zj(215) =  0.92507893290707565236E+00;
    zj(216) =  0.92965485742974005667E+00;
    zj(217) =  0.93406843615772578800E+00;
    zj(218) =  0.93832039777959288365E+00;
    zj(219) =  0.94241156519108305981E+00;
    zj(220) =  0.94634285837340290515E+00;
    zj(221) =  0.95011529752129487656E+00;
    zj(222) =  0.95373000642576113641E+00;
    zj(223) =  0.95718821610986096274E+00;
    zj(224) =  0.96049126870802028342E+00;
    zj(225) =  0.96364062156981213252E+00;
    zj(226) =  0.96663785155841656709E+00;
    zj(227) =  0.96948465950245923177E+00;
    zj(228) =  0.97218287474858179658E+00;
    zj(229) =  0.97473445975240266776E+00;
    zj(230) =  0.97714151463970571416E+00;
    zj(231) =  0.97940628167086268381E+00;
    zj(232) =  0.98153114955374010687E+00;
    zj(233) =  0.98351865757863272876E+00;
    zj(234) =  0.98537149959852037111E+00;
    zj(235) =  0.98709252795403406719E+00;
    zj(236) =  0.98868475754742947994E+00;
    zj(237) =  0.99015137040077015918E+00;
    zj(238) =  0.99149572117810613240E+00;
    zj(239) =  0.99272134428278861533E+00;
    zj(240) =  0.99383196321275502221E+00;
    zj(241) =  0.99483150280062100052E+00;
    zj(242) =  0.99572410469840718851E+00;
    zj(243) =  0.99651414591489027385E+00;
    zj(244) =  0.99720625937222195908E+00;
    zj(245) =  0.99780535449595727456E+00;
    zj(246) =  0.99831663531840739253E+00;
    zj(247) =  0.99874561446809511470E+00;
    zj(248) =  0.99909812496766759766E+00;
    zj(249) =  0.99938033802502358193E+00;
    zj(250) =  0.99959879967191068325E+00;
    zj(251) =  0.99976049092443204733E+00;
    zj(252) =  0.99987288812035761194E+00;
    zj(253) =  0.99994399620705437576E+00;
    zj(254) =  0.99998243035489159858E+00;
    zj(255) =  0.99999759637974846462E+00;
    wj(  1) = 0.69379364324108267170E-05;
    wj(  2) = 0.25157870384280661489E-04;
    wj(  3) = 0.53275293669780613125E-04;
    wj(  4) = 0.90372734658751149261E-04;
    wj(  5) = 0.13575491094922871973E-03;
    wj(  6) = 0.18887326450650491366E-03;
    wj(  7) = 0.24921240048299729402E-03;
    wj(  8) = 0.31630366082226447689E-03;
    wj(  9) = 0.38974528447328229322E-03;
    wj( 10) = 0.46918492424785040975E-03;
    wj( 11) = 0.55429531493037471492E-03;
    wj( 12) = 0.64476204130572477933E-03;
    wj( 13) = 0.74028280424450333046E-03;
    wj( 14) = 0.84057143271072246365E-03;
    wj( 15) = 0.94536151685852538246E-03;
    wj( 16) = 0.10544076228633167722E-02;
    wj( 17) = 0.11674841174299594077E-02;
    wj( 18) = 0.12843824718970101768E-02;
    wj( 19) = 0.14049079956551446427E-02;
    wj( 20) = 0.15288767050877655684E-02;
    wj( 21) = 0.16561127281544526052E-02;
    wj( 22) = 0.17864463917586498247E-02;
    wj( 23) = 0.19197129710138724125E-02;
    wj( 24) = 0.20557519893273465236E-02;
    wj( 25) = 0.21944069253638388388E-02;
    wj( 26) = 0.23355251860571608737E-02;
    wj( 27) = 0.24789582266575679307E-02;
    wj( 28) = 0.26245617274044295626E-02;
    wj( 29) = 0.27721957645934509940E-02;
    wj( 30) = 0.29217249379178197538E-02;
    wj( 31) = 0.30730184347025783234E-02;
    wj( 32) = 0.32259500250878684614E-02;
    wj( 33) = 0.33803979910869203823E-02;
    wj( 34) = 0.35362449977167777340E-02;
    wj( 35) = 0.36933779170256508183E-02;
    wj( 36) = 0.38516876166398709241E-02;
    wj( 37) = 0.40110687240750233989E-02;
    wj( 38) = 0.41714193769840788528E-02;
    wj( 39) = 0.43326409680929828545E-02;
    wj( 40) = 0.44946378920320678616E-02;
    wj( 41) = 0.46573172997568547773E-02;
    wj( 42) = 0.48205888648512683476E-02;
    wj( 43) = 0.49843645647655386012E-02;
    wj( 44) = 0.51485584789781777618E-02;
    wj( 45) = 0.53130866051870565663E-02;
    wj( 46) = 0.54778666939189508240E-02;
    wj( 47) = 0.56428181013844441585E-02;
    wj( 48) = 0.58078616599775673635E-02;
    wj( 49) = 0.59729195655081658049E-02;
    wj( 50) = 0.61379152800413850435E-02;
    wj( 51) = 0.63027734490857587172E-02;
    wj( 52) = 0.64674198318036867274E-02;
    wj( 53) = 0.66317812429018878941E-02;
    wj( 54) = 0.67957855048827733948E-02;
    wj( 55) = 0.69593614093904229394E-02;
    wj( 56) = 0.71224386864583871532E-02;
    wj( 57) = 0.72849479805538070639E-02;
    wj( 58) = 0.74468208324075910174E-02;
    wj( 59) = 0.76079896657190565832E-02;
    wj( 60) = 0.77683877779219912200E-02;
    wj( 61) = 0.79279493342948491103E-02;
    wj( 62) = 0.80866093647888599710E-02;
    wj( 63) = 0.82443037630328680306E-02;
    wj( 64) = 0.84009692870519326354E-02;
    wj( 65) = 0.85565435613076896192E-02;
    wj( 66) = 0.87109650797320868736E-02;
    wj( 67) = 0.88641732094824942641E-02;
    wj( 68) = 0.90161081951956431600E-02;
    wj( 69) = 0.91667111635607884067E-02;
    wj( 70) = 0.93159241280693950932E-02;
    wj( 71) = 0.94636899938300652943E-02;
    wj( 72) = 0.96099525623638830097E-02;
    wj( 73) = 0.97546565363174114611E-02;
    wj( 74) = 0.98977475240487497440E-02;
    wj( 75) = 0.10039172044056840798E-01;
    wj( 76) = 0.10178877529236079733E-01;
    wj( 77) = 0.10316812330947621682E-01;
    wj( 78) = 0.10452925722906011926E-01;
    wj( 79) = 0.10587167904885197931E-01;
    wj( 80) = 0.10719490006251933623E-01;
    wj( 81) = 0.10849844089337314099E-01;
    wj( 82) = 0.10978183152658912470E-01;
    wj( 83) = 0.11104461134006926537E-01;
    wj( 84) = 0.11228632913408049354E-01;
    wj( 85) = 0.11350654315980596602E-01;
    wj( 86) = 0.11470482114693874380E-01;
    wj( 87) = 0.11588074033043952568E-01;
    wj( 88) = 0.11703388747657003101E-01;
    wj( 89) = 0.11816385890830235763E-01;
    wj( 90) = 0.11927026053019270040E-01;
    wj( 91) = 0.12035270785279562630E-01;
    wj( 92) = 0.12141082601668299679E-01;
    wj( 93) = 0.12244424981611985899E-01;
    wj( 94) = 0.12345262372243838455E-01;
    wj( 95) = 0.12443560190714035263E-01;
    wj( 96) = 0.12539284826474884353E-01;
    wj( 97) = 0.12632403643542078765E-01;
    wj( 98) = 0.12722884982732382906E-01;
    wj( 99) = 0.12810698163877361967E-01;
    wj(100) = 0.12895813488012114694E-01;
    wj(101) = 0.12978202239537399286E-01;
    wj(102) = 0.13057836688353048840E-01;
    wj(103) = 0.13134690091960152836E-01;
    wj(104) = 0.13208736697529129966E-01;
    wj(105) = 0.13279951743930530650E-01;
    wj(106) = 0.13348311463725179953E-01;
    wj(107) = 0.13413793085110098513E-01;
    wj(108) = 0.13476374833816515982E-01;
    wj(109) = 0.13536035934956213614E-01;
    wj(110) = 0.13592756614812395910E-01;
    wj(111) = 0.13646518102571291428E-01;
    wj(112) = 0.13697302631990716258E-01;
    wj(113) = 0.13745093443001896632E-01;
    wj(114) = 0.13789874783240936517E-01;
    wj(115) = 0.13831631909506428676E-01;
    wj(116) = 0.13870351089139840997E-01;
    wj(117) = 0.13906019601325461264E-01;
    wj(118) = 0.13938625738306850804E-01;
    wj(119) = 0.13968158806516938516E-01;
    wj(120) = 0.13994609127619079852E-01;
    wj(121) = 0.14017968039456608810E-01;
    wj(122) = 0.14038227896908623303E-01;
    wj(123) = 0.14055382072649964277E-01;
    wj(124) = 0.14069424957813575318E-01;
    wj(125) = 0.14080351962553661325E-01;
    wj(126) = 0.14088159516508301065E-01;
    wj(127) = 0.14092845069160408355E-01;
    wj(128) = 0.14094407090096179347E-01;
    wj(129) = 0.14092845069160408355E-01;
    wj(130) = 0.14088159516508301065E-01;
    wj(131) = 0.14080351962553661325E-01;
    wj(132) = 0.14069424957813575318E-01;
    wj(133) = 0.14055382072649964277E-01;
    wj(134) = 0.14038227896908623303E-01;
    wj(135) = 0.14017968039456608810E-01;
    wj(136) = 0.13994609127619079852E-01;
    wj(137) = 0.13968158806516938516E-01;
    wj(138) = 0.13938625738306850804E-01;
    wj(139) = 0.13906019601325461264E-01;
    wj(140) = 0.13870351089139840997E-01;
    wj(141) = 0.13831631909506428676E-01;
    wj(142) = 0.13789874783240936517E-01;
    wj(143) = 0.13745093443001896632E-01;
    wj(144) = 0.13697302631990716258E-01;
    wj(145) = 0.13646518102571291428E-01;
    wj(146) = 0.13592756614812395910E-01;
    wj(147) = 0.13536035934956213614E-01;
    wj(148) = 0.13476374833816515982E-01;
    wj(149) = 0.13413793085110098513E-01;
    wj(150) = 0.13348311463725179953E-01;
    wj(151) = 0.13279951743930530650E-01;
    wj(152) = 0.13208736697529129966E-01;
    wj(153) = 0.13134690091960152836E-01;
    wj(154) = 0.13057836688353048840E-01;
    wj(155) = 0.12978202239537399286E-01;
    wj(156) = 0.12895813488012114694E-01;
    wj(157) = 0.12810698163877361967E-01;
    wj(158) = 0.12722884982732382906E-01;
    wj(159) = 0.12632403643542078765E-01;
    wj(160) = 0.12539284826474884353E-01;
    wj(161) = 0.12443560190714035263E-01;
    wj(162) = 0.12345262372243838455E-01;
    wj(163) = 0.12244424981611985899E-01;
    wj(164) = 0.12141082601668299679E-01;
    wj(165) = 0.12035270785279562630E-01;
    wj(166) = 0.11927026053019270040E-01;
    wj(167) = 0.11816385890830235763E-01;
    wj(168) = 0.11703388747657003101E-01;
    wj(169) = 0.11588074033043952568E-01;
    wj(170) = 0.11470482114693874380E-01;
    wj(171) = 0.11350654315980596602E-01;
    wj(172) = 0.11228632913408049354E-01;
    wj(173) = 0.11104461134006926537E-01;
    wj(174) = 0.10978183152658912470E-01;
    wj(175) = 0.10849844089337314099E-01;
    wj(176) = 0.10719490006251933623E-01;
    wj(177) = 0.10587167904885197931E-01;
    wj(178) = 0.10452925722906011926E-01;
    wj(179) = 0.10316812330947621682E-01;
    wj(180) = 0.10178877529236079733E-01;
    wj(181) = 0.10039172044056840798E-01;
    wj(182) = 0.98977475240487497440E-02;
    wj(183) = 0.97546565363174114611E-02;
    wj(184) = 0.96099525623638830097E-02;
    wj(185) = 0.94636899938300652943E-02;
    wj(186) = 0.93159241280693950932E-02;
    wj(187) = 0.91667111635607884067E-02;
    wj(188) = 0.90161081951956431600E-02;
    wj(189) = 0.88641732094824942641E-02;
    wj(190) = 0.87109650797320868736E-02;
    wj(191) = 0.85565435613076896192E-02;
    wj(192) = 0.84009692870519326354E-02;
    wj(193) = 0.82443037630328680306E-02;
    wj(194) = 0.80866093647888599710E-02;
    wj(195) = 0.79279493342948491103E-02;
    wj(196) = 0.77683877779219912200E-02;
    wj(197) = 0.76079896657190565832E-02;
    wj(198) = 0.74468208324075910174E-02;
    wj(199) = 0.72849479805538070639E-02;
    wj(200) = 0.71224386864583871532E-02;
    wj(201) = 0.69593614093904229394E-02;
    wj(202) = 0.67957855048827733948E-02;
    wj(203) = 0.66317812429018878941E-02;
    wj(204) = 0.64674198318036867274E-02;
    wj(205) = 0.63027734490857587172E-02;
    wj(206) = 0.61379152800413850435E-02;
    wj(207) = 0.59729195655081658049E-02;
    wj(208) = 0.58078616599775673635E-02;
    wj(209) = 0.56428181013844441585E-02;
    wj(210) = 0.54778666939189508240E-02;
    wj(211) = 0.53130866051870565663E-02;
    wj(212) = 0.51485584789781777618E-02;
    wj(213) = 0.49843645647655386012E-02;
    wj(214) = 0.48205888648512683476E-02;
    wj(215) = 0.46573172997568547773E-02;
    wj(216) = 0.44946378920320678616E-02;
    wj(217) = 0.43326409680929828545E-02;
    wj(218) = 0.41714193769840788528E-02;
    wj(219) = 0.40110687240750233989E-02;
    wj(220) = 0.38516876166398709241E-02;
    wj(221) = 0.36933779170256508183E-02;
    wj(222) = 0.35362449977167777340E-02;
    wj(223) = 0.33803979910869203823E-02;
    wj(224) = 0.32259500250878684614E-02;
    wj(225) = 0.30730184347025783234E-02;
    wj(226) = 0.29217249379178197538E-02;
    wj(227) = 0.27721957645934509940E-02;
    wj(228) = 0.26245617274044295626E-02;
    wj(229) = 0.24789582266575679307E-02;
    wj(230) = 0.23355251860571608737E-02;
    wj(231) = 0.21944069253638388388E-02;
    wj(232) = 0.20557519893273465236E-02;
    wj(233) = 0.19197129710138724125E-02;
    wj(234) = 0.17864463917586498247E-02;
    wj(235) = 0.16561127281544526052E-02;
    wj(236) = 0.15288767050877655684E-02;
    wj(237) = 0.14049079956551446427E-02;
    wj(238) = 0.12843824718970101768E-02;
    wj(239) = 0.11674841174299594077E-02;
    wj(240) = 0.10544076228633167722E-02;
    wj(241) = 0.94536151685852538246E-03;
    wj(242) = 0.84057143271072246365E-03;
    wj(243) = 0.74028280424450333046E-03;
    wj(244) = 0.64476204130572477933E-03;
    wj(245) = 0.55429531493037471492E-03;
    wj(246) = 0.46918492424785040975E-03;
    wj(247) = 0.38974528447328229322E-03;
    wj(248) = 0.31630366082226447689E-03;
    wj(249) = 0.24921240048299729402E-03;
    wj(250) = 0.18887326450650491366E-03;
    wj(251) = 0.13575491094922871973E-03;
    wj(252) = 0.90372734658751149261E-04;
    wj(253) = 0.53275293669780613125E-04;
    wj(254) = 0.25157870384280661489E-04;
    wj(255) = 0.69379364324108267170E-05;
  case 511
    zj(  1) = -0.999999672956734384381;
    zj(  2) = -0.999997596379748464620;
    zj(  3) = -0.999992298136257588028;
    zj(  4) = -0.999982430354891598580;
    zj(  5) = -0.999966730098486276883;
    zj(  6) = -0.999943996207054375764;
    zj(  7) = -0.999913081144678282800;
    zj(  8) = -0.999872888120357611938;
    zj(  9) = -0.999822363679787739196;
    zj( 10) = -0.999760490924432047330;
    zj( 11) = -0.999686286448317731776;
    zj( 12) = -0.999598799671910683252;
    zj( 13) = -0.999497112467187190535;
    zj( 14) = -0.999380338025023581928;
    zj( 15) = -0.999247618943342473599;
    zj( 16) = -0.999098124967667597662;
    zj( 17) = -0.998931050830810562236;
    zj( 18) = -0.998745614468095114704;
    zj( 19) = -0.998541055697167906027;
    zj( 20) = -0.998316635318407392531;
    zj( 21) = -0.998071634524930323302;
    zj( 22) = -0.997805354495957274562;
    zj( 23) = -0.997517116063472399965;
    zj( 24) = -0.997206259372221959076;
    zj( 25) = -0.996872143485260161299;
    zj( 26) = -0.996514145914890273849;
    zj( 27) = -0.996131662079315037786;
    zj( 28) = -0.995724104698407188509;
    zj( 29) = -0.995290903148810302261;
    zj( 30) = -0.994831502800621000519;
    zj( 31) = -0.994345364356723405931;
    zj( 32) = -0.993831963212755022209;
    zj( 33) = -0.993290788851684966211;
    zj( 34) = -0.992721344282788615328;
    zj( 35) = -0.992123145530863117683;
    zj( 36) = -0.991495721178106132399;
    zj( 37) = -0.990838611958294243677;
    zj( 38) = -0.990151370400770159181;
    zj( 39) = -0.989433560520240838716;
    zj( 40) = -0.988684757547429479939;
    zj( 41) = -0.987904547695124280467;
    zj( 42) = -0.987092527954034067190;
    zj( 43) = -0.986248305913007552681;
    zj( 44) = -0.985371499598520371114;
    zj( 45) = -0.984461737328814534596;
    zj( 46) = -0.983518657578632728762;
    zj( 47) = -0.982541908851080604251;
    zj( 48) = -0.981531149553740106867;
    zj( 49) = -0.980486047876721339416;
    zj( 50) = -0.979406281670862683806;
    zj( 51) = -0.978291538324758539526;
    zj( 52) = -0.977141514639705714156;
    zj( 53) = -0.975955916702011753129;
    zj( 54) = -0.974734459752402667761;
    zj( 55) = -0.973476868052506926773;
    zj( 56) = -0.972182874748581796578;
    zj( 57) = -0.970852221732792443256;
    zj( 58) = -0.969484659502459231771;
    zj( 59) = -0.968079947017759947964;
    zj( 60) = -0.966637851558416567092;
    zj( 61) = -0.965158148579915665979;
    zj( 62) = -0.963640621569812132521;
    zj( 63) = -0.962085061904651475741;
    zj( 64) = -0.960491268708020283423;
    zj( 65) = -0.958859048710200221356;
    zj( 66) = -0.957188216109860962736;
    zj( 67) = -0.955478592438183697574;
    zj( 68) = -0.953730006425761136415;
    zj( 69) = -0.951942293872573589498;
    zj( 70) = -0.950115297521294876558;
    zj( 71) = -0.948248866934137357063;
    zj( 72) = -0.946342858373402905148;
    zj( 73) = -0.944397134685866648591;
    zj( 74) = -0.942411565191083059813;
    zj( 75) = -0.940386025573669721370;
    zj( 76) = -0.938320397779592883655;
    zj( 77) = -0.936214569916450806625;
    zj( 78) = -0.934068436157725787999;
    zj( 79) = -0.931881896650953639345;
    zj( 80) = -0.929654857429740056670;
    zj( 81) = -0.927387230329536696843;
    zj( 82) = -0.925078932907075652364;
    zj( 83) = -0.922729888363349241523;
    zj( 84) = -0.920340025470012420730;
    zj( 85) = -0.917909278499077501636;
    zj( 86) = -0.915437587155765040644;
    zj( 87) = -0.912924896514370590080;
    zj( 88) = -0.910371156957004292498;
    zj( 89) = -0.907776324115058903624;
    zj( 90) = -0.905140358813261595189;
    zj( 91) = -0.902463227016165675048;
    zj( 92) = -0.899744899776940036639;
    zj( 93) = -0.896985353188316590376;
    zj( 94) = -0.894184568335559022859;
    zj( 95) = -0.891342531251319871666;
    zj( 96) = -0.888459232872256998890;
    zj( 97) = -0.885534668997285008926;
    zj( 98) = -0.882568840247341906842;
    zj( 99) = -0.879561752026556262568;
    zj(100) = -0.876513414484705269742;
    zj(101) = -0.873423842480859310192;
    zj(102) = -0.870293055548113905851;
    zj(103) = -0.867121077859315215614;
    zj(104) = -0.863907938193690477146;
    zj(105) = -0.860653669904299969802;
    zj(106) = -0.857358310886232156525;
    zj(107) = -0.854021903545468625813;
    zj(108) = -0.850644494768350279758;
    zj(109) = -0.847226135891580884381;
    zj(110) = -0.843766882672708601038;
    zj(111) = -0.840266795261030442350;
    zj(112) = -0.836725938168868735503;
    zj(113) = -0.833144380243172624728;
    zj(114) = -0.829522194637401400178;
    zj(115) = -0.825859458783650001088;
    zj(116) = -0.822156254364980407373;
    zj(117) = -0.818412667287925807395;
    zj(118) = -0.814628787655137413436;
    zj(119) = -0.810804709738146594361;
    zj(120) = -0.806940531950217611856;
    zj(121) = -0.803036356819268687782;
    zj(122) = -0.799092290960841401800;
    zj(123) = -0.795108445051100526780;
    zj(124) = -0.791084933799848361435;
    zj(125) = -0.787021875923539422170;
    zj(126) = -0.782919394118283016385;
    zj(127) = -0.778777615032822744702;
    zj(128) = -0.774596669241483377036;
    zj(129) = -0.770376691217076824278;
    zj(130) = -0.766117819303760090717;
    zj(131) = -0.761820195689839149173;
    zj(132) = -0.757483966380513637926;
    zj(133) = -0.753109281170558142523;
    zj(134) = -0.748696293616936602823;
    zj(135) = -0.744245161011347082309;
    zj(136) = -0.739756044352694758677;
    zj(137) = -0.735229108319491547663;
    zj(138) = -0.730664521242181261329;
    zj(139) = -0.726062455075389632685;
    zj(140) = -0.721423085370098915485;
    zj(141) = -0.716746591245747095767;
    zj(142) = -0.712033155362252034587;
    zj(143) = -0.707282963891961103412;
    zj(144) = -0.702496206491527078610;
    zj(145) = -0.697673076273711232906;
    zj(146) = -0.692813769779114702895;
    zj(147) = -0.687918486947839325756;
    zj(148) = -0.682987431091079228087;
    zj(149) = -0.678020808862644517838;
    zj(150) = -0.673018830230418479199;
    zj(151) = -0.667981708447749702165;
    zj(152) = -0.662909660024780595461;
    zj(153) = -0.657802904699713735422;
    zj(154) = -0.652661665410017496101;
    zj(155) = -0.647486168263572388782;
    zj(156) = -0.642276642509759513774;
    zj(157) = -0.637033320510492495071;
    zj(158) = -0.631756437711194230414;
    zj(159) = -0.626446232611719746542;
    zj(160) = -0.621102946737226402941;
    zj(161) = -0.615726824608992638014;
    zj(162) = -0.610318113715186400156;
    zj(163) = -0.604877064481584353319;
    zj(164) = -0.599403930242242892974;
    zj(165) = -0.593898967210121954393;
    zj(166) = -0.588362434447662541434;
    zj(167) = -0.582794593837318850840;
    zj(168) = -0.577195710052045814844;
    zj(169) = -0.571566050525742833992;
    zj(170) = -0.565905885423654422623;
    zj(171) = -0.560215487612728441818;
    zj(172) = -0.554495132631932548866;
    zj(173) = -0.548745098662529448608;
    zj(174) = -0.542965666498311490492;
    zj(175) = -0.537157119515795115982;
    zj(176) = -0.531319743644375623972;
    zj(177) = -0.525453827336442687395;
    zj(178) = -0.519559661537457021993;
    zj(179) = -0.513637539655988578507;
    zj(180) = -0.507687757533716602155;
    zj(181) = -0.501710613415391878251;
    zj(182) = -0.495706407918761460170;
    zj(183) = -0.489675444004456155436;
    zj(184) = -0.483618026945841027562;
    zj(185) = -0.477534464298829155284;
    zj(186) = -0.471425065871658876934;
    zj(187) = -0.465290143694634735858;
    zj(188) = -0.459130011989832332874;
    zj(189) = -0.452944987140767283784;
    zj(190) = -0.446735387662028473742;
    zj(191) = -0.440501534168875795783;
    zj(192) = -0.434243749346802558002;
    zj(193) = -0.427962357921062742583;
    zj(194) = -0.421657686626163300056;
    zj(195) = -0.415330064175321663764;
    zj(196) = -0.408979821229888672409;
    zj(197) = -0.402607290368737092671;
    zj(198) = -0.396212806057615939183;
    zj(199) = -0.389796704618470795479;
    zj(200) = -0.383359324198730346916;
    zj(201) = -0.376901004740559344802;
    zj(202) = -0.370422087950078230138;
    zj(203) = -0.363922917266549655269;
    zj(204) = -0.357403837831532152376;
    zj(205) = -0.350865196458001209011;
    zj(206) = -0.344307341599438022777;
    zj(207) = -0.337730623318886219621;
    zj(208) = -0.331135393257976833093;
    zj(209) = -0.324522004605921855207;
    zj(210) = -0.317890812068476683182;
    zj(211) = -0.311242171836871800300;
    zj(212) = -0.304576441556714043335;
    zj(213) = -0.297893980296857823437;
    zj(214) = -0.291195148518246681964;
    zj(215) = -0.284480308042725577496;
    zj(216) = -0.277749822021824315065;
    zj(217) = -0.271004054905512543536;
    zj(218) = -0.264243372410926761945;
    zj(219) = -0.257468141491069790481;
    zj(220) = -0.250678730303483176613;
    zj(221) = -0.243875508178893021593;
    zj(222) = -0.237058845589829727213;
    zj(223) = -0.230229114119222177156;
    zj(224) = -0.223386686428966881628;
    zj(225) = -0.216531936228472628081;
    zj(226) = -0.209665238243181194766;
    zj(227) = -0.202786968183064697557;
    zj(228) = -0.195897502711100153915;
    zj(229) = -0.188997219411721861059;
    zj(230) = -0.182086496759252198246;
    zj(231) = -0.175165714086311475707;
    zj(232) = -0.168235251552207464982;
    zj(233) = -0.161295490111305257361;
    zj(234) = -0.154346811481378108692;
    zj(235) = -0.147389598111939940054;
    zj(236) = -0.140424233152560174594;
    zj(237) = -0.133451100421161601344;
    zj(238) = -0.126470584372301966851;
    zj(239) = -0.119483070065440005133;
    zj(240) = -0.112488943133186625746;
    zj(241) = -0.105488589749541988533;
    zj(242) = -0.984823965981192020903E-01;
    zj(243) = -0.914707508403553909095E-01;
    zj(244) = -0.844540400837108837102E-01;
    zj(245) = -0.774326523498572825675E-01;
    zj(246) = -0.704069760428551790633E-01;
    zj(247) = -0.633773999173222898797E-01;
    zj(248) = -0.563443130465927899720E-01;
    zj(249) = -0.493081047908686267156E-01;
    zj(250) = -0.422691647653636032124E-01;
    zj(251) = -0.352278828084410232603E-01;
    zj(252) = -0.281846489497456943394E-01;
    zj(253) = -0.211398533783310883350E-01;
    zj(254) = -0.140938864107824626142E-01;
    zj(255) = -0.704713845933674648514E-02;
    zj(256) = +0.000000000000000000000;
    zj(257) = +0.704713845933674648514E-02;
    zj(258) = +0.140938864107824626142E-01;
    zj(259) = +0.211398533783310883350E-01;
    zj(260) = +0.281846489497456943394E-01;
    zj(261) = +0.352278828084410232603E-01;
    zj(262) = +0.422691647653636032124E-01;
    zj(263) = +0.493081047908686267156E-01;
    zj(264) = +0.563443130465927899720E-01;
    zj(265) = +0.633773999173222898797E-01;
    zj(266) = +0.704069760428551790633E-01;
    zj(267) = +0.774326523498572825675E-01;
    zj(268) = +0.844540400837108837102E-01;
    zj(269) = +0.914707508403553909095E-01;
    zj(270) = +0.984823965981192020903E-01;
    zj(271) = +0.105488589749541988533;
    zj(272) = +0.112488943133186625746;
    zj(273) = +0.119483070065440005133;
    zj(274) = +0.126470584372301966851;
    zj(275) = +0.133451100421161601344;
    zj(276) = +0.140424233152560174594;
    zj(277) = +0.147389598111939940054;
    zj(278) = +0.154346811481378108692;
    zj(279) = +0.161295490111305257361;
    zj(280) = +0.168235251552207464982;
    zj(281) = +0.175165714086311475707;
    zj(282) = +0.182086496759252198246;
    zj(283) = +0.188997219411721861059;
    zj(284) = +0.195897502711100153915;
    zj(285) = +0.202786968183064697557;
    zj(286) = +0.209665238243181194766;
    zj(287) = +0.216531936228472628081;
    zj(288) = +0.223386686428966881628;
    zj(289) = +0.230229114119222177156;
    zj(290) = +0.237058845589829727213;
    zj(291) = +0.243875508178893021593;
    zj(292) = +0.250678730303483176613;
    zj(293) = +0.257468141491069790481;
    zj(294) = +0.264243372410926761945;
    zj(295) = +0.271004054905512543536;
    zj(296) = +0.277749822021824315065;
    zj(297) = +0.284480308042725577496;
    zj(298) = +0.291195148518246681964;
    zj(299) = +0.297893980296857823437;
    zj(300) = +0.304576441556714043335;
    zj(301) = +0.311242171836871800300;
    zj(302) = +0.317890812068476683182;
    zj(303) = +0.324522004605921855207;
    zj(304) = +0.331135393257976833093;
    zj(305) = +0.337730623318886219621;
    zj(306) = +0.344307341599438022777;
    zj(307) = +0.350865196458001209011;
    zj(308) = +0.357403837831532152376;
    zj(309) = +0.363922917266549655269;
    zj(310) = +0.370422087950078230138;
    zj(311) = +0.376901004740559344802;
    zj(312) = +0.383359324198730346916;
    zj(313) = +0.389796704618470795479;
    zj(314) = +0.396212806057615939183;
    zj(315) = +0.402607290368737092671;
    zj(316) = +0.408979821229888672409;
    zj(317) = +0.415330064175321663764;
    zj(318) = +0.421657686626163300056;
    zj(319) = +0.427962357921062742583;
    zj(320) = +0.434243749346802558002;
    zj(321) = +0.440501534168875795783;
    zj(322) = +0.446735387662028473742;
    zj(323) = +0.452944987140767283784;
    zj(324) = +0.459130011989832332874;
    zj(325) = +0.465290143694634735858;
    zj(326) = +0.471425065871658876934;
    zj(327) = +0.477534464298829155284;
    zj(328) = +0.483618026945841027562;
    zj(329) = +0.489675444004456155436;
    zj(330) = +0.495706407918761460170;
    zj(331) = +0.501710613415391878251;
    zj(332) = +0.507687757533716602155;
    zj(333) = +0.513637539655988578507;
    zj(334) = +0.519559661537457021993;
    zj(335) = +0.525453827336442687395;
    zj(336) = +0.531319743644375623972;
    zj(337) = +0.537157119515795115982;
    zj(338) = +0.542965666498311490492;
    zj(339) = +0.548745098662529448608;
    zj(340) = +0.554495132631932548866;
    zj(341) = +0.560215487612728441818;
    zj(342) = +0.565905885423654422623;
    zj(343) = +0.571566050525742833992;
    zj(344) = +0.577195710052045814844;
    zj(345) = +0.582794593837318850840;
    zj(346) = +0.588362434447662541434;
    zj(347) = +0.593898967210121954393;
    zj(348) = +0.599403930242242892974;
    zj(349) = +0.604877064481584353319;
    zj(350) = +0.610318113715186400156;
    zj(351) = +0.615726824608992638014;
    zj(352) = +0.621102946737226402941;
    zj(353) = +0.626446232611719746542;
    zj(354) = +0.631756437711194230414;
    zj(355) = +0.637033320510492495071;
    zj(356) = +0.642276642509759513774;
    zj(357) = +0.647486168263572388782;
    zj(358) = +0.652661665410017496101;
    zj(359) = +0.657802904699713735422;
    zj(360) = +0.662909660024780595461;
    zj(361) = +0.667981708447749702165;
    zj(362) = +0.673018830230418479199;
    zj(363) = +0.678020808862644517838;
    zj(364) = +0.682987431091079228087;
    zj(365) = +0.687918486947839325756;
    zj(366) = +0.692813769779114702895;
    zj(367) = +0.697673076273711232906;
    zj(368) = +0.702496206491527078610;
    zj(369) = +0.707282963891961103412;
    zj(370) = +0.712033155362252034587;
    zj(371) = +0.716746591245747095767;
    zj(372) = +0.721423085370098915485;
    zj(373) = +0.726062455075389632685;
    zj(374) = +0.730664521242181261329;
    zj(375) = +0.735229108319491547663;
    zj(376) = +0.739756044352694758677;
    zj(377) = +0.744245161011347082309;
    zj(378) = +0.748696293616936602823;
    zj(379) = +0.753109281170558142523;
    zj(380) = +0.757483966380513637926;
    zj(381) = +0.761820195689839149173;
    zj(382) = +0.766117819303760090717;
    zj(383) = +0.770376691217076824278;
    zj(384) = +0.774596669241483377036;
    zj(385) = +0.778777615032822744702;
    zj(386) = +0.782919394118283016385;
    zj(387) = +0.787021875923539422170;
    zj(388) = +0.791084933799848361435;
    zj(389) = +0.795108445051100526780;
    zj(390) = +0.799092290960841401800;
    zj(391) = +0.803036356819268687782;
    zj(392) = +0.806940531950217611856;
    zj(393) = +0.810804709738146594361;
    zj(394) = +0.814628787655137413436;
    zj(395) = +0.818412667287925807395;
    zj(396) = +0.822156254364980407373;
    zj(397) = +0.825859458783650001088;
    zj(398) = +0.829522194637401400178;
    zj(399) = +0.833144380243172624728;
    zj(400) = +0.836725938168868735503;
    zj(401) = +0.840266795261030442350;
    zj(402) = +0.843766882672708601038;
    zj(403) = +0.847226135891580884381;
    zj(404) = +0.850644494768350279758;
    zj(405) = +0.854021903545468625813;
    zj(406) = +0.857358310886232156525;
    zj(407) = +0.860653669904299969802;
    zj(408) = +0.863907938193690477146;
    zj(409) = +0.867121077859315215614;
    zj(410) = +0.870293055548113905851;
    zj(411) = +0.873423842480859310192;
    zj(412) = +0.876513414484705269742;
    zj(413) = +0.879561752026556262568;
    zj(414) = +0.882568840247341906842;
    zj(415) = +0.885534668997285008926;
    zj(416) = +0.888459232872256998890;
    zj(417) = +0.891342531251319871666;
    zj(418) = +0.894184568335559022859;
    zj(419) = +0.896985353188316590376;
    zj(420) = +0.899744899776940036639;
    zj(421) = +0.902463227016165675048;
    zj(422) = +0.905140358813261595189;
    zj(423) = +0.907776324115058903624;
    zj(424) = +0.910371156957004292498;
    zj(425) = +0.912924896514370590080;
    zj(426) = +0.915437587155765040644;
    zj(427) = +0.917909278499077501636;
    zj(428) = +0.920340025470012420730;
    zj(429) = +0.922729888363349241523;
    zj(430) = +0.925078932907075652364;
    zj(431) = +0.927387230329536696843;
    zj(432) = +0.929654857429740056670;
    zj(433) = +0.931881896650953639345;
    zj(434) = +0.934068436157725787999;
    zj(435) = +0.936214569916450806625;
    zj(436) = +0.938320397779592883655;
    zj(437) = +0.940386025573669721370;
    zj(438) = +0.942411565191083059813;
    zj(439) = +0.944397134685866648591;
    zj(440) = +0.946342858373402905148;
    zj(441) = +0.948248866934137357063;
    zj(442) = +0.950115297521294876558;
    zj(443) = +0.951942293872573589498;
    zj(444) = +0.953730006425761136415;
    zj(445) = +0.955478592438183697574;
    zj(446) = +0.957188216109860962736;
    zj(447) = +0.958859048710200221356;
    zj(448) = +0.960491268708020283423;
    zj(449) = +0.962085061904651475741;
    zj(450) = +0.963640621569812132521;
    zj(451) = +0.965158148579915665979;
    zj(452) = +0.966637851558416567092;
    zj(453) = +0.968079947017759947964;
    zj(454) = +0.969484659502459231771;
    zj(455) = +0.970852221732792443256;
    zj(456) = +0.972182874748581796578;
    zj(457) = +0.973476868052506926773;
    zj(458) = +0.974734459752402667761;
    zj(459) = +0.975955916702011753129;
    zj(460) = +0.977141514639705714156;
    zj(461) = +0.978291538324758539526;
    zj(462) = +0.979406281670862683806;
    zj(463) = +0.980486047876721339416;
    zj(464) = +0.981531149553740106867;
    zj(465) = +0.982541908851080604251;
    zj(466) = +0.983518657578632728762;
    zj(467) = +0.984461737328814534596;
    zj(468) = +0.985371499598520371114;
    zj(469) = +0.986248305913007552681;
    zj(470) = +0.987092527954034067190;
    zj(471) = +0.987904547695124280467;
    zj(472) = +0.988684757547429479939;
    zj(473) = +0.989433560520240838716;
    zj(474) = +0.990151370400770159181;
    zj(475) = +0.990838611958294243677;
    zj(476) = +0.991495721178106132399;
    zj(477) = +0.992123145530863117683;
    zj(478) = +0.992721344282788615328;
    zj(479) = +0.993290788851684966211;
    zj(480) = +0.993831963212755022209;
    zj(481) = +0.994345364356723405931;
    zj(482) = +0.994831502800621000519;
    zj(483) = +0.995290903148810302261;
    zj(484) = +0.995724104698407188509;
    zj(485) = +0.996131662079315037786;
    zj(486) = +0.996514145914890273849;
    zj(487) = +0.996872143485260161299;
    zj(488) = +0.997206259372221959076;
    zj(489) = +0.997517116063472399965;
    zj(490) = +0.997805354495957274562;
    zj(491) = +0.998071634524930323302;
    zj(492) = +0.998316635318407392531;
    zj(493) = +0.998541055697167906027;
    zj(494) = +0.998745614468095114704;
    zj(495) = +0.998931050830810562236;
    zj(496) = +0.999098124967667597662;
    zj(497) = +0.999247618943342473599;
    zj(498) = +0.999380338025023581928;
    zj(499) = +0.999497112467187190535;
    zj(500) = +0.999598799671910683252;
    zj(501) = +0.999686286448317731776;
    zj(502) = +0.999760490924432047330;
    zj(503) = +0.999822363679787739196;
    zj(504) = +0.999872888120357611938;
    zj(505) = +0.999913081144678282800;
    zj(506) = +0.999943996207054375764;
    zj(507) = +0.999966730098486276883;
    zj(508) = +0.999982430354891598580;
    zj(509) = +0.999992298136257588028;
    zj(510) = +0.999997596379748464620;
    zj(511) = +0.999999672956734384381;
    wj(  1) = 0.945715933950007048827E-06;
    wj(  2) = 0.345456507169149134898E-05;
    wj(  3) = 0.736624069102321668857E-05;
    wj(  4) = 0.125792781889592743525E-04;
    wj(  5) = 0.190213681905875816679E-04;
    wj(  6) = 0.266376412339000901358E-04;
    wj(  7) = 0.353751372055189588628E-04;
    wj(  8) = 0.451863674126296143105E-04;
    wj(  9) = 0.560319507856164252140E-04;
    wj( 10) = 0.678774554733972416227E-04;
    wj( 11) = 0.806899228014035293851E-04;
    wj( 12) = 0.944366322532705527066E-04;
    wj( 13) = 0.109085545645741522051E-03;
    wj( 14) = 0.124606200241498368482E-03;
    wj( 15) = 0.140970302204104791413E-03;
    wj( 16) = 0.158151830411132242924E-03;
    wj( 17) = 0.176126765545083195474E-03;
    wj( 18) = 0.194872642236641146532E-03;
    wj( 19) = 0.214368090034216937149E-03;
    wj( 20) = 0.234592462123925204879E-03;
    wj( 21) = 0.255525589595236862014E-03;
    wj( 22) = 0.277147657465187357459E-03;
    wj( 23) = 0.299439176850911730874E-03;
    wj( 24) = 0.322381020652862389664E-03;
    wj( 25) = 0.345954492129903871350E-03;
    wj( 26) = 0.370141402122251665232E-03;
    wj( 27) = 0.394924138246873704434E-03;
    wj( 28) = 0.420285716355361231823E-03;
    wj( 29) = 0.446209810101403247488E-03;
    wj( 30) = 0.472680758429262691232E-03;
    wj( 31) = 0.499683553312800484519E-03;
    wj( 32) = 0.527203811431658386125E-03;
    wj( 33) = 0.555227733977307579715E-03;
    wj( 34) = 0.583742058714979703847E-03;
    wj( 35) = 0.612734008012225209294E-03;
    wj( 36) = 0.642191235948505088403E-03;
    wj( 37) = 0.672101776960108194646E-03;
    wj( 38) = 0.702453997827572321358E-03;
    wj( 39) = 0.733236554224767912055E-03;
    wj( 40) = 0.764438352543882784191E-03;
    wj( 41) = 0.796048517297550871506E-03;
    wj( 42) = 0.828056364077226302608E-03;
    wj( 43) = 0.860451377808527848128E-03;
    wj( 44) = 0.893223195879324912340E-03;
    wj( 45) = 0.926361595613111283368E-03;
    wj( 46) = 0.959856485506936206261E-03;
    wj( 47) = 0.993697899638760857945E-03;
    wj( 48) = 0.102787599466367326179E-02;
    wj( 49) = 0.106238104885340071375E-02;
    wj( 50) = 0.109720346268191941940E-02;
    wj( 51) = 0.113233376051597664917E-02;
    wj( 52) = 0.116776259302858043685E-02;
    wj( 53) = 0.120348074001265964881E-02;
    wj( 54) = 0.123947911332878396534E-02;
    wj( 55) = 0.127574875977346947345E-02;
    wj( 56) = 0.131228086370221478128E-02;
    wj( 57) = 0.134906674928353113127E-02;
    wj( 58) = 0.138609788229672549700E-02;
    wj( 59) = 0.142336587141720519900E-02;
    wj( 60) = 0.146086246895890987689E-02;
    wj( 61) = 0.149857957106456636214E-02;
    wj( 62) = 0.153650921735128916170E-02;
    wj( 63) = 0.157464359003212166189E-02;
    wj( 64) = 0.161297501254393423070E-02;
    wj( 65) = 0.165149594771914570655E-02;
    wj( 66) = 0.169019899554346019117E-02;
    wj( 67) = 0.172907689054461607168E-02;
    wj( 68) = 0.176812249885838886701E-02;
    wj( 69) = 0.180732881501808930079E-02;
    wj( 70) = 0.184668895851282540913E-02;
    wj( 71) = 0.188619617015808475394E-02;
    wj( 72) = 0.192584380831993546204E-02;
    wj( 73) = 0.196562534503150547732E-02;
    wj( 74) = 0.200553436203751169944E-02;
    wj( 75) = 0.204556454679958293446E-02;
    wj( 76) = 0.208570968849203942640E-02;
    wj( 77) = 0.212596367401472533045E-02;
    wj( 78) = 0.216632048404649142727E-02;
    wj( 79) = 0.220677418916003329194E-02;
    wj( 80) = 0.224731894601603393082E-02;
    wj( 81) = 0.228794899365195972378E-02;
    wj( 82) = 0.232865864987842738864E-02;
    wj( 83) = 0.236944230779380495146E-02;
    wj( 84) = 0.241029443242563417382E-02;
    wj( 85) = 0.245120955750556483923E-02;
    wj( 86) = 0.249218228238276930060E-02;
    wj( 87) = 0.253320726907925325750E-02;
    wj( 88) = 0.257427923948908888092E-02;
    wj( 89) = 0.261539297272236109225E-02;
    wj( 90) = 0.265654330259352828314E-02;
    wj( 91) = 0.269772511525294586667E-02;
    wj( 92) = 0.273893334695947541201E-02;
    wj( 93) = 0.278016298199139435045E-02;
    wj( 94) = 0.282140905069222207923E-02;
    wj( 95) = 0.286266662764757868253E-02;
    wj( 96) = 0.290393082998878368175E-02;
    wj( 97) = 0.294519681581857582284E-02;
    wj( 98) = 0.298645978275408290247E-02;
    wj( 99) = 0.302771496658198544480E-02;
    wj(100) = 0.306895764002069252174E-02;
    wj(101) = 0.311018311158427546158E-02;
    wj(102) = 0.315138672454287935858E-02;
    wj(103) = 0.319256385597434736790E-02;
    wj(104) = 0.323370991590184336368E-02;
    wj(105) = 0.327482034651233969564E-02;
    wj(106) = 0.331589062145094394706E-02;
    wj(107) = 0.335691624518616761342E-02;
    wj(108) = 0.339789275244138669739E-02;
    wj(109) = 0.343881570768790591876E-02;
    wj(110) = 0.347968070469521146972E-02;
    wj(111) = 0.352048336613417922682E-02;
    wj(112) = 0.356121934322919357659E-02;
    wj(113) = 0.360188431545532431869E-02;
    wj(114) = 0.364247399027690353194E-02;
    wj(115) = 0.368298410292403911967E-02;
    wj(116) = 0.372341041620379550870E-02;
    wj(117) = 0.376374872034296338241E-02;
    wj(118) = 0.380399483285952829161E-02;
    wj(119) = 0.384414459846013158917E-02;
    wj(120) = 0.388419388896099560998E-02;
    wj(121) = 0.392413860322995774660E-02;
    wj(122) = 0.396397466714742455513E-02;
    wj(123) = 0.400369803358421688562E-02;
    wj(124) = 0.404330468239442998549E-02;
    wj(125) = 0.408279062042157838350E-02;
    wj(126) = 0.412215188151643401528E-02;
    wj(127) = 0.416138452656509745764E-02;
    wj(128) = 0.420048464352596631772E-02;
    wj(129) = 0.423944834747438184434E-02;
    wj(130) = 0.427827178065384480959E-02;
    wj(131) = 0.431695111253279479928E-02;
    wj(132) = 0.435548253986604343679E-02;
    wj(133) = 0.439386228676004195260E-02;
    wj(134) = 0.443208660474124713206E-02;
    wj(135) = 0.447015177282692726900E-02;
    wj(136) = 0.450805409759782158001E-02;
    wj(137) = 0.454578991327213285488E-02;
    wj(138) = 0.458335558178039420335E-02;
    wj(139) = 0.462074749284080687482E-02;
    wj(140) = 0.465796206403469754658E-02;
    wj(141) = 0.469499574088179046532E-02;
    wj(142) = 0.473184499691503264714E-02;
    wj(143) = 0.476850633375474925263E-02;
    wj(144) = 0.480497628118194150483E-02;
    wj(145) = 0.484125139721057135214E-02;
    wj(146) = 0.487732826815870573054E-02;
    wj(147) = 0.491320350871841897367E-02;
    wj(148) = 0.494887376202437487201E-02;
    wj(149) = 0.498433569972103029914E-02;
    wj(150) = 0.501958602202842039909E-02;
    wj(151) = 0.505462145780650125058E-02;
    wj(152) = 0.508943876461803986674E-02;
    wj(153) = 0.512403472879005351831E-02;
    wj(154) = 0.515840616547381084096E-02;
    wj(155) = 0.519254991870341614863E-02;
    wj(156) = 0.522646286145300596306E-02;
    wj(157) = 0.526014189569259311205E-02;
    wj(158) = 0.529358395244259896547E-02;
    wj(159) = 0.532678599182711857974E-02;
    wj(160) = 0.535974500312596681161E-02;
    wj(161) = 0.539245800482555593606E-02;
    wj(162) = 0.542492204466865704951E-02;
    wj(163) = 0.545713419970309863995E-02;
    wj(164) = 0.548909157632945623482E-02;
    wj(165) = 0.552079131034778706457E-02;
    wj(166) = 0.555223056700346326850E-02;
    wj(167) = 0.558340654103215637610E-02;
    wj(168) = 0.561431645670402467678E-02;
    wj(169) = 0.564495756786715368885E-02;
    wj(170) = 0.567532715799029830087E-02;
    wj(171) = 0.570542254020497332312E-02;
    wj(172) = 0.573524105734693719020E-02;
    wj(173) = 0.576478008199711142954E-02;
    wj(174) = 0.579403701652197628421E-02;
    wj(175) = 0.582300929311348057702E-02;
    wj(176) = 0.585169437382850155033E-02;
    wj(177) = 0.588008975062788803205E-02;
    wj(178) = 0.590819294541511788161E-02;
    wj(179) = 0.593600151007459827614E-02;
    wj(180) = 0.596351302650963502011E-02;
    wj(181) = 0.599072510668009471472E-02;
    wj(182) = 0.601763539263978131522E-02;
    wj(183) = 0.604424155657354634589E-02;
    wj(184) = 0.607054130083414983949E-02;
    wj(185) = 0.609653235797888692923E-02;
    wj(186) = 0.612221249080599294931E-02;
    wj(187) = 0.614757949239083790214E-02;
    wj(188) = 0.617263118612191922727E-02;
    wj(189) = 0.619736542573665996342E-02;
    wj(190) = 0.622178009535701763157E-02;
    wj(191) = 0.624587310952490748541E-02;
    wj(192) = 0.626964241323744217671E-02;
    wj(193) = 0.629308598198198836688E-02;
    wj(194) = 0.631620182177103938227E-02;
    wj(195) = 0.633898796917690165912E-02;
    wj(196) = 0.636144249136619145314E-02;
    wj(197) = 0.638356348613413709795E-02;
    wj(198) = 0.640534908193868098342E-02;
    wj(199) = 0.642679743793437438922E-02;
    wj(200) = 0.644790674400605734710E-02;
    wj(201) = 0.646867522080231481688E-02;
    wj(202) = 0.648910111976869964292E-02;
    wj(203) = 0.650918272318071200827E-02;
    wj(204) = 0.652891834417652442012E-02;
    wj(205) = 0.654830632678944064054E-02;
    wj(206) = 0.656734504598007641819E-02;
    wj(207) = 0.658603290766824937794E-02;
    wj(208) = 0.660436834876456498276E-02;
    wj(209) = 0.662234983720168509457E-02;
    wj(210) = 0.663997587196526532519E-02;
    wj(211) = 0.665724498312454708217E-02;
    wj(212) = 0.667415573186258997654E-02;
    wj(213) = 0.669070671050613006584E-02;
    wj(214) = 0.670689654255504925648E-02;
    wj(215) = 0.672272388271144108036E-02;
    wj(216) = 0.673818741690825799086E-02;
    wj(217) = 0.675328586233752529078E-02;
    wj(218) = 0.676801796747810680683E-02;
    wj(219) = 0.678238251212300746082E-02;
    wj(220) = 0.679637830740619795480E-02;
    wj(221) = 0.681000419582894688374E-02;
    wj(222) = 0.682325905128564571420E-02;
    wj(223) = 0.683614177908911221841E-02;
    wj(224) = 0.684865131599535812903E-02;
    wj(225) = 0.686078663022780697951E-02;
    wj(226) = 0.687254672150094831613E-02;
    wj(227) = 0.688393062104341470995E-02;
    wj(228) = 0.689493739162046825872E-02;
    wj(229) = 0.690556612755588354803E-02;
    wj(230) = 0.691581595475321433825E-02;
    wj(231) = 0.692568603071643155621E-02;
    wj(232) = 0.693517554456992049848E-02;
    wj(233) = 0.694428371707782549438E-02;
    wj(234) = 0.695300980066273063177E-02;
    wj(235) = 0.696135307942366551493E-02;
    wj(236) = 0.696931286915342540213E-02;
    wj(237) = 0.697688851735519545845E-02;
    wj(238) = 0.698407940325846925786E-02;
    wj(239) = 0.699088493783425207545E-02;
    wj(240) = 0.699730456380953992594E-02;
    wj(241) = 0.700333775568106572820E-02;
    wj(242) = 0.700898401972830440494E-02;
    wj(243) = 0.701424289402572916425E-02;
    wj(244) = 0.701911394845431165171E-02;
    wj(245) = 0.702359678471225911031E-02;
    wj(246) = 0.702769103632498213858E-02;
    wj(247) = 0.703139636865428709508E-02;
    wj(248) = 0.703471247890678765907E-02;
    wj(249) = 0.703763909614153052319E-02;
    wj(250) = 0.704017598127683066242E-02;
    wj(251) = 0.704232292709631209597E-02;
    wj(252) = 0.704407975825415053266E-02;
    wj(253) = 0.704544633127951476780E-02;
    wj(254) = 0.704642253458020417748E-02;
    wj(255) = 0.704700828844548013730E-02;
    wj(256) = 0.704720354504808967346E-02;
    wj(257) = 0.704700828844548013730E-02;
    wj(258) = 0.704642253458020417748E-02;
    wj(259) = 0.704544633127951476780E-02;
    wj(260) = 0.704407975825415053266E-02;
    wj(261) = 0.704232292709631209597E-02;
    wj(262) = 0.704017598127683066242E-02;
    wj(263) = 0.703763909614153052319E-02;
    wj(264) = 0.703471247890678765907E-02;
    wj(265) = 0.703139636865428709508E-02;
    wj(266) = 0.702769103632498213858E-02;
    wj(267) = 0.702359678471225911031E-02;
    wj(268) = 0.701911394845431165171E-02;
    wj(269) = 0.701424289402572916425E-02;
    wj(270) = 0.700898401972830440494E-02;
    wj(271) = 0.700333775568106572820E-02;
    wj(272) = 0.699730456380953992594E-02;
    wj(273) = 0.699088493783425207545E-02;
    wj(274) = 0.698407940325846925786E-02;
    wj(275) = 0.697688851735519545845E-02;
    wj(276) = 0.696931286915342540213E-02;
    wj(277) = 0.696135307942366551493E-02;
    wj(278) = 0.695300980066273063177E-02;
    wj(279) = 0.694428371707782549438E-02;
    wj(280) = 0.693517554456992049848E-02;
    wj(281) = 0.692568603071643155621E-02;
    wj(282) = 0.691581595475321433825E-02;
    wj(283) = 0.690556612755588354803E-02;
    wj(284) = 0.689493739162046825872E-02;
    wj(285) = 0.688393062104341470995E-02;
    wj(286) = 0.687254672150094831613E-02;
    wj(287) = 0.686078663022780697951E-02;
    wj(288) = 0.684865131599535812903E-02;
    wj(289) = 0.683614177908911221841E-02;
    wj(290) = 0.682325905128564571420E-02;
    wj(291) = 0.681000419582894688374E-02;
    wj(292) = 0.679637830740619795480E-02;
    wj(293) = 0.678238251212300746082E-02;
    wj(294) = 0.676801796747810680683E-02;
    wj(295) = 0.675328586233752529078E-02;
    wj(296) = 0.673818741690825799086E-02;
    wj(297) = 0.672272388271144108036E-02;
    wj(298) = 0.670689654255504925648E-02;
    wj(299) = 0.669070671050613006584E-02;
    wj(300) = 0.667415573186258997654E-02;
    wj(301) = 0.665724498312454708217E-02;
    wj(302) = 0.663997587196526532519E-02;
    wj(303) = 0.662234983720168509457E-02;
    wj(304) = 0.660436834876456498276E-02;
    wj(305) = 0.658603290766824937794E-02;
    wj(306) = 0.656734504598007641819E-02;
    wj(307) = 0.654830632678944064054E-02;
    wj(308) = 0.652891834417652442012E-02;
    wj(309) = 0.650918272318071200827E-02;
    wj(310) = 0.648910111976869964292E-02;
    wj(311) = 0.646867522080231481688E-02;
    wj(312) = 0.644790674400605734710E-02;
    wj(313) = 0.642679743793437438922E-02;
    wj(314) = 0.640534908193868098342E-02;
    wj(315) = 0.638356348613413709795E-02;
    wj(316) = 0.636144249136619145314E-02;
    wj(317) = 0.633898796917690165912E-02;
    wj(318) = 0.631620182177103938227E-02;
    wj(319) = 0.629308598198198836688E-02;
    wj(320) = 0.626964241323744217671E-02;
    wj(321) = 0.624587310952490748541E-02;
    wj(322) = 0.622178009535701763157E-02;
    wj(323) = 0.619736542573665996342E-02;
    wj(324) = 0.617263118612191922727E-02;
    wj(325) = 0.614757949239083790214E-02;
    wj(326) = 0.612221249080599294931E-02;
    wj(327) = 0.609653235797888692923E-02;
    wj(328) = 0.607054130083414983949E-02;
    wj(329) = 0.604424155657354634589E-02;
    wj(330) = 0.601763539263978131522E-02;
    wj(331) = 0.599072510668009471472E-02;
    wj(332) = 0.596351302650963502011E-02;
    wj(333) = 0.593600151007459827614E-02;
    wj(334) = 0.590819294541511788161E-02;
    wj(335) = 0.588008975062788803205E-02;
    wj(336) = 0.585169437382850155033E-02;
    wj(337) = 0.582300929311348057702E-02;
    wj(338) = 0.579403701652197628421E-02;
    wj(339) = 0.576478008199711142954E-02;
    wj(340) = 0.573524105734693719020E-02;
    wj(341) = 0.570542254020497332312E-02;
    wj(342) = 0.567532715799029830087E-02;
    wj(343) = 0.564495756786715368885E-02;
    wj(344) = 0.561431645670402467678E-02;
    wj(345) = 0.558340654103215637610E-02;
    wj(346) = 0.555223056700346326850E-02;
    wj(347) = 0.552079131034778706457E-02;
    wj(348) = 0.548909157632945623482E-02;
    wj(349) = 0.545713419970309863995E-02;
    wj(350) = 0.542492204466865704951E-02;
    wj(351) = 0.539245800482555593606E-02;
    wj(352) = 0.535974500312596681161E-02;
    wj(353) = 0.532678599182711857974E-02;
    wj(354) = 0.529358395244259896547E-02;
    wj(355) = 0.526014189569259311205E-02;
    wj(356) = 0.522646286145300596306E-02;
    wj(357) = 0.519254991870341614863E-02;
    wj(358) = 0.515840616547381084096E-02;
    wj(359) = 0.512403472879005351831E-02;
    wj(360) = 0.508943876461803986674E-02;
    wj(361) = 0.505462145780650125058E-02;
    wj(362) = 0.501958602202842039909E-02;
    wj(363) = 0.498433569972103029914E-02;
    wj(364) = 0.494887376202437487201E-02;
    wj(365) = 0.491320350871841897367E-02;
    wj(366) = 0.487732826815870573054E-02;
    wj(367) = 0.484125139721057135214E-02;
    wj(368) = 0.480497628118194150483E-02;
    wj(369) = 0.476850633375474925263E-02;
    wj(370) = 0.473184499691503264714E-02;
    wj(371) = 0.469499574088179046532E-02;
    wj(372) = 0.465796206403469754658E-02;
    wj(373) = 0.462074749284080687482E-02;
    wj(374) = 0.458335558178039420335E-02;
    wj(375) = 0.454578991327213285488E-02;
    wj(376) = 0.450805409759782158001E-02;
    wj(377) = 0.447015177282692726900E-02;
    wj(378) = 0.443208660474124713206E-02;
    wj(379) = 0.439386228676004195260E-02;
    wj(380) = 0.435548253986604343679E-02;
    wj(381) = 0.431695111253279479928E-02;
    wj(382) = 0.427827178065384480959E-02;
    wj(383) = 0.423944834747438184434E-02;
    wj(384) = 0.420048464352596631772E-02;
    wj(385) = 0.416138452656509745764E-02;
    wj(386) = 0.412215188151643401528E-02;
    wj(387) = 0.408279062042157838350E-02;
    wj(388) = 0.404330468239442998549E-02;
    wj(389) = 0.400369803358421688562E-02;
    wj(390) = 0.396397466714742455513E-02;
    wj(391) = 0.392413860322995774660E-02;
    wj(392) = 0.388419388896099560998E-02;
    wj(393) = 0.384414459846013158917E-02;
    wj(394) = 0.380399483285952829161E-02;
    wj(395) = 0.376374872034296338241E-02;
    wj(396) = 0.372341041620379550870E-02;
    wj(397) = 0.368298410292403911967E-02;
    wj(398) = 0.364247399027690353194E-02;
    wj(399) = 0.360188431545532431869E-02;
    wj(400) = 0.356121934322919357659E-02;
    wj(401) = 0.352048336613417922682E-02;
    wj(402) = 0.347968070469521146972E-02;
    wj(403) = 0.343881570768790591876E-02;
    wj(404) = 0.339789275244138669739E-02;
    wj(405) = 0.335691624518616761342E-02;
    wj(406) = 0.331589062145094394706E-02;
    wj(407) = 0.327482034651233969564E-02;
    wj(408) = 0.323370991590184336368E-02;
    wj(409) = 0.319256385597434736790E-02;
    wj(410) = 0.315138672454287935858E-02;
    wj(411) = 0.311018311158427546158E-02;
    wj(412) = 0.306895764002069252174E-02;
    wj(413) = 0.302771496658198544480E-02;
    wj(414) = 0.298645978275408290247E-02;
    wj(415) = 0.294519681581857582284E-02;
    wj(416) = 0.290393082998878368175E-02;
    wj(417) = 0.286266662764757868253E-02;
    wj(418) = 0.282140905069222207923E-02;
    wj(419) = 0.278016298199139435045E-02;
    wj(420) = 0.273893334695947541201E-02;
    wj(421) = 0.269772511525294586667E-02;
    wj(422) = 0.265654330259352828314E-02;
    wj(423) = 0.261539297272236109225E-02;
    wj(424) = 0.257427923948908888092E-02;
    wj(425) = 0.253320726907925325750E-02;
    wj(426) = 0.249218228238276930060E-02;
    wj(427) = 0.245120955750556483923E-02;
    wj(428) = 0.241029443242563417382E-02;
    wj(429) = 0.236944230779380495146E-02;
    wj(430) = 0.232865864987842738864E-02;
    wj(431) = 0.228794899365195972378E-02;
    wj(432) = 0.224731894601603393082E-02;
    wj(433) = 0.220677418916003329194E-02;
    wj(434) = 0.216632048404649142727E-02;
    wj(435) = 0.212596367401472533045E-02;
    wj(436) = 0.208570968849203942640E-02;
    wj(437) = 0.204556454679958293446E-02;
    wj(438) = 0.200553436203751169944E-02;
    wj(439) = 0.196562534503150547732E-02;
    wj(440) = 0.192584380831993546204E-02;
    wj(441) = 0.188619617015808475394E-02;
    wj(442) = 0.184668895851282540913E-02;
    wj(443) = 0.180732881501808930079E-02;
    wj(444) = 0.176812249885838886701E-02;
    wj(445) = 0.172907689054461607168E-02;
    wj(446) = 0.169019899554346019117E-02;
    wj(447) = 0.165149594771914570655E-02;
    wj(448) = 0.161297501254393423070E-02;
    wj(449) = 0.157464359003212166189E-02;
    wj(450) = 0.153650921735128916170E-02;
    wj(451) = 0.149857957106456636214E-02;
    wj(452) = 0.146086246895890987689E-02;
    wj(453) = 0.142336587141720519900E-02;
    wj(454) = 0.138609788229672549700E-02;
    wj(455) = 0.134906674928353113127E-02;
    wj(456) = 0.131228086370221478128E-02;
    wj(457) = 0.127574875977346947345E-02;
    wj(458) = 0.123947911332878396534E-02;
    wj(459) = 0.120348074001265964881E-02;
    wj(460) = 0.116776259302858043685E-02;
    wj(461) = 0.113233376051597664917E-02;
    wj(462) = 0.109720346268191941940E-02;
    wj(463) = 0.106238104885340071375E-02;
    wj(464) = 0.102787599466367326179E-02;
    wj(465) = 0.993697899638760857945E-03;
    wj(466) = 0.959856485506936206261E-03;
    wj(467) = 0.926361595613111283368E-03;
    wj(468) = 0.893223195879324912340E-03;
    wj(469) = 0.860451377808527848128E-03;
    wj(470) = 0.828056364077226302608E-03;
    wj(471) = 0.796048517297550871506E-03;
    wj(472) = 0.764438352543882784191E-03;
    wj(473) = 0.733236554224767912055E-03;
    wj(474) = 0.702453997827572321358E-03;
    wj(475) = 0.672101776960108194646E-03;
    wj(476) = 0.642191235948505088403E-03;
    wj(477) = 0.612734008012225209294E-03;
    wj(478) = 0.583742058714979703847E-03;
    wj(479) = 0.555227733977307579715E-03;
    wj(480) = 0.527203811431658386125E-03;
    wj(481) = 0.499683553312800484519E-03;
    wj(482) = 0.472680758429262691232E-03;
    wj(483) = 0.446209810101403247488E-03;
    wj(484) = 0.420285716355361231823E-03;
    wj(485) = 0.394924138246873704434E-03;
    wj(486) = 0.370141402122251665232E-03;
    wj(487) = 0.345954492129903871350E-03;
    wj(488) = 0.322381020652862389664E-03;
    wj(489) = 0.299439176850911730874E-03;
    wj(490) = 0.277147657465187357459E-03;
    wj(491) = 0.255525589595236862014E-03;
    wj(492) = 0.234592462123925204879E-03;
    wj(493) = 0.214368090034216937149E-03;
    wj(494) = 0.194872642236641146532E-03;
    wj(495) = 0.176126765545083195474E-03;
    wj(496) = 0.158151830411132242924E-03;
    wj(497) = 0.140970302204104791413E-03;
    wj(498) = 0.124606200241498368482E-03;
    wj(499) = 0.109085545645741522051E-03;
    wj(500) = 0.944366322532705527066E-04;
    wj(501) = 0.806899228014035293851E-04;
    wj(502) = 0.678774554733972416227E-04;
    wj(503) = 0.560319507856164252140E-04;
    wj(504) = 0.451863674126296143105E-04;
    wj(505) = 0.353751372055189588628E-04;
    wj(506) = 0.266376412339000901358E-04;
    wj(507) = 0.190213681905875816679E-04;
    wj(508) = 0.125792781889592743525E-04;
    wj(509) = 0.736624069102321668857E-05;
    wj(510) = 0.345456507169149134898E-05;
    wj(511) = 0.945715933950007048827E-06;
  otherwise
    disp ('*********************************************');
    disp ('Gauss_Patterson_nested_rule_1d - Fatal error!');
    disp ('Value of alpha must be    ');
    disp ('*********************************************');
    error('Gauss_Patterson_nested_rule_1d - Fatal error!');
end
% shift the Legendre Polynomials in order to fit the intervals    
zj = (0.5*(bjj-ajj)).*zj + (0.5*(ajj+bjj));
wj = (0.5*(bjj-ajj)).*wj;
return
end
