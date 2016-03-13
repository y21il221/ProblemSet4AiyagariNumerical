% PROGRAM NAME: ps4huggett.m
close all
clear, clc
tic

t = cputime;

% PARAMETERS
beta = 0.99; %discount factor 
sigma = 2; % coefficient of risk aversion
alpha = 1/3;
rho = 0.5;
sigma_e = 0.2;
delta = 0.025;
m = 5;
[z, PI] = TAUCHEN(m, rho, sigma_e, 3)    ; % transition matrix
z=exp(z);

[V, D] = eig(PI');
PIstar = V(:, 1)./sum(V(:, 1));

% ASSET VECTOR
a_lo = 0; %lower bound of grid points
a_hi = 80;%upper bound of grid points
num_a = 80;

a_1 = linspace(a_lo, a_hi, num_a); % asset (row) vector


% INITIAL GUESS FOR q
K_min = 30;
K_max = 40;
K_guess = (K_min + K_max) / 2;

% ITERATE OVER ASSET PRICES
aggK = 1 ;
while abs(aggK) >= 0.01 ;

a_1 = linspace(a_lo, a_hi, num_a);
% CURRENT RETURN (UTILITY) FUNCTION
N = PIstar'*z;
r = 1+alpha*(K_guess^(alpha-1))*N^(1-alpha)-delta;
w = (1-alpha)*(K_guess^(alpha))*N^(-alpha);
cons = bsxfun(@minus, r*a_1', a_1);
cons = bsxfun(@plus, cons, permute(w*z', [1 3 2]));
ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
ret(cons<0) = -Inf;

% INITIAL VALUE FUNCTION GUESS
v_guess = zeros(m, num_a);

% VALUE FUNCTION ITERATION
v_tol = 1;
while v_tol >.0001;
   % CONSTRUCT TOTAL RETURN FUNCTION
   v_mat = ret + beta * ...
       repmat(permute(PI * v_guess, [3 2 1]), [num_a 1 1]);
   
   % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
   [vfn, pol_indx] = max(v_mat, [], 2);
   vfn_g = vfn;
   vfn = permute(vfn, [3 1 2]);
   
   v_tol = abs(max(max(v_guess(:) - vfn(:))));
   
   v_guess = vfn; %update value functions
end;


% KEEP DECSISION RULE
pol_indx = permute(pol_indx, [3 1 2]);
pol_fn = a_1(pol_indx);


    %% CUBIC INTERPOLATION ITERATION
pp=pchip(a_1, vfn);
num_a_pr = 80;
a_pr = linspace(a_lo, a_hi, num_a_pr);

cons = bsxfun(@plus, r*a_pr', permute(w*z', [1 3 2]));


v_guess = vfn_g;
vfn = zeros(num_a_pr, 1, 5);
v_mat = zeros(num_a_pr, 1, 5);


v_tol = 1;
while v_tol >.0001;
 for ii = 1:5
     ii;
        for jj = 1:num_a_pr
            
   % CONSTRUCT TOTAL RETURN FUNCTION
   
   f = @(a_opt)-((((cons(jj, 1, ii) - a_opt).^ (1-sigma)) ./ (1 - sigma)) + beta *PI(ii,:) * ppval(pp, a_opt));
       
       if pol_indx(ii, jj) == 1
           a_min = pol_fn(ii,jj);
           a_max = a_pr(pol_indx(ii,jj)+1);
       elseif pol_indx(ii, jj) == size(a_pr, 2)
               a_min = a_pr(pol_indx(ii,jj)-1);
               a_max = pol_fn(ii,jj);
       else
               a_min = a_pr(pol_indx(ii,jj)-1);
               a_max = a_pr(pol_indx(ii,jj)+1);
       end
       %}
       
       a_min=0;
       if a_max>cons(jj, 1, ii)
       a_max=cons(jj, 1, ii);
       end
                      
       GoldenSectionSearch;
   
       vfn(jj, 1, ii) = -f(x_min);
       pol_fn(ii,jj) = x_min;
       
        end
 end
    
   vfn_1 = permute(vfn, [3 1 2]);
   v_tol = abs(max(max(v_guess(:) - vfn(:))))
   pp=pchip(a_pr, vfn_1);
   
   v_guess = vfn; %update value functions
end;
 
 vfn = permute(vfn, [3 1 2]);
 pp=pchip(a_pr, vfn);
 

     %% ITERATION OF THE INTERPOLATION FUNCTION WITH N = 800
a_1=linspace(0, 80, 800);

cons = bsxfun(@minus, r*a_1', a_1);
cons = bsxfun(@plus, cons, permute(w*z', [1 3 2]));
ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
ret(cons<0) = -Inf;

v_mat = zeros(800, 1, 5);


v_guess = zeros(m, 800);
% VALUE FUNCTION ITERATION
v_tol = 1;
while v_tol >.0001;
            
   % CONSTRUCT TOTAL RETURN FUNCTION

   v_mat = ret + beta * ...
      repmat(permute(PI * ppval(pp, a_1), [3 2 1]), [800 1 1]);
   
       [vfn, pol_indx] = max(v_mat, [], 2);
   vfn = permute(vfn, [3 1 2]);
   
   v_tol = abs(max(v_guess(:) - vfn(:)));
   
   v_guess = vfn; %update value functions
end

pol_indx = permute(pol_indx, [3 1 2]);
pol_fn = a_1(pol_indx);


% SET UP INITITAL DISTRIBUTION
Mu = zeros(m,800);
Mu(1, 4) = 1; % initial guess: everyone employed, 0 assets
% Mu = ones(2, num_a); alternative initial guess: same mass in all states
% Mu = Mu_guess / sum(Mu_guess(:)); % normalize total mass to 1

% ITERATE OVER DISTRIBUTIONS
% way 1: loop over all non-zeros states
mu_tol = 1;
while mu_tol > 1e-08
    [emp_ind, a_ind] = find(Mu > 0); % find non-zero indices
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (PI(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
    end

    mu_tol = max(abs(MuNew(:) - Mu(:)));
    
    Mu = MuNew ;
end

% % way 2: use transition matrices
% T_tilde = zeros(num_a, num_a, 2, 2);
% % set up matrices:
% for from_a = 1:num_a
%     for from_s = 1:2
%         T_tilde(from_a, pol_indx(from_s, from_a), from_s, :) = PI(from_s,:);
%     end
% end
% % transition:
% mu_tol = 1;
% while mu_tol > 1e-08
%     MuNew = zeros(size(Mu));
%     for from_s = 1:2
%         for to_s = 1:2
%             MuNew(to_s,:) = MuNew(to_s,:) + Mu(from_s,:) * T_tilde(:,:,from_s,to_s);
%         end
%     end
%     mu_tol = max(abs(Mu(:) - MuNew(:)));
%     Mu = MuNew;
% end


% CHECK AGGREGATE DEMAND
aggK = sum(sum( pol_fn(:) .* Mu(:) )) - K_guess; % Aggregate future assets

if aggK > 0 ;
    K_min = K_guess ;
end ;
if aggK < 0 ;
    K_max = K_guess ;
end ;

display (['K = ', num2str(K_guess)])
display (['Aggregate desired wealth = ', num2str(aggK)]);
display (['New Kmin is ', num2str(K_min), ', new Kmax is ', num2str(K_max)]);
display (['New K is ', num2str((K_max + K_min)/2)]);

K_guess = (K_max + K_min)/2 ;

display (' ') ;

end

    %% PLOT THE POLICY FUNCTION
figure(1)
plot(a_1, pol_fn)
title('Policy function')
xlabel('a_t-1')
ylabel('a_t')
legend('z=0.5', 'z=0.7', 'z=1', 'z=1.4', 'z=2')
    %% PLOT THE WEALTH DISTRIBUTION
figure(2)
plot(a_1, Mu)
title('Wealth distribution')
xlabel('a')
ylabel('share')
legend('z=0.5', 'z=0.7', 'z=1', 'z=1.4', 'z=2')


    %% FIND TOTAL WEALTH DISTRIBUTION AND GINI
agg_wealth = sum(Mu,1) * a_1' ; % wealth is asset holdings plus incomes
wealth_dist = [sum(Mu,1); a_1]';
[A, ordr] = sort(wealth_dist(:,2), 1);
wealth_dist = wealth_dist(ordr,:);

% see formula on wikipedia for computation of gini in discrete
% distributions;
pct_share = cumsum(wealth_dist(:,1));
pct_dist = cumsum( (wealth_dist(:,2) ./ agg_wealth) .* wealth_dist(:,1) );
gini = 1 - sum( ([0; pct_dist(1:end-1)] + pct_dist) .* wealth_dist(:,1) );
display (['Gini coefficient of ', num2str(gini)]);
display (['r is ', num2str(r)]);

    %% PLOT THE LORENZ CURVE
figure(3)
plot(pct_share, pct_dist);
axis([0,1,0,1])
title('Lorenz curve')
xlabel('Cumulative shares of people from lowest to highest wealth')
ylabel('Cumulative shares of wealth')

e = cputime - t;
display (['runtime is ', num2str(e)])

%{
%% CALCULATE CONSUMPTION EQUIVALENTS
cons_FB = PIstar * y_s';
WFB = 1 / (1 - beta) * cons_FB ^ (1 - sigma) ./ (1 - sigma);
display (['Welfare in First Best is ', num2str(WFB)]);

lambda = (WFB ./ vfn) .^ (1 / (1 - sigma)) - 1 ;

gain = lambda(:)' * Mu(:);

display (['Change to FB would be a welfare gain (in pct units of consumption) of ', num2str(gain)]);
%}
