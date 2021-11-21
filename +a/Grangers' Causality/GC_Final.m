%% Returns Granger causality and corresponding p-values

function [F,pval] = GC_Final(X,params)

u.v2struct(params);

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,'F'); % take careful note of arguments!
sig  = significance(pval,0.05,mhtc);
