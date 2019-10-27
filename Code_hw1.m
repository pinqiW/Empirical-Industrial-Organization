
clear

global rdraws p x g sj z w s0 m

temp = xlsread('ps1data','ps1data','A2:J201');
rdraws = load('rdraws');
obs=temp(:,1); m=temp(:,2); p=temp(:,3); x=temp(:,4); g=temp(:,5); sj=temp(:,6); sjg=temp(:,7); 
s0=temp(:,8); w=temp(:,9); z=temp(:,10);

Y = log(sj./s0);
X = [ones(length(Y),1) x g p];

%------------------------------------------%
% QUESTION 1
%------------------------------------------%
%[beta, V, R] = gls(Y,X,0)
beta = inv(X'*X)*X'*Y;
e = Y-X*beta;
var = ((e'*e)/(200-4))*inv(X'*X); 
se = sqrt(diag(var));

% Display Results
[beta se]

%------------------------------------------%
%	QUESTION 2
%------------------------------------------%

% Discussed in tex doc


%------------------------------------------%
%	QUESTION 3
%------------------------------------------%

% Using z as a intsrument, estimate by OLS
W = [ones(length(Y),1) x g w];
betal = inv(X'*W*inv(W'*W)*W'*X)*X'*W*inv(W'*W)*W'*Y;
e = Y-X*betal;
var = ((e'*e)/(200-4))*inv(X'*W*inv(W'*W)*W'*X);
se = sqrt(diag(var));


% Display Results
[betal se]

% Calculate the own and cross price elasticity matrix
% for products with obs 77, 78, 79, 80

	% get subset of matrix
	subsj = sj(77:80, 1);
	
	% Calculate price derivatives
	pricederv = zeros(4,4);
	alpha = betal(4);
	for j=1:4
		for k=1:4
			if (k==j) 	% own price deriv
				pricederv(j,k) = alpha*subsj(j)*( 1-subsj(j));
			else
				pricederv(j,k) = -alpha*subsj(j)*subsj(k);
			end
		end
	end
	

	% Calculate price elasticities
	subp = p(77:80, 1);
	c = [subp'; subp'; subp'; subp'];
	b = [subsj  subsj  subsj  subsj];
	priceratios = c./b;	% find pk/sj
	priceelast = pricederv.*priceratios
	
	% Calculate price markups
	crosspriceE = -diag(priceelast);
	markup = subp./crosspriceE;
	
	% Characteristics
	subz = z(77:80, 1);
	subx = x(77:80, 1);
	subg = g(77:80, 1);
	
	% Display results
	[markup subx subp subz subg]



%------------------------------------------%
%	QUESTION 4
%------------------------------------------%

% Create new charateristic matrix for nested logit IV
lnsjg = log(sjg);
X = [ones(length(Y),1) x g p lnsjg];

% Instrumental variable matrix
Z = [ones(length(Y),1) x g w z];

% 2SLS estimate
betanl = inv(X'*Z*inv(Z'*Z)*Z'*X)*X'*Z*inv(Z'*Z)*Z'*Y;
e = Y-X*betanl;
var = ((e'*e)/(200-5))*inv(X'*Z*inv(Z'*Z)*Z'*X);
se = sqrt(diag(var));

%Display results
[betanl se]




%------------------------------------------%
%	QUESTION 5
%------------------------------------------%

% Calculate the own and cross price elasticity matrix
% for products with obs 77, 78, 79, 80

	% get subset of matrix
	subsj = sj(77:80, 1);
	subsjg = sjg(77:80, 1);
	subg = g(77:80, 1);
	subsg = subsj./subsjg;
	
	pricederv = zeros(4,4);
	
	% Calculate own price derivaties
	rho = betanl(5);
	alpha = betanl(4);
	for j=1:4
		for k=1:4
			if (k==j) 	% own price deriv
				pricederv(j,k) = alpha*subsj(j)*(  (1/(1-rho))*(1-subsjg(j))+subsjg(j)*(1-subsg(j))  );
			
			else
				if (subg(j)==subg(k))	% in the same group
					pricederv(j,k) = -alpha*subsj(j)*((1/(1-rho))*(subsjg(k))+subsjg(k)*(1-subsg(k)));	
				
				else			% in different groups
					pricederv(j,k) = -alpha*subsj(j)*subsjg(k)*subsg(k);				
				
				end
				
			end
		end
	end
	priceelast = pricederv.*priceratios

	
	% Calculate price markups
	crosspriceE = -diag(priceelast);
	markup = subp./crosspriceE;
	

	% Display results
	[markup subx subp subz subg]

%------------------------------------------%
%	QUESTION 6
%------------------------------------------%

%------------------------------------------%
%	QUESTION 7
%------------------------------------------%

% initial value
beta0 = [betanl(1:4); 0]

% Optimizing parameters
options = optimset('TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8)

% Optimization routine
[betahat,fval,exitflag,output] = fminsearch(@mom, beta0, options) 



%------------------------------------------%
%	QUESTION 8
%------------------------------------------%
% Compute standard errors for simulated estimates, betahat


% Find Gamma using numerical derivates
derv = zeros(5,5);
startval = mean(momj(betahat,rdraws))
 

for k=1:5
	newbeta = betahat;
	newbeta(k) = 1.001*newbeta(k);
	newval = mean(momj(newbeta, rdraws));
	derv(:,k) = (newval-startval)./(.001*betahat(k));
end
derv;

% Calculate Gamma
Gamma = derv;

% Calculate V(theta*)
varGj = (1/length(x)^2)*momj(betahat,rdraws)'*momj(betahat,rdraws);

% The final standard errors assuming that A and inv(Vhat) are the same
SE = sqrt(diag(inv(Gamma'*inv(varGj)*Gamma)));
[betahat SE]


%------------------------------------------%
%	QUESTION 9
%------------------------------------------%






