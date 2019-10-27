
function moment1 = momj(param,rdraws)
global p x g sj z w m
beta0=param(1); betax=param(2); beta2g=param(3); sigma=param(5); alpha=param(4);
beta2i_tilt = rdraws';
beta2i = beta2g + sigma*beta2i_tilt;
util_0 = beta0 + betax*x  + alpha *p;
util_1=repmat(util_0,1,50) + beta2i.*g;
sharej0 = zeros(200,50);

for mkt=1:10
    temp = (m==mkt);
    sharej0(temp,:) = exp(util_1(temp,:))./(1+sum(exp(util_1(temp,:)),1));
end
sharej = sum(sharej0,2)/50; 
moment = sj - sharej;
moment0 =[ones(200,1) x g z w];
moment1 = moment.*moment0;

end
