
function tot = mom(param)
global rdraws p x g sj z w m
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
% You don't need A matrix here realy. A=I matrix works
tot0 = zeros(10,5);
Amat = zeros(5,5,10);

AAmat = zeros(5,5,200);
for i=1:200
    AAmat(:,:,i) = moment(i,1)^2*(moment0(i,:)'*moment0(i,:));
end

for mkt=1:10
    temp = (m==mkt);
    Amat(:,:,mkt) = ((1/sum(temp,1))*(sum(AAmat(:,:,temp),3)));
    tot0(mkt,:) = (1/(sum(temp,1)))*sum(moment1(temp,:),1);
   
end
A= eye(5,5);
tot0 = sum(moment1,1)/200;
tot = tot0*A*tot0';
tot=sum(tot);
end
