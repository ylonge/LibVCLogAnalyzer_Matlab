%% common model.
% model based on laplacian distribution.

%% basic parameters.
dIntra = 2/3;
dInter = 5/6;
dNorm = 0.5;

%% distortion model.
% general integral model.
equationIntegral = @(a,b,x) (-0.5.*exp(-a.*x).*((a.*x+1).^2+1)+1+1./(1-exp(-x)).*(-0.5*exp(-(1+a).*x)).*(((1-b).*x+1).^2+1-exp(x).*((b.*x-1).^2+1)));
% simplified integral model.
eqT = @(d,t) 1+t.*exp(-d.*t)./2./(1-exp(-t)).*(t-2*d*t-2);
eqTLow = @(d,t) 1+t.*(t-2*d.*t-2)./2./exp(d.*t);
% library dist model.
eqDLib = @(c1,c2,x) c1.*equationIntegral(dIntra,1-dIntra,x./((c2/2).^0.5));
eqDLibSimpBaseT = @(d,m,q) 2./(m.^2).*eqT(d,m.*q);
eqDLibSimpQuadric = @(d,e,x) -d.*x.^2+e;
eqDLibSimpPower = @(a,k,b,x) a.*x.^k+b;
eqDLibSimpLinear = @(a,s,e,x) a.*s.*e.^0.5.*x;

eqDKey = @(c,s,q,p,n,x) n+p.*(s+eqDLib(c,x)).*eqT(dInter,q./(((s+eqDLib(c,x))/2).^0.5));
eqDKeySimp = @(a1,a2,b,d,c2,c3,f,s,q1) f.*s.*eqDLibSimp(a1,c2+d,q1)+b;

%% rate model.
% general entropy model.
equationEntropy = @(a,x) (exp(-a*x).*(log2((1-exp(-a*x))./(1-exp(-x)))+1+x.*(a*(1-exp(-x))+exp(-x))./(1-exp(-x))/log(2))-log2(1-exp(-a*x)));
% library rate model.
eqRLib = @(u,b,m,q) u.*equationEntropy(dIntra,m.*q)+b;
eqRLibSimpExp = @(c,d,e,g,x) (d.*exp(-e.*x./((c/2).^0.5))+g);
eqRLibSimpPower = @(u,k,b,x) u.*x.^(-k)+b;

eqRKey = @(u1,u2,f1,f2,b1,b2,c1,c2,c3,d,s,a1,q1,q2) s.*equationRLibSimp(u1,f1,b1,c2+d,q1)+b;
eqRKeySimpExp = @(a,b,c,d,e,x) a.*exp(-b.*(c.*x+d).^(-0.5))+e;

