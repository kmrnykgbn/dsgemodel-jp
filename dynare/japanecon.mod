var y c co cr i n k kg b w rk r pi q g tau z y_obs c_obs i_obs g_obs pi_obs r_obs w_obs n_obs b_obs;
varexo  eg ez et er uy uc ui ug upi ur uw un ub;

parameters beta theta varphi alpha nu deltap deltag omega kappa eta psi phipi phiy phib rhog rhoz rhor;
parameters bybar gybar;

beta = 0.996;
theta = 1.5;
varphi = 2;
alpha = 0.33;
deltap = 0.06;
deltag= 0.04;
eta = 0.75;
kappa = 7;
psi = 11;
gybar = 0.08;
bybar = 4;

model (linear);
#rbar = 1/beta;
#rkbar = rbar+deltap-1;
#wbar = (((psi-1)/psi)*((((1-alpha)^(1-alpha))*(alpha^(alpha)))/(rkbar^(alpha))))^(1/(1-alpha));
#nybar = ((((1-alpha)*rkbar)/(alpha*wbar)))^(alpha);
#kybar = ((((1-alpha)*rkbar)/(alpha*wbar)))^(alpha-1);
#iybar = deltap*kybar;
#cybar = 1-iybar-gybar;

co = co(+1)-(1/theta)*(r-pi(+1));
cr = ((wbar*nybar)/cybar)*(w+n)-(1/cybar)*tau;
c = omega*cr+ (1-omega)*co;
n = (1/varphi)*w-(theta/varphi)*c;
i = (1/(1+beta))*i(-1)+(beta/(1+beta))*i(+1)+(kappa/(1+beta))*q;
q = pi(+1)-r+((1-deltap)/(1-deltap+rkbar))*q(+1)+(rkbar/(1-deltap+rkbar))*rk(+1);
pi = beta*pi(+1)+(((1-eta)*(1-beta*eta))/(eta))*((1-alpha)*w+alpha*rk-nu*kg(-1)-z);
k = (1-deltap)*k(-1)+deltap*i;
y = z+alpha*k(-1)+(1-alpha)*n+nu*kg(-1);
y = cybar*c+iybar*i+g;
r = rhor*r(-1)+phipi*pi+phiy*y+er;
n-k(-1) = rk-w;
b = rbar*b(-1)+bybar*rbar*r(-1)-rbar*bybar*pi+g-tau;
kg = (1-deltag)*kg(-1)+deltag*(1/gybar)*g;
tau = phib*b(-1)+et;
g = rhog*g(-1)+eg;
z = rhoz*z(-1)+ez;

y_obs = y+uy;
c_obs = c+uc;
i_obs = i+ui;
g_obs = g+ug;
pi_obs = pi+upi;
r_obs = r+ur;
w_obs = w+uw;
n_obs = n+un;
b_obs = b+ub;
end;

estimated params;
omega, beta_pdf, 0.25, 0.1;
nu, normal_pdf, 0.25, 0.1;
phipi, normal_pdf, 1.05, 0.2;
phiy, normal_pdf, 0.1, 0.1;
phib, normal_pdf, 0.03, 0.01;
rhog, beta_pdf, 0.63, 0.1;
rhoz, beta_pdf, 0.85, 0.1;
rhor, beta_pdf, 0.85, 0.1;
stderr et, inv_gamma_pdf, 0.1, inf;
stderr eg, inv_gamma_pdf, 0.1, inf;
stderr er, inv_gamma_pdf, 0.1, inf;
stderr ez, inv_gamma_pdf, 0.1, inf;
stderr uy, inv_gamma_pdf, 0.1, inf;
stderr uc, inv_gamma_pdf, 0.1, inf;
stderr ui, inv_gamma_pdf, 0.1, inf;
stderr ug, inv_gamma_pdf, 0.1, inf;
stderr upi, inv_gamma_pdf, 0.1, inf;
stderr ur, inv_gamma_pdf, 0.1, inf;
stderr uw, inv_gamma_pdf, 0.1, inf;
stderr un, inv_gamma_pdf, 0.1, inf;
stderr ub, inv_gamma_pdf, 0.1, inf;
end;

varobs y_obs  c_obs i_obs g_obs pi_obs r_obs w_obs n_obs b_obs;

estimation(datafile=jpecon,mode_check,mh_replic=500000,mh_nblocks=2,mh_drop=0.5,mh_jscale=0.5,bayesian_irf);

