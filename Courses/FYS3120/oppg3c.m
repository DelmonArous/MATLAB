clear all;

beta = 0.8;
gamma = 1/sqrt(1-beta^2);

t_tilde = linspace(-10,10,1001);
Deltat_tilde = gamma*sqrt(beta^2*gamma^2*t_tilde.^2 + 1) - beta^2*gamma^2*t_tilde; 

cos_theta = beta*((t_tilde./Deltat_tilde) - 1);
cos_theta_prime = (cos_theta + beta)./(1 + beta.*cos_theta);

frac_nu_nu0 = (gamma.*(1-beta.*cos_theta_prime));

figure;
plot(t_tilde,cos_theta)
xlabel(['$$\tilde{t}$$'],'interpreter','latex','FontSize',12)
ylabel(['$$\cos\theta$$'],'interpreter','latex','FontSize',12)

eqtext = '$$\cos\theta = \beta\left(\frac{\tilde{t}}{\Delta\tilde{t}} - 1 \right)$$';
text(-9, 0.7,eqtext, 'interpreter', 'latex', 'FontSize', 12)

figure;
plot(t_tilde,frac_nu_nu0)
xlabel(['$$\tilde{t}$$'],'interpreter','latex','FontSize',12)
ylabel(['$$\frac{\nu(t)}{\nu_0}$$'],'interpreter','latex','FontSize',12)

eqtext = '$$\frac{\nu(t)}{\nu_0} = \gamma(1-\beta\cos\theta'' )$$';
text(2, 2.5, eqtext, 'interpreter', 'latex', 'FontSize', 12)