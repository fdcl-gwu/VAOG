clear;
close all;

load comp_0


figure(1);
for i=1:length(p_vec)
    plot(t_vec{i}, J_vec{i})
    hold on;
end
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'FontSize',16);

conv_x=[3.625 10; 2.577 10; 1.88 10; 1.864 10]
conv_y=[0.4068 0.02341; 0.07402 2.098e-5;0.06836 2.1e-8; 0.011 2.1e-11]

plot(conv_x',conv_y','k--');

diff(log10(conv_y'))./diff(log10(conv_x'))

legend('$\mathrm{LGVI}_2\;(p=2)$','$\mathrm{LGVI}_4\;(p=4)$','$\mathrm{LGVI}_6\;(p=6)$','$\mathrm{LGVI}_8\;(p=8)$',...
    'Location','SouthWest','interpreter','latex')

xlabel('$t$','interpreter','latex');
ylabel('${f}(R)-{f}(R^*)$','interpreter','latex');

text(5.5,5e-8,'$\mathcal{O}(t^{-11.95})$','interpreter','latex','fontsize',18);
text(5.5,1e-5,'$\mathcal{O}(t^{-8.97})$','interpreter','latex','fontsize',18);
text(5.5,2e-3,'$\mathcal{O}(t^{-6.02})$','interpreter','latex','fontsize',18);
text(5.5,3e-1,'$\mathcal{O}(t^{-2.81})$','interpreter','latex','fontsize',18);



print('comp_0','-depsc');
!cp comp_0.eps ../figs

