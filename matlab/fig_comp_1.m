clear;
close all;

load comp_1


hfig = figure(1);
for i=1:length(p_vec)
    plot(1:N, J_vec{i})
    hold on;
end
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'FontSize',16);

conv_x=[84 1e4; 39 1e4; 39 1e4; 20 1e4];
conv_y=[0.07652 1.278e-6; 0.07402 2.47e-7; 0.02522 8.11e-8; 0.07677 4.894e-8];

plot(conv_x',conv_y','k--');

diff(log10(conv_y'))./diff(log10(conv_x'))

legend('$\mathrm{LGVI}_2\;(p=2)$','$\mathrm{LGVI}_4\;(p=4)$','$\mathrm{LGVI}_6\;(p=6)$','$\mathrm{LGVI}_8\;(p=8)$',...
    'Location','SouthWest','interpreter','latex')

xlabel('$k$','interpreter','latex');
ylabel('${f}(R)-{f}(R^*)$','interpreter','latex');

annotation('textarrow',[0.75 0.6],[0.85 0.68], 'String','$\mathcal{O}(k^{-2.30})$','interpreter','latex','fontsize',18);
annotation('textarrow',[0.75 0.6],[0.78 0.635], 'String','$\mathcal{O}(k^{-2.27})$','interpreter','latex','fontsize',18);
annotation('textarrow',[0.75 0.6],[0.72 0.605], 'String','$\mathcal{O}(k^{-2.28})$','interpreter','latex','fontsize',18);
annotation('textarrow',[0.75 0.6],[0.65 0.595], 'String','$\mathcal{O}(k^{-2.29})$','interpreter','latex','fontsize',18);

print('comp_1','-depsc');
!cp comp_1.eps ../figs

