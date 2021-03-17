clear;
close all;

load comp_3


hfig = figure(1);
for i=1:length(h0_vec)
    plot(t_vec{i}, J_vec{i})
    hold on;
end
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'FontSize',16);

return;



legend('$h_0=0.001$', '$h_0=0.05$', '$h_0=0.01$', '$h_0=0.1$', '$h_0=0.4$',...
    'Location','SouthWest','interpreter','latex')
ylabel('${f}(R)-{f}(R^*)$','interpreter','latex');
xlabel('$t$','interpreter','latex');


return;

figure(1);
print('comp_3','-depsc');


!cp comp_3?.eps ../figs

