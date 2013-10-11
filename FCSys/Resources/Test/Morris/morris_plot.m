function [eff_mean, eff_std] = morris_plot(X, Y)

% we assume that p large is and that therefore
delta = 0.5;

k = size(X,2);
m = k+1;
r = size(X,1)/m;

% compute the elementary effects
xdiff = diff(X);
ydiff = diff(Y);

% eliminate the meaningless datapoints
xdiff(m:m:end,:) = [];
ydiff(m:m:end,:) = [];

for i=1:k
    ind = find(xdiff(:,i)~=0);
    eff = ydiff(ind)./(delta*sign(xdiff(ind,i)));
    eff_mean(i) = mean(eff);
    eff_std(i) = std(eff);
end

%plot elementary effects
plot(eff_mean,eff_std,'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8)
hold on
axisvals = axis;
offset = (axisvals(2)-axisvals(1))*.06;
line_vect = [axisvals(1),0,axisvals(2)];
plot(line_vect, abs(line_vect*sqrt(r)/2), '--b')
axis(axisvals);
for i=1:k
    text(eff_mean(i),eff_std(i)+offset,sprintf('%g',i),'FontSize',15);
end;
hold off
grid on
xlabel('Mean of Elementary Effects')
ylabel('Standard Deviation of Elementary Effects')
title('Method of Morris');

return;