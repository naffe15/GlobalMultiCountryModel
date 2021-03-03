function out=SimpleSeasonalAdjustment(y)





t=1:length(y);

t=t';

x=[ones(length(y),1) t];

b=regre_nan(x,y);

trend=x*b;
out.trend=trend;
residual=y-trend;
out.residual=residual;
q1=mean_nan(residual(1:4:end));
q2=mean_nan(residual(2:4:end));
q3=mean_nan(residual(3:4:end));
q4=mean_nan(residual(4:4:end));

q0=[q1 q2 q3 q4]';

qq=repmat(q0,ceil(length(y)/4),1);

out.sa=y-qq(1:length(y));
