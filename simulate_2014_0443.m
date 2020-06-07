clear all
load data_sim.txt
load def_per.txt
load def_per_ss.txt
load param.txt %_trend_y.txt
load delta.txt

r = 0.04;
coupon = (r+delta)/(1+r);

per_num = param(1);  %HOW MANY PERIODS IN EACH SAMPLE
n = param(2);        %HOW MANY SAMPLES

period_vector = kron(ones(n,1),linspace(1,per_num,per_num)');

y = zeros(per_num, n);
bshort = zeros(per_num, n);
blong = zeros(per_num, n);
q = zeros(per_num, n);
c = zeros(per_num, n);
d = zeros(per_num, n);
excl = zeros(per_num, n);
b_next = zeros(per_num-1, n);
duration = zeros(per_num,n);
market_access = zeros(per_num, n);
%pr_ss = zeros(per_num,n);
z = zeros(per_num, n);

for i=1:n
   y(:,i) = data_sim((i-1)*per_num+1:i*per_num,1); 
   bshort(:,i) = data_sim((i-1)*per_num+1:i*per_num,2); 
   blong(:,i) = data_sim((i-1)*per_num+1:i*per_num,3)*coupon/(1-(1-delta)*exp(-r)); %exp(-r)*(r+delta)/coupon; 
   q(:,i) = data_sim((i-1)*per_num+1:i*per_num,4); 
   c(:,i) = data_sim((i-1)*per_num+1:i*per_num,5); 
   d(:,i) = data_sim((i-1)*per_num+1:i*per_num,7)-1;
   excl(:,i) = data_sim((i-1)*per_num+1:i*per_num,8); 
   %MARKET ACCESS INDICATOR = 1 ==> BAD, HIGH RISK PREMIUM
   %MARKET ACCESS INDICATOR = 0 ==> GOOD, LOW RISK PREMIUM
   market_access(:,i) = data_sim((i-1)*per_num+1:i*per_num,9)-1;
   duration(:,i) = (((1-delta).*q(:,i) +coupon)./(q(:,i)))./(((1-delta).*q(:,i) +coupon)./(q(:,i))-1+delta);
   z(:,i) = data_sim((i-1)*per_num+1:i*per_num,10); 
   indices = find(blong(2:per_num,i)>-0.001);
  
 %  pr_ss(:,i) = data_sim((i-1)*per_num+1:i*per_num,10); 
end
%6.2f


%CREATE OUTPUT SERIES
num = 35; %HOW MANY PERIODS IN EACH SUBSAMPLE
indices = find(sum(d(per_num - num - 25 +1:per_num,:))<1); %LAST DEFAULT: 20 PERIODS BEFORE THE BEGINNING OF EACH SAMPLE.
%b1 = zeros(per_num, n);
%for j=1:n
%    y(1,j) = 1;
 %   for i=2:per_num
 %      y(i,j) = y(i-1,j)*g(i,j);
 %      b1(i,j) = y(i-1,j)*b(i,j);
 %   end
%end

num_observations = length(indices);
last_y = y(per_num - num+1:per_num, indices);  %TRIM THE FIRST 9500 OBSERVATIONS
last_c = c(per_num - num+1:per_num, indices);  %TRIM THE FIRST 9500 OBSERVATIONS
last_q = q(per_num - num+1:per_num, indices);  %TRIM THE FIRST 9500 OBSERVATIONS
last_d = d(per_num - num+1:per_num, indices);
last_short = bshort(per_num - num+1:per_num, indices);
last_long = blong(per_num - num+1:per_num, indices);
last_q = q(per_num - num+1:per_num, indices);
%last_spreadannual = ((coupon./q(per_num - num+1:per_num,indices) +1-delta)/(1+r)).^4 - 1;
%last_spread = ((coupon./q(per_num - num+1:per_num,indices) +1-delta)/(1+r)) - 1;
last_spread = (log(coupon./q(per_num - num+1:per_num, indices) - delta + 1) - r);
last_market_access = market_access(per_num - num+1:per_num, indices);

last_duration = duration(per_num - num+1:per_num, indices);
%MATRIX THAT STORES THE EXTRA SPREAD THAT WOULD BE OBSERVED IF THE GOVERNMENT INCREASED THE DEBT BY 1% OF GDP 
%TO BUY RESERVES. NEED TO REDUCE THE SAMPLE BY 1 OBSERVATION BECAUSE NEED
%THE END OF PERIOD PORTFOLIO INFORMATION.
last_extra_spread = zeros(num-1, num_observations); 


X = [ones(num,1) linspace(1,num,num)'];
y_trend = zeros(num, num_observations);
c_trend = zeros(num, num_observations);

for i=1:num_observations
    y_trend(:,i) = X*(X'*X)^(-1)*X'*last_y(:,i); % hpfilter(last_y,lambda)';            
     c_trend(:,i) = X*(X'*X)^(-1)*X'*last_c(:,i); %hpfilter(last_c, lambda)';           %FILTER log(consumption)
% spread_trend = hpfilter(last_spread, lambda)'; %last_spread * 0; %NO TREND IN SPREADS hpfilter(last_spread, lambda);  %FILTER quaterly spread
%COMPUTE DEVIATIONS FROM TREND  
end

y_dev = last_y - y_trend;
c_dev = last_c - c_trend;
spread_dev = last_spread; % - spread_trend;

fprintf('std of output  = %6.1f \n',100*mean(std(y_dev)))
fprintf('std cons / std y   = %6.1f \n',mean(std(c_dev))/mean(std(y_dev)))
fprintf('std of R_s     = %6.1f \n',100*mean(std(spread_dev)))

matrix_corr = zeros(num_observations,2);
matrix_corr_dres = zeros(num_observations,2);
matrix_corr_ddebt = zeros(num_observations,2);
matrix_corr_dres1 = zeros(num_observations,2);
matrix_corr_ddebt1 = zeros(num_observations,2);
matrix_corr_res = zeros(num_observations,2);
matrix_corr_debt = zeros(num_observations,2);

spread_noss = zeros(num_observations,1);
spread_ss = zeros(num_observations,1);

for i=1:num_observations
    matrix = corrcoef([y_dev(:,i) c_dev(:,i) spread_dev(:,i)  ]);
    matrix_corr(i,:) = matrix(1,2:3);

%     matrix = corrcoef([log(last_short(:,i)) y_dev(2:num,i) spread_dev(2:num,i)]);
%     matrix_corr_dres(i,1) = matrix(1,2);
%     matrix_corr_dres(i,2) = matrix(1,3);
%     
%     matrix = corrcoef([log(last_long(:,i)) y_dev(2:num,i) spread_dev(2:num,i)]);
%     matrix_corr_ddebt(i,1) = matrix(1,2);
%     matrix_corr_ddebt(i,2) = matrix(1,3);
% 
%     matrix = corrcoef([log(last_short(2:num,i)) diff(y_dev(1:num-1,i)) ]);
%     matrix_corr_dres1(i,1) = matrix(1,2);
%     
%     matrix = corrcoef([log(last_long(2:num,i)) diff(y_dev(1:num-1,i)) ]);
%     matrix_corr_ddebt1(i,1) = matrix(1,2);
    
    matrix = corrcoef([last_short(2:num,i) y_dev(1:num-1,i) spread_dev(1:num-1,i)]);
    matrix_corr_res(i,1) = matrix(1,2);
    matrix_corr_res(i,2) = matrix(1,3);
    
    matrix = corrcoef([last_long(2:num,i) y_dev(1:num-1,i) spread_dev(1:num-1,i)]);
    matrix_corr_debt(i,1) = matrix(1,2);
    matrix_corr_debt(i,2) = matrix(1,3);

   %MARKET ACCESS INDICATOR = 1 ==> BAD, HIGH RISK PREMIUM
   %MARKET ACCESS INDICATOR = 0 ==> GOOD, LOW RISK PREMIUM
    
    ind_ss = find(last_market_access(:,i)>0);
    ind_noss = find(last_market_access(:,i)<1);
    spread_noss(i) = mean(last_spread(ind_noss,i));
    spread_ss(i) = mean(last_spread(ind_ss,i));
end

ind_ss = (spread_ss>0);
fprintf('corr(c, y)     = %6.1f \n',mean(matrix_corr(:,1)))
fprintf('corr(R_s,y)    = %6.1f \n',mean(matrix_corr(:,2)))
fprintf('Annual mean default rate (all periods) = %6.2f \n', 100*mean(sum(d)/per_num))
fprintf('Mean debt/y    = %6.1f \n',100*mean(mean(last_long./(exp(last_y)))))
fprintf('Mean (res)  = %6.1f \n',100*mean(mean(last_short./(exp(last_y)))))
fprintf('Mean(res)/mean(gdp)  = %6.1f \n',100*mean(mean(last_short))./mean(mean((exp(last_y)))))
fprintf('E(R_s)         = %6.1f \n',100*mean(mean(last_spread)))
fprintf('Max R_s        = %6.1f \n',100*max(max(last_spread)))

fprintf('corr(res., y)  = %6.1f \n',mean(matrix_corr_res(:,1)))
fprintf('corr(debt.,y) = %6.1f \n',mean(matrix_corr_debt(:,1)))

fprintf('Duration       = %6.1f \n',mean(mean(last_duration)))
%SELECT ONLY SAMPLES IN WHICH THERE ARE SUDDEN STOPS, I.E. THE SPREAD
%DURING SUDDEN STOP TAKES A POSITIVE NUMBER
fprintf('Extra spread during ss = %6.1f \n',100*(mean(spread_ss(ind_ss))-mean(spread_noss(ind_ss))))


