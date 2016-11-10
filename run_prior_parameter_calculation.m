% Copyright (C) 2013-2016 Benjamin Born, Sebastian Merkel, and Johannes Pfeifer
% 
%  This is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  It is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  For a copy of the GNU General Public License,
%  see <http://www.gnu.org/licenses/>.


prior_specification = ...%parameter name, prior type, prior median, prior 5% quant., prior 95% quant.
{'theta', 'gamma_pdf', 3.89, 2.57, 5.81;...
'gamma', 'uniform_pdf', 0.5, 0.05, 0.95;...
'delta', 'inv_gamma_pdf', 0.75, 0.32, 2.45;...
'b', 'beta_pdf', 0.5, 0.17, 0.83;...
};

output_cell = cell(size(prior_specification,1),4);
output_cell(:,1:2) = prior_specification(:,1:2);
test_values = NaN(size(prior_specification,1),1);

x_start=[0.5;.5]; %starting value (always)
for i=1:size(prior_specification,1)
  target=[prior_specification{i,4}, prior_specification{i,3}, prior_specification{i,5}]; %target value
  switch prior_specification{i,2}
    case 'beta_pdf'
      [~,~,~,~,~,bestever] = cmaes('get_beta_ab',x_start,0.2,[],target);
      % test whether results solves the problem
      rand_numbers=betarnd(bestever.x(1),bestever.x(2),10000,1);
      quant_values=quantile(rand_numbers,[0.05,0.5,0.95]);
      test_values(i)=sum((target-quant_values).^2);%sum of squared deviations
      [m,v] = betastat(bestever.x(1),bestever.x(2));
    case 'gamma_pdf'
      [~,~,~,~,~,bestever] = cmaes('get_gamma_ab',x_start,0.2,[],target);
      % test whether results solves the problem
      rand_numbers=gamrnd(bestever.x(1),bestever.x(2),10000,1);
      quant_values=quantile(rand_numbers,[0.05,0.5,0.95]);
      test_values(i)=sum((target-quant_values).^2); %sum of squared deviations
      [m,v] = gamstat(bestever.x(1),bestever.x(2));
    case 'inv_gamma_pdf'
      [~,~,~,~,~,bestever] = cmaes('get_inv_gamma_ab',x_start,0.2,[],target);
      % test whether results solves the problem
      rand_numbers=gamrnd(bestever.x(1),1/bestever.x(2),50000,1); % to get inverse gamma, use transformation of gamma
      quant_values=quantile(1./rand_numbers,[0.05,0.5,0.95]);
      test_values(i)=sum((target-quant_values).^2); %sum of squared deviations
      %no function for mean and variance calculation available. Calculate
      %empirical mean and variance
      m = mean(1./rand_numbers);
      v = var(1./rand_numbers);
    case 'uniform_pdf'
      interval_length = (target(3)-target(1))/0.9;
      lower_bound = target(1)-0.05*interval_length;
      upper_bound = target(3)+0.05*interval_length;
      [m,v] = unifstat(lower_bound,upper_bound);
      test_values(i) = abs(m-target(2));
    otherwise
      error('Distribution not defined')
  end
  output_cell{i,3} = m;
  output_cell{i,4} = sqrt(v);
end

filename = 'prior_specification_Dynare.txt';
fid = fopen(filename,'w');
for i=1:size(output_cell,1)
  fprintf(fid,'%s, %s, %.4f, %.4f;\n',output_cell{i,1},output_cell{i,2},output_cell{i,3},output_cell{i,4});
end
fclose(fid);