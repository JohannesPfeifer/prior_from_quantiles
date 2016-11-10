function outvalue=get_inv_gamma_ab(x,target)
% function outvalue=get_inv_gamma_ab(x,target)
% provides parameters a and b for inverse gamma distribution that matches targets at
% 0.05, 0.5, and 0.95 quantiles
%
% Inputs:
% - x       [2 by 1] vector with starting values for a and b
% - target  [3 by 1] vector with targets for quantiles
% Inputs:
% - outvalue  [2 by 1] vector with values for a and b (uses the Matlab parameterization)
%
% Copyright (C) 2013-2016 Benjamin Born and Johannes Pfeifer
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

randnumbers=gamrnd(x(1),1/x(2),50000,1);
quant_values=quantile(1./randnumbers,[0.05,0.5,0.95]);
outvalue=sum((target-quant_values).^2);