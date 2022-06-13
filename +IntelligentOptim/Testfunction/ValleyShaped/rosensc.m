function [y] = rosensc(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ROSENBROCK FUNCTION, RESCALED
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, x3, x4]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:4
    xxbar(ii) = 15*xx(ii) - 5;
end

sum = 0;
for ii = 1:3
	xibar = xxbar(ii);
	xnextbar = xxbar(ii+1);
	new = 100*(xnextbar-xibar^2)^2 + (1-xibar)^2;
	sum = sum + new;
end

y = (sum - 3.827*10^5) / (3.755*10^5);

end
