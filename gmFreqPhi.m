function B = GMfreq(phi,f);
% function B = GMfreq(omega,f);
% Return the frequency part of the GM spectrum.
%

% $Id$ 

B = (2/pi/f)*cos(phi).*cot(phi);
