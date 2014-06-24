function B = GMfreq(omega,f);
% function B = GMfreq(omega,f);
% Return the frequency part of the GM spectrum.
%

% $Id: gmFreq.m,v 1.1 2005/06/18 23:38:26 jklymak Exp jklymak $ 

B = (2/pi)*f./omega;
B = B.*(omega.^2-f.^2).^(-0.5);

% $Log: gmFreq.m,v $
% Revision 1.1  2005/06/18 23:38:26  jklymak
% Initial revision
%