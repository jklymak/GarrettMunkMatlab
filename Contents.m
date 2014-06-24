% Toolbox for calculating variations of GM spectra.  
% 
% $Id: Contents.m,v 1.1 2008/11/25 22:15:18 jklymak Exp jklymak $
% jklymak@uvic.ca  Jody Klymak
%
% For more information, please see html/index.html
%  
% In general, use the capitalized versions to get the spectra for
% various fields.  Some fields are not fully supported yet.
% 
% GmKx.m
% 	% function S = GmKx(quant,kx,f,N,params);
% 	% Return 1-D horizontal wave-number spectra. 
% GmKxKz.m
% 	% function S = GmKxKz(quant,kx,kz,f,N,params);
% 	% Return 2-D horizontal-vertical wavenumber spectra. 
% GmKz.m
% 	% function S = GmOmKz(quant,om,kz,f,N,params);
% 	% Return frequency-vertical wavenumber spectra. 
% GmOm.m
% 	% function S = GmOm(quant,om,f,N,params);
% 	% Return frequency-vertical wavenumber spectra. 
% GmOmKz.m
% 	% function S = GmOmKz(quant,om,kz,f,N,params);
% 	% Return frequency-vertical wavenumber spectra. 
% Gk91Params.m
% 	% function params = Gk91Params;
% 	% Parameters for the Gregg Kunze '91 vertsion of GM.
% Gm75Params.m
% 	% function params = Gm75Params;
% Gm76Params.m
% 	% function params = Gm76Params;
% gmFreq.m
% 	% function B = GMfreq(omega,f);
% 	% Return the frequency part of the GM spectrum.
% gmFreqPhi.m
% 	% function B = GMfreq(omega,f);
% 	% Return the frequency part of the GM spectrum.
% gmHorizPhi.m
% 	% function A = GMhorizPhi(phi,f,jstar,jp,N,b,N0,z);
% 	% Return the frequency part of the GM spectrum.
% gmVert.m
% 	% function [A] = gMVertKz(kz,f,jstar,jp,N,b,N0,I,s,t,varargin);
% 	% return the vertical wavenumber part of the spectrum.  
