function [fnorm,t,rejratio]=ifestar2(x,t);
%IFESTAR2 Instantaneous frequency estimation using AR2 modelisation. 
%	[FNORM,T2,RATIO]=IFESTAR2(X,T) computes an estimate of the
%       instantaneous frequency of the real signal X at time
%	instant(s) T. The result FNORM lies between 0.0 and 0.5. This
%	estimate is based only on the 4 last signal points, and has
%	therefore an approximate delay of 2.5 points. 
% 
%	X     : real signal to be analyzed.
%	T     : Time instants (must be greater than 4) 
%					(default : 4:length(X)).
%	FNORM : Output (normalized) instantaneous frequency.
%	T2    : Time instants coresponding to FNORM. Since the
%		algorithm can not always give a value, T2 is 
%		different of T. 
%       RATIO : proportion of instants where the algorithm yields
%		an estimation
%
%	Examples : 
%        [x,if]=fmlin(50,0.05,0.3,5); x=real(x); [if2,t]=ifestar2(x);
%        plot(t,if(t),t,if2);
%
%	 N=1100; [deter,if]=fmconst(N,0.05); deter=real(deter);
%	 noise=randn(N,1); NbSNR=101; SNR=linspace(0,100,NbSNR);
%        for iSNR=1:NbSNR,
%         sig=sigmerge(deter,noise,SNR(iSNR));
%	  [if2,t,ratio(iSNR)]=ifestar2(sig); 
%         EQM(iSNR)=norm(if(t)-if2)^2 / length(t) ;
%        end;
%        subplot(211); plot(SNR,-10.0 * log10(EQM)); grid;
%        xlabel('SNR'); ylabel('-10 log10(EQM)');
%        subplot(212); plot(SNR,ratio); grid;
%        xlabel('SNR'); ylabel('ratio');
%
%	 See also  INSTFREQ, KAYTTH, SGRPDLAY.

%	F. Auger, April 1996.
%       This estimator is the causal version of the estimator called
%       "4 points Prony estimator" in the article "Instantaneous
%	frequency estimation using linear prediction with comparisons
%	to the dESAs", IEEE Signal Processing Letters, Vo 3, No 2, p
%	54-56, February 1996. 
%	Copyright (c) 1996 by CNRS (France).
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

if (nargin == 0),
 error('At least one parameter required');
end;
[xrow,xcol] = size(x);
if (xcol~=1),
 error('X must have only one column');
end
x = real(x);

if (nargin == 1),
 t=4:xrow; 
end;

[trow,tcol] = size(t);
if (trow~=1),
 error('T must have only one row'); 
elseif min(t)<4,
 error('The smallest value of T must be greater than 4');
end;


Kappa = x(t-1) .* x(t-2) - x(t  ) .* x(t-3) ;
psi1  = x(t-1) .* x(t-1) - x(t  ) .* x(t-2) ;
psi2  = x(t-2) .* x(t-2) - x(t-1) .* x(t-3) ;
den   = psi1 .* psi2 ;
indices = find(den>0);
arg=0.5*Kappa(indices)./sqrt(den(indices));
indarg=find(abs(arg)>1);
arg(indarg)=sign(arg(indarg));
fnorm = acos(arg)/(2.0*pi);
rejratio = length(indices)/length(t);
t = t(indices);

