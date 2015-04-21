function [margt,margf,E]=margtfr(tfr,t,f) 
%MARGTFR Marginals and energy of a time-frequency representation.
%	[MARGT,MARGF,E]=MARGTFR(TFR,T,F) calculates the time and
%	frequency marginals and the energy of a time-frequency
%	representation. 
%
%	TFR : time-frequency representation (M,N)
%	T   : vector containing the time samples in sec. 
%	       (default : (1:N))
%	F   : vector containing the frequency samples in Hz, not
%	       necessary uniformly sampled. (default : (1:M))
%	MARGT : time marginal
%	MARGF : frequency marginal
%	E     : energy of TFR
%
%	Example :    
%	 S=altes(128,0.05,0.45); TFR=tfrscalo(S,1:128,8,'auto');
%	 [MARGT,MARGF,E] = margtfr(TFR); 
%	 subplot(211); plot(T,MARGT); subplot(212); plot(F,MARGF);
%
%	See also MOMTTFR.

%	P. Goncalves, October 95
%	Copyright (c) 1995 Rice University
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

[M,N] = size(tfr) ;
if nargin == 1,
 t = (1:N);
 f = (1:M)';
elseif nargin == 2,
 f = (1:M)';
end  

E     = real(integ2d(tfr,t,f)/M) ;
margt = real(integ(tfr.',f')/M) ;
margf = real(integ(tfr,t)/N) ;
