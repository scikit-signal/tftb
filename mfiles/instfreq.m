function [fnormhat,t]=instfreq(x,t,L,trace);
%INSTFREQ Instantaneous frequency estimation.
%	[FNORMHAT,T]=INSTFREQ(X,T,L,TRACE) computes the instantaneous 
%	frequency of the analytic signal X at time instant(s) T, using the
%	trapezoidal integration rule.
%	The result FNORMHAT lies between 0.0 and 0.5.
% 
%	X : Analytic signal to be analyzed.
%	T : Time instants	        (default : 2:length(X)-1).
%	L : If L=1, computes the (normalized) instantaneous frequency 
%	    of the signal X defined as angle(X(T+1)*conj(X(T-1)) ;
%	    if L>1, computes a Maximum Likelihood estimation of the
%	    instantaneous frequency of the deterministic part of the signal
%	    blurried in a white gaussian noise.
%	    L must be an integer       	(default : 1).
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                (default : 0).
%	FNORMHAT : Output (normalized) instantaneous frequency.
%	T : Time instants.
%
%	Examples : 
%	 x=fmsin(70,0.05,0.35,25); [instf,t]=instfreq(x); plot(t,instf)
%	 N=64; SNR=10.0; L=4; t=L+1:N-L; x=fmsin(N,0.05,0.35,40);
%	 sig=sigmerge(x,hilbert(randn(N,1)),SNR);
%	 plotifl(t,[instfreq(sig,t,L),instfreq(x,t)]); grid;
%	 title ('theoretical and estimated instantaneous frequencies');
%
%	See also  KAYTTH, SGRPDLAY.

%	F. Auger, March 1994, July 1995.
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

if (nargin == 1),
 t=2:xrow-1; L=1; trace=0.0;
elseif (nargin == 2),
 L = 1; trace=0.0;
elseif (nargin == 3),
 trace=0.0;
end;

if L<1,
 error('L must be >=1');
end
[trow,tcol] = size(t);
if (trow~=1),
 error('T must have only one row'); 
end;

if (L==1),
 if any(t==1)|any(t==xrow),
  error('T can not be equal to 1 neither to the last element of X');
 else
  fnormhat=0.5*(angle(-x(t+1).*conj(x(t-1)))+pi)/(2*pi);
 end;
else
 H=kaytth(L); 
 if any(t<=L)|any(t+L>xrow),
  error('The relation L<T<=length(X)-L must be satisfied');
 else
  for icol=1:tcol,
   if trace, disprog(icol,tcol,10); end;
   ti = t(icol); tau = 0:L;
   R = x(ti+tau).*conj(x(ti-tau));
   M4 = R(2:L+1).*conj(R(1:L));
   
   diff=2e-6;
   tetapred = H * (unwrap(angle(-M4))+pi);
   while tetapred<0.0 , tetapred=tetapred+(2*pi); end;
   while tetapred>2*pi, tetapred=tetapred-(2*pi); end;
   iter = 1;
   while (diff > 1e-6)&(iter<50),
    M4bis=M4 .* exp(-j*2.0*tetapred);
    teta = H * (unwrap(angle(M4bis))+2.0*tetapred);
    while teta<0.0 , teta=(2*pi)+teta; end;
    while teta>2*pi, teta=teta-(2*pi); end;
    diff=abs(teta-tetapred);
    tetapred=teta; iter=iter+1;
   end;
   fnormhat(icol,1)=teta/(2*pi);
  end;
 end;
end;

