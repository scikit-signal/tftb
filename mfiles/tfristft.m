function [x,t] = tfristft(tfr,t,h,trace);
%TFRISTFT Inverse Short time Fourier transform.
%	[X,T]=TFRSTFT(tfr,T,H,TRACE) computes the inverse short-time 
%	Fourier transform of a discrete-time signal X. This function
%	may be used for time-frequency synthesis of signals.
% 
%	X     : signal.
%	T     : time instant(s)          (default : 1:length(X)).
%	H     : frequency smoothing window, H being normalized so as to
%	        be  of unit energy.      (default : Hamming(N/4)). 
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR   : time-frequency decomposition (complex values). The
%	        frequency axis is graduated from -0.5 to 0.5.
%
%	Example :
%        t=200+(-128:127); sig=[fmconst(200,0.2);fmconst(200,0.4)]; 
%        h=hamming(57); tfr=tfrstft(sig,t,256,h,1);
%        sigsyn=tfristft(tfr,t,h,1);
%        plot(t,abs(sigsyn-sig(t)))
% 
%	See also all the time-frequency representations listed in
%	the file CONTENTS (TFR*)

%	F. Auger, November 1996.
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

if (nargin<3),
 error('At least 3 parameters required');
elseif (nargin==3),
 trace=0;
end;

[N,NbPoints]=size(tfr);
[trow,tcol] =size(t);
[hrow,hcol] =size(h); Lh=(hrow-1)/2; 

if (trow~=1),
 error('T must only have one row'); 
elseif (hcol~=1)|(rem(hrow,2)==0),
 error('H must be a smoothing window with odd length');
elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');
elseif (tcol~=NbPoints)
 error('tfr should have as many columns as t has rows.');
end; 

Deltat=t(2:tcol)-t(1:tcol-1); 
Mini=min(Deltat); Maxi=max(Deltat);
if (Mini~=1) & (Maxi~=1),
 error('The tfr must be computed at each time sample.');
end;

h=h/norm(h);

if trace, disp('Inverse Short-time Fourier transform'); end;
tfr=ifft(tfr);

x=zeros(tcol,1);

for icol=1:tcol,
 if trace, disprog(icol,tcol,10); end;
 valuestj=max([1,icol-N/2,icol-Lh]):min([tcol,icol+N/2,icol+Lh]);
 for tj=valuestj,
  tau=icol-tj; indices= rem(N+tau,N)+1; 
  % fprintf('%g %g %g\n',tj,tau,indices);
  x(icol,1)=x(icol,1)+tfr(indices,tj)*h(Lh+1+tau);
 end;
 x(icol,1)=x(icol,1)/sum(abs(h(Lh+1+icol-valuestj)).^2);
end;

if trace, fprintf('\n'); end;
