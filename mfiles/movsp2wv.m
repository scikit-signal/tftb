function M=movsp2wv(N,Np);
%MOVSP2WV Movie illustrating the passage from the spectrogram to the WVD.
%	M=MOVSP2WV(N,Np) generates the movie frames illustrating the passage
%	from the spectrogram to the WVD using different smoothing gaussian
%	windows in the smoothed pseudo-WVD. 
%
%	N : number of points for the signal;
%	Np : number of snapshots (default : 8)
%	M : matrix of movie frames.
%
%	Example : 
%	 M=movsp2wv(128,15); movie(M,5);

%	O. Lemoine - June 1996.
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

if nargin<1,
 error('At least one argument required');
elseif nargin==1,
 Np=8;
end
Np=odd(Np)-1;

M  = moviein(Np);

t1=N/4; t2=3*N/4; t3=t1; t4=t2; 
f1=0.125; f2=f1; f3=0.375; f4=f3; 
t=1:N; 
T=sqrt(2*N);

sig1=amgauss(N,t1,T).*fmconst(N,f1,t1);  
sig2=amgauss(N,t2,T).*fmconst(N,f2,t2);  
sig3=amgauss(N,t3,T).*fmconst(N,f3,t3);  
sig4=amgauss(N,t4,T).*fmconst(N,f4,t4);  

sig=sig1+sig2+sig3+sig4;

h=tftb_window(49,'gauss');
[tfr,t,f]=tfrsp(sig,t,N,h); tfr=tfr(1:N/2,:); f=f(1:N/2);
Max=max(max(tfr)); V=[0.1 0.3 0.5 0.7 0.9]*Max;
contour(t,f,tfr,V);
xlabel('Time'); ylabel('Frequency'); axis('xy')
M(:,1) = getframe;

Ng0=odd(fliplr(linspace(3,N/4,Np/2-1)));
Nh0=odd(linspace(N/2,9*N/10,Np/2-1));

for k=2:Np/2,
 g=tftb_window(Ng0(k-1),'gauss');
 h=tftb_window(Nh0(k-1),'gauss');
 [tfr,t,f]=tfrspwv(sig,t,N,g,h);
 Max=max(max(tfr)); V=[0.1 0.3 0.5 0.7 0.9]*Max;
 contour(t,f,tfr,V);
 xlabel('Time'); ylabel('Frequency'); axis('xy')
 M(:,k) = getframe;
end

[tfr,t,f]=tfrwv(sig);
Max=max(max(tfr)); V=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]*Max;
contour(t,f,tfr,V);
xlabel('Time'); ylabel('Frequency'); axis('xy')
M(:,Np/2+1) = getframe;

for k=Np/2+2:Np,
 M(:,k) =M(:,Np+2-k);
end
