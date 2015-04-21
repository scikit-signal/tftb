function M=movcw4at(N,Np);
%MOVCW4AT Four atoms rotating, analyzed by the Choi-Williams distribution.
%	M=MOVCW4AT(N,Np) generates the movie frames illustrating the passage 
%	from the spectrogram to the WVD using different smoothing gaussian
%	windows in the smoothed pseudo-WVD. 
%
%	N : number of points of the analyzed signal
%	Np : number of snapshots (default : 7)
%	M : matrix of movie frames.
%
%	Example : 
%	 M=movcw4at(128); movie(M,5);

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

if nargin==0,
 error('The number of parameters must be at least one');
elseif nargin==1, 
 Np=7;
end

M   = moviein(Np);
rho = N/3;
T   = sqrt(N);

for k=1:Np,
 theta=(k-1)*pi/(2*(Np-1));
 t1=N/2+rho*cos(theta);
 f1=0.5*(N/2+rho*sin(theta))/N;
 t2=N/2+rho*cos(theta+pi/2);
 f2=0.5*(N/2+rho*sin(theta+pi/2))/N;
 t3=N/2+rho*cos(theta+pi);
 f3=0.5*(N/2+rho*sin(theta+pi))/N;
 t4=N/2+rho*cos(theta+3*pi/2);
 f4=0.5*(N/2+rho*sin(theta+3*pi/2))/N;
 sig1=amgauss(N,t1,T).*fmconst(N,f1,round(t1));  
 sig2=amgauss(N,t2,T).*fmconst(N,f2,round(t2));  
 sig3=amgauss(N,t3,T).*fmconst(N,f3,round(t3));  
 sig4=amgauss(N,t4,T).*fmconst(N,f4,round(t4));  
 sig=sig1+sig2+sig3+sig4;

 [tfr,t,f]=tfrcw(sig,1:N,N,tftb_window(N/2+1),tftb_window(N+1),1);
 Max=max(max(tfr)); V=[0.1 0.3 0.5 0.7 0.9]*Max;
 contour(t,f,tfr,V);pause(.1)
 xlabel('Time'); ylabel('Frequency'); axis('xy')
 title('Choi-Williams distribution');
 M(:,k) = getframe;
end

