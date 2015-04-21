function M=movwv2at(N,Np);
%MOVWV2AT Oscillating structure of the interferences of the WVD.  
%	M=MOVWV2AT(N,Np) generates the movie frames illustrating the 
%	influence of the distance between two components on the oscillating
%	structure of the interferences of the WVD.  
%
%	N : number of points for the signal;
%	Np : number of snapshots (default : 9)
%	M : matrix of movie frames.
%
%	Example : 
%	 M=movwv2at(128,15); movie(M,10);

%	O. Lemoine - May 1996.
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
 Np=9;
end
Np=odd(Np);

M  = moviein(Np);

c  = Np/0.2;

ampl=amgauss(N/2,N/2,N/2); 

for k=1:(Np+1)/2,
 sig=[ampl.*fmconst(N/2,0.25-(k-1)/c);ampl.*fmconst(N/2,0.25+(k-1)/c)];
 [tfr,t,f]=tfrwv(sig); 
 Max=max(max(tfr));V=[0.1 0.3 0.5 0.7 0.9]*Max;
 contour(t,f,tfr,V);xlabel('Time'); ylabel('Frequency'); 
 title('Wigner-Ville distribution'); axis('xy')
 M(:,k) = getframe;
end

for k=(Np+3)/2:Np,
 M(:,k) = M(:,Np+1-k);
end


