function M=movpwdph(N,Np,typesig);
%MOVPWDPH Influence of a phase-shift on the interferences of the pWVD.  
%	M=MOVPWDPH(N,Np,TYPESIG) generates the movie frames illustrating the 
%	influence of a phase-shift between two signals on the interference 
%	terms of the pseudo Wigner-Ville distribution.
%
%	N : number of points for the signal;
%	Np : number of snapshots (default : 8)
%	TYPESIG : type of signal :
%	 'C' : constant frequency modulation (default value) ;
%	 'L' : linear frequency modulation ;
%	 'S' : sinusoidal frequency modulation.
%	M : matrix of movie frames.
%
%	Example : 
%	 M=movpwdph(128,8,'S'); movie(M,10);

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
 typesig='C'; Np=8;
elseif nargin==2,
 typesig='C'; 
end

M  = moviein(Np);

typesig=upper(typesig);

if typesig=='C',

for k=1:Np,
 sig=fmconst(N,0.2)+fmconst(N,0.3)*exp(j*2*pi*k/Np);
 [tfr,t,f]=tfrpwv(sig); 
 Max=max(max(tfr));V=[0.3 0.5 0.7 0.9]*Max;
 contour(t,f,tfr,V);xlabel('Time'); ylabel('Frequency'); 
 title('Pseudo Wigner-Ville distribution');
 M(:,k) = getframe;
end

elseif typesig=='L',

for k=1:Np,
 sig=fmlin(N,0.1,0.5)+fmlin(N,0,0.4)*exp(j*2*pi*k/Np);
 [tfr,t,f]=tfrpwv(sig); 
 Max=max(max(tfr));V=[0.3 0.5 0.7 0.9]*Max;
 contour(t,f,tfr,V);xlabel('Time'); ylabel('Frequency'); 
 title('Pseudo Wigner-Ville distribution');
 M(:,k) = getframe;
end

elseif typesig=='S',

for k=1:Np,
 sig=fmsin(N,0.2,0.4)+fmsin(N,0.1,0.3)*exp(j*2*pi*k/Np);
 [tfr,t,f]=tfrpwv(sig); 
 Max=max(max(tfr));V=[0.3 0.5 0.7 0.9]*Max;
 contour(t,f,tfr,V);xlabel('Time'); ylabel('Frequency'); 
 title('Pseudo Wigner-Ville distribution'); 
 M(:,k) = getframe;
end

else

error('Wrong input parameter for TYPESIG');

end
