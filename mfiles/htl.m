function [ht,rho,theta]=htl(IM,M,N,trace)
% HTL	Hough transform for detection of lines in images.
%	[HT,RHO,THETA]=HTL(IM,M,N,TRACE).
%	From an image IM, computes the integration of the values
%	of the image over all the lines. The lines are parametrized 
%	using polar coordinates. The origin of the coordinates is fixed
%	at the center of the image, and theta is the angle between the
%	VERTICAL axis and the perpendicular (to the line) passing through 
%	the origin. Only the values of IM exceeding 5 % of the maximum 
%	are taken into account (to speed up the algorithm). 
%
%	IM    : image to be analyzed (size Xmax x Ymax).
%	M     : desired number of samples along the radial axis 
%					(default : Xmax).
%	N     : desired number of samples along the azimutal (angle) axis
%					(default : Ymax). 
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                (default : 0).
%	HT    : output matrix (MxN matrix). When called without 
%	        output arguments, HTL displays HT using mesh.
%	RHO   : sequence of samples along the radial axis.
%	THETA : sequence of samples along the azimutal axis.
%
%	Example :
%	 N=64; t=(1:N); y=fmlin(N,0.1,0.3); 
%	 IM=tfrwv(y,t,N); imagesc(IM); pause(1); htl(IM,N,N,1); 

%	O. Lemoine - June 1995.
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
[Xmax,Ymax] = size(IM);

if (nargin == 1),
 M=Xmax; N=Ymax; trace=0;  
elseif (nargin == 2),
 N=Ymax; trace=0;
elseif (nargin == 3),
 trace = 0;
end;

rhomax=sqrt(Xmax^2+Ymax^2)/2;

deltar=rhomax/(M-1);
deltat=2*pi/N;

ht=zeros(M,N);
Max=max(max(IM)); 

if rem(Xmax,2)~=0,
  Xc=(Xmax+1)/2; X0=1-Xc; Xf=Xc-1;
else
  Xc=Xmax/2; X0=1-Xc; Xf=Xc;
end
if rem(Ymax,2)~=0,
  Yc=(Ymax+1)/2; Y0=1-Yc; Yf=Yc-1;
else
  Yc=Ymax/2; Y0=1-Yc; Yf=Yc;
end

if trace, disp('Hough transform - Detection of lines');end
for x=X0:Xf,
 if trace, disprog(x-X0+1,Xmax,10); end
  for y=Y0:Yf,
    if abs(IM(x+Xc,y+Yc))>Max/20,
      for theta=0:deltat:(2*pi-deltat),
        rho=x*cos(theta)-y*sin(theta);
        if ((rho>=0)&(rho<=rhomax)),
          ht(round(rho/deltar)+1,round(theta/deltat)+1)=...
	 	 ht(round(rho/deltar)+1,round(theta/deltat)+1)+...
	 	 IM(x+Xc,y+Yc);
        end
      end
    end
  end
end

rho=0:deltar:rhomax;
theta=0:deltat:(2*pi-deltat);
if trace, disp(' '); end

if nargout==0,
 mesh(theta,rho,ht)
 shading('flat')
 title('Hough transform - Detection of lines');
 Min=min(min(ht));
 Max=max(max(ht));
 axis([min(theta) max(theta) min(rho) max(rho) Min Max]);
 xlabel('Theta');
 ylabel('Rho');
end

