function [xs,ys] = hex_packed(epsilon,q)
% function [xs,ys] = hex_packed_circle(epsilon,q)
%R = R + epsilon/100;  % rounding issue
% hexagonal basis vectors
a1 = [1 0];
a2 = [cos(2*pi/3) sin(2*pi/3)];

xmax = q;
ymax = q;

[M,N] = meshgrid(-xmax:xmax,-ymax:ymax);
xs = a1(1)*M(:) + a2(1)*N(:);
ys = a1(2)*M(:) + a2(2)*N(:);

%rs = sqrt(xs.^2+ys.^2);
%rtax = abs(M(:)+N(:));
%rtax1 = a1(1)*xs + a1(2)*ys;
%rtax2 = a2(1)*xs + a2(2)*ys;
%rtax = max(abs(rtax1),abs(rtax2));
%rtax1=M(:);
%rtax2=N(:);
%bad = (abs(rtax1)>q*epsilon*1.01)|(abs(rtax2)>q*epsilon*1.01);
%rtax(bad) = [];
%rtax1 = rs;
%bad = rs>1.01*q*epsilon;
%
%rtax = abs(xs)+abs(ys);
%bad = rtax>=q*epsilon*(1-1e-8);
%rtax = sqrt(xs.^2+ys.^2);
%bad = [];
bad = (ys > 1.01*tan(pi/3)*(xs+q)) | (ys < 1.01*tan(pi/3)*(xs-q));
%rtax = ys-(q*sqrt(3)/2+sin(2*pi/3)*(xs+q/2));
%bad = (ys>(q*sqrt(3)/2)) | (ys<-(q*sqrt(3)/2)) | (ys>0.99*(q*sqrt(3)/2-sin(pi/3)*(xs-q/2)));

xs(bad) = [];
ys(bad) = [];
%rtax(bad) = [];
xs = xs';
ys = ys';

xs = xs*epsilon;
ys = ys*epsilon;

%xs = XS(:);

