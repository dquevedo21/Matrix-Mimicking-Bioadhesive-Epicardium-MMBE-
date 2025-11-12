%Description: This script generates a txt file that can be used as gcode
%to generate the route of the nozzle to print a sinusoidal fiber patch.
%Written by Diego Quevedo-Moreno (07/08/2020)
clc, clear all, close all
yv=linspace(0,10*pi,500);
Ns=3;%Number of Samples
A=1.5;%Amplitude(mm)
P=.25;%Density 
WL=5;%Wave Length(pi/WL=wavelength)(mm) 
W=.150;%Fiber Width(mm)
SL=50;%Sample Length(mm)
SW=5;%Sample Width(mm)
nL=10;%#Layers
Lt=.025;%Layer thickness
iL=0;%Initial Layer Thickness
dp=.5;%Distance Between points
F= 400;%Printing Speed(mm/min)
nf=(SW/W)*P;%#Fibers
%df=(SW-(2*A)-W)/(nf-1);%Distance between fibers
df=(SW-W)/(nf-1);%Distance between fibers version 2 
Offset=[-20 0];%%Offset if needed [X Y] (mm)
x=0:dp:SL;

%Sinusoidal Generator
for j=1:nf
    for i=1:length(x)
            y(i,j)=A*sin(x(i)*(pi)/WL)+(j-1)*df;
    end
end
%Sinusoidal Pattern Generator
for j=0:nf-1
    if rem(j,2)==0 %Odd direction 
        for i=1:length(x)
            xt(i+(length(x)*j))=x(i);
            yt(i+(length(x)*j))=y(i,(j+1));
        end
    else %Even direction
        for i=1:length(x)
            xt(i+(length(x)*j))= x((SL/dp+1)-i+1);
            yt(i+(length(x)*j))=y(((SL/dp+1)-i+1),(j+1));
        end
    end
end

%%Verical Support Lines
NumVL=0:1:SL/WL*2;
for k=1:length(NumVL)
    if rem(k,2)==1
        vsx(k)=WL/2*(k-1);
        vsx(k+1)=WL/2*(k-1);
        vsy(k)=SW*(k-1);
        vsy(k+1)=SW*(k-1);
    end
end
for k=1:length(vsy)
    if rem((vsy(k)/10),2)==1
        vsy(k)=SW+A+1-4*W;
    else 
        vsy(k)=-A-1-4*W;
    end
end
% vsy(1)=0;
vsy=[0 vsy];
%%Skirt 
svx=[-4 -4 -3 -3 -2 -2 -1 -1];
svy=[0 SW SW 0 0 SW SW 0];
svz=[iL iL iL iL iL iL iL iL];

xtf=flip(xt);
ytf=flip(yt);
%Add Vertical Support lines
xt=[xt vsx];
yt=[yt vsy];
xtf=[xtf flip(vsx)];
ytf=[ytf vsy];
for k=1:nL
    if rem(k,2)==1
        for i=1:length(xt)
            gx(i+length(xt)*(k-1))=xt(i);
            gy(i+length(xt)*(k-1))=yt(i);
            gz(i+length(xt)*(k-1))= ((k-1)*Lt)+iL;
        end
    else
        for i=1:length(xt)
            gx(i+length(xt)*(k-1))=xtf(i);
            gy(i+length(xt)*(k-1))=ytf(i);
            gz(i+length(xt)*(k-1))=((k-1)*Lt)+iL;
        end
    end
end

%Grip Zone
gzvx=0:W/2:10;
for k=1:2:length(gzvx)
gzx(k)=gzvx(k);
gzy(k)=SW+A+1+W;
gzz(k)=Lt*nL;
if k>1
gzx(k-1)=gzvx(k);
gzy(k-1)=-A-1-4*W;
gzz(k-1)=Lt*nL;
end
end

%Horizontal Support Lines
hvy=[0 SW+A+1 SW+A+1 SW+A+1-W SW+A+1-W SW+A+1-2*W SW+A+1-2*W SW+A+1-3*W SW+A+1-3*W SW+A+1-4*W SW+A+1-4*W];
hvy2=[-A-1-4*W -A-1-4*W -A-1-3*W -A-1-3*W -A-1-2*W -A-1-2*W -A-1-W -A-1-W];
hvx=[0 0 SL SL 0 0 SL SL 0 0 SL];
hvx2=[SL 0 0 SL SL 0 0 SL];
hvz=[Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt Lt];

hvx=[hvx hvx2];
hvy=[hvy hvy2];

% Skirt + patch + horizontal support lines + grip zone
gx=[svx gx hvx gzx gzx+SL-10];
gy=[svy gy hvy gzy gzy];
gz=[svz gz hvz gzz gzz];

%skirt+patch+gripzone
% gx=[svx gx gzx gzx+SL-10];
% gy=[svy gy gzy gzy];
% gz=[svz gz gzz gzz];

%%Just grip zone
% gx=[gzx gzx+SL-10];
% gy=[gzy gzy];
% gz=[gzz gzz];

%Skirt + patch 
% gx=[svx gx];
% gy=[svy gy];
% gz=[svz gz];

 
%Offset (if needed)
gx=gx + Offset(1);
gy=gy + Offset(2);



for j=1:Ns
    for i=1:length(gx)
        gtx(i+(length(gx)*(j-1)))=gx(i);
        gty(i+(length(gx)*(j-1)))=(gy(i)+(SW+8)*(j-1));
        gtz(i+(length(gx)*(j-1)))=gz(i);
    end 
end
figure
axis equal 
grid on
hold on
plot3(gtx,gty,gtz,'-');
xlabel('mm');
ylabel('mm');

fid = fopen( 'C:\Users\Diego Quevedo\Dropbox (MIT)\MIT Research Internship\Project Shared Folder\Modeling\3DP\10.27.2020\SampleTest_p50_a0_p5_50x5mm_F400_dp.5_L10_3sam_offset_fw150_Lt25_vl_hl_gz_V.gcode','w');
fprintf(fid,'F%8.6f\r\n',F);

for j = 1:length(gtx)

    fprintf(fid,'LINEAR X%8.6f Y%8.6f Z%8.6f\r\n',gtx(j),gty(j),gtz(j));
    
end

fclose('all');