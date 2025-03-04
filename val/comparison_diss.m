clear all
nx=256;
dx=2*pi/(nx-1);
ns=6;
nu = 6.25e-4;
ut=zeros(nx,nx,nx,ns);
vt=zeros(nx,nx,nx,ns);
wt=zeros(nx,nx,nx,ns);
diss=zeros(ns,1);
dt=0.0002;

ref=importdata('dissref.txt');
cans=importdata('cans.txt');
cudec128e=importdata('cudec_128_eulero.txt');
cudec256e=importdata('cudec_256_eulero.txt');
cudec512e=importdata('cudec_512_eulero.txt');
cudec128ab=importdata('cudec_128_ab.txt');
cudec256ab=importdata('cudec_256_ab.txt');
cudec512ab=importdata('cudec_512_ab.txt');
cudec1024ab=importdata('cudec_1024_ab.txt');
single128ab=importdata('single_128_ab.txt');
single256ab=importdata('single_256_ab.txt');
single512ab=importdata('single_512_ab.txt');


figure(1)
clf
plot(ref(:,1),ref(:,2),'-k','LineWidth',3.0,'DisplayName','Wim M. van Rees (2011) + PS + 256^3')
hold on
plot(cans(:,1),cans(:,2),'-b','LineWidth',3.0,'DisplayName','CaNS (2018) + FD2 + 512^3')
plot(cudec128e(:,1),cudec128e(:,2),'^g','LineWidth',3.0,'DisplayName','MHIT36 + cuDecomp + Euler + 128^3')
plot(cudec256e(:,1),cudec256e(:,2),'sg','LineWidth',10.0,'DisplayName','MHIT36 + cuDecomp + Euler + 256^3')
plot(cudec512e(:,1),cudec512e(:,2),'*g','LineWidth',9.0,'DisplayName','MHIT36 + cuDecomp + Euler + 512^3')
plot(cudec128ab(:,1),cudec128ab(:,2),'^m','LineWidth',3.0,'DisplayName','MHIT36 + cuDecomp + AB + 128^3')
plot(cudec256ab(:,1),cudec256ab(:,2),'sm','LineWidth',7.0,'DisplayName','MHIT36 + cuDecomp + AB + 256^3')
plot(cudec512ab(:,1),cudec512ab(:,2),'*m','LineWidth',7.0,'DisplayName','MHIT36 + cuDecomp + AB + 512^3')
plot(cudec1024ab(:,1),cudec1024ab(:,2),'xm','LineWidth',10.0,'DisplayName','MHIT36 + cuDecomp + AB + 1024^3')
plot(single128ab(:,1),single128ab(:,2),'^r','LineWidth',3.0,'DisplayName','MHIT36 + single + AB + 128^3')
plot(single256ab(:,1),single256ab(:,2),'sr','LineWidth',7.0,'DisplayName','MHIT36 + single + AB + 256^3')
plot(single512ab(:,1),single512ab(:,2),'*r','LineWidth',7.0,'DisplayName','MHIT36 + single + AB + 512^3')
set(gca,'Fontsize',20)
ylabel('Dissipation')
xlabel('Time')
legend show
hold off

