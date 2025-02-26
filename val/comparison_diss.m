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
cudec512e=importdata('cudec_512_eulero.txt');
single128ab=importdata('single_128_ab2.txt');
single256ab=importdata('single_256_ab2.txt');


figure(1)
clf
plot(ref(:,1),ref(:,2),'-k','LineWidth',3.0,'DisplayName','Wim M. van Rees (2011) + PS + 256^3')
hold on
plot(cans(:,1),cans(:,2),'-b','LineWidth',3.0,'DisplayName','CaNS (2018) + FD2 + 512^3')
plot(cudec128e(:,1),cudec128e(:,2),'^g','LineWidth',10.0,'DisplayName','MHIT36 + cuDecomp + Euler + 128^3')
plot(cudec512e(:,1),cudec512e(:,2),'*g','LineWidth',10.0,'DisplayName','MHIT36 + cuDecomp + Euler + 512^3')
plot(single128ab(:,1),single128ab(:,2),'^r','LineWidth',6.0,'DisplayName','MHIT36 + single + AB2 + 128^3')
plot(single256ab(:,1),single256ab(:,2),'sr','LineWidth',10.0,'DisplayName','MHIT36 + single + AB2 + 256^3')
set(gca,'Fontsize',20)
ylabel('Dissipation')
xlabel('Time')
legend show
hold off

