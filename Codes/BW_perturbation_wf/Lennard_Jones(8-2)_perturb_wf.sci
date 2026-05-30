//Posted to GitHub, with Creative Commons, CC BY 4.0 license.
//Code release date: May 30, 2026.
//We welcome all constructive suggestions/feedback.


clear;clc;
close(winsid())

//Parameters in Lennard Jones Model
C0=3.529E+06
C1=1.33E-114
C2=7.7E-23

mp=1.6726*10^-27
mn=1.6749*10^-27
m=(mp*mn)/(mp+mn) //mass of deuteron
e=1.602*10^-19
hbar=1.0545*(10^-34)

N1=1001
N=N1-2
T=zeros(N,N)
V_LJ=zeros(N,N)

a=10*10^-15
h=a/(N+1)
x(1)=h

//For Hamiltonian matrix
for i=1:N
    T(i,i)=-2
    if(i<N)
        x(i+1)=x(i)+h
        T(i,i+1)=1
        T(i+1,i)=1
    end
    V_LJ(i,i)=C0+(C1/(x(i)^8))-(C2/(x(i)^2))
    V1_LJ(i)=V_LJ(i,i) 
end

H_LJ=(((-hbar^2/(2*m*e*h^2))*T)+V_LJ)/(10^6)
[vec_LJ,E_LJ]=spec(H_LJ)

//For bound state
count=0
for i=1:N
    if E_LJ(i,i)<0
        count=count+1
        E1_LJ(count)=E_LJ(i,i)
        phi_LJ(:,count)=vec_LJ(:,i)
    end
end

//For the plot
x_plot=x/10^-15
scf(0)
clf
plot(x_plot,V1_LJ/10^6,'k','LineWidth',2')
for(i=1:count)
    plot(x_plot,E1_LJ(i)/10^6,'r-.','LineWidth',2)
end
aaa=gca
aaa.data_bounds=[0.3,min(V1_LJ)/10^6;3,100]
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4
legend('Lennard Jones potential','Bound state')
xtitle("Lennard Jones Potential",'r in fm','V in MeV')

//results
disp("no. of bound states, for this Lennard Jones potential = "+string(count))
disp("Lennard Jones bound state energy:")
disp(E1_LJ)
disp("Size of the nucleus")
[val,index]=max(phi_LJ(:,1).*phi_LJ(:,1))
disp(x(index))

//Radial Harmonic Oscillator Model - Unperturbed Hamiltonian

//Parameters for RHO
omega=1.745*10^22 //angular frequency
k=m*(omega^2) //force constant
V0_RHO=-19.45D+6//Potential depth

T=zeros(N,N)
V_RHO=zeros(N,N)
hcut=1.05457*10^-34
e=1.602177*10^-19
l=0

r(1)=h 
//For the Hamiltonian matrix
for i=1:N
    T(i,i)=-2
    if(i<N)
        r(i+1)=r(i)+h
        T(i,i+1)=1
        T(i+1,i)=1
    end
    V_RHO(i,i)=V0_RHO+(((k*(r(i)^2)/(2))+(l*(l+1)*hcut^2/(2*m*r(i)^2)))/e)//Modified RHO model
    V1_RHO(i)=V_RHO(i,i) 
end

H_RHO=(((-hcut^2/(2*m*e*h^2))*T)+V_RHO)/(10^6)
[vec_RHO,E_RHO]=spec(H_RHO)
for i=1:N
    E1_RHO(i)=E_RHO(i,i)
end
psi_RHO=vec_RHO

disp("Deuteron as per RHO model--Binding Energy")
disp(E1_RHO(1))

//Wavefunction plot
figure(1)
clf()
plot(r,psi_RHO(:,1),'--','LineWidth',2)
plot(x,phi_LJ(:,1),'r','LineWidth',2)
legend(['RHO';'Lennard-Jones']);
xtitle("RHO vs LJ wavefunctions (un)",'r','u')
aaa=gca
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4

//Potential plot
figure(2)
clf()
plot(x,V1_LJ,'--','LineWidth',2)
plot(r,V1_RHO,'r','LineWidth',2)
legend(['Lennard Jones';'RHO']);
xtitle("RHO vs LJ Potentials",'r','V')
aaa=gca
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4
aaa.data_bounds=[0.3*10^-15,min(V1_LJ);10*10^-15,100*10^6]

//bw perturbation
Pn=psi_RHO(:,1)*(psi_RHO(:,1))'
Qn=eye(N,N)-Pn
En=E1_LJ(1)
EminusH0 = En*eye(N,N)-H_RHO;
Rn1=inv(EminusH0)*Qn
Rn11=Qn*inv(EminusH0)

H1 = H_LJ-H_RHO//Difference Hamiltonian
psi_n=inv(eye(N,N)-(Rn1*H1))*psi_RHO(:,1)

//wavefunctions plot
figure(3)
clf()
psi_RHO_norm=psi_RHO(:,1)./sqrt(inttrap(r,psi_RHO(:,1).*psi_RHO(:,1)))
psi_LJ_norm=phi_LJ(:,1)./sqrt(inttrap(x,phi_LJ(:,1).*phi_LJ(:,1)))
psi_n_norm=psi_n./sqrt(inttrap(r,psi_n.*psi_n))
plot(r/(1D-15),psi_RHO_norm(:,1)/(3.162D+7),'--','LineWidth',2)
plot(r/(1D-15),psi_LJ_norm(:,1)/(3.162D+7),'r','LineWidth',2)
plot(r/(1D-15),psi_n_norm/(3.162D+7),':b','LineWidth',2')
xlabel("$\text{radial coordinate r (fm)}$", "fontsize", 4)
ylabel("$\text{u(r)}\ ( fm )^{(-1/2)}$", "fontsize", 4)
legend("$u_{\text{RHO}}$","$u_{\text{LJ (8-2))}}$","$u_{\text{BW perturbation}}$")
lll = gce(); 
lll.legend_location = "in_upper_right"; 
lll.font_size = 4;        
aaa = gca()
aaa.font_size = 4; 
csvWrite([r psi_RHO_norm(:,1) psi_LJ_norm(:,1) psi_n_norm ],"LJdata.csv")
