//Posted to GitHub, with Creative Commons, CC BY 4.0 license.
//Code release date: May 30, 2026.
//We welcome all constructive suggestions/feedback.


clear 
clc
close(winsid())

//Exponential potential, and unnormalized WF
//Parameters in Exponential
a0=0.93*10^-15
A=111*10^6 //in eV

mp=1.6726*10^-27
mn=1.6749*10^-27
m=(mp*mn)/(mp+mn)//mass of deuteron
e=1.602*10^-19
hbar=1.0545*(10^-34)

N1=1001
N=N1-2
T=zeros(N,N)
V_Fl=zeros(N,N)


a=20*10^-15
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
    V_Fl(i,i)=-A*exp(-x(i)/a0)
    V1_Fl(i)=V_Fl(i,i)
end

H_Fl=(((-hbar^2/(2*m*e*h^2))*T)+V_Fl)/(10^6)
[vec_Fl,E_Fl]=spec(H_Fl)

//For bound states
count=0
for i=1:N
    if E_Fl(i,i)<0
        count=count+1
        E1_Fl(count)=E_Fl(i,i)
        phi_Fl(:,count)=vec_Fl(:,i)
    end
end

//potential
scf(0)
clf(0)
plot(x/10^-15,V1_Fl/10^6,'LineWidth',2)
xtitle('Exponential potential','r(fm)','Energy in MeV')
aaa=gca
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4


//results
disp("no. of bound states, for this Exponential potential = "+string(count))
disp("Exponential Bound state energy:")
disp(E1_Fl)
disp("Size of the nucleus")
[val,index]=max(phi_Fl(:,1).*phi_Fl(:,1))
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
    V_RHO(i,i)=V0_RHO+(((k*(r(i)^2)/(2))+(l*(l+1)*hcut^2/(2*m*r(i)^2)))/e)
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

//Wavefunctions plot
figure(1)
clf()
plot(r,psi_RHO(:,1),'--','LineWidth',2)
plot(x,-phi_Fl(:,1),'r','LineWidth',2)
legend(['RHO';'Exponential']);
xtitle("RHO vs Exponential wavefunctions (u)")
aaa=gca
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4

//Potential plot
figure(2)
clf()
plot(x,V1_Fl,'--','LineWidth',2)
plot(r,V1_RHO,'r','LineWidth',2)
legend(['Exponential';'RHO']);
xtitle("RHO vs Exponential Potentials")
aaa=gca
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4

//bw perturbation
Pn=psi_RHO(:,1)*(psi_RHO(:,1))'
Qn=eye(N,N)-Pn
En=E1_Fl(1)
EminusH0 = En*eye(N,N)-H_RHO;
Rn1=inv(EminusH0)*Qn
Rn11=Qn*inv(EminusH0)

H1 = H_Fl-H_RHO
psi_n=inv(eye(N,N)-(Rn1*H1))*psi_RHO(:,1)

//Wavefunctions
figure(3)
clf()
psi_RHO_norm=psi_RHO(:,1)./sqrt(inttrap(r,psi_RHO(:,1).*psi_RHO(:,1)))
psi_Fl_norm=phi_Fl(:,1)./sqrt(inttrap(x,phi_Fl(:,1).*phi_Fl(:,1)))
psi_n_norm=psi_n./sqrt(inttrap(r,psi_n.*psi_n))
plot(r/(1D-15),psi_RHO_norm(:,1)/(3.162D+7),'--','LineWidth',2)
plot(r/(1D-15),-psi_Fl_norm(:,1)/(3.162D+7),'r','LineWidth',2)
plot(r/(1D-15),psi_n_norm/(3.162D+7),'b:','LineWidth',2)
xlabel("$\text{radial coordinate r (fm)}$", "fontsize", 4)
ylabel("$\text{u(r)}\ (fm )^{(-1/2)}$", "fontsize", 4)
legend("$u_{\text{RHO}}$","$u_{\text{Exponential}}$","$u_{\text{BW perturbation}}$")
lll = gce(); 
lll.legend_location = "in_upper_right"; 
lll.font_size = 4;      
aaa = gca()
aaa.font_size = 4; 
csvWrite([r psi_RHO_norm(:,1) psi_Fl_norm(:,1) psi_n_norm ],"Fluggedata.csv")
