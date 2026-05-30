//Posted to GitHub, with Creative Commons, CC BY 4.0 license.
//Code release date: May 30, 2026.
//We welcome all constructive suggestions/feedback.


clear;clc;
close(winsid())

mp=1.6726*10^-27
mn=1.6749*10^-27
m=(mp*mn)/(mp+mn) //mass of deuteron
omega=1.745*10^22 //angular frequency
hbar=1.0545*(10^-34) //reduced Planck's constant
e=1.602*10^-19
beta1=(sqrt(m*omega/hbar))
V0_RHO=-(19.45D+6)*e//Potential depth

//Parameters in Lennard Jones Model
C0=-3.758E+06*e
C1=15.278E-176*e
C2=2.81E-83*e

h=0.612D-15
r_max=30*10^-15//xmax
r=linspace(h,r_max,1001)
V_RHO=V0_RHO+(m*omega^2*(r.^2)/(2))
f=C0*(ones(1,1001))+(C1./(r.^12))-(C2./(r.^6))-V_RHO

n=10

function Talmi_integral=Talmi(i1,j1)
    integrand=(r.^2).*exp(-(beta1*r).^2).*f.*((beta1*r).^(2*(i1+j1)))
    Talmi_integral=inttrap(r,integrand)
endfunction

function N=Norm_const(k1)
    N=(sqrt(factorial(k1)*2*((beta1)^3)/(gamma(k1+(3/2)))))
endfunction

function ph= pochhammer(a,i)
    if(i==0)
        ph=1
    else
        if(a>=0) 
            if(a==floor(a))
                prod=1
                num=a+i-1
                while(num>=a)
                    prod=prod*num
                    num=num-1
                end
                ph=prod 
            else
                ph=gamma(a+i)/gamma(a)
            end
        else               
            if(i>abs(a))
                ph=0
            else
                ph=((-1)^i)*pochhammer(abs(a)-i+1,i)
            end
        end           
     end  
endfunction

//RHO eigenvalues
for(i=0:n-1)
    epsilon(i+1)=V0_RHO+((2*i+(3/2))*hbar*omega)
end
   
E0=zeros(4)
correction=zeros(n)

//first order correction
E0(1)=(Norm_const(0)*Norm_const(0)*Talmi(0,0))
disp("First order correction: ")
disp(E0(1)/(e*1D+6))
correction(1)=E0(1)

//second order corrections
function ME=Matrix_element(k1,k)
    sumi=0
    for i=0:k
        sumj=0
        for j=0:k1
            sumj=sumj+(((pochhammer(-k,i)*pochhammer(-k1,j))/(pochhammer(3/2,i)*pochhammer(3/2,j)))*Talmi(i,j)/(factorial(i)*factorial(j)))
        end
        sumi=sumi+sumj
    end
    
    ME=Norm_const(k1)*Norm_const(k)*((gamma(k1+3/2)*gamma(k+3/2))/(gamma(k1+1)*gamma(k+1)*gamma(3/2)*gamma(3/2)))*sumi
endfunction
sumk=0
for(k=1:n-1)
    H2k0=Matrix_element(k,0)
    correction(k+1)=((H2k0)^2)/(epsilon(1)-epsilon(k+1))
    sumk=sumk+correction(k+1)
end

E0(2)=sumk
disp("Second order correction")
disp(E0(2)/(e*1D+6))

//Third order correction
sumk1=0
for k1=1:n-1
    sumk=0
    for k=1:n-1
        sumk=sumk+(Matrix_element(0,k1)*Matrix_element(k1,k)*Matrix_element(k,0)/((epsilon(1)-epsilon(k1+1))*(epsilon(1)-epsilon(k+1))))
    end
    sumk1=sumk1+sumk
end
first_term=sumk1
sumk=0
for k=1:n-1
    sumk=sumk+(Matrix_element(k,0)^2/((epsilon(1)-epsilon(k+1))^2))
end
second_term=E0(1)*sumk
E0(3)=first_term-second_term
disp("Third Order correction")
disp(E0(3)/(e*1D+6))


//Fourth order correction

//Fourth order correction

sumk=0
for k=1:n-1
    sumk=sumk+(((Matrix_element(0,k))^2)/(epsilon(k+1)-epsilon(1)))
end
sumk1=0
for k1=1:n-1
    sumk1=sumk1+(Matrix_element(k1,0)^2/((epsilon(k1+1)-epsilon(1))^2))
end
first_term=sumk*sumk1
sumk1=0
for k1=1:n-1
    sumk1=sumk1+(Matrix_element(k1,0)^2/((epsilon(k1+1)-epsilon(1))^3))
end
second_term=(E0(1)^2)*sumk1
sumk1=0
for k1=1:n-1
    sumk=0
    for k=1:n-1
        sumk=sumk+(Matrix_element(0,k1)*Matrix_element(k1,k)*Matrix_element(k,0)*((1/((epsilon(k1+1)-epsilon(1))^2*(epsilon(k+1)-epsilon(1))))+(1/((epsilon(k+1)-epsilon(1))^2*(epsilon(k1+1)-epsilon(1))))))
    end
    sumk1=sumk1+sumk
end
third_term=E0(1)*sumk1
sumk2=0
for k2=1:n-1
     sumk1=0
    for k1=1:n-1
        sumk=0
        for k=1:n-1
            sumk=sumk+(Matrix_element(0,k2)*Matrix_element(k2,k1)*Matrix_element(k1,k)*Matrix_element(k,0)/((epsilon(k2+1)-epsilon(1))*(epsilon(k1+1)-epsilon(1))*(epsilon(k+1)-epsilon(1))))
        end
        sumk1=sumk1+sumk
    end
    sumk2=sumk2+sumk1
end
fourth_term=sumk2   
E0(4)=first_term-second_term+third_term-fourth_term
disp("Fourth Order correction")
disp(E0(4)/(e*1D+6))


disp("retrieved energy for exponential potential")
E=epsilon(1)+E0(1)+E0(2)+E0(3)+E0(4)
disp(E/(e*1D+6))
