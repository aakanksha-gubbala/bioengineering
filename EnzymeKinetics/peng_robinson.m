
syms x positive;

R = 8.31446261815324;

fprintf('Enter everything in SI units\n');
P = input('Enter pressure: ');
T = input('Enter temperature: ');
Pc = input('Enter critical pressure: ');
Tc = input('Enter critical temperature: '); 
w = input('Enter omega: ');

a = 0.45724*R^2*Tc^2/Pc;
b = 0.07780*R*Tc/Pc;
c = (1+(0.37464+1.5422*w-0.26992*w^2)*(1-(T/Tc)^0.5))^2;

x0 = input('Enter initial volume: ');
f = (R*T/(x-b))-(a*c/(x^2+2*x*b-b^2))-P;
df = diff(f);

error = 0.0000000001;
for i=1:1000
     f0=vpa(subs(f,x,x0)); 
     d0=vpa(subs(df,x,x0)); 
     y=x0-f0/d0; 
     err=abs(y-x0);
     if err<error 
        break
     end
     x0=y;
end 

fprintf('The Root is : %e \n',y);
fprintf('No. of Iterations : %d\n',i);

