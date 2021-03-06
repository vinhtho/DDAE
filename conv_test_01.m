E= 1;
A= 1;

B= 1;


xe = @(t)exp(0.1*t);
xed = @(t)0.1*exp(0.1*t);

tau = pi;

f=@(t) E*xed(t)-A*xe(t)-B*xe(t-tau);
E=@(t)E;
A=@(t)A;
B=@(t)B;

phi=xe;

tspan = [0,10];

L=6;

h=2*2.^-(1:L);

X=zeros(1,L);

options.isConst = 1;

for i=1:L
    options.StepSize=h(i);
    [t,x]=solve_lddae(E,A,B,f,tau,phi,tspan,options);
    X(i)=x(end);
end

dX=abs(X-xe(tspan(2)));

plot(log(dX(1:end-1)./dX(2:end))/log(2))
figure
semilogy(dX)