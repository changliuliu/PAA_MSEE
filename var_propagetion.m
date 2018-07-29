% The propogation of variance in PAA 
% 
% Identify the time varying paramters from a linear system and bound the
% prediction error with the MSEE value
% 
% 2014.4.28 
% Changliu Liu


clear;

%% Problem formulation
n=1000;
for i=1:n
    if i<500
        A(i)=1-exp(-i/100);
    else
        A(i)=exp(-(i-500)/400);
    end
    if i>1
        dC(1:2,i-1)=[abs(A(i-1)-A(i));0];
    end
    B(i)=-0.5;
end
Bd=0.005;
DC=[mean(dC(1,:));0];
x(1)=1; u(1)=1;
xhat(1)=1; 
theta(1:2,1)=[0;0];
F={};
F{1}=eye(2);
lambda=0.98; % Forgetting factor


%% Parameter Adaptation
for i=2:n

    x(i)=A(i)*x(i-1)+B(i)*u(i-1)+Bd*randn(1);
    
    phi=[x(i-1);u(i-1)];
    
    xhat(i)=theta(1:2,i-1)'*phi;
    
    F{i}=(F{i-1}-F{i-1}*phi*phi'*F{i-1}/(lambda+phi'*F{i-1}*phi))/lambda;
    
    theta(1:2,i)=theta(1:2,i-1)+F{i}*[x(i-1);u(i-1)]*(x(i)-xhat(i));
    
    u(i)=randn(1)/10-0.8*u(i-1);

end

figure(1)
subplot(211)
plot(theta');
hold on
plot([A;B]','--');
hold off
legend('Estimated A','Estimated B','True A','True B')
title('Parameter Adaptation Result')
subplot(212)
plot(x-xhat);

%% Mean squared estimation error
varx(1)=0;vartheta{1}=eye(2);

thetatilde(1:2,1)=dC(:,1);

for i=2:n
    
    varx(i)=[x(i-1);u(i-1)]'*vartheta{i-1}*[x(i-1);u(i-1)]+Bd'*Bd;
    
    vartheta{i}=vartheta{i-1}+F{i}*[x(i-1);u(i-1)]*varx(i)*[x(i-1);u(i-1)]'*F{i}-vartheta{i-1}*[x(i-1);u(i-1)]*[x(i-1);u(i-1)]'*F{i}-(vartheta{i-1}*[x(i-1);u(i-1)]*[x(i-1);u(i-1)]'*F{i})';
    
    dC(:,i-1)=DC;
    thetatilde(1:2,i)=(eye(2)-[x(i-1);u(i-1)]*[x(i-1);u(i-1)]'*F{i-1})*thetatilde(1:2,i-1)+dC(:,i-1);
    
    vartheta{i}=vartheta{i}+2.*thetatilde(1:2,i)*dC(:,i-1)'-dC(:,i-1)*dC(:,i-1)';
    
end

hold on
%plot([sqrt(abs(varx));-sqrt(abs(varx))]','r'); % Show sigma bound
plot([3*sqrt(abs(varx));-3*sqrt(abs(varx))]','g'); % Show 3sigma bound
plot([std(x-xhat)*ones(1,n);-std(x-xhat)*ones(1,n)]','--k') % Show standard deviation
axis([0 1000 -0.15 0.15])
hold off
title('Prediction Error Profile')
legend('Prediction error','3\sigma(k|k-1)','-3\sigma(k|k-1)','Std of the Data')
xlabel('Time Step k')

    
    
    