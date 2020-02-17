//Cálculo da matriz de correlação
function Rx = comat(x,N)
    px = xcorr(x,N,'biased');
    for i = 1:N+1
        Rx(:,i) = px((N+2-i):(2*(N+1)-i));
    end
endfunction

function [y,e,w]=filterLMS(x,d,mi,N)
    xn=zeros(N+1,1);
    w=zeros(N+1,1);
    M=length(d);
    for i=1:M
        xn=[x(i);xn(1:N)];
        y(i)= w(:,i)'*xn;
        e(i)=d(i)-y(i);
        w = [w w(:,i)+2*mi*e(i)*xn];
    end 
endfunction

function [y,e,w]=filterNLMS(x,d,mi,gama,N)
    xn=zeros(N+1,1);
    w=zeros(N+1,1);
    M=length(d);
    for i=1:M
        xn=[x(i);xn(1:N)];
        y(i)= w(:,i)'*xn;
        e(i)=d(i)-y(i);
        w= [w w(:,i)+mi*(1/(gama+xn'*xn))*e(i)*xn];
    end 
endfunction    

function [y,e,w]=filterRLS(x,d,lambda,delta,N)
    xn=zeros(N+1,1);
    w = zeros(N+1,1);
    Sd = delta*eye(N+1,N+1);
    M=length(d);
    for i=1:M
        xn=[x(i);xn(1:N)];
        eta=d(i)-xn'*w(:,i);
        psi=Sd*xn;
        Sd=((1/lambda)*(Sd-(psi*psi'/(lambda+xn'*psi))));
        w(:,i+1)= w(:,i)+eta*Sd*xn;
        y(i)=w(:,i+1)'*xn;
        e(i)=d(i)-y(i);        
    end
endfunction
