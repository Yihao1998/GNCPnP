function [R,t] = GNCPnP(XX,xx)

P=XX;
Qp=xx;

minerror2 = 0.02; %for the synthetic experiments (f = 800)

% [M, Cw, Alph] = PrepareData(P,Qp);

prev_sv = Inf;
weights1=ones(1,length(P));

for kk=1:100
        
    [M, ~, ~] = PrepareData(P,Qp);
    m=length(P);
    
    for ii=1:length(weights1)
        N(2*ii-1,:) = M(2*ii-1,:)*weights1(ii);
        N(2*ii,:) = M(2*ii,:)*weights1(ii);
    end
    
    [~,~,v] = svd(N'*N);

    error21    = M(1:2:end,:) * v(:,end);
    error22    = M(2:2:end,:) * v(:,end);
    error2     = sqrt(error21.^2 + error22.^2);

    [sv, ~] = sort(error2);        
    med = sv(floor(m/2)); 
    
    minerror=max(med,minerror2);
    
    if kk == 1
          mu = 0.05;
    end

    ub = (mu+1)/mu * minerror;
    lb = (mu)/(mu+1) * minerror; 

    weights1=[];
    
    for k = 1:length(error2)
        if error2(k) - ub >= 0
            weights1(k) = 0;

        elseif error2(k) - lb <= 0
            weights1(k) = 1;
        else                
            weights1(k) = sqrt( minerror*mu*(mu+1)/error2(k)) - mu;
        end
    end
    
    mu=mu*3;
    
    
    if (med >= prev_sv)
        break;
    else
        prev_sv = med;
    end
end


[Mfinal, ~, ~] = PrepareData(P,Qp);

error21f    = Mfinal(1:2:end,:) * v(:,end);
error22f    = Mfinal(2:2:end,:) * v(:,end);
error2f     = sqrt(error21f.^2 + error22f.^2);

for k = 1:length(error2f)
    if error2f(k) - ub >= 0
        weights1(k) = 0;

    elseif error2f(k) - lb <= 0
        weights1(k) = 1;
    else                
        weights1(k) = sqrt( minerror*mu*(mu+1)/error2f(k)) - mu;
    end
end


fff=1;
for ff=1:length(weights1)
    
   if weights1(ff)==1
       fP(:,fff)=P(:,ff);
       fQ(:,fff)=Qp(:,ff);
       fff=fff+1;
   end
end

[~, cw, ~] = PrepareData(fP,fQ);
K = v(:,end-4+1:end);   

[Ri,~, ~] = KernelPnP(cw, K, 4, 1);

[R,t, ~,~,~]=objpose2(fP,fQ,Ri);














