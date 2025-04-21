function [res_x,f_sample] = realMultichannel_GAMP(yy,f_sample,Fs,L)

N=length(f_sample);
for tt=1:L
    [T(tt),~]=size(yy(tt).cluster);
    Q_T=real_trans(T(tt));
    yy(tt).cluster=Q_T*yy(tt).cluster;
    yy(tt).cluster=[real(yy(tt).cluster), imag(yy(tt).cluster)];
    [~,M(tt)]=size(yy(tt).cluster);
    norm_y(tt)=norm(yy(tt).cluster,'fro')/sqrt(T(tt)*M(tt));
    yy(tt).cluster=yy(tt).cluster/norm_y(tt);
    
    T_all=[T(tt)-1:-2:1-T(tt)]';
    A_temp= exp(-1i*pi*T_all/Fs*f_sample)/sqrt(N);
    B_temp=(-1i*pi*T_all/Fs).*A_temp;
    A=Q_T*A_temp;
    B=Q_T*B_temp;
    AB(tt)=struct('basis',A,'deriv',B);
    S=  zeros(T(tt),M(tt));
    S_hat(tt)=struct('meanS',S);
    V_x(:,:,tt)   =  ones(N,M(tt));
    X_hat(:,:,tt)=zeros(N,M(tt));
    beta0(tt)=1;
end
reslu=f_sample(2)-f_sample(1);

%% GAMP initialization
maxiter=100;
rho=0.8;
hold_max=10^10;
hold_min=10^(-10);
delta=ones(N,1);
converged = false;
truncated_active=false;
iter = 0;
tol=1e-5;
etc=20;

while ~converged && iter<maxiter
    %% update X_hat and V_x
    for tt=1:L                   
        V_p=ones(T(tt),1)*(sum(V_x(:,:,tt),1)/N);                            % (5a)  from next sub-block
        V_p=min(max(V_p,hold_min),hold_max);
        P_hat=(AB(tt).basis*X_hat(:,:,tt)) - (V_p.*S_hat(tt).meanS);        % (5b)  from next sub-block
              
        tao_z(tt).varZ=(V_p.*beta0(tt))./(V_p+beta0(tt));
        z_hat(tt).meanZ=tao_z(tt).varZ.*((yy(tt).cluster)./beta0(tt)+P_hat./V_p);
        V_s=(1-tao_z(tt).varZ./V_p)./V_p;
        V_s=min(max(V_s,hold_min),hold_max);
        S_hat(tt).meanS=(z_hat(tt).meanZ-P_hat)./V_p;
              
        % variable update part
        V_r= N./( sum(V_s) );
        V_r=min(max(V_r,hold_min),hold_max);
        R_hat= X_hat(:,:,tt) + (V_r.*( AB(tt).basis'* S_hat(tt).meanS));    % (7b)
        % step 4
        X_hat_old=X_hat(:,:,tt);
        V_x_old=V_x(:,:,tt);
        V_x(:,:,tt)=1./(  1./delta  + 1./V_r   );           % (8a)
        V_x(:,:,tt)=min(max(V_x(:,:,tt),hold_min),hold_max);
        X_hat(:,:,tt)= V_x(:,:,tt).*( 0 +  R_hat./V_r   );                    % (8b) remove  0.*delta
        X_hat(:,:,tt)= (1-rho)*X_hat(:,:,tt) + (rho)*X_hat_old;
        V_x(:,:,tt)= (1-rho)*V_x(:,:,tt) + (rho)*V_x_old;
    end
    
    
   %% update beta0 and delta
    mu=X_hat;
    Sigma=V_x;
    Exx =  sum(sum(mu.*conj(mu),2),3) + sum(sum(Sigma,2),3);

    beta_old=beta0;  
    for ll=1:L
        resid=(yy(ll).cluster-z_hat(ll).meanZ);
        term=sum(sum(tao_z(ll).varZ,2));
        beta0(ll)=( norm(resid(:), 'fro')^2+  real(term)   )/(T(ll)*2);
        beta0(ll)=  0.02*beta0(ll) + 0.98*beta_old(ll);
    end

 
    delta_last=delta;
    sum_temp1=sum(Exx,2);
    delta=sum_temp1;
    

       %% off-grid
      if ~truncated_active  
        Pm=sum(sum( mu.*conj(mu), 2),3);
        [~,sort_ind]=sort(Pm, 'descend');
        idx=sort_ind(1:etc);
        P=[];v=[];
        for t=1:L
            BB=AB(t).deriv(:,idx);
            BHB =BB' * BB;
%             P2=diag(Sigma(idx,1,t)+Sigma(idx,2,t)); 
            P2=diag(sum(Sigma(idx,:,t),2));
            P(:,:,t) = real( conj(BHB) .* ((mu(idx,:,t) * mu(idx,:,t)') +   P2   )  );
            Tes=diag(sum(Sigma(:,:,t),2));
            v2= real(diag(BB' * AB(t).basis * Tes(:,idx)   ));
            v(:,t) = sum( real(conj(mu(idx,:,t)) .* (BB' * (yy(t).cluster - AB(t).basis * mu(:,:,t)))),2) -   v2;
        end
        Sum_P=sum(P,3);Sum_v=sum(v,2);
        temp_grid=Sum_v./diag(Sum_P);
        temp_grid=temp_grid';
        
        theld=reslu/20*0.95^(iter);
        ind_small=find(abs(temp_grid)<theld);
        temp_grid(ind_small)=sign(temp_grid(ind_small))*theld;
        ind_unchang=find (abs(temp_grid)>reslu);
        temp_grid(ind_unchang)=sign(temp_grid(ind_unchang)) * reslu/20;
        f_sample(idx)=f_sample(idx) + temp_grid;
        for tt=1:L
            T_all=[T(tt)-1:-2:1-T(tt)]';
            AB(tt).basis(:,idx)=real_trans(T(tt)) * (exp(-1i*pi*T_all/Fs*f_sample(idx))/sqrt(N));
            AB(tt).deriv(:,idx)=real_trans(T(tt)) *((-1i*pi*T_all/Fs).*AB(tt).basis(:,idx));
        end
        
     else
         
         PP=[];vv=[];
         for t=1:L
             BHB = AB(t).deriv' * AB(t).deriv;
             P2=  diag(sum(Sigma(:,:,t),2));
             PP(:,:,t) = real( conj(BHB) .* ((mu(:,:,t) * mu(:,:,t)') +   P2   )  );
             v2= real(diag(AB(t).deriv' * AB(t).basis * diag(sum(Sigma(:,:,t),2)) ));
             vv(:,t) = sum( real(conj(mu(:,:,t)) .* (AB(t).deriv' * (yy(t).cluster - AB(t).basis * mu(:,:,t)))),2) -   v2;
         end
         sum_PP=sum(PP,3);sum_vv=sum(vv,2);
         vect1=[1:N]';
         sum_PP=vect1'*sum_PP*vect1;
         temp_grid=(vect1'*sum_vv)/sum_PP;
         theld=reslu/20*0.95^(iter);
         ind_small=find(abs(temp_grid)<theld);
         temp_grid(ind_small)=sign(temp_grid(ind_small))*theld;
         ind_unchang=find (abs(temp_grid)>reslu);
         temp_grid(ind_unchang)=sign(temp_grid(ind_unchang)) * reslu/20;
         f_sample=f_sample + temp_grid*vect1';
         for t=1:L
             T_all=[T(t)-1:-2:1-T(t)]';
             AB(t).basis=real_trans(T(t)) *( exp(-1i*pi*T_all/Fs*f_sample)/sqrt(N) );
             AB(t).deriv=real_trans(T(t)) *( (-1i*pi*T_all/Fs).*AB(t).basis );
         end
     end
     
     if iter==50
         Pm2=delta;
         X_hat=[];V_x=[];
         fn=search_Pm(Pm2,f_sample);
         fn_all= fn*[1:1:10];
         truncated_active=true;
         f_sample=fn_all;
         N=length(f_sample);
         for t=1:L
             T_all=[T(t)-1:-2:1-T(t)]';
             AB(t).basis=real_trans(T(t)) *( exp(-1i*pi*T_all/Fs*f_sample)/sqrt(N) );
             AB(t).deriv=real_trans(T(t)) *( (-1i*pi*T_all/Fs).*AB(t).basis );
             X_hat(:,:,t)=zeros(N,M(t));
             V_x(:,:,t)=ones(N,M(t));
         end        
         delta=ones(N,1)*1;
         etc=min(etc,N);
         delta_last=100;
     end
        
        
        if iter<50
            theld2=50;
            ind_remove= find (delta<(max(delta)/theld2));
            for tt=1:L
                AB(tt).basis(:,ind_remove)=[];
                AB(tt).deriv(:,ind_remove)=[];
            end
            delta(ind_remove)=[];
            delta_last(ind_remove)=[];
            f_sample(ind_remove)=[];
            X_hat(ind_remove,:,:)=[];
            V_x(ind_remove,:,:)=[];
            N=length(f_sample);
            etc=min(etc,N);
        end
      
    %% stopping criteria
    erro=  max(max(abs(delta - delta_last)));  
    if erro < tol || iter >= maxiter
        converged = true;
    end
    iter = iter + 1;

end

res_x=sqrt(abs(delta))*norm_y;  
