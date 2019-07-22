T = 1;
r = 0.1;
sigma = 0.5;
K = 5;
S0 = [4.5;5;5.5];
N_step = 1000;
n_path = 100000;
dti = T/N_step;

% 2.1 Exact Solution
Sti = zeros(N_step+1,n_path,length(S0));
disPo = zeros(n_path,length(S0));
for l = 1:length(S0)
for j = 1: n_path
dWti = normrnd(0,sqrt(dti),N_step,1);
Sti(1,j,l) = S0(l,1);
for i = 1:N_step
    Sti(i+1,j,l) = Sti(i,j,l)*exp((r-0.5*sigma^2)*dti+sigma*dWti(i,1));
end
disPo(j,l) = exp(-r*T)*subplus(mean(Sti(1:N_step,j,l))-K);
end
[mean_es(:,l),sigma_es(:,l),muci_es(:,l),sigmaci_es(:,l)] = normfit(disPo(:,l));
std_error_es(:,l) = sigma_es(:,l)/sqrt(n_path);
end

% 2.2 Approximate Dynamic
Sti = zeros(N_step+1,n_path,length(S0));
dSti = zeros(N_step,n_path,length(S0));
disPo_AD = zeros(n_path,length(S0));
for l = 1:length(S0)
for j = 1: n_path
dWti = normrnd(0,sqrt(dti),N_step,1);
Sti(1,j,l) = S0(l,1);
for i = 1:N_step
    dSti(i,j,l) = Sti(i,j,l)*(r*dti+sigma*dWti(i,1));
    Sti(i+1,j,l) = Sti(i,j,l) + dSti(i,j,l);
end
disPo_AD(j,l) = exp(-r*T)*subplus(mean(Sti(1:N_step,j,l))-K);
end
[mean_ad(:,l),sigma_ad(:,l),muci_ad(:,l),sigmaci_ad(:,l)] = normfit(disPo_AD(:,l));
std_error_ad(:,l) = sigma_ad(:,l)/sqrt(n_path);
end    


