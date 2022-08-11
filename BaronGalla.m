%% please run the code in sections
%% Question 1



N=8; %% N represents the total number of sepcies
clear zeros
a1 = zeros(N,N); % this implements a NxN matrix full of zeroes
zerosx = zeros(N/2,N/2); % this implements an N/2xN/2 matrix full of zeroes
ve = ones(1,N/2);  % This implements a vector of N/2 elements full of ones
vn = 5*ones(1,N/2); % This implements a vector of N/2 elements that are 5s
Du = diag(ve); % This turns the ve variable into a N/2xN/2 matrix where the diagonals are 1s
Dv = diag(vn); % This turns the vn variable into a N/2xN/2 matrix where the the diagonals are 5s
clear zeros
du=diag(ve); % turns the ve vector of 1s into a diagonal vector with 1s as the diagonal numbers
dv=diag(ve);% turns the ve vector of 1s into a diagonal vector with 1s as the diagonal numbers
q=1; % variable q set as 1
sd=0.1; % Standard deviation of the random variables is set as 0.1
D = [Du zerosx; zerosx Dv]; % This is used to combine the submatrices into 1 matrix
d1 = [du zerosx; zerosx dv]; %This is used to combine the submatrices into 1 matrix, the matrix is called
                             %d1 as to not cause confusion with the global
                             %variable d
rng(4)
A = normrnd(0,sd,size(a1)); % this generates a random numbers
                           % from a guassian distribution with standard
                           % deviation 0.1 and mean 0. Using the size
                           % funciton allows for the matrix to be outputted
                           % based on the size of the original matrix.

                       
Mq = -q^(2)*D-d1+A;  %this is the function to find the problem matrix Mq                         

%% question 2

N=100; %%N represents the total number of species
clear zeros
a1 = zeros(N,N); % this implements a NxN matrix full of zeroes
zerosx = zeros(N/2,N/2); % this implements an N/2xN/2 matrix full of zeroes
ve = ones(1,N/2);  % This implements a vector of N/2 elements full of ones
vn = 5*ones(1,N/2); % This implements a vector of N/2 elements that are 5s
Du = diag(ve); % This turns the ve variable into a N/2xN/2 matrix where the diagonals are 1s
Dv = diag(vn); % This turns the vn variable into a N/2xN/2 matrix where the the diagonals are 5s
clear zeros
du=diag(ve); % turns the ve vector of 1s into a diagonal vector with 1s as the diagonal numbers
dv=diag(ve);% turns the ve vector of 1s into a diagonal vector with 1s as the diagonal numbers
q=0; % variable q set as 0
sd=0.1; %standard deviation of the random variables is set as 0.1
D = [Du zerosx; zerosx Dv]; % This is used to combine the submatrices into 1 matrix
d1 = [du zerosx; zerosx dv]; %This is used to combine the submatrices into 1 matrix, the matrix is called
                             %d1 as to not cause confusion with the global
                             %variable d
rng(4)
A = normrnd(0,sd,size(a1)); % this generates a random numbers
                           % from a guassian distribution with standard
                           % deviation 0.1 and mean 0

                       
Mq = -q^(2)*D-d1+A; %this is the function to find the problem matrix Mq

w = eig(Mq); % This finds the eigenvalues of the Mq matrix
wr= real(w); % We take only the real values of the eigenvalues and place it in a vector
wi = imag(w); % We take only the imaginary of the eigenvalues and place them in a vector

q2=0.2; %this creates a new q variable for the next matrix example
Mq2 = -q2^(2)*D-d1+A; %this is the function to find the problem matrix Mq

w2 = eig(Mq2); % This finds the eigenvalues of the Mq matrix
wr2= real(w2); % We take only the real values of the eigenvalues and place it in a vector
wi2 = imag(w2); % We take only the imaginary of the eigenvalues and place them in a vector

q3=0.5; %this creates a new q variable for the next matrix example
Mq3 = -q3^(2)*D-d1+A; %this is the function to find the problem matrix Mq

w3 = eig(Mq3); % This finds the eigenvalues of the Mq matrix
wr3= real(w3); % We take only the real values of the eigenvalues and place it in a vector
wi3 = imag(w3); % We take only the imaginary of the eigenvalues and place them in a vector



tiledlayout(2,2) %Sets a 2x2 layout of plots
nexttile
plot(wr,wi,'.'); % This plots the real parts of the eigenvalues against the imaginary parts of the eigenvalues.
title("Mq eigenvalues real plotted against imaginary where q=0")
xlabel("Re[\omega]","FontSize",10)
ylabel("IM(\omega)","FontSize",10,"Rotation",0)
hold on
xline(0, '--black');
hold off
nexttile
plot(wr2,wi2,'.'); % This plots the real parts of the eigenvalues against the imaginary parts of the eigenvalues.
title("Mq eigenvalues real plotted against imaginary where q=0.2")
xlabel("Re[\omega]","FontSize",10)
ylabel("IM(\omega)","FontSize",10,"Rotation",0)
hold on
xline(0, '--black');
hold off
nexttile
plot(wr3,wi3,'.'); % This plots the real parts of the eigenvalues against the imaginary parts of the eigenvalues.
title("Mq eigenvalues real plotted against imaginary where q=0.5")
xlabel("Re[\omega]","FontSize",10)
ylabel("IM(\omega)","FontSize",10,"Rotation",0)
hold on
xline(0, '--black');
hold off

%% question 3

N=100; %%N represents the total number of species
clear zeros
a1 = zeros(N,N); % this implements a NxN matrix full of zeroes
zerosx = zeros(N/2,N/2); % this implements an N/2xN/2 matrix full of zeroes
ve = ones(1,N/2);  % This implements a vector of N/2 elements full of ones
vn = 5*ones(1,N/2); % This implements a vector of N/2 elements that are 5s
Du = diag(ve); % This turns the ve variable into a N/2xN/2 matrix where the diagonals are 1s
Dv = diag(vn); % This turns the vn variable into a N/2xN/2 matrix where the the diagonals are 5s
clear zeros
du=diag(ve); % turns the ve vector of 1s into a diagonal vector with 1s as the diagonal numbers
dv=diag(ve);% turns the ve vector of 1s into a diagonal vector with 1s as the diagonal numbers
sd=0.11; % sets the standard deviation of the random deviations in this case to be 0.11
D = [Du zerosx; zerosx Dv]; % This is used to combine the submatrices into 1 matrix
d1 = [du zerosx; zerosx dv]; %This is used to combine the submatrices into 1 matrix, the matrix is called
                             %d1 as to not cause confusion with the global
                             %variable d
rng(4)
A = normrnd(0,sd,size(a1)); % this generates random numbers
                           % from a guassian distribution with standard
                           % deviation 0.11 and mean 0
wmax= zeros(11,1); %sets a zero vector of size 11

for q2=0:10
  Mq = -(q2/10)^(2)*D-d1+A;

  w = eig(Mq); %finds all the eigenvalues
  wr= real(w); %takes only the real parts of the eigenvalues
  wi = imag(w); %takes only the imaginary parts of the eigenvalues

  wmax(q2+1)=max(wr); %takes the highest real eigenvalue for each q value
end

sd=0.19; 
A = normrnd(0,sd,size(a1)); %this generates random numbers from a gaussian distribution
                            %with standard deviation 0.19 and mean 0
wmax2= zeros(11,1);

for q2=0:10
  Mq = -(q2/10)^(2)*D-d1+A;

  w = eig(Mq); %finds all the eigenvalues
  wr= real(w); %takes only the real parts of the eigenvalues
  wi = imag(w); %takes only the imaginary parts of the eigenvalues

  wmax2(q2+1)=max(wr); %takes the highest real eigenvalue for each q value
end

sd=0.24;
A = normrnd(0,sd,size(a1)); %this generates random numbers from a gaussian distribution
wmax3= zeros(11,1);         %with standard deviation 0.24 and mean 0

for q2=0:10
  Mq = -(q2/10)^(2)*D-d1+A;

  w = eig(Mq); %finds all the eigenvalues
  wr= real(w); %takes only the real parts of the eigenvalues
  wi = imag(w); %takes only the imaginary parts of the eigenvalues

  wmax3(q2+1)=max(wr); %takes the highest real eigenvalue for each q value
end

q = 0:0.1:1;
line = zeros(11,1);
tiledlayout(2,2)
nexttile
plot(q,wmax)
title("q plotted against Mq matrix maximum real eigenvalues, standard deviation 0.11")
xlabel("q","FontSize",10)
ylabel("Re[\omega_{max}]","FontSize",10,"Rotation",0)
hold on
plot(q,line,'--','Color','black')
hold off
nexttile
plot(q,wmax2)
title("q plotted against Mq matrix maximum real eigenvalues, standard deviation 0.19")
xlabel("q","FontSize",10)
ylabel("Re[\omega_{max}]","FontSize",10,"Rotation",0)
hold on
plot(q,line,'--','Color','black')
hold off
nexttile
plot(q,wmax3)
title("q plotted against Mq matrix maximum real eigenvalues, standard deviation 0.24")
xlabel("q","FontSize",10)
ylabel("Re[\omega_{max}]","FontSize",10,"Rotation",0)
hold on
plot(q,line,'--','Color','black')
hold off
%% question 4

global a b c d e
%set parameters a,b,c,d,e as global so they can be used across functions
a=1;
b=1;
c=1;
d=1;
e=0.1;
%load the derivative functions from the loadderivs file
Load = loadderivs


%this is used to find the equilibria
temp=fsolve(@Load.loadderivs1,1+rand(1,2));
for i = 1:100
    temp=[temp;fsolve(@Load.loadderivs1,1+rand(1,2))];
end
eq=unique(round(temp,8),'rows');
tiledlayout(2,3)

%this is used to show 5 graphs showing the system reaching the stable fixed point
%over time
for i= 1:5
nexttile
y0=rand(1,2); %sets a random starting point to show they always reach the same fixed point
tspan = [0 30];
[t,y]=ode45(@Load.loadderivs2,tspan,y0);
plot(t,y)

end

%% question 5

global a b c d e
%set parameters a,b,c,d,e as global so they can be used across functions
a=1;
b=1;
c=1;
d=1;
e=0.1;
%Load the derivative functions from the loadderivs file
Load = loadderivs;
tiledlayout(3,3)
% adds a delay to the system to see if an instability is caused
for delay=0.5:0.1:1.3
nexttile
y0=rand(1,2);
tspan = [0 100];
sol=dde23(@Load.loadderivs3,delay,y0,tspan);
plot(sol.x,sol.y)
title("System points changing overtime with delay " + delay)
end