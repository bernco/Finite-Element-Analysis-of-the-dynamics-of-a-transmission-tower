clear all
close all
%% PARAMETERS TO CHOOSE:

Nelement=9; %Nb of elements for each beams

NB_EIGENVALUES=15; %Nb of eigenvalue taking into consideration
%% 
%% I-Modelling and modal analysis
%% 
%% 1/PARAMETERS TOWER


rho=8*10^3;        % [kg/m3] material density
E=210e9;        % [Pa] Young's modulus
nu=0.3;%Poisson coefficient
Jx= 2.25*(0.07/2)^4;%Moment of torsion
Iy=0.07*0.07^3/12; %Moment of inertia
Iz=0.07*0.07^3/12;%Moment of inertia
A=0.07*0.07; % [m2] Cross section area
rt=sqrt(Jx/A);%Gyroscopic radius
G=E/(2*(1+nu)); %IN Pa shear modulus of steel

%% 2/PARAMETERS FE
% List of coordinates of all the beams: [NODES1 0 NODES2]
COORD_BEAMS=10^-3*[0 0 0 0 750 750 5000;10000 0 0 0 9250 750 5000;10000 10000 0 0 9250 9250 5000;0 10000 0 0 750 9250 5000;750 750 5000 0 9250 750 5000;9250 750 5000 0 9250 9250 5000;9250 9250 5000 0 750 9250 5000;750 9250 5000 0 750 750 5000;0 0 0 0 9250 750 5000;10000 0 0 0 9250 9250 5000;10000 10000 0 0 750 9250 5000;0 10000 0 0 750 750 5000;1500 1500 10000 0 8500 1500 10000;8500 1500 10000 0 8500 8500 10000;8500 8500 10000 0 1500 8500 10000;1500 8500 10000 0 1500 1500 10000;750 750 5000 0 1500 1500 10000;9250 750 5000 0 8500 1500 10000;9250 9250 5000 0 8500 8500 10000;750 9250 5000 0 1500 8500 10000;9250 750 5000 0 1500 1500 10000;9250 9250 5000 0 8500 1500 10000;750 9250 5000 0 8500 8500 10000;750 750 5000 0 1500 8500 10000;2250 2250 15000 0 7750 2250 15000;7750 2250 15000 0 7750 7750 15000;7750 7750 15000 0 2250 7750 15000;2250 7750 15000 0 2250 2250 15000;1500 1500 10000 0 2250 2250 15000;8500 1500 10000 0 7750 2250 15000;8500 8500 10000 0 7750 7750 15000;1500 8500 10000 0 2250 7750 15000;1500 1500 10000 0 7750 2250 15000;8500 1500 10000 0 7750 7750 15000;8500 8500 10000 0 2250 7750 15000;1500 8500 10000 0 2250 2250 15000;3000 3000 20000 0 7000 3000 20000;7000 3000 20000 0 7000 7000 20000;7000 7000 20000 0 3000 7000 20000;3000 7000 20000 0 3000 3000 20000;2250 2250 15000 0 3000 3000 20000;7750 2250 15000 0 7000 3000 20000;7750 7750 15000 0 7000 7000 20000;2250 7750 15000 0 3000 7000 20000;7750 2250 15000 0 3000 3000 20000;7750 7750 15000 0 7000 3000 20000;2250 7750 15000 0 7000 7000 20000;2250 2250 15000 0 3000 7000 20000;3000 3000 20000 0 3000 3000 23000;7000 3000 20000 0 7000 3000 23000;7000 7000 20000 0 7000 7000 23000;3000 7000 20000 0 3000 7000 23000;3000 3000 23000 0 7000 3000 23000;7000 3000 23000 0 7000 7000 23000;7000 7000 23000 0 3000 7000 23000;3000 7000 23000 0 3000 3000 23000;3000 3000 20000 0 5000 -3000 20000;7000 3000 20000 0 5000 -3000 20000;3000 3000 23000 0 5000 -3000 20000;7000 3000 23000 0 5000 -3000 20000;7000 7000 20000 0 5000 13000 20000;3000 7000 20000 0 5000 13000 20000;7000 7000 23000 0 5000 13000 20000;3000 7000 23000 0 5000 13000 20000];

%List of coordinates of all NODES: [NODES]

NODES=[COORD_BEAMS(:,1:3);COORD_BEAMS(:,5:7)];
NODES=unique(NODES,'rows'); %To keep only one time one node (because one node can appear several times because it can belongs to 1,2 or 3 beams)
nb_nodes_tower=length(NODES);
nb_beams_tower=length(COORD_BEAMS);

%% 3/Discretisation of the tower

% * Parameters of the discretization: We have to test differents values of Nelement 
% in order to discuss on this parameter
% 
% What I think:
% 
% -If Nelement is to small (2,3,4) , the discrezation is too large and so 
% the results will not be correct (far from the reality)
% 
% -If Nelement too big: it will require a long time to compute the matix 
% and it is too much --> problem of convergence
% 
% 
% 
% To check in running NX for differents values of discretization.
% 
% * Discritization through x-axis in the local axis
% 
% |-----------------------|-----------------------|-----------------------|-----------------------|
% 
% N1                     N1,1                    N1,2                   N1,3                      N2             
%                     
% 
% N1,1, N1,2, N1,3 : secoundary nodes

for i=1:nb_beams_tower
    R=ThVib_Group8_Rotation_function(COORD_BEAMS(i,1:3),COORD_BEAMS(i,5:7));%Rotation matrix to go from global to local
    COORD_NODE1=R*COORD_BEAMS(i,1:3)';%Coordinate of node 1 in local axis
    COORD_NODE2=R*COORD_BEAMS(i,5:7)'; %Coordinate of node 1 in local axis
    X_DISCRETISATION=linspace(COORD_NODE1(1), COORD_NODE2(1), Nelement); %Discretization throug x-axis
    
    counter_secondary_nodes=2;%We start in taking the second node because the first one is a N1 a primary nodes
    counter_colomn=1;
    while counter_secondary_nodes<Nelement
        NODES_FEM_BIS(i,counter_colomn:counter_colomn+2)=inv(R)*[X_DISCRETISATION(counter_secondary_nodes),COORD_NODE1(2),COORD_NODE1(3)]';
        %[X_DISCRETISATION(counter_colomn),COORD_NODE1(2),COORD_NODE1(3)]=[x-coordinate given by linspace, y-coordinate: same than the primary nodes, z-coordinate: same than the primary nodes]
        %y and z are constant in the local axis.
        %*inv(R) --> to go back into the global axis
        counter_colomn=counter_colomn+3;%Increase by three to fill the next three coloms
        counter_secondary_nodes=counter_secondary_nodes+1;%Increase by one to take the next x-coordinate of the next secondary node
    end 

end
%%
% Above, in NODES_FEM_BIS, all the secondary nodes are refered: in one row, 
% all the nodes belonging to a same beam are refered.
%  
% - Get the list of coordinates of each nodes after applying FEM
% 
% Below, the purpose of the code is to reference all the node in the list 
% called NODES_FEM: one row=one node

counter_beam=1;%we start in considering the first beam
counter_row=1;%We start to fill the first row
while counter_beam<=length(NODES_FEM_BIS) %We want to go through all the beam
    for k=1:3:(Nelement-2)*3 %k is the index of the colomn
        NODES_FEM(counter_row,1:3)=NODES_FEM_BIS(counter_beam,k:k+2);%We reference each secondary nodes in one big matrix
        counter_row=counter_row+1;%One you fill the row with one node, we increase by one to fill the next row
    end
    counter_beam=counter_beam+1;%Once we reference all the secondary nodes of the beam we go on the next one

end

NODES_FEM=[NODES;NODES_FEM];%We add to the secondary nodes, the primary nodes
%NODES_FEM=[PRIMARY NODES;SECONDARY NODES]

%%
% * Get the list of coordinate for each elements after applying FEM


NODES_FEM_TER=[COORD_BEAMS(:,1:3),NODES_FEM_BIS,COORD_BEAMS(:,5:7)];
%NODES_FEM_TER=[Node 1, secondary nodes, Node 2]
counter_row=1;%We start to fill the first colomn
counter_beam=1;% We start taking the first beam
while counter_beam<=length(NODES_FEM_TER) %We want to go through all the beams (64)

    for k=1:3:Nelement*3-3 %k represents the index of the colomn
        %We know that two consecutives nodes in the matrix are linked and so they form one element
        COORD_FEM(counter_row,1:3)=NODES_FEM_TER(counter_beam,k:k+2);%Nodes 1
        
        COORD_FEM(counter_row,4:6)=NODES_FEM_TER(counter_beam,k+3:k+5);%Nodes 2
        
         counter_row=counter_row+1;%One you have [Nodes 1 Nodes2], we increase the counter_row to fill the next row with the next element
    end
    counter_beam=counter_beam+1;%One you went through all the elements of one beam, you can go through the next beam
end

        
%% 
%  Plot after applying FEM

figure(1)
%Plot the segment betwenn two nodes
for i=1:length(COORD_BEAMS)
    hold on
 
    plot3([COORD_BEAMS(i,1),COORD_BEAMS(i,5)],[COORD_BEAMS(i,2),COORD_BEAMS(i,6)],[COORD_BEAMS(i,3),COORD_BEAMS(i,7)],'r')
    
end

%Plot the nodes
for i=1:length(NODES_FEM)
    hold on
    
    scatter3(NODES_FEM(i,1),NODES_FEM(i,2),NODES_FEM(i,3),'o')
    
    end
     

%% Localization matrix
% * Degrees of freedom

NB_NODES_FEM=length(NODES_FEM);
NB_DOF=6*NB_NODES_FEM;%6 DOF for each nodes

% List of degrees of freedom (DoFs related to each node)
dofList = [(1:6:NB_DOF-5)' (2:6:NB_DOF-4)' (3:6:NB_DOF-3)' (4:6:NB_DOF-2)' (5:6:NB_DOF-1)' (6:6:NB_DOF)'];

%  Localisation matrix: [degrees of freedom of node 1, degrees of freedom of node 2]= 12 colomns=6*2nodes
% 
% For each element, we are going to try to find the degrees of freedom of 
% both nodes and then we can put together those DOF using locel(i,:)= [dofs1,dofs2].
% 
% To do that, for each element, we are going to:
% 
%     -Take the coordinate of node 1 and node2
% 
%     -Then we went to determine the degrees of freedom of those nodes. So 
% we are going to try to indentify those nodes in the list of node: NODES_FEM.
% 
%     -One you have the row index of Node 1 and Node 2, you can have the 
% DOF of those nodes going in dofList(index1,:) and dofList(index2,:).

for i=1:length(COORD_FEM)
    NODE1=COORD_FEM(i,1:3);
    NODE2=COORD_FEM(i,4:6);
    index1=1;%We start at the first row of NODES_FEM
    while (NODES_FEM(index1,1)~=NODE1(1)) || (NODES_FEM(index1,2)~=NODE1(2)) || (NODES_FEM(index1,3)~=NODE1(3))
       %While you do not have identify your node in the list of nodes, you will increanse the index, to try the next values
       % | : or --> While x or y or z are not equal we increase the index to go further
        index1=index1+1;      
    
    end
    %Exactly the same but for the seconde node
    index2=1;
    while (NODES_FEM(index2,1)~=NODE2(1)) ||(NODES_FEM(index2,2)~=NODE2(2))||(NODES_FEM(index2,3)~=NODE2(3))
        index2=index2+1;      
    end
    INDEX_COORD(i,:)=[index1 index2];
    dofs1 = dofList(index1, :); %Degrees of freedom of the first node
    dofs2 = dofList(index2, :); %Degrees of freedom of the seconde node
    locel(i,:) = [dofs1,dofs2];
end

    
%% Stiffeness and mass matrix

 Ks=zeros(NB_DOF, NB_DOF);%Initial matrix
 Ms=zeros(NB_DOF, NB_DOF);%Initial matrix

for i=1:length(COORD_FEM) %Going through all the elements
    
    L=sqrt((COORD_FEM(i,1)-COORD_FEM(i,4))^2+(COORD_FEM(i,2)-COORD_FEM(i,5))^2+(COORD_FEM(i,3)-COORD_FEM(i,6))^2);
    l=L; %Length of your element
    
    %Elementary stiffeness matrix in LOCAL AXIS !!
    Kel =[E*A/l 0 0 0 0 0 -E*A/l 0 0 0 0 0;
    0 12*E*Iz/l^3 0 0 0 6*E*Iz/l^2 0 -12*E*Iz/l^3 0 0 0 6*E*Iz/l^2;
    0 0 12*E*Iy/l^3 0 -6*E*Iy/l^2 0 0 0 -12*E*Iy/l^3 0 -6*E*Iy/l^2 0;
    0 0 0 G*Jx/l 0 0 0 0 0 -G*Jx/l 0 0;
    0 0 -6*E*Iy/l^2 0 4*E*Iy/l 0 0 0 6*E*Iy/l^2 0 2*E*Iy/l 0;
    0 6*E*Iz/l^2 0 0 0 4*E*Iz/l 0 -6*E*Iz/l^2 0 0 0 2*E*Iz/l;
    -E*A/l 0 0 0 0 0 E*A/l 0 0 0 0 0;
    0 -12*E*Iz/l^3 0 0 0 -6*E*Iz/l^2 0 12*E*Iz/l^3 0 0 0 -6*E*Iz/l^2;
    0 0 -12*E*Iy/l^3 0 6*E*Iy/l^2 0 0 0 12*E*Iy/l^3 0 6*E*Iy/l^2 0;
    0 0 0 -G*Jx/l 0 0 0 0 0 G*Jx/l 0 0;
    0 0 -6*E*Iy/l^2 0 2*E*Iy/l 0 0 0 6*E*Iy/l^2 0 4*E*Iy/l 0;
    0 6*E*Iz/l^2 0 0 0 2*E*Iz/l 0 -6*E*Iz/l^2 0 0 0 4*E*Iz/l];
    

 %Elementary mass matrix in LOCAL AXIS !!
    Mel = rho*A*l*...
    [1/3 0 0 0 0 0 1/6 0 0 0 0 0;
    0 13/35 0 0 0 11*l/210 0 9/70 0 0 0 -13*l/420;
    0 0 13/35 0 -11*l/210 0 0 0 9/70 0 13*l/420 0;
    0 0 0 rt^2/3 0 0 0 0 0 rt^2/6 0 0;
    0 0 -11*l/210 0 l^2/105 0 0 0 -13*l/420 0 -l^2/140 0;
    0 11*l/210 0 0 0 l^2/105 0 13*l/420 0 0 0 -l^2/140;
    1/6 0 0 0 0 0 1/3 0 0 0 0 0;
    0 9/70 0 0 0 13*l/420 0 13/35 0 0 0 -11*l/210;
    0 0 9/70 0 -13*l/420 0 0 0 13/35 0 11*l/210 0;
    0 0 0 rt^2/6 0 0 0 0 0 rt^2/3 0 0;
    0 0 13*l/420 0 -l^2/140 0 0 0 11*l/210 0 l^2/105 0;
    0 -13*l/420 0 0 0 -l^2/140 0 -11*l/210 0 0 0 l^2/105];

    % ROtation matrix
    R=ThVib_Group8_Rotation_function(COORD_FEM(i,1:3),COORD_FEM(i,4:6));
    
    T=blkdiag(R,R,R,R);%Block diagonal matrix --> 12*12 matrix
    %Here, we need a 12*12 rotation matrix to do the matricial product with Kel or Mel
    %One block R is for 1 block (x,y,z) either for the translation or the rotation
%% Stiffness and mass matrices of the element expressed in structural axes
    KeS=T'*Kel*T;
    MeS=T'*Mel*T;
    
    %Here we use the locel matrix to insert the elementary matrix at the good place in the big matrix
    Ks(locel(i,:),locel(i,:)) = Ks(locel(i,:),locel(i,:)) + KeS;
    
    Ms(locel(i,:),locel(i,:)) = Ms(locel(i,:),locel(i,:)) + MeS;
 
end

%% Lumped masses

%ADD the lumped mass at both extremity of the tower: 13,14 are the index of the nodes at both extremities
Ms(dofList(13,1:3),dofList(13,1:3))=Ms(dofList(13,1:3),dofList(13,1:3))+300;
Ms(dofList(14,1:3),dofList(14,1:3))=Ms(dofList(14,1:3),dofList(14,1:3))+300;



%% Boundary conditions
% Below, the purpose of the code is to find the nodes wich are clamped in order 
% to delete those elements in the main matrix Ks and Ms

i=1;
index=1;%We are going through all the nodes to find the clamped nodes
NODES_CLAMPED=zeros(4);
while index<5
    if NODES_FEM(i,3)==0 %One node is clamped if the z-coordinate is null
        NODES_CLAMPED(index)=i;%We add the index of each clamped nodes
        index=index+1;
    end
    i=i+1;
end

%For each clamped nodes, we have to remove the rows and colomns related to the degrees of freedom o those nodes
for k=4:-1:1 %Step=-1 : beacause we start from the last clamped nodes and finish with the first one 
    %PB: if you start by removing the first row, the second row will become the first one so you will have a problem 
    for i=6:-1:1 %For each degrees of freedom
        x=NODES_CLAMPED(k);
        Ks(6*(x-1)+i,:)=[];
        Ks(:,6*(x-1)+i)=[];
        Ms(6*(x-1)+i,:)=[];
        Ms(:,6*(x-1)+i)=[];

    end
end
   
NDOF=length(Ks);%Number of degrees of freedom after removing the clamped nodes


%% Solve the problem
%%
[MODE, D]=eigs(Ks, Ms,NB_EIGENVALUES,'sm'); %Solve the problem
%eigs(Ks, Ms,8,'sm'): 8: ask for the 8 first eigenvalues and 'sm': to get a rank list
w = sqrt(diag(D));  % [rad/s]
f = w /2/pi;        % [Hz]

%% Graphic representation of the deflection
if NB_EIGENVALUES>3
    no_plot=3;
else
    no_plot=NB_EIGENVALUES;
end
for mod=1:no_plot
    i=1;
    for k=1:length(NODES_FEM)
        if NODES_FEM(k,3)~=0
           %To get the deflection, you just have to add to your coordinate, the values of the eigenvector 
           NODES_FEM_DEFL(k,1)=5*MODE(i,mod)+NODES_FEM(k,1);
           NODES_FEM_DEFL(k,2)=5*MODE(i+1,mod)+NODES_FEM(k,2);
           NODES_FEM_DEFL(k,3)=5*MODE(i+2,mod)+NODES_FEM(k,3);
           i=6+i;
        end
        if NODES_FEM(k,3)==0 %If it is a clamped node, there are no deflections so the coordinates are still the same
            
           NODES_FEM_DEFL(k,1)=NODES_FEM(k,1);
           NODES_FEM_DEFL(k,2)=NODES_FEM(k,2);
           NODES_FEM_DEFL(k,3)=NODES_FEM(k,3);
           
        end
    
    end

    for k=1:length(COORD_FEM)
    
        COORD_FEM_DEFL(k,1:3)=NODES_FEM_DEFL(INDEX_COORD(k,1),:);
        COORD_FEM_DEFL(k,4:6)=NODES_FEM_DEFL(INDEX_COORD(k,2),:);
    end
    

    figure(mod+1)
    for i=1:length(COORD_FEM_DEFL)

        hold on
        plot3([COORD_FEM_DEFL(i,1),COORD_FEM_DEFL(i,4)],[COORD_FEM_DEFL(i,2),COORD_FEM_DEFL(i,5)],[COORD_FEM_DEFL(i,3),COORD_FEM_DEFL(i,6)],'b','LineWidth',2)
    end



   title(['Mode No.',' ',int2str(mod),' ','for',' ',int2str(Nelement-1),' ', 'elements'])
end
    

   
%% Damping matrix
%%
DAMPING_COEFFICIENT=[0.01,0.01]; %Damping for the two first modes
W=2*pi*f; %Eight first pulsations (undamped problem)

%%

% * Calculate the coefficents a and b

a=2*DAMPING_COEFFICIENT(1)/(W(1)+W(2));
b=a*(W(1)*W(2));
% a=8.51*10^-4;
% b=0.106;
C=a*Ks+b*Ms;
%%  Other damping ratio

for i=3:length(W)
    DAMPING_COEFFICIENT(i)=0.5*(a*W(i)+b/W(i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% II-Modal superposition methods %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
LAMBDA=W';%Vector of the eigenvalues
W_DISSIPATION=LAMBDA.*sqrt(1-DAMPING_COEFFICIENT'.^2); %Vec tor with the eigenvalues for the damped system
%% Solve the pb
%  Time discretization

temporal_step=0.01; %Time step
T=30;%Final time
t=0:temporal_step:T; %t will be a list with the time going from 0 to T with a step of 0.01

%%  Load force
% 
% P will be a matrix: NDOF*length(t)
% 
% We are going to apply only two loads (one positive and one negative at 
% each extremity)
% 
% So in P you will lot of zeros and only two rows not null
% 
% Fz=force*sin(2*pi/period*t): harmonic force applied at the node 13
% 
% Fz=-force*sin(2*pi/period*t): harmonic force applied at the node 14

period=0.4;
force=0.4 * 10^3; %Force 0.4 kN

P=zeros(NDOF,length(t));

P(dofList(13-2,3),:)=force*sin(2*pi/period*t);%-2 because the index changed after having removed the DOF of the clamped nodes
P(dofList(14-2,3),:)=-force*sin(2*pi/period*t);
%% 


q0=zeros(NDOF,1); %Initialization off the displacement at 0
dq0=zeros(NDOF,1);%Initialization off the velocity at 0


for i=1:NB_EIGENVALUES
    A=0;%Initialization of A at 0
    B=0;%Initialization of B at 0

 
    mu=MODE(:,i)'*Ms*MODE(:,i);
    PHI(i,:)=MODE(:,i)'*P/mu;
    
    h(i,:)=temporal_step/W_DISSIPATION(i)*exp(-DAMPING_COEFFICIENT(i)*LAMBDA(i)*t).*sin(W_DISSIPATION(i)*t);
    
    n_1(i,:)=exp(-DAMPING_COEFFICIENT(i)*LAMBDA(i)*t).*(A*cos(W_DISSIPATION(i)*t)+B*sin(W_DISSIPATION(i)*t));

    n_2=conv(PHI(i,:),h(i,:));
    n(i,:)=n_1(i,:)+n_2(1:length(t));
end
    
%%  Modal displacement method
%%


q_displacement=zeros(NDOF,length(t));
q_displacement=(n'*MODE')';

%% Modal acceleration method
%% 

q_22=zeros(NDOF,length(t));
q_acceleration=zeros(NDOF,length(t));

for i=1:NB_EIGENVALUES
    mu=MODE(:,i)'*Ms*MODE(:,i);
    q_22=q_22+(PHI(i,:)'/LAMBDA(i)^2*MODE(:,i)')';%Same than below
    %q_22=q_22+MODE(i,:)*MODE(i,:)'*P/(LAMBDA(i)^2*mu);
end 
q_acceleration=q_displacement+inv(Ks)*P-q_22;

%% Direct integration
%%
% 

gamma=0.5;
beta=0.25;
d2q0=inv(Ms)*(P(:,1)-C*dq0-Ks*q0);%Initialization of q dot dot 0

d2q(1,:)=d2q0;
dq(1,:)=dq0;
q(1,:)=q0;

d2q_star(1,:)=d2q0;
dq_star(1,:)=dq0;
q_star(1,:)=q0;

S=Ms+gamma*temporal_step*C+beta*temporal_step^2*Ks;
S=inv(S);
for temp=1:length(t)-1
    
    %Prediction
    dq_star(temp+1,:)=dq(temp,:)+(1-gamma)*temporal_step*d2q(temp,:);
    q_star(temp+1,:)=q(temp,:)+temporal_step*dq(temp,:)+(0.5-beta)*temporal_step^2*d2q(temp,:);
    %COmputation of accelerations
    d2q(temp+1,:)=S*(P(:,temp+1)-C*dq_star(temp+1,:)'-Ks*q_star(temp+1,:)');
    %Correction
    dq(temp+1,:)=dq_star(temp+1,:)+temporal_step*gamma*d2q(temp+1,:);
    q(temp+1,:)=q_star(temp+1,:)+temporal_step^2*beta*d2q(temp+1,:);
end

  
%% FFT
%%


Y_disp=fft(q_displacement(dofList(13-2,1),:));
Y_acc = fft(q_acceleration(dofList(13-2,1),:));
Y_newmark=fft(q(:,dofList(13-2,1)));
Y_disp(1) = [];
Y_acc(1) = [];
Y_newmark(1) = [];

n = length(Y_disp);
power_disp = abs(Y_disp(1:floor(n/2))).^2;
power_acc = abs(Y_acc(1:floor(n/2))).^2;
power_newmark = abs(Y_newmark(1:floor(n/2))).^2;
nyquist = 1/2;
freq_disp = (1:n/2)/(n/2)*nyquist;
freq_acc = (1:n/2)/(n/2)*nyquist;
freq_newmark = (1:n/2)/(n/2)*nyquist;
freq_disp=freq_disp*100;
freq_acc=freq_acc*100;
freq_newmark=freq_newmark*100;

figure()
subplot(3,1,1)
plot(freq_disp(1:500),power_disp(1:500),'k')
title('Modal displacements: X direction')
xlabel('frequency')

subplot(3,1,2)
plot(freq_acc(1:500),power_acc(1:500),'k')
title('Modal accelerations: X direction')
xlabel('frequency')

subplot(3,1,3)
plot(freq_newmark(1:500),power_newmark(1:500),'k')
title('Direct time integration: X direction')
xlabel('frequency')


Y_disp=fft(q_displacement(dofList(13-2,2),:));
Y_acc = fft(q_acceleration(dofList(13-2,2),:));
Y_newmark=fft(q(:,dofList(13-2,2)));
Y_disp(1) = [];
Y_acc(1) = [];
Y_newmark(1) = [];

n = length(Y_disp);
power_disp = abs(Y_disp(1:floor(n/2))).^2;
power_acc = abs(Y_acc(1:floor(n/2))).^2;
power_newmark = abs(Y_newmark(1:floor(n/2))).^2;
nyquist = 1/2;
freq_disp = (1:n/2)/(n/2)*nyquist;
freq_acc = (1:n/2)/(n/2)*nyquist;
freq_newmark = (1:n/2)/(n/2)*nyquist;
freq_disp=freq_disp*100;
freq_acc=freq_acc*100;
freq_newmark=freq_newmark*100;

figure()
subplot(3,1,1)
plot(freq_disp(1:500),power_disp(1:500),'k')
title('Modal displacements: Y direction')
xlabel('frequency')
subplot(3,1,2)
plot(freq_acc(1:500),power_acc(1:500),'k')
title('Modal accelerations: Y direction')
xlabel('frequency')
subplot(3,1,3)
plot(freq_newmark(1:500),power_newmark(1:500),'k')
title('Direct time integration: Y direction')
xlabel('frequency')


Y_disp=fft(q_displacement(dofList(13-2,3),:));
Y_acc = fft(q_acceleration(dofList(13-2,3),:));
Y_newmark=fft(q(:,dofList(13-2,3)));
Y_disp(1) = [];
Y_acc(1) = [];
Y_newmark(1) = [];

n = length(Y_disp);
power_disp = abs(Y_disp(1:floor(n/2))).^2;
power_acc = abs(Y_acc(1:floor(n/2))).^2;
power_newmark = abs(Y_newmark(1:floor(n/2))).^2;
nyquist = 1/2;
freq_disp = (1:n/2)/(n/2)*nyquist;
freq_acc = (1:n/2)/(n/2)*nyquist;
freq_newmark = (1:n/2)/(n/2)*nyquist;
freq_disp=freq_disp*100;
freq_acc=freq_acc*100;
freq_newmark=freq_newmark*100;





figure()
subplot(3,1,1)
plot(freq_disp(1:500),power_disp(1:500),'k')
title('Modal displacements: Z direction')
xlabel('frequency')

subplot(3,1,2)
plot(freq_acc(1:500),power_acc(1:500),'k')
title('Modal accelerations: Z direction')
xlabel('frequency')

subplot(3,1,3)
plot(freq_newmark(1:500),power_newmark(1:500),'k')
title('Direct time integration: Z direction')
xlabel('frequency')



%% Plot modal displacement
%%
%%Dynamic response -- X axis
figure
hold on
plot(t,q_displacement(dofList(13-2,1),:),'b')
hold on
plot(t,q_acceleration(dofList(13-2,1),:),'r')
hold on
plot(t,q(:,dofList(13-2,1)),'g')
legend('Modal displacement','Modal acceleration','Time integration')
figure

subplot(3,1,1)
plot(t,q_displacement(dofList(13-2,1),:),'k')
title('Modal displacements:X direction')
xlabel('Time')
subplot(3,1,2)
plot(t,q_acceleration(dofList(13-2,1),:),'k')
title('Modal Accelerations:X direction')
xlabel('Time')
subplot(3,1,3)
plot(t,q(:,dofList(13-2,1)),'k')
title('Modal displacements:X direction')
xlabel('Time')

%%Dynamic response -- Y axis
figure
hold on
plot(t,q_displacement(dofList(13-2,2),:),'b')
hold on
plot(t,q_acceleration(dofList(13-2,2),:),'r')
hold on
plot(t,q(:,dofList(13-2,2)),'g')
legend('Modal displacement','Modal acceleration','Time integration')
title('Dynamic response Y-axis')

figure

subplot(3,1,1)
plot(t,q_displacement(dofList(13-2,2),:),'k')
title('Modal displacements:Y direction')
xlabel('Time')
subplot(3,1,2)
plot(t,q_acceleration(dofList(13-2,2),:),'k')
title('Modal Accelerations:Y direction')
xlabel('Time')
subplot(3,1,3)
plot(t,q(:,dofList(13-2,2)),'k')
title('Modal displacements:Y direction')
xlabel('Time')

%%Dynamic response -- Z axis
figure
hold on
plot(t,q_displacement(dofList(13-2,3),:),'b')
hold on
plot(t,q_acceleration(dofList(13-2,3),:),'r')
hold on
plot(t,q(:,dofList(13-2,3)),'g')
legend('Modal displacement','Modal acceleration','Time integration')
title('Dynamic response Z-axis')

figure

subplot(3,1,1)
plot(t,q_displacement(dofList(13-2,3),:),'k')
title('Modal displacements:Z direction')
xlabel('Time')
subplot(3,1,2)
plot(t,q_acceleration(dofList(13-2,3),:),'k')
title('Modal Accelerations:Z direction')
xlabel('Time')
subplot(3,1,3)
plot(t,q(:,dofList(13-2,3)),'k')
title('Modal displacements:Z direction')
xlabel('Time')

