% Data   X=dlmread('sensors_data.txt');
% goi lenh [V,U1, J] = FCMChuanWUSN(X1, 3, 2, 0.01, 1, BS, 0.002, 5);
function output = FCMChuanWUSN(X,C,m,Eps,maxTest, BS, alpha, kU)
%clc
output = {};
[N_init,r]=size(X);
efs = 10; % pJ/bit/m^2
emp = 0.0013;   
ground = 20;
BK = 50;
Eelec= 50; % nJ/bit = 10^-9 J
Eda =5 ; % pJ /bit = 10^-12 J
initE = 10;
up = 700;
%khoi tao random V
for i = 1:C
    for j = 1:r
        V(i,j) =  min(X(:,j))   + rand() *  (max(X(:,j)) - min(X(:,j)));
 % V(i,j) = randi(10);
 %            V(i,j) = min(X(:,j)); 
     end       
end
%Vong lap
dem = 0;
flag = 0;
r3 = 1;


matrixEnery = zeros(N_init, 3, 'single');
matrixEnery(:, 1) = initE;
matrixEnery(:,3) = 1;
matrixEnery(:,4) = 0;
EJ = [];
colors=hsv(C); 
J_min = initE *N_init;
lostNode = [];
aliveNode = [];
%stopNode = 2;
%stopNode = round(N *0.1)
numberDeadnodes = 0;
numDeadUn = 0;
repeatRounds = 2;
roundsStop = 0;
lostUnderNode = [];
X_next = X;
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
while (flag~=1)
    % Theo thuat toan FCM-Chuan de tim U,V
    [N,r]=size(X_next);
    while (1>0)
        for i = 1:N
           for j = 1:C
              %Tinh tu so
              TuBS = 0;
              TuXi = 0;
              TuXiU = 0;
              for k = 1:r
                  TuBS = TuBS + power(V(j,k)-BS(1,k),2);
                  TuXi = TuXi + power(X_next(i,k)-V(j,k),2);
              end    
              TuBS = sqrt(TuBS);
              TuXi = sqrt(TuXi);
              if X(i, 3) < ground
                  for k = 1:r
                   TuXiU = TuXi + power(X_next(i,k)-V(j,k),2);
                  end 
                   TuXBaiU = sqrt(TuXiU);
              end 

              Tu = emp * power(TuBS,4) + (N - C)*efs*power(TuXi, 2) + kU*((1/log(10)) + 8.69 * alpha)* TuXiU ;

              %Tinh mau so

              tong = 0;
              for l =1:C
                  MauBS = 0;
                  MauXk = 0;
                  MauXkU = 0;
                  for k = 1:r
                      MauBS = MauBS + power(V(l,k)-BS(1,k),2);
                      MauXk = MauXk + power(X_next(i,k)-V(l,k),2);
                      if(X(i, 3) < ground )
                          MauXkU = MauXkU + power(X_next(i,k)-V(l,k),2);
                      end
                  end    
                  MauBS = sqrt(MauBS);
                  MauXk = sqrt(MauXk);
                  MauXkU = sqrt(MauXkU);     

                  Mau = emp * power(MauBS,4) + (N - C)*efs*power(MauXk, 2) + kU*((1/log(10)) + 8.69 * alpha)* MauXkU ;

                  tong = tong + power(Tu/Mau, 1/(m-1)); 
              end
              U(j,i) = 1/tong;
           end    
        end
            %tinh khoang cach tongX(i) den BS)
            tongKC = [];
           for k1 = 1:r
                tempKC = 0;
               for i1 = 1:N
                      tempKC = tempKC +  (BS(1,k1) - X_next(i1,k1));
                end    
                tongKC(1,k1) = tempKC;
           end
            %Tinh V(t+1) tu U ra W
            for j = 1:C
               tongU = 0;
               for k = 1:N
                   tongU = tongU + power(U(j,k),m);
               end
               R = [-100 100];
               xi(j) = rand(1)*range(R) + min(R);
               Rz = [-10 10];
               xz(j) = rand(1)*range(Rz) + min(Rz);
               for t = 1:r
                   A = 4*C*emp*tongU;
                   B = 0;
                   CT = 2*(N-C)*efs*power(10,-12)*tongU;
                   D1 =  CT*tongKC(1,t);
                   D2 = kU*((1/log(10)) + 8.69*alpha)*tongU;
                   D = D1 - D2;
                   %CT1 = -k*((1/log(10)) + 8.69*alpha)*tongU;
                   W1 = power ((nthroot((27*A * power(D, 2) + 4 * power(CT,3)), 2)/ (2 * power(3,1.5) * A)) - D/(2*A), 1/3);
                   W2 = CT/ (3*A * power ((nthroot((27*A * power(D, 2) + 4 * power(CT,3)), 2)/ (2 * power(3,1.5) * A)) - D/(2*A), 1/3));

                   %W1_1 = power((nthroot((27*A*power(CT,2)*power(D1,2) + 54*A*CT*D1*(-D2)+ 27*A*power(D2,2) + 4 * power(CT,3)) ,2)/(2 * power(3,3/2) * A)) - ((D1*CT - D2)/(2*A)),1/3);
                   %W1_2 = B/(3*A*W1_1);

                   %kiem tra phan 
                   %y =  roots([A 0 B CT]);
                   %y(imag(y) ~= 0) = [];
                   if t ~= 3
                      W(j,t) = W1 - W2 + BS(1,t) + xi(j);
                   end
                   if t == 3
                      W(j,t) = W1 - W2 + BS(1,t) + xz(j);
                   end

               end
            end   
            zt = '00000000000000000000000000 lan'

            %So sanh W va V        
            saiso = 0;
            for i = 1:C
               for j = 1:r              
                  saiso = saiso + power(W(i,j)-V(i,j),2); 
               end
            end
            saiso = sqrt(saiso);

            %Kiem tra sai so voi Eps
            if (saiso <= (kU*Eps))
                break;
            else          
                %Lap tiep: Gan V = W
                for i7 = 1:C
                    for j = 1:r
                        V(i,j)= W(i,j);
                    end
                end
            end
             %Kiem tra voi so lan lap lon nhat
             if(dem>=maxTest)
                break;
             end
             %Tang so vong lap len
               dem = dem + 1;

     % saiso
    end
    
    MaxU=max(U)
    XTempt = X_next;
    %Dinh vi cum bang MaxU trong moi sensor
    U
    for k=1:N
        for j=1:C
            if U(j,k)==MaxU(k)
                XTempt(k,4) = j;
            end
        end
    end
    cluster = {};
  
    for j1 = 1:C
        clusterMembers = [];
        t1 = 1;
        for k1 = 1:N
            if XTempt(k1, 4) == j1
                clusterMembers(t1, :) = X_next(k1, :) ;
                t1 = t1 + 1;
            end
        end
        cluster{j1,1} = clusterMembers;
    end
    centers = zeros(C,3);
    clusterTempCH = cluster;
    %cai dat lai clusterTempCH
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%');
    for j1 = 1:C
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(isempty(cluster{j1,1}) == 0)
           [L1,r1]=size(cluster{j1,1});
          
           index = 0;
           flagStop = false;
           clusterTem = cluster{j1,1};
           clusterTem(:,4) = 8;
           numC = 0;
           while flagStop == false
              minTemp = 0;
              indexTemp = 1;
              
              for cl4 = 1: L1
                  if clusterTem(cl4, 4) ~=5
                      dis = calcSumDistDataPoint3X(V(j1, :), clusterTem(cl4, :));
                      if cl4 == 1
                          minTemp = dis;
                      else 
                          if dis < minTemp
                             indexTemp = cl4;
                          end
                      end 
                  end
              end
               
              if(clusterTem(indexTemp, 3) > 20)
                  %find index in cluster => set value
                  clusterTempCH{j1,1}(indexTemp, 4) = 1;
                  centers(j1,:) = cluster{j1,1}(indexTemp, :);  
                  flagStop = true;
                  index = indexTemp;
                  
                  indexC = find(ismember(X, cluster{j1,1}(indexTemp, :), 'row')) ;
                  matrixEnery(indexC, 4) = matrixEnery(indexC, 4) + 1;
              else 
                  for cl5 = 1:L1
                      if cl5 == indexTemp
                         clusterTem(indexTemp,4) = 5
                         numC = numC + 1;
                         break;
                      end
                  end
                  %%%%%% xet truong hop neu tat ca cac nut nam duoi long
                  %%%%%% dat
                  if numC == L1
                       for cl4 = 1: L1
                           dis = calcSumDistDataPoint3X(V(j1, :), clusterTem(cl4, :));
                           if cl4 == 1
                              minTemp = dis;
                           else 
                              if dis < minTemp
                                 indexTemp = cl4;
                              end
                           end 
                       end
                      clusterTempCH{j1,1}(indexTemp, 4) = 1;
                      indexC = find(ismember(X, cluster{j1,1}(indexTemp, :), 'row')) 
                      matrixEnery(indexC, 4) = matrixEnery(indexC, 4) + 1;
                      centers(j1,:) = cluster{j1,1}(indexTemp, :);  
                      flagStop = true;
                      index = indexTemp;
                  end 
              end 
              
          end 
          for cl1 = 1: L1
              if cl1 ~= index
                  disIn = calcSumDistDataPoint3X(cluster{j1,1}(index, :), cluster{j1,1}(cl1, :));
                  if (disIn/4) <= (2* BK)
                     clusterTempCH{j1,1}(cl1, 4) = 2;       
                    
                  else 
                     clusterTempCH{j1,1}(cl1, 4) = 0;
                  end
              end
           end 
       end
    end
    totalEnergyConsume = 0;
    for j2 = 1:C
        if(isempty(cluster{j2,1}) == 0)
               fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%% Cluster');
             [L2,r2]=size(cluster{j2,1});
             enerCluster = 0;
             numSenEachCluster = 0;
             indexCH = 0;
              for cl2 = 1: L2
                  if clusterTempCH{j2,1}(cl2, 4) == 2
                      enerU = 0;
                      disXToCH = calcSumDistDataPoint3X(centers(j2,:), cluster{j2,1}(cl2, :))
                      if(cluster{j2,1}(cl2, 3) < ground)
                           enerU = ((1/log(10)) + 8.69*alpha) *up * disXToCH;
                      end
                      enerCH = up*8*Eelec*(10^-9) + up*8*efs*10^-12*(disXToCH^2) + enerU;
                      indexX = find(ismember(X, cluster{j2,1}(cl2, :), 'row'));
                        if (matrixEnery(indexX, 1) >=enerCH)
                            matrixEnery(indexX, 1) = matrixEnery(indexX, 1) - enerCH;
                            matrixEnery(indexX, 2) = matrixEnery(indexX, 2) + enerCH;
                            %fprintf('vt = %d\n',vt);
                            enerCluster = enerCluster + enerCH;
                            numSenEachCluster = numSenEachCluster +1;
                        else
                            %matrixEnery(indexX, 1) = 0;
                            if  matrixEnery(indexX, 3) == 1
                                fprintf('vt dead node = %d\n',indexX)
                                numberDeadnodes = numberDeadnodes + 1
                                
                               if(clusterTempCH{j2,1}(cl2, 3) == 20)
                                   numDeadUn = numDeadUn + 1
                               end
                               matrixEnery(indexX, 3)  = 0
                            end 
%                              if numberDeadnodes >= stopNode
%                                   
%                                   if roundsStop >= repeatRounds
%                                     flag = 1;
%                                     fprintf('vt = %d\n',cl2);
%                                     
%                                   else 
%                                      roundsStop = roundsStop + 1;
%                                   end    
%                                     
%                              end
                        end
                       
                  end 
                  if clusterTempCH{j2,1}(cl2, 4) == 1
                      indexCH = cl2;
                  end
              end
              numSenEachCluster
              enerCluster
              enerBS = up*8*numSenEachCluster*Eelec*10^-9 + up*8*numSenEachCluster*Eda*10^-12
              dBS = calcSumDistDataPoint3X(cluster{j2,1}(indexCH, :),BS)
              enerBS =  enerBS +  up*8*emp*10^-12*dBS^4 %Nang luong ket noi den BS
              fprintf('enerBS = %f\n',enerBS);
              indexC = find(ismember(X, cluster{j2,1}(indexCH, :), 'row'));
            
              if ( matrixEnery(indexC, 1) >=enerBS)
                 matrixEnery(indexC, 1) = matrixEnery(indexC, 1) - enerBS;
                 matrixEnery(indexC, 2) =  matrixEnery(indexC, 2) + enerBS;
                 enerCluster = enerCluster + enerBS;      
                 totalEnergyConsume = totalEnergyConsume + enerCluster;
              else
                 if matrixEnery(indexC, 3)  == 1
                     %matrixEnery(indexC, 1) = 0;
                     fprintf('vt dead node for CH = %d\n',indexC);
                     numberDeadnodes = numberDeadnodes + 1;
                     if(cluster{j2,1}(indexCH, 3) == 20)
                         numDeadUn = numDeadUn + 1;
                     end
                     matrixEnery(indexC, 3) = 0;
                 end
%                   if numberDeadnodes >= stopNode
%                        
%                         if roundsStop >= repeatRounds
%                           flag = 1;
%                           fprintf('vt = %d\n',cl2);
%                         else 
%                           roundsStop = roundsStop + 1;
%                         end 
%                         
%                    end 
                
              end   
               
                 
        end
    end
   
   totalConnectNodes = 0; 
   for j2 = 1:C 
     if(isempty(cluster{j2,1}) == 0)
        [L2,r2]=size(cluster{j2,1});
        for cl2 = 1: L2
             if clusterTempCH{j2,1}(cl2, 4) ~= 0
                 totalConnectNodes = totalConnectNodes + 1;
             end
        end 
     end 
   end
   
 %  figure, hold on
    %for i5=1:C
    %  if(isempty(cluster{i5,1}) == 0)
    %    scatter3(cluster{i5,1}(:,1),cluster{i5,1}(:,2), cluster{i5,1}(:,3),[],colors(i5,:),'filled'); 
    %  end
    %end 
 
      
    %scatter3(cluster(:,1),centers(:,2), centers(:,3),[],'Marker','o'); 
   % plot3(centers(:,1), centers(:,2), centers(:,3), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
   % hold off
    %title(['Cluster Result Round ' num2str(r3)]);
    %xlabel('X');
    %ylabel('Y');
    %zlabel('Z');
    
%     matrixEnery
%     J_min_truoc= J_min
    
%     roundsStop
%     stopNode
%     numberDeadnodes
    J_min = round(J_min - totalEnergyConsume,10)
    0.5 *initE * N_init
    EJ (r3) = J_min;
    lostNode(r3) = N_init- numberDeadnodes;
    if (J_min < (0.03 *initE *N_init) ||  numberDeadnodes > (0.8* N_init))
       flag = 1;
    else
        if totalEnergyConsume == 0
            if roundsStop >= repeatRounds 
                flag = 1;
                fprintf('stop because here');
            else 
                roundsStop = roundsStop + 1;
            end 
        end 
    end
    
    r3 = r3 + 1;
    lostUnderNode(r3)=  kU - numDeadUn
    
    X_next = [];
    indexNext = 1;
    for ik = 1: N_init
        if  matrixEnery(ik, 2) < 9.99999 &&  matrixEnery(ik, 3) == 1
            X_next(indexNext, :) = X(ik, :);
         
            indexNext = indexNext + 1;
        end 
    end
   
end

  %scatter3(cluster{2,1}(:,1),cluster{2,1}(:,2), cluster{2,1}(:,3),[],colors(1,:),'filled'); 
  %[SN, EJ1] = LEACH_CC();
  [EJFCM, lostNumberFCM, lostNumberUnderGroundFCM] = FCMChuanRouting(X,C,m,Eps,maxTest, BS, alpha, kU);

%   figure(1)
%   plot(EJ, 'b')
%   hold on
%   plot(EJFCM, 'r') 
%   title({'Total energy consumption of FCM-WUSN & FCM '});
%   legend({'FCM-WUSN','FCM' },'Location','Best')
%   xlabel('Rounds');
%   ylabel('Energy Consumption');
%   
%   hold off;
%   
%   figure(3)
%   yyaxis left
%   xlabel('Rounds');
%   ylabel('Alived Nodes ')
%   hold on
%   plot(lostNode, 'b')
%   hold on
%   plot(lostNumberFCM, 'g') 
%   title({'Total alived nodes '});
%   
%   yyaxis right
%   ylabel('Alived underground nodes')
%   hold on
%   plot(lostUnderNode, 'y')
%   hold on
%   plot(lostNumberUnderGroundFCM, 'r') 
%   legend({'Total alived nodes(FCM-WUSN)','Total alived nodes(FCM)', 'Alived un-nodes(FCM-WUSN) ', 'Alived un-nodes(FCM) '},'Location','Best')
%   %legend({'Total alived nodes(FCM-WUSN)'},'Location','Best')
%   hold off;
  
  output{1,1} = EJ;
  %csvwrite(fileName_FCMWUSN_EC,EJ);
  
  output{2,1} = EJFCM;
  %csvwrite(fileNameFCM_EC,EJFCM);
  
  output{3,1} = lostNode;
  %csvwrite(fileNameULN,lostNode);
  
  output{4,1} = lostNumberFCM;
  %csvwrite(fileName_FCM_LN,lostNumberFCM);
  
  output{5,1} = lostUnderNode;
  %csvwrite(fileName_FCMWUSN_LUN,lostUnderNode);
  
  output{6,1} = lostNumberUnderGroundFCM;
  %csvwrite(fileName_FCM_LUN,lostNumberUnderGroundFCM);
  
end


function ch = findCHInCluster (cluster, V, j1 )
   if(isempty(cluster) == 0)
           [L1,r1]=size(cluster);
          
           index = 0;
           flagStop = false;
           clusterTem = cluster
           clusterTem(:,4) = 8;
           while flagStop == false
              minTemp = 0;
              indexTemp = 1;
              
              for cl4 = 1: L1
                  if clusterTem(cl4, 4) ~=5
                      dis = calcSumDistDataPoint3X(V(j1, :), clusterTem(cl4, :))
                      if cl4 == 1
                          minTemp = dis;
                      else 
                          if dis < minTemp
                             indexTemp = cl4;
                          end
                      end 
                  end
              end
               
              if(clusterTem(indexTemp, 3) > ground)
                  %find index in cluster => set value
                  cluster{j1,1}(index, 4) = 1;
                  ch = cluster{j1,1}(index, :);  
                  flagStop = true;
                 
              else 
                  clusterTem = [];
                  for cl5 = 1:L1
                      if cl5 == indexTemp
                         clusterTem(:,4) = 5;
                      end
                  end
 
              end 
           end 
          
       end
end
function dt = calcSumDistDataPoint3X(data, X)
dt = sqrt((data(1) - X(1))^2 + (data(2)- X(2))^2 + (data(3) - X(3))^2);

end

function t = calcSumDistDataPoint2X(data, X)
temp = data - X(ones(size(data, 1), 1), :);
temp = temp.^2;
temp = sum(temp, 2);
temp = sqrt(temp);
t = sum(temp);
end

function PBM_value = PBM(numClust, data, center, U)
E_1 = calcSumDistDataPoint2X(data, mean(data));

maxU = max(U);
E_k = 0;
for i = 1 : numClust
    index = find(U(i,:) == maxU);
    clustData = data(index, :);
    E_k = E_k +  calcSumDistDataPoint2X(clustData, center(i, :));
end

D_k = 0;
for i = 1:numClust-1
    for j = i+1:numClust
        D_k = max(D_k, norm(center(i, :) - center(j, :)));
    end
end
    
PBM_value = (E_1 * D_k / (numClust * E_k)) ^ 2;
end


function SWC_value = SWC(numClust, data, U)
maxU = max(U);
SWC_value = 0;

wb = waitbar(0);
dem = 0;

for i = 1 : numClust
    index = find(U(i,:) == maxU);
    if size(index, 2) == 1
        continue;
    end
    clustData = data(index, :);
    for j = index
        a_i_j = calcSumDistDataPoint2X(clustData, data(j, :)) / size(index, 2);
        b_i_j = inf;
        for k = 1 : numClust
            if k ~= i
                index_k = find(U(k,:) == maxU);
                clustData_k = data(index_k, :);
                d_k_j = calcSumDistDataPoint2X(clustData_k, data(j, :)) / size(index_k, 2);
                b_i_j = min(b_i_j, d_k_j);
            end
        end
        SWC_value = SWC_value + (b_i_j - a_i_j) / max(a_i_j, b_i_j);
        dem = dem + 1;
        
%         waitbar(dem / size(data, 1), wb);
    end
end
SWC_value = SWC_value / size(data, 1);
end


        
function DB_value = DB(numClust, data, center, U)
    maxU = max(U);
    
    for i = 1 : numClust
        index = find(U(i,:) == maxU);
        clustData = data(index, :);
        stdData = std(clustData, 1);
        S(i) = sqrt(stdData(1)^2 + stdData(2)^2 + stdData(3)^2);
    end
    
    DB_value = 0;
    
    for i = 1:numClust
        maxSM = 0;
        for j = 1:numClust
           if j ~= i
               temp = (S(i) + S(j))/norm(center(i, :) - center(j, :));
               maxSM = max(maxSM, temp);
           end
        end
        DB_value = DB_value + maxSM;
    end
    
    DB_value = DB_value/numClust;
end

    
function IFV_value = IFV(numClust, data, sizeData, center, U)
    sigmaD = 0;
    sum = 0;
    
    for i = 1:numClust
        tg1 = 0;
        tg2 = 0;
        for j = 1:sizeData
            if U(j,i) == 0 
                U(j, i) = eps;
            end
            if U(j, i) == 1 
                U(j, i) = 1 - eps;
            end
            
            tg1 = tg1 + log(U(j, i))/log(2);
            tg2 = tg2 + U(j, i)^2;
            sigmaD = sigmaD + norm(data(j, :) - center(i, :))^2;
        end
        
        tg = (log(numClust)/log(2) - tg1/sizeData)^2;
        tg2 = tg2/sizeData;
        
        sum = sum + tg * tg2;
    end
    
    sigmaD = sigmaD/(numClust * sizeData);
    
    calcSDmax = 0;
    for i = 1:numClust-1
        for j = i+1:numClust
            calcSDmax = max(calcSDmax, norm(center(i, :) - center(j, :))^2);
        end
    end
    
    IFV_value  = (sum * calcSDmax) / (sigmaD * numClust);
end