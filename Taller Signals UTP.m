%%%%%%%% Taller Señales UTP 05.02.18 %%%%%%%%%

%%% creado por Sisi Ma
%%% actualizado silvia 03.15.16

%% 1 EEG ANALYSIS WITH PCA
%% 1.1 PCA is sensitive to outliers

% 1. Download data from
% https://archive.ics.uci.edu/ml/datasets/EEG+Eye+State

% 2. Import the data. Compute the sampling frequency of the EEG data.
Tab = readtable('EEGEyeState2.txt','Delimiter',',','ReadVariableNames',false);
Tab.Properties.VariableNames = {'AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8',...
    'FC6','F4','F8','AF4','eyeDetection'};
samplingRate = height(Tab)/117;
sprintf('The sampling rate is: %g Hz',samplingRate)

% 3. Plot the data collected on each EEG probe with box plot.
D = table2array(Tab(:,1:14));
figure(1)
boxplot(D)
xlabel('Electrode','Fontweight','bold','FontSize',16)
ylabel('Voltage (mV)','Fontweight','bold','FontSize',16)
title('Boxplot of EEG data','Fontweight','bold','FontSize',18)
box off
set(gca,'TickDir','out')
set(gcf,'Color','w')

% to find linear indices of outliers execute
% gname
% place your cursor on outlier and click on mouse, press return when done.

% 4. PCA analysis of EEG data
[COEFF, SCORE, LATENT, ~, EXPLAINED] = pca(D);
comp = 1:1:14;
figure(2)
plot(comp,EXPLAINED,'b-o')
xlabel('Principal Component','Fontweight','bold','FontSize',16)
ylabel('% Variance Explained','Fontweight','bold','FontSize',16)
title('Scree plot of PCA','Fontweight','bold','FontSize',18)
box off
set(gca,'TickDir','out')
set(gcf,'Color','w')

% 5. Now delete the outliers from you dataset. 
Do = D;
[i,~] = find(D>10^5);
i = unique(i);
Do(i,:) = nan;
[i,~] = find(D<10^2);
i = unique(i);
Do(i,:) = nan;
Do = Do(~any(isnan(Do),2),:);

% 6. Conduct principal component analysis on the EEG data. Select the same
% number of principle component as you did before (so that we can visually
% compare the biplots). plot the biplot.
[COEFFo, SCOREo, LATENTo, ~, EXPLAINEDo] = pca(Do);
% For the biplot, I will compute the scores (the representation of the data
% in the new space and the loadings (eigenvectors).

figure(3)
subplot(1,2,1)
b1 = biplot(COEFFo(:,1:3),'scores',SCOREo(:,1:3),'varlabels',...
    Tab.Properties.VariableNames(1:14));
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
title('3 first PC without outliers')
subplot(1,2,2)
b2 = biplot(COEFF(:,1:3),'scores',SCORE(:,1:3),'varlabels',...
    Tab.Properties.VariableNames(1:14));
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
title('3 first PC with outliers')
hold off

% Compare the biplot generate from data with or without outliers removed.
% Is PCA sensitive to outliers? Find two probes that are close on the
% biplot generated from data with outliers but are far on the biplot
% generated from data without outliers. Explain the discrepancy.

%% 1.2 PCA interpretation

% 1. From now on, use the data with outliers removed. Determine number of
% components to keep from the scree plot.

figure(4)
plot(comp,EXPLAINEDo,'m-o')
xlabel('Principal Component','Fontweight','bold','FontSize',16)
ylabel('% Variance Explained','Fontweight','bold','FontSize',16)
title('Scree plot of PCA without outliers','Fontweight','bold','FontSize',18)
box off
set(gca,'TickDir','out')
set(gcf,'Color','w')

% I will chose three PCs for visualization purposes.

% 2. Project the original data to the direction of the first k principle
% components respectively (k is the number of principle component you
% selected in the previous step). plot the distribution of projected data
% (one histogram per direction).

% I will compute these projections by doing the dot product of the data
% matrix (de-meaned so distributions will be centered around 0) with each
% of the first principal components
projections = nan(length(Do),3);

figure(5)
for k = 1:3;
    projections(:,k) = Do*COEFFo(:,k);
    subplot(2,2,k)
    hist(projections(:,k),1000000)
    xlabel('Principal Component','Fontweight','bold','FontSize',14)
    ylabel('# Data points','Fontweight','bold','FontSize',14)
    title(sprintf('Histogram of projection on PC%g',k),'Fontweight',...
        'bold','FontSize',14)
end

% 3. Compute the variance of the data after the projecting to the direction
% of the first k principle components respectively. Verify that the
% variance of the projected data decreases.
var(projections(:,1))
var(projections(:,2))
var(projections(:,3))

% 4. Reconstruct the recordings on probe F7 using the selected principle
% components. Which principle component contribute the most to the
% reconstruction of probe F7, why? Quantify the Reconstruction Accuracy
% (you may need to take care of the effect of centering depends on the
% metric you chose).
X1 = SCOREo(:,1) * COEFFo(:,1)';
X1 = bsxfun(@plus,mean(Do,1),X1);

X2 = SCOREo(:,2) * COEFFo(:,2)';
X2 = bsxfun(@plus,mean(Do,1),X2);

X3 = SCOREo(:,3) * COEFFo(:,3)';
X3 = bsxfun(@plus,mean(Do,1),X3);

error1 = norm(Do(:,2) - X1(:,2))
error2 = norm(Do(:,2) - X2(:,2))
error3 = norm(Do(:,2) - X3(:,2))

% the third principal component contributes the most to the reconstruction
% of probe F7 because it's error is the smallest.

Xapprox = SCOREo(:,1:3) * COEFFo(:,1:3)';
Xapprox = bsxfun(@plus,mean(Do,1),Xapprox);

% 5. Compare the original recording on probe F7 and the reconstructed data
% with a plot of your choice.
figure(6)
s = scatter(Do(:,2),Xapprox(:,2),80);
hold on
s.MarkerFaceColor = 'm';
s.MarkerEdgeColor = 'w';
l = refline(1,0);
l.LineStyle = '--';
l.Color = 'k';
l.LineWidth = 3;
xlabel('Original F7 recordings','Fontweight','bold','FontSize',14)
ylabel('Reconstructed data','Fontweight','bold','FontSize',14)
title('Data Reconstruction using 3 PCs','Fontweight','bold','FontSize',18)


%% 1.3. PCA and Data Visualization

% 1. Project all data points onto the k (if your k is > 3, set k = 3 for
% visualization purpose) dimensional principle component space.
% had already done this above, now for 3D:
projections(:,1:3);

% 2. Plot the projection of the data in the k dimensional space defined by
% the principle components (k dimensional scatter plot). Use different
% color or symbol to represent the data collected at different eye states
% (column 15). Do you see a separation between the data collected during
% different eye states?
eyes = Tab.eyeDetection;
eyes(899,:) = nan;
eyes(10387,:) = nan;
eyes(13180,:) = nan;
eyes(11510,:) = nan;
eyes = eyes(~any(isnan(eyes),2),:);

idC = find(eyes==1); % eyes closed is coded as 1
idO = find(eyes==0); % eyes open is coded as 0

figure(7)
s1 = scatter3(projections(idC,1),projections(idC,2),projections(idC,3),100);
hold on
s1.MarkerFaceColor = 'b';
s1.MarkerEdgeColor = 'b';
s2 = scatter3(projections(idO,1),projections(idO,2),projections(idO,3),100);
s2.MarkerFaceColor = 'r';
s2.MarkerEdgeColor = 'r';
xlabel('PC 1','Fontweight','bold','FontSize',14)
ylabel('PC 2','Fontweight','bold','FontSize',14)
zlabel('PC 3','Fontweight','bold','FontSize',14)
title('Data Projected on first 3 PCs','Fontweight','bold','FontSize',18)
legend('Eyes Closed','Eyes Open','Location','Best');

% 3. Now, let's incorporate the information that the data are recorded
% sequentially in time. plot the projection of the data to the first k
% principle components respectively, with time on the x-axis and voltage on
% the y-axis. Use different color or symbol to represent the data collected
% at different eye states (column 15). What do you observe?

figure(8)
for k = 1:3
subplot(2,2,k)
s1 = scatter(idC,projections(idC,k),4);
hold on
s1.MarkerFaceColor = 'b';
s1.MarkerEdgeColor = 'b';
s2 = scatter(idO,projections(idO,k),4);
s2.MarkerFaceColor = 'r';
s2.MarkerEdgeColor = 'r';
xlabel('time','Fontweight','bold','FontSize',14)
ylabel('voltage (mV)','Fontweight','bold','FontSize',14)
legend('Eyes Closed','Eyes Open','Location','Best');
title(sprintf('Data Projected on PC%g',k),'Fontweight','bold','FontSize',14)
end

% 4. How can you predict different eye states using information from EEG
% recordings? State the rule in words.
% See report.


% 5. Implement your rule and graph the result. Quantify the accuracy of
% your prediction.
                               
ix0 = unique([idC(1) idC(diff([0 idC'])>1)']);                                       
ix1 = unique([idO(1) idO(diff([0 idO'])>1)']);        
ixv = sort([ix0 ix1 length(eyes)]);                    
for k1 = 1:length(ixv)-1
    section{k1} = (ixv(k1):ixv(k1+1)-1)';
end
celldisp(section)

% we will train the classifier on 20 of those sections and then test the
% classfier on 4.

for x = 1:20
    if mod(x,2) == 0
    closed = [closed; projections(section{x},:)];
    else
    open = [open; projections(section{x},:)];
    end
end

Group1 = cell(length(closed),1);
Group1(:) = {'Closed'};
Group2 = cell(length(open),1);
Group2(:) = {'Open'};

Train = [closed;open];
Class = [Group1;Group2];

obj = fitcdiscr(Train,Class);

[label1,~] = predict(obj,projections(section{21},:));
[label2,~] = predict(obj,projections(section{22},:));
[label3,~] = predict(obj,projections(section{23},:));
[label4,~] = predict(obj,projections(section{24},:));

sum(strcmp(label1,'Open'))/length(label1)*100
sum(strcmp(label2,'Closed'))/length(label2)*100
sum(strcmp(label3,'Open'))/length(label3)*100
sum(strcmp(label4,'Closed'))/length(label4)*100

%% 2. PCA vs. ICA

% 1. What is the key difference between PCA and ICA?

% 2. Create a dataset such that PCA and ICA recovers similar structure from
% the data. Use proper graphs to illustrate.
mu=[0,0,0];

sigma_x1=1;
sigma_x2=2;
sigma_x3=3;
rho = 0.6;
rho2 = 0.4;
rho3 = 0.7;

sigma=[sigma_x1^2,rho*sigma_x1*sigma_x2,rho2*sigma_x1*sigma_x3;...
    rho*sigma_x1*sigma_x2,sigma_x2^2,rho3*sigma_x2*sigma_x3;...
    rho2*sigma_x3*sigma_x1,rho3*sigma_x3*sigma_x2,sigma_x3^2];
rng default  % For reproducibility
S=mvnrnd(mu,sigma,1000);

[COEFFs, SCOREs, LATENTs, ~, EXPLAINEDs] = pca(S);
[Zica,A,T,Mu]  = myICA(S',2);
Zr = T \ pinv(A) * Zica + repmat(Mu,1,1000);

figure(9)
g = scatter3(S(:,1),S(:,2),S(:,3),4);
g.MarkerFaceColor = 'k';
g.MarkerEdgeColor = 'k';
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
title('Data in original space')
set(gca,'fontsize',14);
set(gcf,'color','w');
axis square;
s_mean=mean(S);
s_cov=cov(S);
disp(['means of sample : ' num2str(s_mean)])
disp('covariance of sample : ')
disp(num2str(s_cov));
hold on;
p = plot3([0 LATENTs(1,1)*COEFFs(1,1)],[0 LATENTs(1,1)*COEFFs(2,1)],...
    [0 LATENTs(1,1)*COEFFs(3,1)],'r-');
p2 = plot3([0 LATENTs(2,1)*COEFFs(1,2)],[0 LATENTs(2,1)*COEFFs(2,2)],...
    [0 LATENTs(2,1)*COEFFs(3,2)],'b-');
p3 = plot3([0 7*A(1,1)],[0 7*A(1,2)],[0 7*A(1,3)],'m-');
p4 = plot3([0 7*A(2,1)],[0 7*A(2,2)],[0 7*A(2,3)],'g-');
p.LineWidth = 4;p2.LineWidth = 4;p3.LineWidth = 4;p4.LineWidth = 4;
t1 = text(LATENTs(1,1)*COEFFs(1,1),LATENTs(1,1)*COEFFs(2,1),LATENTs(1,1)...
    *COEFFs(3,1),'PC 1');
t2 = text(LATENTs(2,1)*COEFFs(1,2),LATENTs(2,1)*COEFFs(2,2),LATENTs(2,1)...
    *COEFFs(3,2),'PC 2');
t3 = text(7*A(1,1),7*A(1,2),7*A(1,3),'IC 1');
t4 = text(7*A(2,1),7*A(2,2),7*A(2,3),'IC 2');
t1.Color = 'r'; t1.FontWeight = 'bold';t2.Color = 'b'; t2.FontWeight = 'bold';
t3.Color = 'm'; t3.FontWeight = 'bold';t4.Color = 'g'; t4.FontWeight = 'bold';

% 3. Create a dataset such that PCA and ICA recovers distinct structure
% from the data. Use proper graphs to illustrate.
m = 3;
n= 1000;
mu = 0;
sigma = 2;
skew = 1;
kurt = 4;

[r,type] = pearsrnd(mu,sigma,skew,kurt,m,n);
R = r';
[COEFFr, SCOREr, LATENTr, ~, EXPLAINEDr] = pca(R);
[Zicar,Ar,Tr,Mur]  = myICA(R',2);

figure(10)
g1 = scatter3(R(:,1),R(:,2),R(:,3),4);
g1.MarkerFaceColor = 'k';
g1.MarkerEdgeColor = 'k';
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
title('Data in original space')
set(gca,'fontsize',14);
set(gcf,'color','w');
axis square;
s_mean=mean(S);
s_cov=cov(S);
disp(['means of sample : ' num2str(s_mean)])
disp('covariance of sample : ')
disp(num2str(s_cov));
hold on;
p = plot3([0 LATENTr(1,1)*COEFFr(1,1)],[0 LATENTr(1,1)*COEFFr(2,1)],...
    [0 LATENTr(1,1)*COEFFr(3,1)],'r-');
p2 = plot3([0 LATENTr(2,1)*COEFFr(1,2)],[0 LATENTr(2,1)*COEFFs(2,2)],...
    [0 LATENTr(2,1)*COEFFr(3,2)],'b-');
p3 = plot3([0 7*Ar(1,1)],[0 7*Ar(1,2)],[0 7*Ar(1,3)],'m-');
p4 = plot3([0 7*Ar(2,1)],[0 7*Ar(2,2)],[0 7*Ar(2,3)],'g-');
p.LineWidth = 4;p2.LineWidth = 4;p3.LineWidth = 4;p4.LineWidth = 4;
t1 = text(LATENTr(1,1)*COEFFr(1,1),LATENTr(1,1)*COEFFr(2,1),LATENTr(1,1)...
    *COEFFr(3,1),'PC 1');
t2 = text(LATENTr(2,1)*COEFFr(1,2),LATENTr(2,1)*COEFFr(2,2),LATENTr(2,1)...
    *COEFFr(3,2),'PC 2');
t3 = text(7*Ar(1,1),7*Ar(1,2),7*Ar(1,3),'IC 1');
t4 = text(7*Ar(2,1),7*Ar(2,2),7*Ar(2,3),'IC 2');
t1.Color = 'r'; t1.FontWeight = 'bold';t2.Color = 'b'; t2.FontWeight = 'bold';
t3.Color = 'm'; t3.FontWeight = 'bold';t4.Color = 'g'; t4.FontWeight = 'bold';











