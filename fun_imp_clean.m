%% load data positions and features
clear 
close all
tic
filename='Aggregation.txt';%Trips_Sync Spiral R15 Pathbased Jain Flame D31 Compound Aggregation;
path=['./data/',filename];
[lons,lats,id]=textread(path, '%f,%f,%d');
dataPts=[lons,lats];
%dataPts=dataPts(1:50:end,:);
%normalize
dataPts = normalize(dataPts);
locations = dataPts';
ND=size(locations,2);
fprintf('number of input data: %d\n',ND);
%% projection into grids
unit = 0.001; %projection size for normalized data
minx = min(locations,[],2);
maxx = max(locations,[],2);
bound=maxx-minx;
nbins=uint16(bound./unit)';
edges1 = linspace(minx(1), maxx(1), nbins(1)+1);
edges1 = [-Inf edges1(2:end-1) Inf];
edges2 = linspace(minx(2), maxx(2), nbins(2)+1);
edges2 = [-Inf edges2(2:end-1) Inf];
p = size(locations,2);
bin = zeros(p,2);
[~,bin(:,2)] = histc(locations(1,:),edges1);
[~,bin(:,1)] = histc(locations(2,:),edges2);
H= accumarray(bin,1,nbins([2 1]));
index_all=H>0;
%% middle filter to get rid of isolated peaks
%H = medfilt2(H_ori,[11,11]);
%% find neighborhood distance with integral image
Percent = 0.02; % 2% the super-parameter
DC = mex_dc_evaluate(H,ND,Percent);
dc = DC*unit;
%dc=4;
%% find representation points in H
threshold = 0;
Flag = H>threshold; % find pixels contain points
ND_rep = sum(Flag(index_all)); % number of representation points
% ND_rep = ND;
fprintf('number of remain points: %d\n',ND_rep);
data_rep = zeros(2,ND_rep);
loc_rep = zeros(2,ND_rep);
den_rep = zeros(1,ND_rep);
k=1;
for i=1:ND
    x = bin(i,1);
    y = bin(i,2);
    if Flag(x,y)
        data_rep(:,k)=locations(:,i);
        loc_rep(:,k)=[x;y];
        den_rep(k) = H(x,y);
        Flag(x,y) = 0;
        k=k+1;
    end
end
%% calculate densities with permutohedral algorithm
features = data_rep; %loc_rep;
filter_in = den_rep;
sigma_vec = dc;  %dc*2 dc*2 dc dc dc/2 dc dc/2.0 dc dc
sigma_dim = 2;
filter_out = mex_permutohedral(features, filter_in, sigma_vec, sigma_dim);
rho = filter_out; 
%% distance matrix used by distance between different points
x = 1:size(H,2);
y = 1:size(H,1);
[X,Y] = meshgrid(x,y);
grid_dist = sqrt(X.^2+Y.^2); %the distance matrix
%% splits the whole area into sub-regions
delta = zeros(1,ND_rep);
nneigh = zeros(1,ND_rep);
[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;
cnt=1;
for ii=2:ND_rep
    delta(ordrho(ii))= grid_dist(end,end);
    for jj=1:ii-1
        pos_diff = abs(loc_rep(:,ordrho(ii)) - loc_rep(:,ordrho(jj)))+1;
        tmp_dist = grid_dist(pos_diff(1),pos_diff(2));
        
        %xx(cnt,1) = tmp_dist/dist(ordrho(ii),ordrho(jj));
        cnt = cnt + 1;
        
        if(tmp_dist <delta(ordrho(ii)))
            delta(ordrho(ii))=tmp_dist;
            nneigh(ordrho(ii))=ordrho(jj);
        end
    end
end
delta(ordrho(1))=max(delta(:));
disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')
name1='DECISION_GRAPH1';
fid = fopen(name1, 'w');
for i=1:ND_rep
    fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
end
fclose(fid);
%% choose peaks
disp('Select a rectangle enclosing cluster centers')
scrsz = get(0,'ScreenSize');
figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);
ind = 1:ND_rep;
gamma =rho.*delta;
[gamma_sorted,gamma_indx]=sort(gamma,'descend'); %Ö¸±êÅÅÐò

subplot(2,1,1)
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')
toc
% get the cluster by ploting a rectange on the figure
subplot(2,1,1)
rect = getrect(1);

tic % restart the timer
rhomin=rect(1);
deltamin=rect(2);
NCLUST=0; % number of cluster centers
cl = zeros(1,ND_rep)-1;
icl = zeros(1,ND_rep); %inverse id of cluster result
for i=1:ND_rep
    if ( (rho(i)>rhomin) && (delta(i)>deltamin))
        NCLUST=NCLUST+1;
        cl(i)=NCLUST;
        icl(NCLUST)=i;
    end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
%% assignation
disp('Performing assignation')
for i=1:ND_rep
    if (cl(ordrho(i))==-1)
        cl(ordrho(i)) = cl(nneigh(ordrho(i)));
    end
end
%% no halo
halo = cl;
%% with halo
%  bord_rho = zeros(1,NCLUST);
%  if (NCLUST>1)
%     for i=1:ND_rep-1
%         for j=i+1:ND_rep
%             pos_diff = abs(loc_rep(:,i) - loc_rep(:,j))+1;
%             tmp_dist = grid_dist(pos_diff(1),pos_diff(2));
%             if ((cl(i)~=cl(j))&&(tmp_dist<=DC))
%                 rho_aver=(rho(i)+rho(j))/2.;
%                 if (rho_aver>bord_rho(cl(i)))
%                     bord_rho(cl(i))=rho_aver;
%                 end
%                 if (rho_aver>bord_rho(cl(j)))
%                     bord_rho(cl(j))=rho_aver;
%                 end
%             end
%         end
%     end
%     for i=1:ND_rep
%         if (rho(i)<bord_rho(cl(i)))
%             halo(i)=0;
%         end
%     end
% end
 
% %print the statistical result
%  for i=1:NCLUST
%    nc=0;
%    nh=0;
%    for j=1:ND_rep
%      if (cl(j)==i) 
%        nc=nc+1;
%      end
%      if (halo(j)==i) 
%        nh=nh+1;
%      end
%    end
%    fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
%  end
%% plot the result
cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   subplot(2,1,1)
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
subplot(2,1,2)
Y1 = data_rep';
plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Improved Cluster result:Trips_Sync','FontSize',12.0)
xlabel ('X')
ylabel ('Y')
A=zeros(ND_rep,2);
for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND_rep
    if (halo(j)==i)
      nn=nn+1;
      A(nn,1)=Y1(j,1);
      A(nn,2)=Y1(j,2);
    end
  end
  hold on
  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
%% plot cluster result
figure,hold on
set (gcf,'Position',[300,300,500,500], 'color','w')
set(gca,'Position',[.08 .12 .72 .76]);
cmap=[1,0,0;0,1,0;0,0,1];
tchar=['*','x','^','.','d','p','h','v','+','<'];
tag=cell(1,NCLUST+1);
%tag=cell(1,NCLUST+2);
%plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Clustering of Aggregation by GC-DPC','FontSize',12.0)
A=zeros(ND,2);
for i=1:NCLUST
  nn=0;
  %ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND_rep
    if (halo(j)==i)
      nn=nn+1;
      A(nn,1)=Y1(j,1);
      A(nn,2)=Y1(j,2);
    end
  end
  plot(A(1:nn,1),A(1:nn,2),tchar(ceil(i/3)),'MarkerSize',3,'MarkerFaceColor',cmap(mod(i,3)+1,:),'MarkerEdgeColor',cmap(mod(i,3)+1,:));
  tag(i)={strcat(num2str(i),'-prototy')};
end
index=halo==0;
A=Y1(index,:);
plot(Y1(icl(1:NCLUST),1),Y1(icl(1:NCLUST),2),'s','MarkerSize',10,'MarkerFaceColor',[1,1,0],'MarkerEdgeColor',[1,0,1]);
tag(end)={'-center'};
% if sum(index)>0
%   plot(A(:,1),A(:,2),tchar(ceil((NCLUST+1)/3)),'MarkerSize',3,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0]);
%   tag(end-1)={'-halo'};
% end
legend(tag);
%% project the represent points back
halo_all = zeros(1,ND); % labels for all
cl_all =  zeros(1,ND);
icl_all = zeros(1,ND);
for i = 1:ND_rep
    x = loc_rep(1,i);
    y = loc_rep(2,i);
    halo_all(bin(:,1)==x & bin(:,2)==y) = halo(i);
    cl_all(bin(:,1)==x & bin(:,2)==y) = cl(i);
    if icl(i)
        x = loc_rep(1,icl(i));
        y = loc_rep(2,icl(i));
        indx = find(bin(:,1)==x & bin(:,2)==y);
        icl_all(i) = indx(1);
        
    end
end
%% write the result into file
faa = fopen(strcat('result_',filename(1:end-3)), 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster centers')
disp('column 3:cluster assignation without halo control')
disp('column 4:cluster assignation with halo control')
disp('column 5:ground true id');
for i=1:ND
   fprintf(faa, '%i %i %i %i %i\n',i,icl_all(i),cl_all(i),halo_all(i),id(i));
end
fclose(faa);
toc