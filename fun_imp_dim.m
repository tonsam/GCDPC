%% improved DPC for any dimention input
clear 
close all
tic
%filename='./data/dim6.txt';
%input = textread(filename);
%dim = size(input,2);
filename='R15.txt';%S1 Spiral R15 Pathbased Jain Flame D31 Compound Aggregation;
path=['./data/',filename];
[lons,lats,id_cl]=textread(path, '%f,%f,%d');
dataPts=[lons,lats];
dataPts = normalize(dataPts);
locations = dataPts;
ND=size(dataPts,1);
dim = size(dataPts,2);
%% project to grids
unit = 0.001; %projection size for normalized data
minx = min(locations,[],1);
maxx = max(locations,[],1);
bound=maxx-minx;
bin = zeros(size(locations));
nbins=uint16(bound./unit)';
for i = 1:dim
    edges = linspace(minx(i), maxx(i), nbins(i)+1);
    edges = [-Inf edges(2:end-1) Inf];
    [~,bin(:,i)] = histc(locations(:,i),edges);
end
%% accumulate arrays 
% H = accumarray(bin,1,nbins');
density_map = containers.Map;
key_format = ['%0',num2str(numel(num2str(unit))-2),'d '];
index=zeros(size(dataPts,1),1);
for i=1: size(dataPts,1)
    id = num2str(bin(i,:),key_format);
    if ~density_map.isKey(id)
        density_map(id) = 1;
        index(i,1)=1;
    else
        density_map(id) = density_map(id)+1;
    end
end
%% get sorted densities
sorted_keys = keys(density_map);
sorted_den = zeros(size(sorted_keys,2),dim+1);
for i=1: size(sorted_keys,2)
    id = sorted_keys{i};
    sorted_den(i,1:dim) = str2num(id);
    sorted_den(i,dim+1) = density_map(id);
end
board=max(sorted_den(:,1:dim));  
R=min(board);
DC=integral(sorted_den,ND,0.021);%Compound 0.21
dc = DC*unit;
%% find representation points 
index=logical(index);
data_rep=locations(index,:);
loc_rep=bin(index,:);
ND_rep=size(loc_rep,1);
den_rep=zeros(1,ND_rep);
for i=1:ND_rep
    id = num2str(loc_rep(i,:),key_format);
    den_rep(i)=density_map(id);
end
%% calculate densities with permutohedral algorithm
features = data_rep'; %loc_rep;
filter_in = den_rep;
sigma_vec = dc;  %dc*2 dc*2 dc dc/2 dc/2 dc dc/2.0 dc*1.93 dc
sigma_dim = dim;
filter_out = mex_permutohedral(features, filter_in, sigma_vec, sigma_dim);
rho = filter_out; 
%% splits the whole area into sub-regions
delta = zeros(1,ND_rep);
nneigh = zeros(1,ND_rep);
[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;
cnt=1;
for ii=2:ND_rep
    delta(ordrho(ii))=norm(max(bin));
    for jj=1:ii-1
        tmp_dist = norm(loc_rep(ordrho(ii),:)-loc_rep(ordrho(jj),:));       
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
%% 选择峰值
disp('Select a rectangle enclosing cluster centers')
scrsz = get(0,'ScreenSize');
figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);
ind = 1:ND_rep;
gamma =rho.*delta;
[gamma_sorted,gamma_indx]=sort(gamma,'descend'); %指标排序

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
%% plot the result
cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   subplot(2,1,1)
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
subplot(2,1,2)
Y1 = data_rep;
plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Improved Cluster result:Trips_Sync','FontSize',12.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND_rep
 A(i,1)=0.;
 A(i,2)=0.;
end
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
tchar=['*','^','.','d','p','h','x','v','+','<','x','>'];
tag=cell(1,NCLUST+1);
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
  plot(A(1:nn,1),A(1:nn,2),tchar(ceil(i/3)),'MarkerSize',5,'MarkerFaceColor',cmap(mod(i,3)+1,:),'MarkerEdgeColor',cmap(mod(i,3)+1,:));
  tag(i)={strcat(num2str(i),'-prototy')};
end
plot(data_rep(icl(1:NCLUST),1),data_rep(icl(1:NCLUST),2),'s','MarkerSize',10,'MarkerFaceColor',[1,1,0],'MarkerEdgeColor',[1,0,1]);
tag(end)={'-center'};
legend(tag);
%% plot the ground truth
 figure,hold on;
 title ('Ground truth','FontSize',12.0)
 for i = 1:size(id_cl,1)
     if(id_cl(i)<1)
         plot(dataPts(i,1),dataPts(i,2),tchar(ceil(id_cl/3)),'MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
     end
     %ic=int8((id_cl(i)*64.)/(max(id_cl(:))*1.));
     plot(dataPts(i,1),dataPts(i,2),tchar(ceil(id_cl(i)/3)),'MarkerSize',4,'MarkerFaceColor',cmap(mod(id_cl(i),3)+1,:),'MarkerEdgeColor',cmap(mod(id_cl(i),3)+1,:));
 end
 hold off;
%% project the represent points back
halo_all = zeros(1,ND); % labels for all
cl_all =  zeros(1,ND);
icl_all = zeros(1,ND);
for i = 1:ND_rep
    index_last=ones(ND,1);
    for j=1:dim
        index_last=index_last&(bin(:,j)==loc_rep(i,j));
    end
    halo_all(:,index_last) = halo(i);
    cl_all(:,index_last) = cl(i);
    index_last=ones(ND,1);
    if icl(i)
        for m=1:dim
        index_last=index_last&(bin(:,m)==loc_rep(icl(i),m));  
        end
        indx=find(index_last==1);
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
   fprintf(faa, '%i %i %i %i %i\n',i,icl_all(i),cl_all(i),halo_all(i),id_cl(i));
end
fclose(faa);
toc