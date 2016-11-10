clear
close all
tic
filename='./data/Aggregation.txt'; %Trips_Sync Spiral R15 Pathbased Jain Flame D31 Compound Aggregation;
[lons,lats,id]=textread(filename, '%f,%f,%d');
% lons = (lons-103.8).*10000;
% lats = (lats-1.2).*10000;
locations =[lons,lats]';
ND=size(locations,2);
fprintf('number of input data: %d\n',ND);
% dataPts=[lons(1:100:end),lats(1:100:end)];
dataPts=[lons,lats];
% data normalize
% dataPts = normalizeData(dataPts);
ND=size(dataPts,1);
dist = zeros(ND);
N = ND*(ND-1)/2;
xx = zeros(N,1);
%% similarity matrix calculation?
disp('Calculate distance matrix between each pair of points')
cnt = 1;
for i=1:ND
    for j=i+1:ND
        dist(i,j) = norm(dataPts(i,:)-dataPts(j,:));
%         dist(i,j) = (dataPts(i,1)-dataPts(j,1))^2+(dataPts(i,2)-dataPts(j,2))^2;
        dist(j,i) = dist(i,j);
        xx(cnt,1) = dist(i,j);
        cnt = cnt + 1;
    end
end
%% difine the neighourhood distance?
disp('Define neighbourhood')
percent=2; % the only superparameter
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);
position=round(N*percent/100);
sda=sort(xx); 
dc=sda(position); % top 2% smallest threshold
% dc = 2.5*3/sqrt(2);
fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
maxd=sda(end); % the maximum distance
%% interpolate the density
rho = zeros(1,ND); 
% with global Gaussian kernel
for i=1:ND-1
  for j=i+1:ND
     theta = exp(-1.0*double(dist(i,j)/dc)^2);
     rho(i)=rho(i)+theta;
     rho(j)=rho(j)+theta;
  end
end
%%  find density peaks by search
disp('Calculate distance of different desity peaks')
delta = zeros(1,ND); % distance with neighbour point
nneigh = zeros(1,ND); % id of the neighbour
[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.; % distance for point with the largest density
nneigh(ordrho(1))=0;  % id for point with the largest density

for ii=2:ND % find distance and id by search via densities and distances?
    delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
delta(ordrho(1))=max(delta(:));
disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')
% wirte the decision graph into file
fid = fopen('DECISION_GRAPH0', 'w');
for i=1:ND
   fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
end
fclose(fid);
%% select cluster centers with decision graph?
disp('Select a rectangle enclosing cluster centers')
scrsz = get(0,'ScreenSize');
figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);
% ind = 1:ND;
% gamma =rho.*delta;
% [gamma_sorted,gamma_indx]=sort(gamma,'descend');

subplot(2,1,1)
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')
% get center points manually
subplot(2,1,1)
toc 

rect = getrect(1);
tic %restart the timer

rhomin=rect(1);
deltamin=rect(2);
NCLUST=0; % cluster number
cl = zeros(1,ND)-1;
icl = zeros(1,ND); 
for i=1:ND
    if ( (rho(i)>rhomin) && (delta(i)>deltamin))
        NCLUST=NCLUST+1;
        cl(i)=NCLUST;
        icl(NCLUST)=i;
    end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST); 
%% assignation
disp('Performing assignation')
for i=1:ND
    if (cl(ordrho(i))==-1)
        cl(ordrho(i)) = cl(nneigh(ordrho(i)));
    end
end
%% no halo
halo = cl;
%% halo
% bord_rho = zeros(1,NCLUST);
% if (NCLUST>1)
% %     bord_rho = 0;
%     for i=1:ND-1
%         for j=i+1:ND
%             if ((cl(i)~=cl(j))&&(dist(i,j)<=dc))
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
%     for i=1:ND
%         if (rho(i)<bord_rho(cl(i)))
%             halo(i)=0;
%         end
%     end
% end
% 
% %find statistics of cluster result
% for i=1:NCLUST
%   nc=0;
%   nh=0;
%   for j=1:ND
%     if (cl(j)==i) 
%       nc=nc+1;
%     end
%     if (halo(j)==i) 
%       nh=nh+1;
%     end
%   end
%   fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
% end
%% plot the result?
cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   subplot(2,1,1)
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
subplot(2,1,2)
Y1 = dataPts';
plot(Y1(1,:),Y1(2,:),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND
 A(i,1)=0.;
 A(i,2)=0.;
end
for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND
    if (halo(j)==i)
      nn=nn+1;
      A(nn,1)=Y1(1,j);
      A(nn,2)=Y1(2,j);
    end
  end
  hold on
  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end

%% plot with ground truth
% figure,hold on;
% cmap=colormap;
% for i = 1:size(id,1)
%     ic=int8((id(i)*64.)/(max(id(:))*1.));
%     plot(lons(i),lats(i),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
% end
% hold off;
%% save the result?
faa = fopen('CLUSTER_ASSIGNATION0', 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster centers')
disp('column 3:cluster assignation without halo control')
disp('column 4:cluster assignation with halo control')
for i=1:ND
   fprintf(faa, '%i %i %i %i\n',i,icl(i),cl(i),halo(i));
end
fclose(faa);
toc