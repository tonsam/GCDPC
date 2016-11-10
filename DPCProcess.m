function [cl,halo]=DPCProcess(loc,den,bin,H,mode)
tic
    close(figure(1));
%     addpath('./permutohedral-double/x64/Debug');
%     addpath('./dc_evaluate/x64/Debug');
    ND=size(loc,2);
    fprintf('number of the data:%d\n',ND);
    Percent = 0.02; % 2% the super-parameter
    fprintf('average percentage of neighbours (hard coded): %5.6f\n', Percent);
    DC = mex_dc_evaluate(H,ND,Percent);
    unit = 0.001; %projection size for normalized data
    dc = DC*unit;
    features =loc; %loc_rep;
    filter_in = den;
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
    delta = zeros(1,ND);
    nneigh = zeros(1,ND);
   [~,ordrho]=sort(rho,'descend');
    delta(ordrho(1))=-1.;
    nneigh(ordrho(1))=0;
    cnt=1;
    for ii=2:ND
       delta(ordrho(ii))= grid_dist(end,end);
       for jj=1:ii-1
          pos_diff = abs(bin(:,ordrho(ii)) - bin(:,ordrho(jj)))+1;
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
%% Ñ¡Ôñ·åÖµ
   disp('Select a rectangle enclosing cluster centers')
   scrsz = get(0,'ScreenSize');
   figure('Position',[6 100 scrsz(3)/2. scrsz(4)/1.3]);
   plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
   title ('Decision Graph','FontSize',15.0)
   xlabel ('\rho')
   ylabel ('\delta')
   toc
   rect = getrect(1);

   rhomin=rect(1);
   deltamin=rect(2);
   NCLUST=0; % number of cluster centers
   cl = zeros(1,ND)-1;
   icl = zeros(1,ND); %inverse id of cluster result
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
%% assignation
   if nargin==5&&strcmp(mode,'remove_halo');
        halo = cl;
        bord_rho = zeros(1,NCLUST);
        if (NCLUST>1)
            for i=1:ND-1
                for j=i+1:ND
                   pos_diff = abs(bin(:,i) - bin(:,j))+1;
                   tmp_dist = grid_dist(pos_diff(1),pos_diff(2));
                   if ((cl(i)~=cl(j))&&(tmp_dist<=DC))
                        rho_aver=(rho(i)+rho(j))/2.;
                       if (rho_aver>bord_rho(cl(i)))
                           bord_rho(cl(i))=rho_aver;
                       end
                       if (rho_aver>bord_rho(cl(j)))
                          bord_rho(cl(j))=rho_aver;
                       end
                   end
                end
            end
            for i=1:ND
               if (rho(i)<bord_rho(cl(i)))
                   halo(i)=0;
               end
            end
        end
 
% %print the statistical result
        for i=1:NCLUST
            nc=0;
            nh=0;
            for j=1:ND
               if (cl(j)==i) 
                  nc=nc+1;
               end
               if (halo(j)==i) 
                  nh=nh+1;
               end
            end
            fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
       end
 %% plot cluster result
      cmap=colormap;
      for i=1:NCLUST
         ic=int8((i*64.)/(NCLUST*1.));
         hold on
         plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
      end
      figure(2);
      hold on
      Y1 = loc';
      plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
      title ('Improved Cluster result:Trips_Sync','FontSize',12.0)
      xlabel ('X')
      ylabel ('Y')
      A=zeros(ND,2);
      for i=0:NCLUST
         nn=0;
         ic=randi([1,size(cmap,1)],1);
         for j=1:ND
            if (halo(j)==i)
               nn=nn+1;
               A(nn,1)=Y1(j,1);
               A(nn,2)=Y1(j,2);
            end
         end
         hold on
         if(i>0)
            plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
         else
            plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0]);
         end
      end
      plot(Y1(icl(1:NCLUST),1),Y1(icl(1:NCLUST),2),'s','MarkerSize',8,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,1]);
      hold off
   else 
       cmap=colormap;
       for i=1:NCLUST
          ic=int8((i*64.)/(NCLUST*1.));
          hold on
          plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
       end
       figure(2);
       hold on
       Y1 = loc';
       plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
       title ('Improved Cluster result:Trips_Sync','FontSize',12.0)
       xlabel ('X')
       ylabel ('Y')
       A=zeros(ND,2);
       for i=1:NCLUST
          nn=0;
          ic=randi([1,size(cmap,1)],1);
          for j=1:ND
             if (cl(j)==i)
              nn=nn+1;
              A(nn,1)=Y1(j,1);
              A(nn,2)=Y1(j,2);
             end
          end
          hold on
          plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
       end
       plot(Y1(icl(1:NCLUST),1),Y1(icl(1:NCLUST),2),'s','MarkerSize',8,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,1]);
       hold off
    end
 toc
end  