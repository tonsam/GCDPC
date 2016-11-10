function fun_imp_hierarchy(parameter_name,mode)
tic
close all
%[~,~,cl_all,~]=textread(filename_CLUSTER,'%d %d %d %d');%change path
load(parameter_name);
dataRep=data_rep';
NDrep=size(loc_rep,2);
NCLUSTER=max(cl);
cl_new=0;
locations=0;
Nm=0;
for n=1:NCLUSTER
     ind=find(cl==n);
     loc=dataRep(ind,:);
     bin=loc_rep(:,ind)';
     den=den_rep(ind);
     H_new=zeros(size(H));
     ND=size(loc,1);
     for i=1:ND
        H_new(bin(i,1),bin(i,2))=H(bin(i,1),bin(i,2));
     end
     row=1-all(H_new==0, 1);
     col=1-all(H_new==0, 2);
     row=find(row~=0);
     col=find(col~=0);
     H_new=H_new(col(1):col(end),row(1):row(end));
     bin(:,1)=bin(:,1)-col(1)+1;
     bin(:,2)=bin(:,2)-row(1)+1;
     if nargin==2&&strcmp(mode,'remove_halo')
        [cl_each,halo]=DPCProcess(loc',den,bin',H_new,mode); 
     end
     if nargin==1
         [cl_each]=DPCProcess(loc',den,bin',H_new); 
     end
     NCLUST=max(cl_each);
     if n==1
         locations=loc;
         cl_new=cl_each;
      else
         locations=[locations;loc];
         cl_each=cl_each+Nm;
         cl_new=[cl_new,cl_each];
    end
    if nargin==2&&strcmp(mode,'remove_halo')
        if n==1
          cl_halo=halo;
        else
          indd= halo~=0;
          temp=zeros(1,ND);
          temp(indd)=Nm;
          halo=halo+temp;
          cl_halo=[cl_halo,halo];
        end
    end
    Nm=NCLUST+Nm;
    fprintf('number 0f label:%d\n',Nm);
end
%% save
name1='CLUSTER_ASSIGNATION_0901';
% th=filename_CLUSTER(end);
% name1=strcat(name1,th);
if nargin==2&&strcmp(mode,'remove_halo')
   faa = fopen(name1, 'w');
   for i=1:NDrep
       fprintf(faa, '%i %i %i\n',i,cl_new(i),cl_halo(i));
   end
   fclose(faa);
else
   faa = fopen(name1, 'w');
   for i=1:NDrep
       fprintf(faa, '%i %i \n',i,cl_new(i));
   end
   fclose(faa);
end  
name2='locationsof_0901_present';
fac = fopen(name2, 'w');
for i=1:NDrep
   fprintf(fac, '%f %f\n',locations(i,1),locations(i,2));
end
fclose(fac);
toc
end