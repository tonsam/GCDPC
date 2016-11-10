function dc=integral(dens,ND,percent) 
  %clear;
  %load('dens.mat')
  %dens=sorted_den;
   %percent=0.02;
   N=size(dens,1);
   dims=size(dens,2)-1;
   coordinate=dens(:,1:dims);
   Threshold_num=round(ND*(ND-1)*percent);
   candidate=zeros(2^dims,dims);
   for i=1:dims
       atom=[ones(2^(dims-i),1);-ones(2^(dims-i),1)];
       candidate(:,i)=repmat(atom,2^(i-1),1);
   end
   record=zeros(2,2^dims);
   board=max(dens(:,1:dims));  
   R=min(board);
%    for r=1:R
%       tmp_num=0;
%       for i=1:N
%           cube_num=0;
%           for j=1:2^dims
%               cor_alt=max(0,coordinate(i,:)+r*candidate(j,:));
%               record(2,j)=(-1)^(sum(candidate(j,:)==-1));
%               %cor_alt=min(board,coordinate(i,:)+r*candidate(j,:));
%               Flags_last=ones(N,1);
%               for m=1:dims
%                   Flags=coordinate(:,m)<=cor_alt(1,m);
%                   Flags_last=Flags_last&Flags;
%               end
%               record(1,j)=sum(dens(Flags_last,end));
%               cube_num=record(1,j)*record(2,j)+cube_num;
%           end
%           cube_num=(cube_num-1)*dens(i,end);
%           tmp_num=tmp_num+cube_num;
%       end

   for r=1:R
      tmp_num=0;
      for i=1:N
          cube_num=0;
          for j=1:2^dims
              cor_alt=max(0,coordinate(i,:)+r*candidate(j,:));
              record(2,j)=(-1)^(sum(candidate(j,:)==-1));
              temp_den = dens;
              for m=1:dims
                  Flags = temp_den(:,m)<=cor_alt(1,m);
                  temp_den = temp_den(Flags,:);
              end
              record(1,j)=sum(temp_den(:,end));
              cube_num=record(1,j)*record(2,j)+cube_num;
          end
          cube_num=(cube_num-1)*dens(i,end);
          tmp_num=tmp_num+cube_num;
      end
      if(tmp_num>Threshold_num)
             break;
       end
   end
   dc=r;
end