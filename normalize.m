function patches_norm = normalize(patches)

% patches = patches(:,1);
% Squash data to [0.1, 0.9] since we use sigmoid as the activation
% function in the output layer

% Remove DC (mean of images). 
patches = bsxfun(@minus, patches, mean(patches));
% B1 = padarray(max(patches),[size(patches,1)-1,0],'replicate','post');
% B2 = padarray(min(patches),[size(patches,1)-1,0],'replicate','post');
% patches_norm = (patches - B2)./ (B1 - B2);

% Truncate to +/-3 standard deviations and scale to -1 to 1
% pstd = 3 * std(patches);
% B = padarray(pstd,[size(patches,1)-1,0],'replicate','post');
% % patches = max(min(patches, pstd), -pstd) ./ pstd;
%  patches = max(min(patches, B),-B)./B;
% % Rescale from [-1,1] to [0.1,0.9]
% patches_norm = (patches + 1) * 0.4 + 0.1;
 pstd = 3 * std(patches(:));
 patches = max(min(patches, pstd), -pstd) ./ pstd;
 patches_norm = (patches + 1) * 0.4 + 0.1;
 
end
