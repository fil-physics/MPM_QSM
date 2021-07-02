function [ im ] = changeImageSize( im,new_size )
% [ im ] = changeImageSize( im,new_size )
% Emma Dixon's script
% im - image matrix
% new_size - can be smaller or larger matrix size
if mod(size(im,1),2) ~=0
    even_im = zeros(size(im,1)+1,size(im,2),size(im,3));
    even_im = im(2:end,:,:);
    im = even_im;
end

if mod(size(im,2),2) ~=0
    even_im = zeros(size(im,1),size(im,2)+1,size(im,3));
    even_im = im(:,2:end,:);
    im = even_im;
end

if mod(size(im,3),2) ~=0
    even_im = zeros(size(im,1),size(im,2),size(im,3)+1);
    even_im = im(:,:,2:end);
    im = even_im;
end


sizeDiff = new_size - size(im);
sizeDiff = sizeDiff + mod(sizeDiff,2) ;

for dim = 1:3
    if sizeDiff(dim)<0
        if dim==1
            im = im(-1*sizeDiff(1)/2 +1 : -1*sizeDiff(1)/2 + new_size(1),:,:);
        elseif dim==2
            im = im(:,-1*sizeDiff(2)/2 +1 : -1*sizeDiff(2)/2 + new_size(2),:);
        elseif dim ==3
            im = im(:,:,-1*sizeDiff(3)/2 +1 : -1*sizeDiff(3)/2 + new_size(3));
        end
    elseif sizeDiff(dim)>0
        if dim==1
            im = padarray(im, [sizeDiff(dim)/2 0 0]);
        elseif dim==2
            im = padarray(im, [0 sizeDiff(dim)/2 0]);
        elseif dim ==3
            im = padarray(im, [0 0 sizeDiff(dim)/2]);
        end
        
    elseif sizeDiff(dim) == 0
        im = im;
    end
end
end