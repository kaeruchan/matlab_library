clc;
clear all;
filename='test.jpg';%input file name
img=imread(filename);%input image
R=img(:,:,1);
G=img(:,:,2);
B=img(:,:,3);
[x,y,z]=size(img);
for i=1:x
    for j=1:y
        if ((R(i,j)>=100)&&(R(i,j)<=255)&&(G(i,j)<110)&&(B(i,j)<110))%red
            R(i,j)=255;%new 
            G(i,j)=70;
            B(i,j)=70;
        else%others
            %Gray = R*0.299 + G*0.587 + B*0.114    
            gray=R(i,j)*0.299+G(i,j)*0.587+B(i,j)*0.114;
            if gray>120  %if gray bigger than
                R(i,j)=255;
                G(i,j)=255;
                B(i,j)=255;
            end
        end
    end
end
disp('Transform completed');
for i=1:x
    disp(['Recoloring ',num2str(i/x*100),'%']); %Progress
    for j=1:y
        res(i,j,1) = R(i,j);
        res(i,j,2) = G(i,j);
        res(i,j,3) = B(i,j);
    end
end
imwrite(res,'stripes2.png'); %save image
disp('image saved');