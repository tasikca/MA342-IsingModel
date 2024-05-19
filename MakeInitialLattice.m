%%
clear all
close all hidden
%% Square Lattice
n = 30;

Square_Lattice = -1 + 2*randi([0 1], n);
[rows, cols] = size(Square_Lattice);

%Sqaure
count = 1;
for i = 1:1:n
    for j = 1:1:n
        if i == n
            A(count+1) = Square_Lattice(1,j)*Square_Lattice(i,j);
        else
            A(count+1) = Square_Lattice(i+1,j)*Square_Lattice(i,j);
        end
        if j == n
            A(count) = Square_Lattice(i,1)*Square_Lattice(i,j);
        else
            A(count) = Square_Lattice(i,j+1)*Square_Lattice(i,j);
        end
        count = count + 2;
    end
end

% %Triangle
% count = 1;
% for i = 1:1:n
%     for j = 1:1:n
%         if i == n && j == n
%             A(count+2) = Square_Lattice(1,1)*Square_Lattice(i,j);
%         elseif i == n
%             A(count+2) = Square_Lattice(1,j+1)*Square_Lattice(i,j);
%         elseif j ==n
%             A(count+2) = Square_Lattice(i+1,1)*Square_Lattice(i,j);
%         else
%             A(count+2) = Square_Lattice(i+1,j+1)*Square_Lattice(i,j);
%         end
%         if i == n
%             A(count+1) = Square_Lattice(1,j)*Square_Lattice(i,j);
%         else
%             A(count+1) = Square_Lattice(i+1,j)*Square_Lattice(i,j);
%         end
%         if j == n
%             A(count) = Square_Lattice(i,1)*Square_Lattice(i,j);
%         else
%             A(count) = Square_Lattice(i,j+1)*Square_Lattice(i,j);
%         end
%         count = count + 3;
%     end
% end

% %Hexagon
% count = 1;
% for i = 1:2:n
%     for j = 1:1:n
%         if mod(j,2) ~= 0 %odd cols
%             if i == n
%                 A(count) = Square_Lattice(1,j+1)*Square_Lattice(i,j);
%             else
%                 A(count) = Square_Lattice(i+1,j+1)*Square_Lattice(i,j);
%             end
%             if j == n
%                 A(count+1) = Square_Lattice(i,1)*Square_Lattice(i,j);
%                 A(count+2) = Square_Lattice(i,2)*Square_Lattice(i,j);
%             elseif j == n-1
%                 A(count+1) = Square_Lattice(i,j+1)*Square_Lattice(i,j);
%                 A(count+2) = Square_Lattice(i,1)*Square_Lattice(i,j);
%             else
%                 A(count+1) = Square_Lattice(i,j+1)*Square_Lattice(i,j);
%                 A(count+2) = Square_Lattice(i,j+2)*Square_Lattice(i,j);
%             end
%         else
%             if i == n
%                 A(count) = Square_Lattice(2,j)*Square_Lattice(i,j);
%             else
%                 A(count) = Square_Lattice(i+1,j)*Square_Lattice(i,j);
%             end
%             if j == n
%                 A(count+1) = Square_Lattice(i,1)*Square_Lattice(i,j);
%                 A(count+2) = Square_Lattice(i,2)*Square_Lattice(i,j);
%             elseif j == n-1
%                 A(count+1) = Square_Lattice(i,j+1)*Square_Lattice(i,j);
%                 A(count+2) = Square_Lattice(i,1)*Square_Lattice(i,j);
%             else
%                 A(count+1) = Square_Lattice(i,j+1)*Square_Lattice(i,j);
%                 A(count+2) = Square_Lattice(i,j+2)*Square_Lattice(i,j);
%             end
%         end
%         
%         count = count + 3;
%     end
% end


J = 1;

E = -J*sum(A);

N = n*n;

P = zeros(n,n);
sq = Square_Lattice;
sig(:,:,1) = sq;
E_N(1) = energy(Square_Lattice,n);
jj = 2;
ii = 2;
for T = 1:0.03:4
    for k = 1:1:500 % run 1 temperature this many times
        for i = 1:rows*cols
            above = mod(i+-2,rows)+1+floor((i-1)/rows)*rows;
            below = mod(i,rows)+1+floor((i-1)/rows)*rows;
            left = mod(i-rows-1,rows*cols)+1;
            right = mod(i+rows-1,rows*cols)+1;
            deltaE = -2*Square_Lattice(i)*(Square_Lattice(above) + Square_Lattice(below) + Square_Lattice(left) + Square_Lattice(right));
            if deltaE<=0
                P(i) = exp(deltaE/T); % T in denominator
            else
                P(i) = 1;
            end
            y = rand(1);
            if P(i) >= y && rand(1)<=0.6
                sq(i) = -1*Square_Lattice(i);
            end
        end
        Square_Lattice = sq;
        sig(:,:,jj) = Square_Lattice;
        E_N(ii) = energy(Square_Lattice,n)/N;
        ii = ii + 1;
    end
    ET_N(jj) = mean(Square_Lattice,"all");
    jj = jj + 1;
    ii = 1;
    E_N(ii) = energy(Square_Lattice,n)/N;
end

% vertLine = linspace(200,200,50);
% y_vec = linspace(-2,0.1,50);
% 
% figure (1)
% plot(E_N,'o')
% hold on 
% plot(vertLine,y_vec,'-')
% xlabel('Simulation Iteration')
% ylabel('Energy Per Lattice Point (E/N)')
% axis([0 500 -2 0.1])
% hold off

% for i=1:size(sig,3)
% Z = sig(:,:,i);
% pcolor(Z)
% drawnow
% end

filename = '30x30sqaure.gif';
for i=1:size(sig,3)
    Z = sig(:,:,i);
    pcolor(Z);
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if i==1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end

T_vec = 1:0.03:4;
T_vec(length(T_vec)+1) = T_vec(length(T_vec)) + 0.03;
x_vec = linspace(2.2,2.2,50);
y_vec = linspace(-1,0.4,50);

figure (1)
plot(T_vec,ET_N,'-o')
hold on
plot(x_vec,y_vec,'r-')
xlabel('Temperature (T)')
ylabel('Sample Energy Per N')
hold off

figure (2)
Z1 = sig(:,:,1);
pcolor(Z1)
drawnow

figure (3)
Z2 = sig(:,:,6);
pcolor(Z2)
drawnow

figure (4)
Z3 = sig(:,:,15);
pcolor(Z3)
drawnow

figure (5)
Z4 = sig(:,:,102);
pcolor(Z4)
drawnow

T_vec = 1:0.003:4;

