clc;clear
obj = GA_Optimzer( );
obj.LowBound = -1 * ones(1,30);
obj.HighBound = 1 * ones(1,30);
obj.obj_function = @ackley;
epoch = 1000;
obj.N_Itr = epoch;
tic
[x,xset]= obj.Optim();
toc

yset = zeros(1,epoch);
for j = 1 : epoch
    yset(j) = ackley(xset(:,j));
end
 plot(yset)