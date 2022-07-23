classdef GA_Optimzer < handle
    %   Genetic Algorithm

    properties
        N_Pop = 50;     %  Number of Population  :: Positive Int
        N_Itr = 1000;   %  Number of Iteration   :: Positive Int 

        LowBound  = []; %  LowBound of Variables :: Double (1,N)  N--Dimension of Variables
        HighBound = []; %  HighBound of Variables :: Double (1,N) N--Dimension of Variables

        Mutate_Pe = 0.01; % Mutate Pecentage of Gene :: Double -- from 0 to 1
        N_Pc = 0.9;       % Offspring Percentage of parnet :: Double -- from 0 to 1

        obj_function; % Objective Function :: Handle 
    end

    methods
        %%  
        function obj = GA_Optimzer()
        end
        %%  选择算子  roulette wheel selection  N-1
        function Single_Pop = RWS_Select(obj,Pop)
            AccFit = zeros(1,obj.N_Pop);
            Propertion = zeros(1,obj.N_Pop);
            Delta = zeros(1,obj.N_Pop);
            Randomindex = rand();
            for i = 1 : obj.N_Pop
                AccFit(i) = AccFit(i) + Pop(i).Fit;
            end
            for i = 1 : obj.N_Pop
                Propertion(i) =  AccFit(i) ./ AccFit(obj.N_Pop);
            end
            for i = 1 : obj.N_Pop
                Delta(i) = abs(Randomindex-Propertion(i));
            end
            [~,Single_Pop_index] = min(Delta);
            Single_Pop = Pop(Single_Pop_index);
        end
        %%  算术杂交算子 arithmetical crossover 2-2
        function results = Ari_Cross(obj,GeneA,GeneB)
            lamdaA = rand();
            lamdaB = 1 - lamdaA;
            resulta = lamdaA .* GeneA +lamdaB .* GeneB;
            resultb = lamdaA .* GeneB +lamdaB .* GeneA;
            results = [resulta;resultb];
        end
        %%  仿射杂交算子 affine crossover 2-2
        function results = Aff_Cross(obj,GeneA,GeneB)
            lamdaA = -1 + 2*rand();
            lamdaB = 1 - lamdaA;
            resulta = lamdaA .* GeneA +lamdaB .* GeneB;
            resultb = lamdaA .* GeneB +lamdaB .* GeneA;
            results = [resulta;resultb];
        end
        %%  边界杂交算子 boundary crossover 2-1
        function results = Bou_Cross(obj,GeneA,GeneB)
            vec = length (obj.LowBound);
            results = zeros(1,length(obj.LowBound));
            for i = 1 : vec
                results(i) = sqrt(rand().*(GeneA(i)^2) + (1-rand())*(GeneB(i)^2));
            end
        end
        %%  定向杂交算子 direction-based crossover 2-1
        function results = Dir_Cross(obj,GeneA,GeneB)
            results = rand().*(GeneB-GeneA) + GeneB;
        end
        %%  均匀变异算子 normal mutate 1-1
        function GeneT = Nor_Mutate(obj,GeneT)
            if rand() < obj.Mutate_Pe
                Randomindex = randi([1, length(GeneT)]); %随机选择变异位置
                GeneT(Randomindex) = min(GeneT) + (max(GeneT)-min(GeneT)) * rand(); % 染色体中值变异
            end
        end
        %%  正态变异算子  NDX mutata 1-1
        function GeneT = NDX_Mutate(obj,GeneT)
            sigma = std(GeneT);
            mu = 0;
            if rand() < obj.Mutate_Pe
                for i = 1 : length(GeneT)
                    GeneT(i) =  GeneT(i) + normrnd(mu,sigma);
                end
            end
        end
        %%  非一致变异 non-uniform mutation
        function GeneT = Non_Mutate(obj,GeneT,itr)% itr当前进化代数
            b = randi([2,5]);
            LB = min(GeneT);
            UB = max(GeneT);
            if rand() > 0.5
                if rand() < obj.Mutate_Pe
                    for i = 1 : length(GeneT)
                        GeneT(i) = GeneT(i) + (UB-GeneT(i))*(1-((rand())^((1-itr/obj.N_Itr)^b)));
                    end
                end
            else
                if rand() < obj.Mutate_Pe
                    for i = 1 : length(GeneT)
                        GeneT(i) = GeneT(i) - (GeneT(i)-LB)*(1-((rand())^((1-itr/obj.N_Itr)^b)));
                    end
                end
            end
        end
        %%  优化入口
        function [SolutionX,SolutionXset] = Normal_Optim(obj)
            % Initialize
            Pop = [];
            for i = 1 : obj.N_Pop
                Pop(i).Gene = unifrnd(obj.LowBound,obj.HighBound,1,length(obj.LowBound));
                Pop(i).Fit = obj.obj_function(Pop(i).Gene);
            end
            %
            nC = round(obj.N_Pop * obj.N_Pc / 2) * 2;
            %
            Template.Gene = [ ];
            Template.Fit = [ ];
            SolutionXset = zeros(obj.N_Itr,length(Pop(1).Gene));
            Parent = Pop;
            %
            for itr = 1 : obj.N_Itr
                %
                Offspring = repmat(Template, nC/2, 2);
                for j = 1 : nC/2
                    %
                    objA = obj.RWS_Select(Parent);
                    objB = obj.RWS_Select(Parent);
                    %
                    results = obj.Ari_Cross(objA.Gene,objB.Gene);
                    Offspring(j, 1).Gene = results(1,:);
                    Offspring(j, 2).Gene = results(2,:);
                end
                Offspring = Offspring(:);
                %
                for k = 1 :nC
                    Offspring(k).Gene = obj.Non_Mutate(Offspring(k).Gene,itr);
                    Offspring(k).Fit = obj.obj_function(Offspring(k).Gene);
                end
                %
                Parent = Parent(:);
                newPop = [Parent;Offspring];
                %
                [~, so] = sort([newPop.Fit],'ascend');
                newPop = newPop(so);
                %
                Parent = newPop(1:obj.N_Pop);
                SolutionXset(itr,:) = Parent(1).Gene;
            end
            SolutionX = Parent(1).Gene;
        end
        %% Elitist Selection GA
        function [SolutionX,SolutionXset] = Elitist_Optim(obj)
        end
        %% Father-offspring Combined Selection GA
        function [SolutionX,SolutionXset] = Combined_Optim(obj)
        end
        %% Sharing GA
        function [SolutionX,SolutionXset] = Sharing_Optim(obj)
        end
    end
end



