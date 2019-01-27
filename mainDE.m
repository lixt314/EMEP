clear all
clc

%addpath('C:\Users\WYF\Desktop\SIMLR-SIMLR\MATLAB\src')

rand('state', 0);   

for F=0.4:0.1:0.4
    for CR=0.1:0.1:0.1
RRR=[];

problemSet = [1 : 6];
for problemIndex = 1:6

    problem = problemSet(problemIndex)
   switch problem
        case 1
            Data= importdata('Data_Buettner.mat') 
            X=Data.in_X;         
        case 2
           Data= importdata('Data_Deng.mat') 
           X=Data.in_X;    
       case 3
           Data= importdata('Data_Ginhoux.mat') 
           X=Data.in_X;    
       case 4
            Data= importdata('Data_Pollen.mat')
            X=Data.in_X;          
       case 5
           Data= importdata('Data_Ting.mat') 
           X=Data.in_X;

       case 6
           Data= importdata('Data_Treutlin.mat') 
           X=Data.in_X;
   
   end
X1=normalizeData(X);
K = length(unique(Data.true_labs)) % the number of clusters in the final clustering (using in consensus functions)
truelabels = Data.true_labs; 

vPredictClass=[];
     for i = 2: 20
        option.algorithm='nmfrule';
        [A,Y]=nmf(X',i,option); 
        Cl = kmeans(Y',K);
        while length(unique(Cl)) ~= K;
        Cl = kmeans(Y',K);
        end
        vPredictClass = [vPredictClass Cl];
     end
 
    D=size(vPredictClass,2);
    M=3;
    p1 = [99 13  7  5  4  3  3  2  3];
    p2 = [ 0  0  0  0  1  2  2  2  2];
    p1 = p1(M-1);
    p2 = p2(M-1);
    [N,W] = F_weight(p1,p2,M);

    W(W==0) = 0.000001;
    T = 4;
    B = zeros(N);
    for i = 1 : N
        for j = i : N
            B(i,j) = norm(W(i,:)-W(j,:));
            B(j,i) = B(i,j);
        end
    end
    [~,B] = sort(B,2);
    B = B(:,1:T);
    MaxFES = 1000;
    Population = rand(N,D);
    Population1 = Population>0.5;
    Dim=D+2;
    Para=[ randi([1,3],1,N)' randi([1,3],1,N)'];
  
    for i=1:N
     Index=find(Population1(i,:)==1);
     if size(Index,2)==0
               Index=1; 
     end
      [FunctionValue(i,1:2), POP_label(i,:)] =objectivefunction(vPredictClass(:,Index),Para(i,:),K,truelabels,X1);
      FunctionValue(i,3) =sum(Population1(i,:));
    end


    Z = min(FunctionValue);
    Coding='Binary';
    Boundary=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    Generations=10;
    for Gene = 1 : Generations
            %Fmax = max(FunctionValue);
            %Fmin = Z;
            %FunctionValue = (FunctionValue-repmat(Fmin,N,1))./repmat(Fmax-Fmin,N,1);
        for i = 1 : N
           P = 1:N;

            k = randperm(length(P));
      
            Offspring = F_generator(Population(P(k(1)),:),Population(P(k(2)),:),Population(P(k(3)),:),Boundary, F, CR);
            Offspring1=Offspring>0.5;
            Index=find(Offspring1==1);
              if size(Index,2)==0
               Index=1; 
            end
            [OffFunValue(1:2),Off_label] =objectivefunction(vPredictClass(:,Index),Para(i,:),K,truelabels,X1);
            OffFunValue(3) =sum(Offspring1);
            Z = min(Z,OffFunValue);
            %OffFunValue = (OffFunValue-Fmin)./(Fmax-Fmin);
            for j = 1 : T
                    g_old = max(abs(FunctionValue(B(i,j),:)-Z).*W(B(i,j),:));
                    g_new = max(abs(OffFunValue-Z).*W(B(i,j),:));
                if g_new < g_old
                    %更新当前向量的个体
                    Population(B(i,j),:) = Offspring;
                    FunctionValue(B(i,j),:) = OffFunValue;
                    POP_label(B(i,j),:)=Off_label;
                else
                    Para(B(i,j),:)=[ randi([1,3],1,1)' randi([1,3],1,1)'];
                end
            end

            %反归一化
            %FunctionValue = FunctionValue.*repmat(Fmax-Fmin,N,1)+repmat(Fmin,N,1);

        end
        
        
    end
    for i=1:1:N
        NNMI(i,:) = Cal_NMI(POP_label(i,:),truelabels);
        [ARAR(i,:),RI,MI,HI]=RandIndex(POP_label(i,:),truelabels);
    end
    Results(problemIndex,:)=max(NNMI)
end

end
    
end
