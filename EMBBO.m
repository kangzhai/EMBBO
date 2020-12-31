function [Alpha_score,Cbest,Alpha_pos]=EMBBO(f,a,b,d,M,N)%��Ч���ںϵ��������ѧ�Ż��㷨
% Efficient and Merged Biogeography based optimization(EMBBO)
%
lambdaLower = 0.0; % lower bound for immigration probabilty per gene
lambdaUpper = 1; % upper bound for immigration probabilty per gene
Ib = 1; % max immigration rate for each island
%Eb = 1; % max emigration rate, for each island
P = N; % max species count, for each island
dn=d+1;
%��ʼ��
X=rand(N,d+1);
aa=repmat(a,N,1);
bb=repmat(b,N,1);
X(:,1:d)=aa+(bb-aa).*X(:,1:d);
X(:,d+1)=feval(f, X(:,1:d));
X = PopSort(X);

SpeciesCount=zeros(1);
lambda=zeros(1);
%mu=zeros(1);%��Ϊ���ð���ѧϰѡ�񷽰������Բ���Ҫ�˲���
for j = 1 : N
    if X(j,d+1) < inf
        SpeciesCount(j) = P - j;
    else
        SpeciesCount(j) = 0;
    end
    lambda(j) = Ib * (1 - SpeciesCount(j) / P);
    %mu(j) = Eb * SpeciesCount(j) / P;
end
lambdaMin = min(lambda);
lambdaMax = max(lambda);
lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda - lambdaMin) / (lambdaMax - lambdaMin);

Cbest=zeros(1,M);
i=0;
while i<M
    ps=1-(1-0.65)*i/M;
    p=i/M;
    Island=X(:,1:d);
    Num=ceil(N*rand);
    for k = 1 : N
        if k==Num
            Island(k,:)=a+(b-X(1,1:d));%����ѧϰ
        else
            alfa=rand;
            %�����ǵ�ά������ȫά�����ں�
            if rand<=ps
                j1=ceil(rand*d);
            else
                j1=1:d;
            end
            for j = j1
                if rand < lambdaScale(k)
                    % Pick a habitat from which to obtain a feature
                    SelectIndex=ceil((k-1)*rand);%����ѧϰѡ��
                    if rand<0.5
                        if rand<p
                            Island(k,j) =X(SelectIndex,j);
                        else
                            cnum=ceil(rand*d);Island(k,j) =X(1,cnum);%��������
                        end
                    else
                        Island(k,j) =X(SelectIndex,j)+(rand-0.5)*(X(SelectIndex,j)-X(k,j));%����ʽ��������
                    end
                else
                    rnum=ceil(N*rand);
                    while rnum==k,rnum=ceil(N*rand);end
                    rnum1=ceil(N*rand);
                    while rnum1==k ||rnum1==rnum,rnum1=ceil(N*rand);end
                    if rand<0.5
                        rnum2=ceil(N*rand);
                        while rnum2==k ||rnum2==rnum ||rnum2==rnum1,rnum2=ceil(N*rand);end
                        rnum3=ceil(N*rand);
                        while rnum3==k ||rnum3==rnum ||rnum3==rnum1 || rnum3==rnum2,rnum3=ceil(N*rand);end
                        Island(k,j) = X(k,j)+alfa*(X(rnum2,j)-X(rnum3,j)+X(rnum1,j)-X(rnum,j));%�����ֱ���1
                    else
                        Island(k,j) = X(k,j)+alfa*(X(1,j)-X(k,j)+X(rnum1,j)-X(rnum,j));%�����ֱ���2
                    end
                end
            end
        end
        cc=Island(k,:)<a;Island(k,cc)=a(cc);cc=Island(k,:)>b;Island(k,cc)=b(cc);%Խ�紦��
    end
    fit=feval(f,Island);
    for j=1:N%̰��ѡ��
        value=fit(j);
        if X(j,dn)>value
            X(j,dn)=value;
            X(j,1:d)=Island(j,:);
        end
    end
    X = PopSort(X);%��Ⱥ����
    i=i+1;
    Cbest(i)=X(1,dn);
end
Alpha_score=X(1,dn);%����ֵ
Alpha_pos=X(1,1:d);%���Ž�
