clc
close 

it=50;
num=2;
for select= 1:4
  
diff_feature(select,it) % plotta le differenze
pca(select,num,it);
plot_rocc(select)
end


function pca(select,num,it)
[A,name] = select_data(select);
Z = zscore(A);
R = corrcoef(Z);

eigv=eig(R);
eigv =eigv/sum(eigv);
eigv = sort(eigv,'descend')*100;
eigv = eigv(eigv>=0.02);
close
bar(eigv);
ax= gca;
ax.XTick = 1: length(eigv);
xlabel("Dimensione")
ylabel("Percentuale di varianza spiegata")
matlab2tikz([char(name) '_autovalori.tex'],'showInfo', false)
saveas(gcf,[char(name) '_autovalori'],'jpeg')

[U,~]=eigs(R,num);
y = Z* U; %nuove variabili

if num==2
name_ax = {"Dim 1 (" + num2str(eigv(1),3) + "%)","Dim 2 (" + num2str(eigv(2),3) + "%)"};
end

%% cluster grado del tumore a 3
load("Dati","deg","pdl1","patient_name");
C =KMEANS(y,3,it); 
[~,C]=valuta(C,deg);
plot_cluster(y,C,patient_name,name_ax, [name  '_degree'],{'degree 1','degreee 2','degree 3'});
plot_cluster(y,deg,patient_name,name_ax,[name '_degree_real'],{'degree 1','degreee 2','degree 3'});

close all
img=confusionchart(deg,C);
saveas(img,[name '_confusional_degree' num],'jpeg')

%% cluster grado del tumore condensando 2 e 3 
C =KMEANS(y,2,it); 
deg2 = deg;
deg2(deg==3)=2;
[~,C]=valuta(C,deg2);
plot_cluster(y,C,patient_name,name_ax,[name +'_degree_cond'],{'degree 1','degreee >1'});
plot_cluster(y,deg2,patient_name,name_ax,[name '_degree_cond_real'],{'degree 1','degree >1'});

close all
img=confusionchart(deg2,C);
saveas(img,[name '_confusional_degree_cond' num],'jpeg')


%% cluster PDL1 (cluster 2 = positivo)
PDL1 = (pdl1>=1) + 1;
C = KMEANS(y,2,it); 
[~,C]=valuta(C,PDL1);
plot_cluster(y,C,patient_name,name_ax,[name '_PDL'],{'PDL1 -','PDL1 +'})
plot_cluster(y,PDL1,patient_name,name_ax,[name '_PDL_real'],{'PDL1 -','PDL1 +'})

close;
img=confusionchart(PDL1,C);
saveas(img,[name '_confusional_PDL' num],'jpeg')

end

function C =KMEANS(A,k,IT)
% C =KMEANS(A,k,IT) itera IT volte llyoid(A,k)
% C vettore di numeri da 1 a k, se C(i) = j allora l'oggetto i e' assegnato
% al j-esimo cluster 
[C,cost] = llyoid(A,k);
for i = 1: IT
    [Cn,costn]=llyoid(A,k);
    if costn<cost 
        C= Cn;
        cost = costn;
    end
end
end

function [C,cost]=llyoid(A,k)
% [C,cost]=llyoid(A,k) esegue l'algoritmo di llyoid per clusterizzare le
% righe di A in k cluster
N =size(A,1);
C =randi(k,N,1);
Cpred  = zeros(size(C));
D = zeros(N,k);
while ~isequal(C,Cpred)
    z =calcola_centroidi(A,k,C);
    for i=1:k
        B = A-z(i,:);
        D(:,i) = sum(B.^2,2);
    end
    [~,new] = min(D,[],2);
    Cpred = C;
    C = new;
end
cost = sum(min(D,[],2));

end

function z=calcola_centroidi(A,k,C)
% z=calcola_centroidi(A,k,C) restituisce il centroide dei k cluster
% definiti da C
z = zeros(k,size(A,2));
for i = 1:k
    if isempty(A(C==i,:))
        z(i,:)=0;
    else
        z(i,:) = mean(A(C==i,:));
    end
 
end
end


function [err,D] = valuta(C,H)
% confronta il cluster ottenuto C con il cluster desiderato, cercando quale
% permutazione degli indici riduce il numero di errori
% La funzione restituisce il numero di errori minimi e il cluster che
% minimizza gli errori
k = max(C);
perm  = perms(1:k);
v = factorial(k);
err = zeros(v,1);
for j = 1: v 
    D = zeros(size(C));
    for i = 1:k
        D(C==i) = perm(j,i);
    end
   
    err(j)= sum(D~=H);
end
[err,index] = min(err);
for i = 1:k
    D(C==i) = perm(index,i);
end
end


function plot_cluster(y,C,point_name,name_ax,filename,legends)
% plotta i punti colorando ogni cluster con un colore diverso (massimo 6)
close all
col = 'krgbcmy';
k = max(C);
for i =1:k
    if size(y,2)==2
        x1 = y(C==i,1);
        x2= y(C==i,2);
        plot(x1,x2, [col(i) 'o'],'MarkerFaceColor' ,col(i), 'MarkerSize',6);
        a=text(x1,x2 , point_name(C==i),'VerticalAlignment','top','HorizontalAlignment','left');
        set(a,"Color",col(i));
        hold on
    else
        % 
    end
end

set(gca,'yticklabels',[])
set(gca,'xticklabels',[])

grid minor;

xlabel(name_ax(1)) ; ylabel(name_ax(2)); 
plot(gca().XLim,[0 0],'k')
plot([0,0],gca().YLim,'k')
lg = legend(legends);
lg.String = legends;
lg.Location='best';


cleanfigure
matlab2tikz([filename '.tex'],'showInfo', false)
saveas(gcf,filename,'jpeg')
end

function plot_rocc(select)
[A,name]= select_data(select);
load("Dati.mat",'feature_name','pdl1');
N =length(feature_name);
PDL1 = (pdl1>=1) + 1;

for i = 1: N
    x = rocc(A(:,i),PDL1);
    plot(x(:,1),x(:,2),'k','LineWidth',1);
    hold on 
    plot([1 0 ],[0 1],'k','LineWidth',0.01,'LineStyle','-.')
    xlabel('specificita (%)')
    set(gca,'XDir','reverse')
    ylabel('sensibilita (%)')
    title(feature_name(i));
    matlab2tikz([name '_' char(feature_name(i)) '.tex'])
    saveas(gcf,[name '_' char(feature_name(i))],'jpeg')
    hold off
end

end


function x= rocc(y,H)
% positivi quelli con valore maggiore della soglia
% in H i positivi hanno valore 2 
m= min(y);
M =max(y);
N = 1000;
x = zeros(N+1,2);
tax = linspace(m,M,N);
%sensibilita veri positivi/malati -> y 
% specificita veri negativi/ sani -> x
% x i falsi positivi 
tax(end+1) =M+1;
% y i veri positivi 
for i = 1: N+1 
    tasso = tax(i);
    pos_test = y>tasso;
    x(i,1)= sum((pos_test==0) & (H==1))/sum(H==1);
    x(i,2)= sum((pos_test ==1) & (H==2))/sum(H==2);
end
end

function diff_feature(select,it)
[A,name]= select_data(select);
C = KMEANS(A,2,it);
Z= calcola_centroidi(A,2,C);

load("Dati.mat",'feature_name')
val1 = find(sum(Z>50)==2);
if ~isempty(val1)
name1  = feature_name(val1);
plot(Z(1,val1),'ro-')
hold on
plot(Z(2,val1),'kx-')
ax = gca;
ax.XTick = 1:length(val1);
ax.XLim(1)=1;
xticklabels(name1)
ax.XTickLabelRotation=90;
legend('group 1','group 2')
xlabel('variable')
ylabel('value')
matlab2tikz([char(name) '_diff_feature1.tex'])
saveas(gcf,[char(name)  '_diff_feature1'],'jpeg')
end

close all

val2 = 1:length(Z);
if ~isempty(val2)
val2(val1)=[];
name2 = feature_name(val2);
plot(Z(1,val2),'ro-')
hold on
plot(Z(2,val2),'kx-')
ax = gca;
ax.XTick = 1:length(val2);
ax.XLim(1)=1;
xticklabels(name2)
ax.XTickLabelRotation=90;
legend('group 1','group 2')
xlabel('variable')
ylabel('value')
matlab2tikz([char(name)  '_diff_feature2.tex'])
saveas(gcf,[char(name) '_diff_feature2'],'jpeg')
end
end

function [A,name] = select_data(select)

switch select
    case 1 
        name = "Equal";
        load("Dati.mat",'Equal')
        A = Equal;
    case 2
       name = "Llyod";
        load("Dati.mat",'Lloyd')
        A = Lloyd;
    case 3
        name = "Noquant";
        load("Dati.mat",'Noquant')
        A = Noquant;
    case 4
        name = "Uniform";
        load("Dati.mat",'Uniform')
        A = Uniform;
end
name =char(name);
end
