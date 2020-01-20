%%%%%%%% ------------------- simulador_viga.m ------------------- %%%%%%%%%
%%%
%%% Autor: Guilherme Nishino Fontebasso
%%% GitHub: guifonte
%%% Data: Setembro 2018
%%% Unicamp - IM349 A - Teoria Técnica e Mecânica dos Sólidos 1
%%% Função: Programa que executa o cálculo simbólico de vigas genéricas
%%%         e retorna os gráficos de esforços internos.
%%%         Permite a adição de forças de forças e momentos concentrados,
%%%         rótulas e forças distribuídas de qualquer ordem de grandeza.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vynum, Mznum, EIzzthetaznum, EIzzvnum] = simulador_viga(table,comprL,VinculacaoE,VinculacaoD,MzE,FyE,MzD,FyD)

%%%%% ---------- Variáveis padrão ---------- %%%%%
syms x EIzzv(x) v(x) EIzzthetaz(x) thetaz(x) E Izz Mz(x) Vy(x);
syms q(x) q1(x) q2(x) q3(x) q4(x);
C = sym('c_',[1 4]);
syms L;
assume(L > 0)



%%%%% ---------- Vinculação esquerda ---------- %%%%% 
syms FE ME;
if VinculacaoE == "Livre"
    VincE = Vinculacao.Livre;
elseif VinculacaoE == "Engastada"
    VincE = Vinculacao.Engastada;
else
    VincE = Vinculacao.Apoiada;
end
%Padrão
FE = 0;
ME = 0;


%%%%% ---------- Vinculação direita ---------- %%%%% 
syms FD MD;
if VinculacaoD == "Livre"
    VincD = Vinculacao.Livre;
elseif VinculacaoD == "Engastada"
    VincD = Vinculacao.Engastada;
else
    VincD = Vinculacao.Apoiada;
end
%Padrão
FD = 0;
MD = 0;


%%%%% ---------- Lista de apoios (posições) ---------- %%%%% 

%Pega números de Apoios
numApoios = 0;
for i = 1:size(table,1)
    res = table{i,1} == "Apoio";
    if res == 1
        numApoios = numApoios + 1;
    end
end

Apoios = sym(zeros(numApoios,1)); %Posições

%Preenche o vetor de Apoios
countApoios = 1;
for i = 1:size(table,1)
    res = table{i,1} == "Apoio";
    if res == 1
        Apoios(countApoios) = table{i,3}/comprL*L;
        countApoios = countApoios + 1;
    end
end
Rys = sym('R_y',[1,length(Apoios)]);



%%%%% ---------- Lista de rótulas (posições) ---------- %%%%% 

%Pega números de Apoios
numRotulas = 0;
for i = 1:size(table,1)
    res = table{i,1} == "Rótula";
    if res == 1
        numRotulas = numRotulas + 1;
    end
end

Rotulas = sym(zeros(numRotulas,1)); %Posições

%Preenche o vetor de Apoios
countRotulas = 1;
for i = 1:size(table,1)
    res = table{i,1} == "Rótula";
    if res == 1
        Rotulas(countRotulas) = table{i,3}/comprL*L;
        countRotulas = countRotulas + 1;
    end
end
Thetas = sym('theta_z',[1,length(Rotulas)]);



%%%%% ---------- Lista de Forças Concentradas ---------- %%%%% 

%Pega o número de Forças Concentradas
numForcasCon = 0;
for i = 1:size(table,1)
    res = table{i,1} == "Força Concentrada";
    if res == 1
        numForcasCon = numForcasCon + 1;
    end
end
%[Símbolo, Posição, Força]
Flista = sym(zeros(numForcasCon,3));

%Preenche o vetor de Forças Concentradas
countForcasCon = 1;
for i = 1:size(table,1)
    res = table{i,1} == "Força Concentrada";
    if res == 1
        Flista(countForcasCon,2) = table{i,3}/comprL*L;
        Flista(countForcasCon,3) = table{i,2};
        countForcasCon = countForcasCon + 1;
    end
end

%Cria uma lista de Forças Concentradas de valores únicos,
%para denominação de símbolos
FConUnicas = sym(zeros(numForcasCon+2,2));
startInd = 1;
if FyE ~= 0
    FConUnicas(startInd,1) = FyE;
    startInd = startInd + 1;
end

if FyD ~= 0
    if FyD ~= FyE
        FConUnicas(startInd,1) = FyD;
        startInd = startInd + 1; 
    end
end
for i = 1:numForcasCon
    found = 0;
    for j = 1:size(FConUnicas,1)
        if FConUnicas(j,1) == Flista(i,3)
           found = 1; 
        break
        end
    end
    if found == 0
        FConUnicas(startInd,1) = Flista(i,3);
        startInd = startInd + 1;
    end
end

Fs = sym('F_',[1,startInd-1]);

if length(Fs) == 1
    syms F;
    Fs(1) = F;
end

for i = 1:startInd-1
    FConUnicas(i,2) = Fs(i);
end

for i = 1:startInd - 1
    if FyE == FConUnicas(i,1)
        FE = FConUnicas(i,2);
    end
    if FyD == FConUnicas(i,1)
        FD = FConUnicas(i,2);
    end
end
for i = 1:numForcasCon
    for j = 1:startInd - 1
        if Flista(i,3) == FConUnicas(j,1)
            Flista(i,1) = FConUnicas(j,2);
        end
    end
end
numFConUnicas = startInd - 1;

%%%%% ---------- Lista de Forças distribuidas ---------- %%%%%
numForcasDist = 0;
for i = 1:size(table,1)
    res = table{i,1} == "Força Distribuída";
    if res == 1
        numForcasDist = numForcasDist + 1;
    end
end
%[Símbolo, Começo, Fim, Força, Ordem]
Fdist = sym(zeros(numForcasDist,5));

countForcasDist = 1;
for i = 1:size(table,1)
    res = table{i,1} == "Força Distribuída";
    if res == 1
        Fdist(countForcasDist,2) = table{i,3}/comprL*L;
        Fdist(countForcasDist,3) = table{i,4}/comprL*L;
        Fdist(countForcasDist,4) = table{i,2};
        Fdist(countForcasDist,5) = table{i,5};
        countForcasDist = countForcasDist + 1;
    end
end

FDistUnicas = sym(zeros(numForcasDist,2));
startInd = 1;

for i = 1:numForcasDist
    found = 0;
    for j = 1:size(FDistUnicas,1)
        if FDistUnicas(j,1) == Fdist(i,4)
           found = 1; 
        break
        end
    end
    if found == 0
        FDistUnicas(startInd,1) = Fdist(i,4);
        startInd = startInd + 1;
    end
end
qs = sym('q_',[1,startInd-1]);
if length(qs) == 1
    syms q_0
    qs(1) = q_0;
end

for i = 1:startInd-1
    FDistUnicas(i,2) = qs(i);
end

for i = 1:numForcasDist
    for j = 1:startInd - 1
        if Fdist(i,4) == FDistUnicas(j,1)
            Fdist(i,1) = FDistUnicas(j,2);
        end
    end
end
numFDistUnicas = startInd - 1;

%%%%% ---------- Lista de Momentos Concentrados ---------- %%%%% 

%Pega o número de Momentos Concentrados
numMomCon = 0;
for i = 1:size(table,1)
    res = table{i,1} == "Momento Concentrado";
    if res == 1
        numMomCon = numMomCon + 1;
    end
end
%[Símbolo, Posição, Valor]
Mlista = sym(zeros(numMomCon,3));

%Preenche o vetor de Forças Concentradas
countMomCon = 1;
for i = 1:size(table,1)
    res = table{i,1} == "Momento Concentrado";
    if res == 1
        Mlista(countMomCon,2) = table{i,3}/comprL*L;
        Mlista(countMomCon,3) = table{i,2};
        countMomCon = countMomCon + 1;
    end
end

%Cria uma lista de Momentos Concentrados de valores únicos,
%para denominação de símbolos
MConUnicas = sym(zeros(numMomCon+2,2));
startInd = 1;
if MzE ~= 0
    MConUnicas(startInd,1) = MzE;
    startInd = startInd + 1;
end

if MzD ~= 0
    if MzD ~= MzE
        MConUnicas(startInd,1) = MzD;
        startInd = startInd + 1; 
    end
end
for i = 1:numMomCon
    found = 0;
    for j = 1:size(MConUnicas,1)
        if MConUnicas(j,1) == Mlista(i,3)
           found = 1; 
        break
        end
    end
    if found == 0
        MConUnicas(startInd,1) = Mlista(i,3);
        startInd = startInd + 1;
    end
end

Ms = sym('M_',[1,startInd-1]);

if length(Ms) == 1
    syms M;
    Ms(1) = M;
end

for i = 1:startInd-1
    MConUnicas(i,2) = Ms(i);
end

for i = 1:startInd - 1
    if MzE == MConUnicas(i,1)
        ME = MConUnicas(i,2);
    end
    if MzD == MConUnicas(i,1)
        MD = MConUnicas(i,2);
    end
end
for i = 1:numMomCon
    for j = 1:startInd - 1
        if Mlista(i,3) == MConUnicas(j,1)
            Mlista(i,1) = MConUnicas(j,2);
        end
    end
end
numMConUnicas = startInd - 1;



%%%%% ---------- Equação do carregamento ---------- %%%%%

q1(x) = 0;
for i = 1:length(Apoios)
    q1(x) = piecewise(x<=Apoios(i),q1 + 0, x>Apoios(i), q1 + Rys(i));
end
for i = 1:size(Flista,1)
    q1(x) = piecewise(x<=Flista(i,2),q1 + 0, x>Flista(i,2), q1 + Flista(i,1));%[Símbolo, Posição, Força]
end
for i = 1:size(Fdist,1)
    q1(x) = piecewise(x<=Fdist(i,2),q1 + 0, x>Fdist(i,2), q1 + (1/factorial(Fdist(i,5)+1))*Fdist(i,1)*(x-Fdist(i,2))^(Fdist(i,5)+1));%[Símbolo, Começo, Fim, Força, Ordem]
    if simplify(Fdist(i,3) < L)
        q1(x) = piecewise(x<=Fdist(i,3),q1 + 0, x>Fdist(i,3), q1 - (1/factorial(Fdist(i,5)+1))*Fdist(i,1)*(x-Fdist(i,3))^(Fdist(i,5)+1));
    end
end

q2(x) = 0;
for i = 1:length(Apoios)
    q2(x) = piecewise(x<=Apoios(i),q2 + 0, x>Apoios(i), q2 + Rys(i)*(x-Apoios(i)));
end
for i = 1:size(Flista,1)
    q2(x) = piecewise(x<=Flista(i,2),q2 + 0, x>Flista(i,2), q2 + Flista(i,1)*(x-Flista(i,2)));
end
for i = 1:size(Mlista,1)
    q2(x) = piecewise(x<=Mlista(i,2),q2 + 0, x>Mlista(i,2), q2 - Mlista(i,1));
end
for i = 1:size(Fdist,1)
    q2(x) = piecewise(x<=Fdist(i,2),q2 + 0, x>Fdist(i,2), q2 + (1/factorial(Fdist(i,5)+2))*Fdist(i,1)*(x-Fdist(i,2))^(Fdist(i,5)+2));
    if simplify(Fdist(i,3) < L)
        q2(x) = piecewise(x<=Fdist(i,3),q2 + 0, x>Fdist(i,3), q2 - (1/factorial(Fdist(i,5)+2))*Fdist(i,1)*(x-Fdist(i,3))^(Fdist(i,5)+2));
    end
end

q3(x) = 0;
for i = 1:length(Apoios)
    q3(x) = piecewise(x<=Apoios(i),q3 + 0, x>Apoios(i), q3 + (1/2)*Rys(i)*(x-Apoios(i))^2);
end
for i = 1:length(Rotulas)
    q3(x) = piecewise(x<=Rotulas(i),q3 + 0, x>Rotulas(i), q3 + Thetas(i));
end
for i = 1:size(Flista,1)
    q3(x) = piecewise(x<=Flista(i,2),q3 + 0, x>Flista(i,2), q3 + (1/2)*Flista(i,1)*(x-Flista(i,2))^2);
end
for i = 1:size(Mlista,1)
    q3(x) = piecewise(x<=Mlista(i,2),q3 + 0, x>Mlista(i,2), q3 - Mlista(i,1)*(x-Mlista(i,2)));
end
for i = 1:size(Fdist,1)
    q3(x) = piecewise(x<=Fdist(i,2),q3 + 0, x>Fdist(i,2), q3 + (1/factorial(Fdist(i,5)+3))*Fdist(i,1)*(x-Fdist(i,2))^(Fdist(i,5)+3));
    if simplify(Fdist(i,3) < L)
        q3(x) = piecewise(x<=Fdist(i,3),q3 + 0, x>Fdist(i,3), q3 - (1/factorial(Fdist(i,5)+3))*Fdist(i,1)*(x-Fdist(i,3))^(Fdist(i,5)+3));
    end
end

q4(x) = 0;
for i = 1:length(Apoios)
    q4(x) = piecewise(x<=Apoios(i),q4 + 0, x>Apoios(i), q4 + (1/6)*Rys(i)*(x-Apoios(i))^3);
end
for i = 1:length(Rotulas)
    q4(x) = piecewise(x<=Rotulas(i),q4 + 0, x>Rotulas(i), q4 + Thetas(i)*(x-Rotulas(i)));
end
for i = 1:size(Flista,1)
    q4(x) = piecewise(x<=Flista(i,2),q4 + 0, x>Flista(i,2), q4 + (1/6)*Flista(i,1)*(x-Flista(i,2))^3);
end
for i = 1:size(Mlista,1)
    q4(x) = piecewise(x<=Mlista(i,2),q4 + 0, x>Mlista(i,2), q4 - (1/2)*Mlista(i,1)*(x-Mlista(i,2))^2);
end
for i = 1:size(Fdist,1)
    q4(x) = piecewise(x<=Fdist(i,2),q4 + 0, x>Fdist(i,2), q4 + (1/factorial(Fdist(i,5)+4))*Fdist(i,1)*(x-Fdist(i,2))^(Fdist(i,5)+4));
    if simplify(Fdist(i,3) < L)
        q4(x) = piecewise(x<=Fdist(i,3),q4 + 0, x>Fdist(i,3), q4 - (1/factorial(Fdist(i,5)+4))*Fdist(i,1)*(x-Fdist(i,3))^(Fdist(i,5)+4));
    end
end



%%%%% ---------- Condicoes de contorno ---------- %%%%%

eqCounter = 4 + length(Apoios);
Eq = sym('eq', [1 eqCounter]);
EqVy = [];
EqMz = [];
EqThetaz = [];
Eqvz = [];

%Vinculações
%Extremidade Esquerda
if(VincE == Vinculacao.Apoiada)
    Eqvz = [Eqvz; [0,0]];
    EqMz = [EqMz; [0,ME]];
elseif(VincE == Vinculacao.Engastada)
    Eqvz = [Eqvz; [0,0]];
    EqThetaz = [EqThetaz; [0,0]];
else %Vinculacao.Livre
    EqVy = [EqVy; [0,FE]];
    EqMz = [EqMz; [0,ME]];
end
%Extremidade Direita
if(VincD == Vinculacao.Apoiada)
    Eqvz = [Eqvz; [L,0]];
    EqMz = [EqMz; [L,MD]];
elseif(VincD == Vinculacao.Engastada)
    Eqvz = [Eqvz; [L,0]];
    EqThetaz = [EqThetaz; [L,0]];
else %Vinculacao.Livre
    EqVy = [EqVy; [L,FD]];
    EqMz = [EqMz; [L,MD]];
end

%Apoios
for i = 1:length(Apoios)
    Eqvz = [Eqvz; [Apoios(i),0]];
end

%Rótulas
for i = 1:length(Rotulas)
    EqMz = [EqMz; [Rotulas(i),0]];
end



%%%%% ---------- Integração da equação diferencial ---------- %%%%%

Vy(x) = C(1);
Mz(x) = int(Vy,x) + C(2);
EIzzthetaz(x) = int(Mz,x) + C(3);
EIzzv(x) = int(EIzzthetaz,x) + C(4);
Vy = Vy + q1;
Mz = Mz + q2;
EIzzthetaz = EIzzthetaz + q3;
EIzzv = EIzzv + q4;
%latex(EIzzv)



%%%%% ---------- Constantes de integração ---------- %%%%%

count = 1;
for i = count:size(Eqvz,1)
    Eq(count) = EIzzv(Eqvz(i,1)) == Eqvz(i,2);
    count = count + 1;
end
for i = 1:size(EqThetaz,1)
    Eq(count) = EIzzthetaz(EqThetaz(i,1)) == EqThetaz(i,2);
    count = count + 1;
end
for i = 1:size(EqMz,1)
    Eq(count) = Mz(EqMz(i,1)) == EqMz(i,2);
    count = count + 1;
end
for i = 1:size(EqVy,1)
    Eq(count) = Vy(EqVy(i,1)) == EqVy(i,2);
    count = count + 1;
end

Incog = [C,Rys,Thetas];
[A,B] = equationsToMatrix(Eq,Incog);
X = simplify(linsolve(A,B));



%%%%% ---------- Equações Finais (simbólicas) ---------- %%%%%

Vyf = Vy;
Mzf = Mz;
EIzzthetazf = EIzzthetaz;
EIzzvf = EIzzv;

for i = 1:length(Incog)
   Vyf(x) = subs(Vyf(x),Incog(i),X(i));
   Mzf(x) = subs(Mzf(x),Incog(i),X(i));
   EIzzthetazf(x) = subs(EIzzthetazf(x),Incog(i),X(i));
   EIzzvf(x) = subs(EIzzvf(x),Incog(i),X(i));
end



%%%%% ---------- Substituição para valores numéricos ---------- %%%%%

Vynum = Vyf;
Mznum = Mzf;
EIzzthetaznum = EIzzthetazf;
EIzzvnum = EIzzvf;

Vynum(x) = subs(Vynum(x),L,comprL);
Mznum(x) = subs(Mznum(x),L,comprL);
EIzzthetaznum(x) = subs(EIzzthetaznum(x),L,comprL);
EIzzvnum(x) = subs(EIzzvnum(x),L,comprL);

for i = 1:numFConUnicas
    Vynum(x) = subs(Vynum(x),FConUnicas(i,2),FConUnicas(i,1));
    Mznum(x) = subs(Mznum(x),FConUnicas(i,2),FConUnicas(i,1));
    EIzzthetaznum(x) = subs(EIzzthetaznum(x),FConUnicas(i,2),FConUnicas(i,1));
    EIzzvnum(x) = subs(EIzzvnum(x),FConUnicas(i,2),FConUnicas(i,1));
end

for i = 1:numFDistUnicas
    Vynum(x) = subs(Vynum(x),FDistUnicas(i,2),FDistUnicas(i,1));
    Mznum(x) = subs(Mznum(x),FDistUnicas(i,2),FDistUnicas(i,1));
    EIzzthetaznum(x) = subs(EIzzthetaznum(x),FDistUnicas(i,2),FDistUnicas(i,1));
    EIzzvnum(x) = subs(EIzzvnum(x),FDistUnicas(i,2),FDistUnicas(i,1));
end

for i = 1:numMConUnicas
    Vynum(x) = subs(Vynum(x),MConUnicas(i,2),MConUnicas(i,1));
    Mznum(x) = subs(Mznum(x),MConUnicas(i,2),MConUnicas(i,1));
    EIzzthetaznum(x) = subs(EIzzthetaznum(x),MConUnicas(i,2),MConUnicas(i,1));
    EIzzvnum(x) = subs(EIzzvnum(x),MConUnicas(i,2),MConUnicas(i,1));
end
% Vynum
% Mznum
% EIzzthetaznum
% EIzzvnum
%latex(EIzzv)

%Vy(100)
end



%%%%% ---------- Diagrama de esforços internos e deflexões ---------- %%%%%
% interval = [0,L];
% 
% subplot(4,1,1);
% fplot(Vynum,interval);
% title('V_y(x)')
% subplot(4,1,2);
% fplot(Mznum,interval);
% title('M_z(x)')
% subplot(4,1,3);
% fplot(EIzzthetaznum,interval);
% title('EI_{zz}\theta_z(x)')
% subplot(4,1,4);
% fplot(EIzzvnum,interval);
% title('EI_{zz}v(x)')
