% AUTHOR : ILHAM AZKA RAMADHAN 
% NIM    : 15020131

%{
  TODO
  - check equation in every ldtype
  - display displacement,  
  - wrap display code in a function
  - wrap long code in a seperate function
  - refactor code
%}

fprintf("\nINPUT NOMOR SOAL ANTARA 1-2\n")
filename = input("No Soal: ") + ".txt";
data = dlmread(filename);

% PART I
d = 1;
NJ = data(d);
COORD = zeros(NJ, 3);
for i=1:NJ
    COORD(i, 1) = data(d+1, 1);
    COORD(i, 2) = data(d+1, 2);
    COORD(i, 3) = data(d+1, 3);
    d = d + 1;
end

% PART II
d = d + 1;
NS = data(d);
NCJT = 6;
MSUP = zeros(NS, NCJT+1);
for i=1:NS
    for j=1:NCJT+1
        MSUP(i, j) = data(d+1, j);
    end
    d = d + 1;
end

% PART III
d = d + 1;
NMP = data(d);
EM = zeros(NMP, 2);
for i=1:NMP
    EM(i, 1) = data(d+1, 1);
    EM(i, 2) = data(d+1, 2);
    d = d + 1;
end

% PART IV
d = d + 1;
NCP = data(d);
CP = zeros(NCP, 4);
for i=1:NCP
    CP(i, 1) = data(d+1, 1);
    CP(i, 2) = data(d+1, 2);
    CP(i, 3) = data(d+1, 3);
    CP(i, 4) = data(d+1, 4);
  end
    d = d + 1;

% PART V
d = d + 1;
NM = data(d);
MPRP = zeros(NM, 5);
for i=1:NM
    MPRP(i, 1) = data(d+1, 1);
    MPRP(i, 2) = data(d+1, 2);
    MPRP(i, 3) = data(d+1, 3);
    MPRP(i, 4) = data(d+1, 4);
    MPRP(i, 5) = data(d+1, 5);
    d = d + 1;
end

% PART VIa
d = d + 1;
NJL = data(d);
if NJL > 0
    JP = zeros(NJL, 1);
    PJ = zeros(NJL, NCJT);
    for i=1:NJL
        JP(i) = data(d+1);
        for j=1:NCJT
            PJ(i, j) = data(d+1, j+1);
        end
        d = d + 1;
    end
end

d = d + 1;
% PART VIb
NML = data(d);
if NML > 0
    MP = zeros(NML, 3);
    PM = zeros(NML, 4);
    for i=1:NML
        MP(i, 1) = data(d+1, 1);
        MP(i, 2) = data(d+1, 2);
        MP(i, 3) = data(d+1, 3);
        if MP(i, 2) == 1 || MP(i, 2) == 2 || MP(i, 2) == 5
            PM(i, 1) = data(d+1, 4);
            PM(i, 3) = data(d+1, 6);
        elseif MP(i, 2) == 3 || MP(i, 2) == 6
            PM(i, 1) = data(d+1, 4);
            PM(i, 3) = data(d+1, 6);
            PM(i, 4) = data(d+1, 7);
        elseif MP(i, 2) == 4
            PM(i, 1) = data(d+1, 4);
            PM(i, 2) = data(d+1, 5);
            PM(i, 3) = data(d+1, 6);
            PM(i, 4) = data(d+1, 7);
        end
        d = d + 1;
    end
end
        
% DISPLAYING GENERAL STRUCTURAL INFORMATION
fprintf("GENERAL STRUCTURAL\n");
fprintf("Structure Type: Truss\n");
fprintf("Number of Joints: %d\n", NJ);
fprintf("Number of Members: %d\n", NM);
fprintf("Number of Material Property Sets (E): %d\n", NMP);
fprintf("Number of Cross-Sectional Property Sets: %d\n\n", NCP);

fprintf("\nJOINT COORDINATES\n\n");
NJTJ = zeros(NJ, 1);
for i=1:NJ
    NJTJ(i, 1) = i;
end

NJT_label = {'Joint No' 
             'X Coordinate' 
             'Y Coordinate' 
             'Z Coordinate'};

NJT = table(NJTJ, COORD(:, 1), COORD(:, 2), COORD(:, 3), 'VariableNames', NJT_label);
disp(NJT)

fprintf("\nSUPPORT DATA\n\n");
MSUPT = string(MSUP);
for i=1:NS
    for j=2:NCJT+1
        if MSUPT(i, j) == "1"
            MSUPT(i, j) = "Yes";
        elseif MSUPT(i, j) == "0"
            MSUPT(i, j) = "No";
        end
    end
end

MSUPT_label = {'Joint No' 
               'X Restraint'
               'Y Restraint'
               'Z Restraint'
               'X Rot. Restraint'
               'Y Rot. Restraint'
               'Z Rot. Restraint'};

MSUPTJ = table(double(MSUPT(:, 1)), MSUPT(:, 2), MSUPT(:, 3), MSUPT(:, 4), MSUPT(:, 5), MSUPT(:, 6), MSUPT(:, 7), 'VariableNames', MSUPT_label);
disp((MSUPTJ))

fprintf("\nMATERIAL PROPERTIES\n\n");
EMTJ = zeros(NMP, 1);
for i=1:NMP
    EMTJ(i, 1) = i;
end
EMT_label = {'Material No' 
             'Elasticity Modulus'
             'Shear Modulus'};

EMT = table(EMTJ, EM(:, 1), EM(:, 2), 'VariableNames', EMT_label);
disp(EMT)

fprintf("\nMEMBER DATA\n\n");
NMTJ = zeros(NM, 1);
for i=1:NM
    NMTJ(i, 1) = i;
end
NMT_label = {'No' 
             'Joint Beginning' 
             'Joint End' 
             'No. Material' 
             'No. Prop CrossSectional' 
             'Angel Of Roll'};

NMT = table(NMTJ, MPRP(:, 1), MPRP(:, 2), MPRP(:, 3), MPRP(:, 4), MPRP(:, 5), 'VariableNames', NMT_label);
disp(NMT)

fprintf("\nCROSS-SECTIONAL PROPERTIES\n\n");
CPTJ = zeros(NCP, 1);
for i=1:NCP
    CPTJ(i, 1) = i;
end
CPT_label = {'Property No' 
             'Area'
             'Inertia-Z'
             'Inertia-Y'
             'Torsion Constant'};

CPT = table(CPTJ, CP(:, 1), CP(:, 2), CP(:, 3), CP(:, 4), 'VariableNames', CPT_label);
disp(CPT)

fprintf("\nJOINT LOAD\n\n");
PJTJ = zeros(NJL, 1);
for i=1:NJL
    PJTJ(i, 1) = i;
end
PJT_label = {'No' 
             'X Forces' 
             'Y Forces' 
             'Z Forces' 
             'X Moment' 
             'Y Moment' 
             'Z Moment'};

PJT = table(JP(:, 1), PJ(:, 1), PJ(:, 2), PJ(:, 3), PJ(:, 4), PJ(:, 5), PJ(:, 6), 'VariableNames', PJT_label);
disp(PJT)

fprintf("\nMEMBER LOAD\n\n")
PMT_label = {'Member No' 
             'Load Type' 
             'Axis' 
             'Load Magnitude (W or M)' 
             'Load Intensity (W2)' 
             'Distance L1' 
             'Distance L2'};

PMT = table(MP(:, 1), MP(:, 2), MP(:, 3), PM(:, 1), PM(:, 2), PM(:, 3), PM(:, 4), 'VariableNames', PMT_label);
disp(PMT)
% PART VII
NR = 0;
% numbering reaction
for i=1:NS
    for i1=2:NCJT+1
        if MSUP(i, i1) == 1
            NR = NR + 1;
        end
    end
end
NDOF = NCJT * NJ - NR;
fprintf("NDOF: %d\n\n", NDOF);
fprintf("NR: %d\n\n", NR);

% PART VIII
NSC = zeros(NCJT * NJ, 1);
j = 0;
k = NDOF;
for i=1:NJ
    icount = 0;
    for i1=1:NS
        % if joint numbering == joint support number
        if MSUP(i1, 1) == i 
            icount = 1;
            for i2=1:NCJT
                i3 = (i-1) * NCJT + i2;
                if MSUP(i1, i2+1) == 1
                    k = k + 1;
                    NSC(i3) = k;
                else
                    j = j + 1;
                    NSC(i3) = j;
                end
            end
        end
    end
    if icount == 0
        for i2=1:NCJT
            i3 = (i-1) * NCJT + i2;
            j = j + 1;
            NSC(i3) = j;
        end
    end
end

% PART IX
S  = zeros(NDOF, NDOF);
P  = zeros(NDOF, 1);
GK = zeros(2*NCJT, 2*NCJT);
BK = zeros(2*NCJT, 2*NCJT);
T  = zeros(2*NCJT, 2*NCJT);
FF = zeros(2*NCJT, 1);
QF = zeros(2*NCJT, 1);

for im=1:NM
    JB = MPRP(im, 1); 
    JE = MPRP(im, 2);
    AOR = MPRP(im, 5);

    % MATERIAL PROPERTIES
    I = MPRP(im, 3); E = EM(I, 1); G = EM(I, 2);

    % CROSS SECTIONAL PROPERTIES
    I = MPRP(im, 4); A = CP(I, 1); ZI = CP(I, 2); YI = CP(I, 3); J = CP(I, 4);

    XB = COORD(JB, 1); 
    YB = COORD(JB, 2); 
    ZB = COORD(JB, 3);

    XE = COORD(JE, 1); 
    YE = COORD(JE, 2); 
    ZE = COORD(JE, 3);

    BL = sqrt((XE-XB)^2 + (YE-YB)^2 + (ZE - ZB)^2);

    RXX = (XE-XB)/BL;
    RXY = (YE-YB)/BL;
    RXZ = (ZE-ZB)/BL;

    RYX = (-RXX * RXY * cosd(AOR) - RXZ * sind(AOR)) / (sqrt(RXX^2 + RXZ^2));
    RYY = sqrt(RXX^2 + RXZ^2) * cosd(AOR);
    RYZ = (-RXY * RXZ * cosd(AOR) + RXX * sind(AOR)) / (sqrt(RXX^2 + RXZ^2));

    RZX = (RXX * RXY * sind(AOR) - RXZ * cosd(AOR)) / (sqrt(RXX^2 + RXZ^2));
    RZY = -sqrt(RXX^2 + RXZ^2) * sind(AOR);
    RZZ = (RXY * RXZ * sind(AOR) + RXX * cosd(AOR)) / (sqrt(RXX^2 + RXZ^2));

    [BK] = MSTIFFL(E, G, A, ZI, YI, J, BL, NCJT, BK);
    [T] = MTRANS(AOR, RXX, RXY, RXZ, RYX, RYY, RYZ, RZX, RZY, RZZ, NCJT, T);
    [GK] = MSTIFFG(NCJT, BK, T, GK);
    [S] = STORES(JB, JE, NCJT, NDOF, NSC, GK, S);

    if NML > 0
        QF = zeros(2*NCJT, 1);
        for IML=1:NML
            if im == MP(IML, 1)
                [QF] = MFEFLL(IML, BL, MP, PM, QF);
                break;
            end
        end
        [FF] = MFEFG(NCJT, T, QF, FF);
        [P] = STOREPF(JB, JE, NCJT, NDOF, NSC, FF, P);
    end
end

fprintf("\n\n==STRUCTURE STIFFNESS MATRIX==\n");
format compact
disp(num2str(S))

% PART X
for i=1:NJL
    i1 = JP(i);
    i2 = (i1 - 1) * NCJT;
    for j=1:NCJT
        i2 = i2 + 1;
        N = NSC(i2);
        if N <= NDOF
            P(N) = P(N) + PJ(i, j);
        end
    end
end

PN = zeros(NDOF);
for i=1:NDOF
    PN(i) = i;
end

fprintf("\n\n==P-Pf==\n\n")
PT = table(PN(:, 1), P(:, 1), 'VariableNames', {'GNN' 'Joint_Load'});
disp(PT)

% PART XI
% P = inv(S) * P;
% Manual way
for i=1:NDOF
    Z1 = S(i, i);
    for j=1:NDOF
        S(i, j) = S(i, j)/Z1;
    end
    P(i) = P(i)/Z1;
    for k=1:NDOF
        if k ~= i
            Z = S(k, i);
            for m=1:NDOF
                S(k, m) = S(k, m) - S(i, m) * Z;
            end
            P(k) = P(k) - P(i) * Z;
        end
    end
end

% PART XII
BK = zeros(2*NCJT, 2*NCJT);
T  = zeros(2*NCJT, 2*NCJT);

V  = zeros(2*NCJT, 1);
U  = zeros(2*NCJT, 1);
Q  = zeros(2*NCJT, 1);
F  = zeros(2*NCJT, 1);
QF = zeros(2*NCJT, 1);

R  = zeros(NR, 1);

MAF  = zeros(NM, 1);
NMAF = zeros(NM, 1);

for im=1:NM
    JB  = MPRP(im, 1); 
    JE  = MPRP(im, 2);
    AOR = MPRP(im, 5);

    % MATERIAL PROPERTIES
    I = MPRP(im, 3); E = EM(I, 1); G = EM(I, 2);

    % CROSS SECTIONAL PROPERTIES
    I = MPRP(im, 4); A = CP(I, 1); ZI = CP(I, 2); YI = CP(I, 3); J = CP(I, 4);

    % BEGINNING MEMBER COORDINATE
    XB = COORD(JB, 1); 
    YB = COORD(JB, 2); 
    ZB = COORD(JB, 3);
    % ENDING MEMBER COORDINATE
    XE = COORD(JE, 1); 
    YE = COORD(JE, 2); 
    ZE = COORD(JE, 3);
    % MEMBER LENGTH
    BL = sqrt((XE-XB)^2 + (YE-YB)^2 + (ZE-ZB)^2);

    RXX = (XE-XB)/BL;
    RXY = (YE-YB)/BL;
    RXZ = (ZE-ZB)/BL;

    RYX = (-RXX * RXY * cosd(AOR) - RXZ * sind(AOR)) / (sqrt(RXX^2 + RXZ^2));
    RYY = sqrt(RXX^2 + RXZ^2) * cosd(AOR);
    RYZ = (-RXY * RXZ * cosd(AOR) + RXX * sind(AOR)) / (sqrt(RXX^2 + RXZ^2));

    RZX = (RXX * RXY * sind(AOR) - RXZ * cosd(AOR)) / (sqrt(RXX^2 + RXZ^2));
    RZY = -sqrt(RXX^2 + RXZ^2) * sind(AOR);
    RZZ = (RXY * RXZ * sind(AOR) + RXX * cosd(AOR)) / (sqrt(RXX^2 + RXZ^2));

    % fprintf("\nMEMBER %d\n", im);
    % fprintf("%d\t%d\t%d\n%d\t%d\t%d\n%d\t%d\t%d\n", RXX, RXY, RXZ, RYX, RYY, RYZ, RZX, RZY, RZZ);

    [V]  = MDISPG(JB, JE, NCJT, NDOF, NSC, P, V);
    [T]  = MTRANS(AOR, RXX, RXY, RXZ, RYX, RYY, RYZ, RZX, RZY, RZZ, NCJT, T);
    [U]  = MDISPL(NCJT, V, T, U);
    [BK] = MSTIFFL(E, G, A, ZI, YI, J, BL, NCJT, BK);

    % INIT ALL QF ELEMENT TO ZERO
    QF = zeros(2*NCJT, 1);

    if NML > 0
        for IML=1:NML
            if im == MP(IML, 1);
                [QF] = MFEFLL(IML, BL, MP, PM, QF);
                break;
            end
        end
    end
    [Q] = MFORCEL(NCJT, BK, U, Q, QF);
    [F] = MFORCEG(NCJT, T, Q, F);
    [R] = STORER(JB, JE, NCJT, NDOF, NSC, F, R);
    % fprintf("\n\nLOCAL MEMBER FORCE %d\n\n", im);
    % disp(num2str(Q));
    % fprintf("\n\nGLOBAL MEMBER FORCE %d\n\n", im);
    % disp(num2str(F));
end

% ===FORMATTED DISPLAY===

DISPLACEMENT_DISPLAY(NJ, NCJT, NS, NSC, P, NDOF);
REACTION_FORCE_DISPLAY(NJ, NCJT, NR, NS, NSC, R, MSUP);

% fprintf("\n\nSTRUCTURE NUMBERING (NSC)\n\n");
% disp(NSC);

% ====== FUNCTION DECLARATION ======
function DISPLACEMENT_DISPLAY(NJ, NCJT, NS, NSC, P, NDOF)
    n = NJ * NCJT;
    temp_table_d = zeros(n, 2);
    b = 1;
    % PREPARING DISPLACEMENT DISPLAY
    for i=1:n
        temp_table_d(i, 1) = NSC(i, 1);
        if NSC(i, 1) == b && b <= NDOF
            temp_table_d(i, 2) = P(b, 1);
            b = b + 1;
        end
    end

    % SPLITTING DISPLACEMENT INTO X AND Y AXIS
    u1     = zeros(NJ, 1);
    u2     = zeros(NJ, 1);
    u3     = zeros(NJ, 1);
    u4     = zeros(NJ, 1);
    u5     = zeros(NJ, 1);
    u6     = zeros(NJ, 1);
    n_disp = zeros(NJ, 1);
    b = 1;
    for i=1:NJ
        u1(i) = temp_table_d(b  , 2);
        u2(i) = temp_table_d(b+1, 2);
        u3(i) = temp_table_d(b+2, 2);
        u4(i) = temp_table_d(b+3, 2);
        u5(i) = temp_table_d(b+4, 2);
        u6(i) = temp_table_d(b+5, 2);
        n_disp(i) = i;
        b = b + NCJT;
    end
    fprintf("\n\nDISPLACEMENT\n\n");
    disp_label = {'No Joint' 
                  'Translation_X' 
                  'Translation_Y' 
                  'Translation_Z' 
                  'Rotation_X' 
                  'Rotation_Y' 
                  'Rotation_Z'};

    disp(table(n_disp, u1, u2, u3, u4, u5, u6, 'VariableNames', disp_label))
end

function REACTION_FORCE_DISPLAY(NJ, NCJT, NR, NS, NSC, R, MSUP)
    temp_table_r = zeros(NS * 2, 2);
    e = NJ * NCJT - NR + 1;
    DOF = NJ * NCJT;
    temp_r = zeros(NR, 2);
    for i=e:DOF
        temp_r(i - e + 1, 1) = i;
        temp_r(i - e + 1, 2) = R(i - e + 1);
    end
    support_gnn = zeros(NS * 2, 1);
    for i=1:NS
        temp_table_r(i * NCJT, 1) = NSC(MSUP(i, 1) * NCJT);
        temp_table_r(i * NCJT - 1, 1) = NSC(MSUP(i, 1) * NCJT - 1);
        temp_table_r(i * NCJT - 2, 1) = NSC(MSUP(i, 1) * NCJT - 2);
        temp_table_r(i * NCJT - 3, 1) = NSC(MSUP(i, 1) * NCJT - 3);
        temp_table_r(i * NCJT - 4, 1) = NSC(MSUP(i, 1) * NCJT - 4);
        temp_table_r(i * NCJT - 5, 1) = NSC(MSUP(i, 1) * NCJT - 5);
    end
    % MATCH TEMP_R AND SUPPORT_GNN
    % disp(num2str(temp_r))
    for i=1:NS*NCJT
        for j=1:NR
            if temp_table_r(i, 1) == temp_r(j, 1)
                temp_table_r(i, 2) = temp_r(j, 2);
            end
        end
    end
    % fprintf("\n\nTEMP TABLE R\n\n");
    % disp(num2str(temp_table_r))
    % SPLIT TO X AND Y
    r1 = zeros(NS, 1);
    r2 = zeros(NS, 1);
    r3 = zeros(NS, 1);
    r4 = zeros(NS, 1);
    r5 = zeros(NS, 1);
    r6 = zeros(NS, 1);
    b = 1;
    for i=1:NS
        r1(i) = temp_table_r(b  , 2);
        r2(i) = temp_table_r(b+1, 2);
        r3(i) = temp_table_r(b+2, 2);
        r4(i) = temp_table_r(b+3, 2);
        r5(i) = temp_table_r(b+4, 2);
        r6(i) = temp_table_r(b+5, 2);
        b = b + NCJT;
    end
    fprintf("\n\nREACTION FORCE\n\n");
    react_label = {'No Joint' 
                   'Forces_X' 
                   'Forces_Y' 
                   'Forces_Z' 
                   'Moment_X' 
                   'Moment_Y' 
                   'Moment_Z'};

    disp(table(MSUP(:, 1), r1, r2, r3, r4, r5, r6, 'VariableNames', react_label))
end

function [GK] = MSTIFFG(NCJT, BK, T, GK)
    TS = zeros(2*NCJT, 2*NCJT);
    GK = zeros(2*NCJT, 2*NCJT);
    GK = transpose(T) * BK * T;
    % for i=1:2*NCJT
    %     for j=1:2*NCJT
    %         for k=1:2*NCJT
    %             TS(i, j) = TS(i, j) + BK(i, k) * T(k, j);
    %         end
    %     end
    % end

    % for i=1:2*NCJT
    %     for j=1:2*NCJT
    %         for k=1:2*NCJT
    %             GK(i, j) = GK(i, j) + T(k, i) * TS(k, j);
    %         end
    %     end
    % end

end

function [S] = STORES(JB, JE, NCJT, NDOF, NSC, GK, S)
    for i=1:2*NCJT
        i1 = (JE - 1) * NCJT + (i - NCJT);
        if i <= NCJT
            i1 = (JB - 1) * NCJT + i;
        end
        N1 = NSC(i1);
        if N1 <= NDOF
            for j=1:2*NCJT
                i1 = (JE - 1) * NCJT + (j - NCJT);
                if j <= NCJT
                    i1 = (JB - 1) * NCJT + j;
                end
                N2 = NSC(i1);
                if N2 <= NDOF
                    S(N1, N2) = S(N1, N2) + GK(i, j);
                end
            end
        end
    end
end

function [FF] = MFEFG(NCJT, T, QF, FF)
    FF = zeros(2*NCJT, 1);
    FF = FF + transpose(T) * QF;
    % for i=1:2*NCJT
    %     for j=1:2*NCJT
    %         FF(i) = FF(i) + T(j, i) * QF(j);
    %     end
    % end
    % fprintf("\nT\n");
    % disp(T)
    % fprintf("\nFF\n");
    % disp(FF)
end

function [QF] = MFEFLL(IML, BL, MP, PM, QF)
    % INIT ALL FEM ELEMENT TO ZERO
    FAB = 0; FSBY = 0; FSBZ = 0; FTB = 0; FMBY = 0; FMBZ = 0;
    FAE = 0; FSEY = 0; FSEZ = 0; FTE = 0; FMEY = 0; FMEZ = 0;
    LDTYPE = MP(IML, 2);
    LDDIR = MP(IML, 3);
    switch LDTYPE
    case 1
        BW  = PM(IML, 1);
        BL1 = PM(IML, 3);
        BL2 = BL - BL1;
        % CALCULATE FSB, FMB, FSE, FME
        if LDDIR == 2
            FSBY = BW * BL2^2 * (3*BL1 + BL2) / BL^3;
            FMBZ = BW * BL1 * BL2^2 / BL^2;
            FSEY = BW * BL2^2 * (BL1 + 3*BL2) / BL^3;
            FMEZ = -BW * BL1^2 * BL2 / BL^2;
        elseif LDDIR == 3
            FSBZ = BW * BL2^2 * (3*BL1 + BL2) / BL^3;
            FMBY = BW * BL1 * BL2^2 / BL^2;
            FSEZ = BW * BL2^2 * (BL1 + 3*BL2) / BL^3;
            FMEY = -BW * BL1^2 * BL2 / BL^2;
        end
    case 2
        BM  = PM(IML, 1);
        BL1 = PM(IML, 3); 
        BL2 = BL - BL1;
        if LDDIR == 2
            FSBZ = -6 * BM * BL1 * BL2 / BL^3;
            FMBY = BM * BL2 * (BL2 - 2*BL1) / BL^2;
            FSEZ = 6 * BM * BL1 * BL2 / BL^3;
            FMEY = BM * BL2 * (BL1 - 2*BL2) / BL^2;
        elseif LDDIR == 3
            FSBY = -6 * BM * BL1 * BL2 / BL^3;
            FMBZ = BM * BL2 * (BL2 - 2*BL1) / BL^2;
            FSEY = 6 * BM * BL1 * BL2 / BL^3;
            FMEZ = BM * BL2 * (BL1 - 2*BL2) / BL^2;
        end
    case 3
        W   = PM(IML, 1);
        BL1 = PM(IML, 3); BL2 = PM(IML, 4);
        if LDDIR == 2
            FSBY = (W * BL) / 2 * (1 - BL1 * (2*BL^3 - 2*BL*BL1^2 + BL1^3) / BL^4 - BL2^3 * (2*BL - BL2) / BL^4);
            FMBZ = (W * BL^2) / 12 * (1 - BL1^2 * (6*BL^2 - 8*BL1*BL + 3*BL1^2) / BL^4 - BL2^3 * (4*BL - 3*BL2) / BL^4); 
            FSEY = (W * BL) / 2 * (1 - BL2 * (2*BL^3 - 2*BL*BL2^2 + BL2^3) / BL^4 - BL1^3 * (2*BL - BL1) / BL^4);
            FMEZ = -(W * BL^2) / 12 * (1 - BL2^2 * (6*BL^2 - 8*BL2*BL + 3*BL2^2) / BL^4 - BL1^3 * (4*BL - 3*BL1) / BL^4); 
        elseif LDDIR == 3
            FSBZ = (W * BL) / 2 * (1 - BL1 * (2*BL^3 - 2*BL*BL1^2 + BL1^3) / BL^4 - BL2^3 * (2*BL - BL2) / BL^4);
            FMBY = (W * BL^2) / 12 * (1 - BL1^2 * (6*BL^2 - 8*BL1*BL + 3*BL1^2) / BL^4 - BL2^3 * (4*BL - 3*BL2) / BL^4); 
            FSEZ = (W * BL) / 2 * (1 - BL2 * (2*BL^3 - 2*BL*BL2^2 + BL2^3) / BL^4 - BL1^3 * (2*BL - BL1) / BL^4);
            FMEY = -(W * BL^2) / 12 * (1 - BL2^2 * (6*BL^2 - 8*BL2*BL + 3*BL2^2) / BL^4 - BL1^3 * (4*BL - 3*BL1) / BL^4); 
        end
    case 4
        W1  = PM(IML, 1); W2  = PM(IML, 2);
        BL1 = PM(IML, 3); BL2 = PM(IML, 4);
        if LDDIR == 2
            L = W1 * (BL - BL1)^3 / (20*BL^3);
            A = 7*BL + 8*BL1;
            B = BL2 * (3*BL + 2*BL1) / (BL - BL1);
            B1 = 1 + BL2 / (BL - BL1) + BL2^2 / (BL - BL1)^2;
            C = 2*BL2^4 / (BL - BL1)^3;
            R = W2 * (BL - BL1)^3 / (20*BL^3);
            D = 3*BL + 2*BL1;
            D1 = 1 + BL2 / (BL - BL1) + BL^2 / (BL - BL1)^2;
            E = BL2^3 / (BL - BL1)^2;
            E1 = 2 + (15*BL - 8*BL2) / (BL - BL1);
            FSBY = L * (A - B * B1 + C) + R * (D * D1 - E * E1);

            L = W1 * (BL - BL1)^3 / (60*BL^2);
            A = 3 * (BL + 4*BL1);
            B = BL2 * (2*BL + 3*BL1) / (BL - BL1);
            B1 = 1 + BL2 / (BL - BL1) + BL2^2 / (BL - BL1)^2;
            C = 3*BL2^4 / (BL - BL1)^3;
            R = W2 * (BL - BL1)^3 / (60*BL^2);
            D = 2*BL + 3*BL1;
            D1 = 1 + BL2 / (BL - BL1) + BL^2 / (BL - BL1)^2;
            E = 3*BL2^3 / (BL - BL1)^2;
            E1 = 1 + (5*BL - 4*BL2) / (BL - BL1);
            FMBZ = L * (A - B * B1 + C) + R * (D * D1 - E * E1);

            FSEY = (W1 + W2) / 2 * (BL - BL1 - BL2) - FSBY;
            FMEZ = (BL - BL1 - BL2) * (W1 * (-2*BL + 2*BL1 - BL2) - W2 * (BL - BL1 + 2*BL2)) / 6 + FSBY * BL - FMBZ;
        elseif LDDIR == 3
            L = W1 * (BL - BL1)^3 / (20*BL^3);
            A = 7*BL + 8*BL1;
            B = BL2 * (3*BL + 2*BL1) / (BL - BL1);
            B1 = 1 + BL2 / (BL - BL1) + BL2^2 / (BL - BL1)^2;
            C = 2*BL2^4 / (BL - BL1)^3;
            R = W2 * (BL - BL1)^3 / (20*BL^3);
            D = 3*BL + 2*BL1;
            D1 = 1 + BL2 / (BL - BL1) + BL^2 / (BL - BL1)^2;
            E = BL2^3 / (BL - BL1)^2;
            E1 = 2 + (15*BL - 8*BL2) / (BL - BL1);
            FSBZ = L * (A - B * B1 + C) + R * (D * D1 - E * E1);

            L = W1 * (BL - BL1)^3 / (60*BL^2);
            A = 3 * (BL + 4*BL1);
            B = BL2 * (2*BL + 3*BL1) / (BL - BL1);
            B1 = 1 + BL2 / (BL - BL1) + BL2^2 / (BL - BL1)^2;
            C = 3*BL2^4 / (BL - BL1)^3;
            R = W2 * (BL - BL1)^3 / (60*BL^2);
            D = 2*BL + 3*BL1;
            D1 = 1 + BL2 / (BL - BL1) + BL^2 / (BL - BL1)^2;
            E = 3*BL2^3 / (BL - BL1)^2;
            E1 = 1 + (5*BL - 4*BL2) / (BL - BL1);
            FMBY = L * (A - B * B1 + C) + R * (D * D1 - E * E1);

            FSEZ = (W1 + W2) / 2 * (BL - BL1 - BL2) - FSBY;
            FMEY= (BL - BL1 - BL2) * (W1 * (-2*BL + 2*BL1 - BL2) - W2 * (BL - BL1 + 2*BL2)) / 6 + FSBY * BL - FMBZ;
        end
    case 5
        BW  = PM(IML, 1);
        BL1 = PM(IML, 3);
        BL2 = BL - BL1;
        FAB = BW * BL2 / BL;
        FAE = BW * BL1 / BL;
    case 6
        W   = PM(IML, 1);
        BL1 = PM(IML, 3); BL2 = PM(IML, 4);
        FAB = W * (BL - BL1 - BL2) * (BL - BL1 + BL2) / (2*BL);
        FAE = W * (BL - BL1 - BL2) * (BL + BL1 - BL2) / (2*BL);
    end

    QF(1)  = QF(1)  + FAB;
    QF(2)  = QF(2)  + FSBY;
    QF(3)  = QF(3)  + FSBZ;
    QF(4)  = QF(4)  + FTB;
    QF(5)  = QF(5)  + FMBY;
    QF(6)  = QF(6)  + FMBZ;
    QF(7)  = QF(7)  + FAE;
    QF(8)  = QF(8)  + FSEY;
    QF(9)  = QF(9)  + FSEZ;
    QF(10) = QF(10) + FTE;
    QF(11) = QF(11) + FMEY;
    QF(12) = QF(12) + FMEZ;

end

function [V] = MDISPG(JB, JE, NCJT, NDOF, NSC, P, V)
    V = zeros(2*NCJT);
    J = (JB - 1) * NCJT;
    for i=1:NCJT
        J = J + 1;
        N = NSC(J);
        if N <= NDOF
            V(i) = P(N);
        end
    end
    J = (JE - 1) * NCJT;
    for i=(NCJT+1):2*NCJT
        J = J + 1;
        N = NSC(J);
        if N <= NDOF
            V(i) = P(N);
        end
    end
end

function [M] = INIT_COLUMN_TRANS(AOR, RXX, RXY, R_T)
    M(1, 1) = 0;
    M(1, 2) = RXY;
    M(1, 3) = 0; 
    M(2, 1) = -RXY * cosd(AOR);
    M(2, 2) = 0;
    M(2, 3) = sind(AOR);
    M(3, 1) = RXY * sind(AOR);
    M(3, 2) = 0;
    M(3, 3) = cosd(AOR);
end

function [M] = INIT_NORMAL_TRANS(RXX, RXY, RXZ, RYX, RYY, RYZ, RZX, RZY, RZZ, M)
    M(1, 1) = RXX;
    M(1, 2) = RXY;
    M(1, 3) = RXZ; 
    M(2, 1) = RYX;
    M(2, 2) = RYY;
    M(2, 3) = RYZ;
    M(3, 1) = RZX;
    M(3, 2) = RZY;
    M(3, 3) = RZZ;
end

function [T] = MTRANS(AOR, RXX, RXY, RXZ, RYX, RYY, RYZ, RZX, RZY, RZZ, NCJT, T)
    T = zeros(2*NCJT, 2*NCJT);
    R_T = zeros(3, 3);
    if RXX == 0 && RXZ == 0
        % fprintf("\nCOLUMN MEMBER\n")
        R_T = INIT_COLUMN_TRANS(AOR, RXX, RXY, R_T);
    else
        % fprintf("\nBEAM MEMBER\n")
        R_T = INIT_NORMAL_TRANS(RXX, RXY, RXZ, RYX, RYY, RYZ, RZX, RZY, RZZ, R_T);
    end

    T(1:3, 1:3)     = R_T(1:3, 1:3);
    T(4:6, 4:6)     = R_T(1:3, 1:3);
    T(7:9, 7:9)     = R_T(1:3, 1:3);
    T(10:12, 10:12) = R_T(1:3, 1:3);

end

function [U] = MDISPL(NCJT, V, T, U)
    U = zeros(2*NCJT);
    U = U + T * V;
    % for i=1:2*NCJT
    %     for j=1:2*NCJT
    %         U(i) = U(i) + T(i, j) * V(j);
    %     end
    % end
end

function [M] = INIT_TOP_LEFT_STIFF(E, G, A, ZI, YI, J, BL, NCJT)
    M = zeros(NCJT, NCJT);
    M(1, 1) = E * A / BL;
    M(2, 2) = 12 * E * ZI / (BL^3);
    M(3, 3) = 12 * E * YI / (BL^3);
    M(4, 4) = G * J / BL;
    M(5, 5) = 4 * E * YI / BL;
    M(6, 6) = 4 * E * ZI / BL;
    % LEFT DIAGONAL
    M(6, 2) = 6 * E * ZI / (BL^2);
    M(5, 3) = -6 * E * YI / (BL^2);
    M(3, 5) = -6 * E * YI / (BL^2);
    M(2, 6) = 6 * E * ZI / (BL^2);
end

function [M] = INIT_TOP_RIGHT_STIFF(E, G, A, ZI, YI, J, BL, NCJT)
    M = zeros(NCJT, NCJT);
    M(1, 1) = -E * A / BL;
    M(2, 2) = -12 * E * ZI / (BL^3);
    M(3, 3) = -12 * E * YI / (BL^3);
    M(4, 4) = -G * J / BL;
    M(5, 5) = 2 * E * YI / BL;
    M(6, 6) = 2 * E * ZI / BL;
    % LEFT DIAGONAL
    M(6, 2) = -6 * E * ZI / (BL^2);
    M(5, 3) = 6 * E * YI / (BL^2);
    M(3, 5) = -6 * E * YI / (BL^2);
    M(2, 6) = 6 * E * ZI / (BL^2);
end

function [M] = INIT_BOTTOM_RIGHT_STIFF(E, G, A, ZI, YI, J, BL, NCJT)
    M = zeros(NCJT, NCJT);
    M(1, 1) = E * A / BL;
    M(2, 2) = 12 * E * ZI / (BL^3);
    M(3, 3) = 12 * E * YI / (BL^3);
    M(4, 4) = G * J / BL;
    M(5, 5) = 4 * E * YI / BL;
    M(6, 6) = 4 * E * ZI / BL;
    % LEFT DIAGONAL
    M(6, 2) = -6 * E * ZI / (BL^2);
    M(5, 3) = 6 * E * YI / (BL^2);
    M(3, 5) = 6 * E * YI / (BL^2);
    M(2, 6) = -6 * E * ZI / (BL^2);
end

function [BK] = MSTIFFL(E, G, A, ZI, YI, J, BL, NCJT, BK)
    BK = zeros(2*NCJT, 2*NCJT);
    TOP_LEFT     = zeros(NCJT, NCJT);
    TOP_RIGHT    = zeros(NCJT, NCJT);
    BOTTOM_LEFT  = zeros(NCJT, NCJT);
    BOTTOM_RIGHT = zeros(NCJT, NCJT);

    TOP_LEFT     = INIT_TOP_LEFT_STIFF(E, G, A, ZI, YI, J, BL, NCJT);
    TOP_RIGHT    = INIT_TOP_RIGHT_STIFF(E, G, A, ZI, YI, J, BL, NCJT);
    BOTTOM_LEFT  = transpose(TOP_RIGHT);
    BOTTOM_RIGHT = INIT_BOTTOM_RIGHT_STIFF(E, G, A, ZI, YI, J, BL, NCJT);

    BK(1:6, 1:6)   = BK(1:6, 1:6)   + TOP_LEFT;
    BK(1:6, 7:12)  = BK(1:6, 7:12)  + TOP_RIGHT;
    BK(7:12, 7:12) = BK(7:12, 7:12) + BOTTOM_RIGHT;
    BK(7:12, 1:6)  = BK(7:12, 1:6)  + BOTTOM_LEFT;

    % fprintf("\nBK\n");
    % disp(BK);

end

function [Q] = MFORCEL(NCJT, BK, U, Q, QF)
    Q = QF;
    % Q = Q + BK * U;
    for i=1:2*NCJT
        for j=1:2*NCJT
            Q(i) = Q(i) + BK(i, j) * U(j);
        end
    end
end

function [F] = MFORCEG(NCJT, T, Q, F)
    F = zeros(2*NCJT, 1);
    % F = F + T * Q;
    for i=1:2*NCJT
        for j=1:2*NCJT
            F(i) = F(i) + T(j, i) * Q(j);
        end
    end
end

function [R] = STORER(JB, JE, NCJT, NDOF, NSC, F, R)
    for i=1:2*NCJT
        i1 = (JE - 1) * NCJT + (i - NCJT);
        if i <= NCJT
            i1 = (JB - 1) * NCJT + i;
        end
        N = NSC(i1);
        if N > NDOF
            R(N - NDOF) = R(N - NDOF) + F(i);
        end
    end
end

function [P] = STOREPF(JB, JE, NCJT, NDOF, NSC, FF, P)
    for i=1:2*NCJT
        i1 = (JE - 1) * NCJT + (i - NCJT);
        if i <= NCJT
            i1 = (JB - 1) * NCJT + i;
        end
        N1 = NSC(i1);
        if N1 <= NDOF
            P(N1) = P(N1) - FF(i);
        end
    end
end

