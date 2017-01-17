

fbmodel.NB = 5;

fbmodel.jtype = {'fb', 'Rx', 'Rz', 'Rx', 'Rz'};

fbmodel.parent = [0 1 2 1 4];

X2 = plux(eye(3), [0 -.2 0]) ;
X3 = plux(eye(3), [0  -1 0]) ;
X4 = plux(eye(3), [0  .2 0]) ;
X5 = plux(eye(3), [0   1 0]) ;

fbmodel.Xtree = {eye(6), X2, X3, X4, X5};
Ib = [eye(3)*.2 zeros(3,3) ; zeros(3,3) 2*eye(3)];
hl = [0 -.5 0];
Il = [eye(3)*.1+skew(hl)*skew(hl)' skew(hl) ; skew(hl)' eye(3)];
fbmodel.I = {Ib, Il, Il, Il, Il};
qt = rand(1,4);
qt = qt/norm(qt);
q = { [qt 0 0 0]', -3, 2, 1, 4};
qd= { rand(6,1)*50, 50 , 20 , -4.430, 30.22 };
dt = 1e-9;
for i = 1:length(q)
    qnew{i} = q{i} + jdcalc(fbmodel.jtype{i}, q{i}, qd{i}) * dt;
end
[C1 H1 Hd] = Cmat(fbmodel, q, qd)
C1 = cell2mat(C1)
H1 = cell2mat(H1)
Hd = cell2mat(Hd);

[C2 H2] = Cmat(fbmodel, qnew, qd)
H2 = cell2mat(H2);
Hd2 = (H2 - H1) / dt;
Hd - C1 -C1'
Hd - Hd2
