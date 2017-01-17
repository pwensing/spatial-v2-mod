
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
fbmodel.a_grav = [ 0 0 0 0 0 -9.81]';

qt = rand(1,4);
qt = qt/norm(qt);
q = { [qt 0 0 0]', -3, 2, 1, 4};
qd= { rand(6,1)*5, 5 , 2 , -4.430, 3.22 };
qdr = { rand(6,1)*5, rand(1) , rand(1) , rand(1), rand(1) };
qddr= { rand(6,1)*5, rand(1) , rand(1) , rand(1), rand(1) };

dt = 1e-9;

a = cell(fbmodel.NB,1);
for i = 1:fbmodel.NB
    a{i} = inertiaMatToVec( fbmodel.I{i} );
end
a = cell2mat(a);



tau = SlotineLiID(fbmodel, qddr, qdr, qd, q);
tau = cell2mat(tau);

Y   = SlotineLiY(fbmodel, qddr, qdr, qd,q);
Y = cell2mat(Y);

tau2 = Y*a
tau - tau2