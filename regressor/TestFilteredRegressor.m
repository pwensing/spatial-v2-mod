
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
qd= { rand(6,1)*5, 5 , 2 , -4.430, 3.22 }';
qdr = { rand(6,1)*5, rand(1) , rand(1) , rand(1), rand(1) };
qddr= { rand(6,1)*5, rand(1) , rand(1) , rand(1), rand(1) };


a = cell(fbmodel.NB,1);
for i = 1:fbmodel.NB
    a{i} = inertiaMatToVec( fbmodel.I{i} );
end
a = cell2mat(a);

[C1 H1 Hd] = Cmat(fbmodel, q, qd);
C1 = cell2mat(C1);
H1 = cell2mat(H1);

[Mqd, CTqd] = MandCTqd(fbmodel,q,qd);


qdv = cell2mat(qd)

eMqd = norm( Mqd - H1*qdv)
eCTqd = C1'*qdv - CTqd;

lambda = .2;

Y =  MandCTY(fbmodel, qd, q, lambda) ;


Yb = Y(1:6,:);
Yb2 =  MandCTYb(fbmodel, qd, q, lambda) ;

eYb = norm( Yb2 + Yb )


eY = norm( -lambda*Mqd-CTqd-Y*a)




