
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


fbmodel_nograv = fbmodel;
fbmodel_nograv.a_grav = [0 0 0 0 0 0]';


a = cell(fbmodel.NB,1);
for i = 1:fbmodel.NB
    a{i} = inertiaMatToVec( fbmodel.I{i} );
end
a = cell2mat(a);





qt = rand(1,4);
qt = qt/norm(qt);
q = { [qt 0 0 0]', -3, 2, 1, 4};
qd= { rand(6,1)*5, 5 , 2 , -4.430, 3.22 }';
qdr = { rand(6,1)*5, rand(1) , rand(1) , rand(1), rand(1) };
qddr= { rand(6,1)*5, rand(1) , rand(1) , rand(1), rand(1) };

zqd = {zeros(6,1) ; 0 ; 0; 0; 0};

dt = 1e-5;
t = 0;
lambda = 1e-3;

W1 = zeros( length(cell2mat(qd)),fbmodel.NB*10); % Compute explicitly from filtered Y
W23 = zeros( length(cell2mat(qd)), fbmodel.NB*10);
qdd = { rand(6,1), rand(1), rand(1), rand(1), rand(1)}';

tauf = zeros( length(cell2mat(qd)) ,1 );
tauf2 = tauf;


[W0] = SlotineLiY(fbmodel_nograv,qd, zqd, zqd, q);
[p0, ~] = MandCTqd(fbmodel,q,qd);


kk = 0;
while t < 1000*dt
    Y = SlotineLiY(fbmodel, qdd, qd, qd, q);
    Yf= MandCTY(fbmodel, qd,q,lambda);
    Yg =SlotineLiY(fbmodel, zqd, zqd, zqd, q);
    
    [Mqd, CTqd] = MandCTqd(fbmodel,q,qd);
    g = SlotineLiID(fbmodel, zqd, zqd, zqd, q);
    
    
    err = -lambda*Mqd - CTqd +g - (Yf+Yg)*a;
    
    tau = SlotineLiID(fbmodel, qdd, qd, qd, q);
    
    tauf = tauf + dt*lambda*(tau - tauf);    
    tauf2= tauf2+ dt*lambda*(-lambda*Mqd - CTqd +g - tauf2);
    
    W1 = W1 + dt*lambda * (Y - W1) ;
    W23 = W23 + dt*lambda*(Yf + Yg - W23);    
    %return
    for i = 1:length(q)
        qnew{i} = q{i} + jdcalc(fbmodel.jtype{i}, q{i}, qd{i}) * dt;
        qd{i}   = qd{i} + dt * qdd{i};
    end    
    q = qnew;
    q{1}(1:4) = ( q{1}(1:4) ) / norm( q{1}(1:4) );
    t = t+dt;
    t
end

[pT, ~] = MandCTqd(fbmodel,q,qd);
[WT] = SlotineLiY(fbmodel_nograv,qd, zqd, zqd, q);

wt = lambda * exp(-lambda* t);
w0 = lambda;

W2 = WT*w0 - W0*wt + W23;


errPt = norm( pT - WT*a)
errP0 = norm( p0 - W0*a)
errW23 = norm( tauf2 - W23*a)
errW1  = norm( tauf  - W1*a)
errW2  = norm( tauf  - W2*a)





    
    
    





