
%% Set up a simple model
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

% Create a corresponding model without any gravity
fbmodel_nograv = fbmodel;
fbmodel_nograv.a_grav = [0 0 0 0 0 0]';

a = cell(fbmodel.NB,1);
for i = 1:fbmodel.NB
    a{i} = inertiaMatToVec( fbmodel.I{i} );
end
a = cell2mat(a);

%% Set initial state
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

% Inialize regressors
W = zeros( length(cell2mat(qd)),fbmodel.NB*10);          % Compute explicitly from filtered Y
U = zeros( 6, fbmodel.NB*10);                            % Computed based on new equations
V = zeros( 6, fbmodel.NB*10);

qdd = { rand(6,1), rand(1), rand(1), rand(1), rand(1)}';

tauf = zeros( length(cell2mat(qd)) ,1 );
tauf2 = tauf;

% get initial matrix with entries p_ik
[p0, ~] = pwTerms(fbmodel,qd,q, lambda);

kk = 0;
while t < 100*dt
    % Compute Y directly and filter
    Y = SlotineLiY(fbmodel, qdd, qd, qd, q);
    W = W + dt*lambda * (Y - W) ;
    
    % Alternately, compute p_ik and w_ik terms
    [p w] = pwTerms(fbmodel, qd,q,lambda);

    % w * a = (the first six rows of) "lambda*Mqd + CTqd -g"
    [Mqd, CTqd] = MandCTqd(fbmodel,q,qd);
    g = SlotineLiID(fbmodel, zqd, zqd, zqd, q);
    
    % form temp quantity
    tt = lambda*Mqd + CTqd -g;
    
    % Check error
    % err = tt(1:6) - w*a
    
    % Compute full torque 
    tau = SlotineLiID(fbmodel, qdd, qd, qd, q);
    
    % Filter it
    tauf = tauf + dt*lambda*(tau - tauf);    
    
    % Filter this other quanitity as well
    tauf2= tauf2+ dt*lambda*(lambda*Mqd + CTqd -g - tauf2);
    
    % Filter w_ik to get V_ik
    V = V + dt*lambda*(w - V);   
    
    % integrate
    for i = 1:length(q)
        qnew{i} = q{i} + jdcalc(fbmodel.jtype{i}, q{i}, qd{i}) * dt;
        qd{i}   = qd{i} + dt * qdd{i};
    end    
    q = qnew;
    q{1}(1:4) = ( q{1}(1:4) ) / norm( q{1}(1:4) );
    t = t+dt;
end

% Compute final matrix with entries p_ik
[pT, ~] = pwTerms(fbmodel, qd, q, lambda);

% Compute U from impulse response of filter
imp_response_t = lambda * exp(-lambda* t);
imp_response_0 = lambda;
U = imp_response_0*pT - imp_response_t*p0;

% Check Our computation
Wbtest = U - V;
filteredBaseRegressorErr = max( max( abs ( Wbtest - W(1:6,:) ) ) )
