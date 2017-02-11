function  [p, w]  = pwTerms( model, qd,q , lambda)
% The entries of p, and w, denoted p_ik and w_ik are computed as in the 
%   algorithm of Wensing, Ugurlu, and Slotine IROS 2017 (hopefully)

a_grav = get_gravity(model);
p = cell(1, model.NB);
w = cell(1, model.NB);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    X1{i} = eye(6);
    ag1 = Xup{i} * a_grav;
  else
    X1{i}= Xup{i} * X1{ model.parent(i) };
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end

  p{i} = zeros(6,10);
  w{i} = zeros(6,10);
  
  for k = 1:10
     ak = zeros(10,1); ak(k) = 1;
     Ik = inertiaVecToMat(ak);
     p{i}(:,k) = X1{i}'* Ik * v{i};
     w{i}(:,k) = (lambda*eye(6) - crf(v{1}) ) * p{i}(:,k) + X1{i}' * Ik * X1{i} * ag1;
  end  
end
p = cell2mat(p);
w = cell2mat(w);

