function  Y  = MandCTYb( model, qd,q , lambda, factorFunction)
% Indirect regressor term Yb s.t. Yb a = Sb^T ( lambda M qd + C^T qd )
%

if nargin == 4
    factorFunction = @(I,v)(factorFunctions(I,v, 1));
end

Y = cell(1, model.NB);
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    X1{i} = eye(6);
  else
    X1{i}= Xup{i} * X1{ model.parent(i) };
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end
  Sd{i} = crm(v{i})*S{i};

  h{i} = zeros(6,10);
  b{i} = zeros(6,10);
  for k = 1:10
     ak = zeros(10,1); ak(k) = 1;
     Ik = inertiaVecToMat(ak);
     h{i}(:,k) = Ik * v{i};
  end  
  Y{1,i} =  (lambda * eye(6) - crf(v{1}) ) *X1{i}' * h{i};
end
Y = cell2mat( Y );
