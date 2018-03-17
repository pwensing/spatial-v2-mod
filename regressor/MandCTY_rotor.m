function  Y  = MCTY_rotor( model, qd,q , lambda, factorFunction)
% Indirect regressor term Y s.t. Y a = -lambda M qd - C^T qd
%

if nargin == 4
    factorFunction = @(I,v)(factorFunctions(I,v, 1));
end

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    vd{i} = crm(v{i}) * vJ;
    
    vrot{i} = vJ * model.gr{i};
    vdrot{i} = crm(v{i}) * vJ * model.gr{i};
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    vd{i}= Xup{i}*vd{model.parent(i)} + crm(v{i})*vJ;
    
    vrot{i} = Xup{i}*v{model.parent(i)} + vJ*model.gr{i};
    vdrot{i}= Xup{i}*vd{model.parent(i)} + crm(v{i})*vJ*model.gr{i};
  end
  Sd{i} = crm(v{i})*S{i};

  h{i} = zeros(6,10);
  b{i} = zeros(6,10);
  for k = 1:10
     ak = zeros(10,1); ak(k) = 1;
     Ik = inertiaVecToMat(ak);
     h{i}(:,k) = Ik * v{i};
     b{i}(:,k) = factorFunction(Ik, v{i})' * v{i};
  end  
end
Y = cell(model.NB, model.NB);
for i = model.NB:-1:1
  Y{i,i} = -lambda* S{i}' * h{i} -  Sd{i}'* h{i} - S{i}' * b{i};
  hji = h{i};
  bji = b{i};
  j = i;
  while model.parent(j) ~= 0
      hji = Xup{j}' * hji;
      bji = Xup{j}' * bji;
      j = model.parent(j);
      Y{j,i} = -lambda* S{j}' * hji -  Sd{j}'* hji - S{j}' * bji;
  end  
end

for i = 1:model.NB
    for j = 1:model.NB
        if size(Y{i,j},1) == 0
           Y{i,j} = zeros(length(qd{i}), 10);
        end
    end
end
Y = cell2mat(Y);
