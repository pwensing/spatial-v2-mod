function  [Mqd, CTqd]  = MandCTqd( model, q , qd, factorFunction)

if nargin == 3
    factorFunction = @(I,v)(factorFunctions(I,v, 3));
end

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    vd{i} = crm(v{i}) * vJ;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    vd{i}= Xup{i}*vd{model.parent(i)} + crm(v{i})*vJ;
  end
  Sd{i} = crm(v{i})*S{i};
  B{i} = factorFunction(model.I{i}, v{i});
  h{i} = model.I{i}*v{i};
  b{i} = B{i}' * v{i};
end

Mqd = cell(model.NB, 1);
CTqd= cell(model.NB, 1);

for i = model.NB:-1:1    
  if model.parent(i) ~= 0
     p = model.parent(i);
     h{p} = h{p} + Xup{i}'*h{i};
     b{p} = b{p} + Xup{i}'*b{i};
  end
  Mqd{i} = S{i}' * h{i};
  CTqd{i}= Sd{i}'* h{i} + S{i}' * b{i};
end

Mqd = cell2mat(Mqd);
CTqd = cell2mat(CTqd);
