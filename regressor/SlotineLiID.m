function  tau = SlotineLiID( model, qddr, qdr,qd,q)

if nargin == 5
    factorFunction = @(I,v)(factorFunctions(I,v, 3));
end

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  vJr= S{i}*qdr{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    w{i} = vJr;
    wd{i}= -Xup{i}*model.a_grav + S{i} * qddr{i} + crm(v{i})*vJr;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    w{i} = Xup{i}*w{model.parent(i)} + vJr;
    wd{i}= Xup{i}*wd{model.parent(i)} + S{i} * qddr{i} + crm(v{i})*vJr;
  end
  f{i} = model.I{i}*wd{i} + factorFunction(model.I{i}, v{i})*w{i};
end

tau = cell(model.NB,1);
for i = model.NB:-1:1    
  if model.parent(i) ~= 0
     p = model.parent(i);
     f{p} = f{p} + Xup{i}'*f{i};
  end
  tau{i} = S{i}'*f{i};
end

tau = cell2mat(tau);
