function  [Y, Yrot, data] = SlotineLiY_rotor( model, qddr, qdr,qd,q, factorFunction)

if nargin == 5
    factorFunction = @(I,v)(factorFunctions(I,v, 3));
end
Xup = repmat({zeros(6,6)},model.NB,1);
Xrot = repmat({zeros(6,6)},model.NB,1);
X0 = repmat({zeros(6,6)},model.NB,1);
X0rot = repmat({zeros(6,6)},model.NB,1);

v   = repmat({zeros(6,1)},model.NB,1);
w   = repmat({zeros(6,1)},model.NB,1);
wd  = repmat({zeros(6,1)},model.NB,1);
vrot   = repmat({zeros(6,1)},model.NB,1);
wrot   = repmat({zeros(6,1)},model.NB,1);
wdrot  = repmat({zeros(6,1)},model.NB,1);
F  = repmat({zeros(6,10)},model.NB,1);

Frot  = repmat({zeros(6,10)},model.NB,1);

a_grav = get_gravity(model);
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  [ XJrot, ~] = jcalc( model.jtype{i}, q{i}*model.gr{i} );
  
  vJ = S{i}*qd{i};
  vJr= S{i}*qdr{i};
  
  Srot{i} = S{i} * model.gr{i};
 
  vJrrot= Srot{i}*qdr{i};
  vJrot = Srot{i}*qd{i};
  
  
  Xup{i} = XJ * model.Xtree{i};
  Xrot{i} = XJrot * model.Xrotor{i};
  
  if model.parent(i) == 0
    X0{i}    = Xup{i};
    X0rot{i} = Xrot{i};
    
    
    
    
    vrot{i}  = vJrot;
    wrot{i}  = vJrrot;
    wdrot{i} = -Xrot{i}*a_grav + Srot{i} * qddr{i} + crm(vrot{i})*vJrrot;
    
    v{i} = vJ;
    w{i} = vJr;
    wd{i}= -Xup{i}*a_grav + S{i} * qddr{i} + crm(v{i})*vJr;
  else
      
    X0{i}    = Xup{i} *  X0{model.parent(i)};
    X0rot{i} = Xrot{i} * X0{model.parent(i)};
    
    vrot{i}  = Xrot{i}*v{model.parent(i)} + vJrot;
    wrot{i}  = Xrot{i}*w{model.parent(i)} + vJrrot;
    wdrot{i} = Xrot{i}*wd{model.parent(i)} + Srot{i} * qddr{i}+ crm(vrot{i})*vJrrot;
    
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    w{i} = Xup{i}*w{model.parent(i)} + vJr;
    wd{i}= Xup{i}*wd{model.parent(i)} + S{i} * qddr{i} + crm(v{i})*vJr;
  end
  for k = 1:10
     ak = zeros(10,1); ak(k) = 1;
     Ik = inertiaVecToMat(ak);
     F{i}(:,k)    = Ik * wd{i}    + factorFunction(Ik, v{i})    * w{i};
     Frot{i}(:,k) = Ik * wdrot{i} + factorFunction(Ik, vrot{i}) * wrot{i};
  end
end

Y    = repmat( {zeros(1,10)} , model.NB, model.NB);
Yrot = repmat( {zeros(1,10)} , model.NB, model.NB);

for i = model.NB:-1:1
  Y{i,i} = S{i}' * F{i};
  Yrot{i,i} = Srot{i}' * Frot{i};
  
  Fji = F{i};
  Frotji = Frot{i};
  
  Fji    = Xup{i}' * Fji;
  Frotji = Xrot{i}' * Frotji;
  
  j = i;
  j = model.parent(j);
  while j ~= 0
      Y{j,i} = S{j}' * Fji;
      Yrot{j,i} = S{j}' * Frotji;
      
      Fji = Xup{j}' * Fji;
      Frotji = Xup{j}' * Frotji;
      j = model.parent(j);
  end
end

for i = 1:model.NB
    for j = 1:model.NB
        if size(Y{i,j},1) == 0
           Y{i,j} = zeros(length(qd{i}), 10);
           Yrot{i,j} = zeros(length(qd{i}), 10);
        end
    end
end
Yrot{1,1} = zeros(6,10);

data.X0 = X0;
data.X0rot = X0rot;

Y    = cell2mat(Y);
Yrot = cell2mat(Yrot);
