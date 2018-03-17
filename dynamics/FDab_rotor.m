function  [qdd, v] = FDab_rotor( model, q, qd, tau, f_ext )

% FDab  Forward Dynamics via Articulated-Body Algorithm
% FDab(model,q,qd,tau,f_ext,grav_accn)  calculates the forward dynamics of
% a kinematic tree via the articulated-body algorithm.  q, qd and tau are
% vectors of joint position, velocity and force variables; and the return
% value is a vector of joint acceleration variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.
a_grav = get_gravity(model);
NJ = 0;
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  
  [ XJrot, ~ ] = jcalc( model.jtype{i}, q{i}*model.gr{i} );
  
  Srot{i} = S{i} * model.gr{i};
  
  vJ = S{i}*qd{i};
  vJrot = Srot{i}*qd{i} ;  
  
  Xup{i} = XJ * model.Xtree{i};
  Xuprot{i} = XJrot * model.Xrotor{i};
  
  
  
  
  if model.parent(i) == 0
    v{i} = vJ;
    vrot{i} = vJrot;
    
    c{i} = zeros(size(a_grav));		% spatial or planar zero vector
    crot{i} = zeros(size(a_grav));
    
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    vrot{i} = Xuprot{i}*v{model.parent(i)} + vJrot;
    c{i} = crm(v{i}) * vJ;
    crot{i} = crm(vrot{i}) * vJrot;
  end
  
  IA{i} = model.I{i};
  pA{i} = crf(v{i}) * model.I{i} * v{i};
  pArot{i} = crf(vrot{i}) * model.Irot{i} * vrot{i};
  
  NJ = NJ + length(qd{i});
end

if nargin == 5
  pA = apply_external_forces( model.parent, Xup, pA, f_ext );
end

for i = model.NB:-1:1
  U{i}  = IA{i} * S{i};
  Ur{i} = model.Irot{i} * Srot{i};
  d{i}  = S{i}' * U{i} + Srot{i}'*Ur{i};
  
  u{i} = tau{i} - S{i}'*pA{i} - Srot{i}'*pArot{i} - U{i}'*c{i} - Ur{i}'*crot{i};
  Utot{i} = Xup{i}' * U{i} + Xuprot{i}' * Ur{i};
  
  if model.parent(i) ~= 0
      
    Ia = Xup{i}' * IA{i} * Xup{i} + Xuprot{i}'*model.Irot{i}*Xuprot{i};  
    Ia = Ia - Utot{i}/d{i}*Utot{i}';
    
    pa = Xup{i}'*(pA{i} + IA{i}*c{i}) + Xuprot{i}'*(pArot{i} + model.Irot{i} * crot{i});    
    pa = pa + Utot{i} * (d{i}\u{i});
    
    IA{model.parent(i)} = IA{model.parent(i)} + Ia;
    pA{model.parent(i)} = pA{model.parent(i)} + pa;
  end
end

qdd = cell(model.NB,1);
for i = 1:model.NB
  if model.parent(i) == 0
     ap= -a_grav; 
  else
     ap = a{model.parent(i)};
  end
  qdd{i} = d{i}\(u{i} - Utot{i}'*ap);
  a{i} = Xup{i} * ap + S{i}*qdd{i} + c{i};
end
