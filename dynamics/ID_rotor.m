function  tau = ID_rotor( model, q, qd, qdd, f_ext )

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

a_grav = get_gravity(model);

% Xuprot   = repmat({zeros(6,6)},model.NB,1);
% Xup   = repmat({zeros(6,6)},model.NB,1);
% 
% f      = repmat({zeros(6,1)},model.NB,1);
% frot   = repmat({zeros(6,1)},model.NB,1);
% 
% v   = repmat({zeros(6,1)},model.NB,1);
% a   = repmat({zeros(6,1)},model.NB,1);
% vrot   = repmat({zeros(6,1)},model.NB,1);
% arot   = repmat({zeros(6,1)},model.NB,1);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  [ XJrot, ~ ] = jcalc( model.jtype{i}, q{i}*model.gr{i} );
  if strcmp('fb',model.jtype{i})
      XJrot = eye(6);
  end
  
  Srot{i} = S{i} * model.gr{i};
  
  vJ = S{i}*qd{i};
  vJrot = vJ * model.gr{i};
  
  Xup{i} = XJ * model.Xtree{i};
  Xuprot{i} = XJrot * model.Xrotor{i};
  
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd{i};
    vrot{i} = vJrot;
    arot{i} = Xuprot{i}*(-a_grav) + Srot{i}*qdd{i};
    
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd{i} + crm(v{i})*vJ;
    vrot{i} = Xuprot{i}*v{model.parent(i)} + vJrot;
    arot{i} = Xuprot{i}*a{model.parent(i)} + Srot{i}*qdd{i} + crm(vrot{i})*vJrot;
    
  end
  f{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
  frot{i} = model.Irot{i}*arot{i} + crf(vrot{i})*model.Irot{i}*vrot{i};
end

% if nargin == 5
%   f = apply_external_forces( model.parent, Xup, f, f_ext );
% end
tau = cell(length(qd),1);
for i = model.NB:-1:1
  tau{i} = S{i}' * f{i} + Srot{i}' * frot{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i} + Xuprot{i}'*frot{i};
  end
end
